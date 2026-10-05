"""Independent replay and declared-window diagnostics; no model calls."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann,electron_volt
from scipy.spatial.distance import pdist
from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_thermal_identity import canonical_temperature_terms
from solver_v1.silicon_thermal_research import block_statistics,split_rhat
from results.silicon_wafer_feasibility.audit_force_thermal import replay

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def arrays(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def checkpoint(path):
    z=arrays(path)
    return dict(u=z['u'],phi=float(z['phi']),gradient=z['gradient'],observation=z['observation'],
        completed=int(z['completed']),retained_evaluation=int(z['retained_evaluation']),
        rng_state=json.loads(str(z['rng_json']))),json.loads(str(z['transaction_json']))

def replay_extension(folder,parent_folder,source,target,domain,protocol):
    files=sorted(folder.glob('evaluation_*.npz'));values=[arrays(p) for p in files]
    checkpoints=sorted(folder.glob('checkpoint_*.npz'))
    draws=protocol['additional_draws_per_chain']
    if len(checkpoints)!=draws+1:raise ValueError('complete checkpoint sequence required')
    kT=Boltzmann/electron_volt*protocol['temperature_K'];root_error=dict(coordinate=0.,phi=0.,score=0.,path=0.,delta_h=0.)
    for i,(path,v) in enumerate(zip(files,values)):
        if path.name!=f'evaluation_{i:04d}.npz':raise ValueError('missing evaluation')
        r=source['positions']+(source['basis']@(target.reference.mean+
            np.sqrt(kT)*np.linalg.solve(target.reference.lower.T,v['u']))).reshape(source['positions'].shape)
        root_error['coordinate']=max(root_error['coordinate'],float(abs(r-v['positions']).max()))
        if target.inside(r)!=bool(v['inside']):raise ValueError('domain changed')
        if bool(v['inside']):
            phi=(float(v['energy'])-float(source['energy']))/kT-.5*v['u']@v['u']
            score=np.linalg.solve(target.reference.lower,-source['basis'].T@v['forces'].ravel())/np.sqrt(kT)-v['u']
            observation=np.array([float(v['energy']),float((r-domain['center']).ravel()@domain['direction']),
                float(np.sqrt(np.mean(np.sum((r[target.free]-domain['center'][target.free])**2,axis=1))))])
            root_error['phi']=max(root_error['phi'],abs(phi-float(v['correction'])))
            root_error['score']=max(root_error['score'],float(abs(score-v['gradient']).max()))
            np.testing.assert_allclose(v['observation'],observation,rtol=0,atol=1e-11)
        elif not np.isposinf(v['correction']):raise ValueError('outside finite target')
    old=arrays(parent_folder/'chain.npz');parent_summary=read(parent_folder/'summary.json')
    parent_index=int(old['retained_evaluation'][-1]);parent_record=arrays(parent_folder/f'evaluation_{parent_index:04d}.npz')
    current,tx0=checkpoint(checkpoints[0])
    if tx0 is not None or current['completed']!=0 or current['retained_evaluation']!=0:
        raise ValueError('invalid initial checkpoint')
    for key,left,right in [('u',current['u'],old['samples'][-1]),('score',current['gradient'],parent_record['gradient']),
        ('phi',current['phi'],old['potentials'][-1]),('observation',current['observation'],old['observations'][-1])]:
        np.testing.assert_array_equal(left,right)
    if current['rng_state']!=parent_summary['final_rng_state']:raise ValueError('parent RNG not restored')
    np.testing.assert_allclose(values[0]['energy'],parent_record['energy'],rtol=0,atol=1e-8)
    np.testing.assert_allclose(values[0]['forces'],parent_record['forces'],rtol=0,atol=1e-8)
    rng=np.random.default_rng();rng.bit_generator.state=current['rng_state']
    cursor=1;retained=[];accepted=[];outside_count=0;deltas=[]
    cosine,sine=np.cos(protocol['step']),np.sin(protocol['step'])
    for i,path in enumerate(checkpoints[1:]):
        next_state,tx=checkpoint(path)
        if path.name!=f'checkpoint_{i+1:04d}.npz' or next_state['completed']!=i+1:
            raise ValueError('checkpoint ordering')
        p=rng.standard_normal(len(current['u']));length=int(rng.integers(protocol['steps'][0],protocol['steps'][1]+1))
        if (tx['steps']!=length or tx['iteration']!=i or tx['global_proposal']!=56+i
                or tx['start_evaluation']!=current['retained_evaluation'] or tx['first_evaluation']!=cursor):
            raise ValueError('transaction or RNG changed')
        u=current['u'].copy();gradient=current['gradient'].copy()
        h0=.5*(u@u+p@p)+current['phi'];outside=False
        for j in range(length):
            p-=.5*protocol['step']*gradient
            u,p=cosine*u+sine*p,-sine*u+cosine*p
            value=values[cursor]
            if str(value['context'])!='sampling' or int(value['iteration'])!=i:raise ValueError('wrong path attribution')
            root_error['path']=max(root_error['path'],float(abs(u-value['u']).max()))
            cursor+=1
            if not bool(value['inside']):outside=True;break
            gradient=value['gradient'];p-=.5*protocol['step']*gradient
        if tx['last_evaluation']!=cursor-1:raise ValueError('path endpoint')
        if outside:
            accept=False;difference=None;outside_count+=1
            if tx['uniform'] is not None or tx['delta_h'] is not None:raise ValueError('outside consumed RNG')
        else:
            difference=float(.5*(u@u+p@p)+float(value['correction'])-h0)
            uniform=float(rng.random());accept=bool(np.log(uniform)<min(0.,-difference))
            if uniform!=tx['uniform']:raise ValueError('uniform RNG changed')
            root_error['delta_h']=max(root_error['delta_h'],abs(difference-tx['delta_h']))
        if accept!=tx['accepted'] or outside!=tx['domain_rejected']:raise ValueError('Metropolis decision changed')
        expected_index=cursor-1 if accept else current['retained_evaluation']
        if expected_index!=tx['retained_evaluation'] or expected_index!=next_state['retained_evaluation']:
            raise ValueError('retained rejection repeat lost')
        expected=values[expected_index] if accept else current
        for key,recordkey in [('u','u'),('phi','correction' if accept else 'phi'),('gradient','gradient'),('observation','observation')]:
            np.testing.assert_allclose(next_state[key],expected[recordkey],rtol=0,atol=1e-11)
        if rng.bit_generator.state!=next_state['rng_state']:raise ValueError('committed RNG mismatch')
        retained.append(expected_index);accepted.append(accept);deltas.append(difference);current=next_state
    if cursor!=len(values):raise ValueError('unused evaluations in a complete run')
    if max(root_error.values())>1e-9:raise ValueError('independent replay tolerance exceeded')
    result=dict(evaluations=len(values),inside_model_calls=sum(bool(v['inside']) for v in values),
        acceptance_fraction=float(np.mean(accepted)),distinct_retained_evaluations=len(set(retained)),
        retained_evaluation_indices=retained,domain_rejected_trajectories=outside_count,
        maximum_errors=root_error,all_rng_states_replayed=True,
        median_absolute_delta_h=float(np.median([abs(d) for d in deltas if d is not None])))
    return result,[values[i] for i in retained]

def terms(records,source,target,domain):
    thermal=[];components=[];features=[];bounds=[];common_thermal=[];common_components=[];atom_components=[]
    for value in records:
        r=value['positions'];f=value['forces']
        t=canonical_temperature_terms(r,f,target.free,domain['center'],source['positions'],temperature_K=300.,
            halfwidth_A=3.,minimum_pair_A=1.,taper_width_A=.4)
        residual=t['numerator_components_eV']-target.reference.thermal_energy*t['divergence_components']
        if abs(float(residual.sum())-t['residual_eV'])>1e-11:raise ValueError('thermal component mismatch')
        thermal.append([t['numerator_eV'],t['denominator'],t['residual_eV']]);components.append(residual.sum(axis=0))
        atom_components.append(residual)
        shared=canonical_temperature_terms(r,f,target.free,domain['center'],domain['center'],temperature_K=300.,
            halfwidth_A=3.,minimum_pair_A=1.,taper_width_A=.4)
        shared_residual=shared['numerator_components_eV']-target.reference.thermal_energy*shared['divergence_components']
        if abs(float(shared_residual.sum())-shared['residual_eV'])>1e-11:raise ValueError('common-reference decomposition')
        common_thermal.append([shared['numerator_eV'],shared['denominator'],shared['residual_eV']])
        common_components.append(shared_residual.sum(axis=0))
        features.append(value['observation']);bounds.append([float(abs(r[target.free]-domain['center'][target.free]).max()),float(pdist(r).min())])
    return tuple(np.asarray(v) for v in (thermal,components,features,bounds,common_thermal,common_components,atom_components))

def describe(thermal,components,features,bounds):
    block=block_statistics(thermal,blocks=8);mean=thermal.mean(axis=0)
    return dict(correlated_draws=len(thermal),mean_features=features.mean(axis=0).tolist(),
        mean_identity_residual_eV=float(mean[2]),mean_numerator_eV=float(mean[0]),mean_divergence=float(mean[1]),
        residual_block_SE_diagnostic_eV=float(block['block_standard_error'][2]),
        residual_first_half_eV=float(block['first_half_mean'][2]),residual_second_half_eV=float(block['second_half_mean'][2]),
        Cartesian_xyz_mean_residual_eV=components.mean(axis=0).tolist(),
        maximum_box_displacement_A=float(bounds[:,0].max()),minimum_pair_distance_A=float(bounds[:,1].min()))

def main(a):
    if a.output.exists():raise ValueError('fresh audit output required')
    p=read(a.pilot/'protocol.json');completion=read(a.pilot/'summary.json')
    if not completion['complete']:raise ValueError('retain partial result; not a completed extension')
    parent=a.root/p['parent_package'];oldpilot=parent/'atomistic_pilot';oldprotocol=read(oldpilot/'protocol.json')
    if sha(parent/'package_manifest.json')!=p['parent_manifest_sha256']:raise ValueError('parent fingerprint changed')
    for relative,digest in read(parent/'package_manifest.json')['files'].items():
        if sha(a.root/relative)!=digest:raise ValueError('parent modified')
    for relative,digest in p['inherited_source_sha256'].items():
        if sha(a.root/relative)!=digest:raise ValueError('original science source changed')
    if sha(a.root/'solver_v1/silicon_chain_continuation.py')!=p['continuation_core_sha256']:raise ValueError('continuation core changed')
    if sha(a.root/'results/silicon_wafer_feasibility/continue_force_thermal.py')!=p['runner_sha256']:raise ValueError('runner changed')
    domain=arrays(oldpilot/'common_domain.npz');results=[];allthermal=[];allcomponents=[];allfeatures=[];allbounds=[]
    allcommonthermal=[];allcommoncomponents=[];allatomcomponents=[]
    for i,name in enumerate(('loading8','return8')):
        source=arrays(a.root/f'results/silicon_initiation_v13/dense_{name}/raw_hessian.npz')
        target=FixedGripTarget(source['positions'],source['basis'],source['hessian'],source['gradient'],300.,domain['center'])
        # Revalidate the entire parent chain with the pre-existing independent replay.
        parent_row,_,_=replay(a.root,oldpilot,i,source,target,oldprotocol,domain['center'])
        row,newrecords=replay_extension(a.pilot/f'chain_{i}',oldpilot/f'chain_{i}',source,target,domain,p)
        oldrecords=[arrays(oldpilot/f'chain_{i}/evaluation_{j:04d}.npz') for j in parent_row['retained_evaluation_indices']]
        t,c,f,b,ct,cc,ac=terms(oldrecords+newrecords,source,target,domain)
        row['parent_evaluations_replayed']=parent_row['evaluation_count'];results.append(row)
        allthermal.append(t);allcomponents.append(c);allfeatures.append(f);allbounds.append(b)
        allcommonthermal.append(ct);allcommoncomponents.append(cc);allatomcomponents.append(ac)
    thermal,components,features,bounds=[np.asarray(v) for v in (allthermal,allcomponents,allfeatures,allbounds)]
    commonthermal,commoncomponents,atomcomponents=[np.asarray(v) for v in (allcommonthermal,allcommoncomponents,allatomcomponents)]
    if sum(r['inside_model_calls'] for r in results)!=completion['new_model_calls']:raise ValueError('model call count')
    windows={}
    for label,select in [('parent_32',slice(0,32)),('extension_64',slice(32,None)),('last_32',slice(-32,None)),('combined_96',slice(None))]:
        windows[label]=dict(chains=[describe(thermal[i,select],components[i,select],features[i,select],bounds[i,select]) for i in range(2)],
            classical_split_rhat=split_rhat(features[:,select]).tolist(),
            common_reference_identity=[describe(commonthermal[i,select],commoncomponents[i,select],features[i,select],bounds[i,select]) for i in range(2)])
    local=[]
    for i in range(2):
        average=atomcomponents[i,-32:].mean(axis=0)
        local.append(dict(last_32_sum_absolute_atom_axis_mean_residuals_eV=float(abs(average).sum()),
            last_32_max_absolute_atom_axis_mean_residual_eV=float(abs(average).max()),
            atom_axis_residual_is_not_damage_or_probability=True))
    summary=dict(completed_utc=datetime.now(timezone.utc).isoformat(),chains=results,windows=windows,
        feature_names=['total_energy_eV','common_endpoint_projection_A','rms_from_common_center_A'],
        new_actual_model_calls=completion['new_model_calls'],actual_calculation_seconds=completion['elapsed_seconds'],
        parent_calls_replayed=sum(r['parent_evaluations_replayed'] for r in results),new_calls_in_audit=0,
        local_component_diagnostics=local,
        identity_test_fields=['g_i=(r_i-own_source_i) * box_factor_i * smooth_pair_taper',
                              'g_i=(r_i-common_center_i) * box_factor_i * smooth_pair_taper'],
        identity_probe_reference_does_not_change_target_temperature=True,
        windows_declared_in_source_not_selected_for_agreement=True,additional_burn_discarded=0,
        all_samples_correlated=True,block_diagnostics_not_confidence_intervals=True,
        classical_split_rhat_not_rank_normalized=True,Cartesian_axes_not_verified_crack_coordinates=True,
        kBT_fitted=False,equilibrium_certified=False,free_energy_estimated=False,material_approved=False,
        first_crack_states_verified=False,crack_probability_estimated=False,physical_clock=None,
        new_MD=0,new_DFT=0,production_changed=False,audit_source_sha256=sha(Path(__file__)))
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'diagnostics.npz',thermal_terms=thermal,Cartesian_components=components,features=features,domain_bounds=bounds,
        common_reference_thermal_terms=commonthermal,common_reference_Cartesian_components=commoncomponents,
        atom_axis_residuals_eV=atomcomponents)
    dump(a.output/'summary.json',summary)
    dump(a.output/'input_manifest.json',{str(path.relative_to(a.root)).replace('\\','/'):sha(path) for path in sorted(a.pilot.rglob('*')) if path.is_file()})
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.4),layout='constrained')
    for i,label in enumerate(('Initial 8% start','Returned 8% start')):
        for axis,values in zip(axes,(features[i,:,0],features[i,:,1],thermal[i,:,2])):
            axis.plot(values,label=label);axis.axvline(31.5,color='grey',ls='--',lw=.8)
    for axis,title,unit in zip(axes,('Actual energy','Common projection','Boundary-safe thermal identity'),('eV','Angstrom','Residual (eV)')):
        axis.set(title=title,ylabel=unit,xlabel='Retained HMC index (NOT physical time)');axis.grid(alpha=.2)
    axes[0].legend(fontsize=8);axes[2].axhline(0,color='black',lw=.7)
    fig.savefig(a.output/'continuation_diagnostics.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
