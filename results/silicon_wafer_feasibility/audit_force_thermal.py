"""Independent cached replay: no new model calls, MD, DFT or clock inference."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann,electron_volt
from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_thermal_identity import canonical_temperature_terms
from solver_v1.silicon_thermal_research import split_rhat,block_statistics

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def arrays(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}
def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def replay(root,pilot,chain,source,target,protocol,common):
    folder=pilot/f'chain_{chain}'
    records=json.loads((folder/'transactions.json').read_text())
    metadata=json.loads((folder/'summary.json').read_text())
    if not metadata['complete'] or len(records)!=protocol['warmup']+protocol['draws']:
        raise ValueError('chain incomplete; retain partial data instead of certifying replay')
    files=sorted(folder.glob('evaluation_*.npz'))
    values=[arrays(path) for path in files]
    thermal_energy=Boltzmann/electron_volt*protocol['temperature_K']
    direction=arrays(pilot/'common_domain.npz')['direction']
    max_coordinate=max_phi=max_gradient=max_path=max_hamiltonian=0.
    for index,(path,value) in enumerate(zip(files,values)):
        if path.name!=f'evaluation_{index:04d}.npz':raise ValueError('missing evaluation')
        u=value['u']
        position=source['positions']+(source['basis']@(target.reference.mean+
            np.sqrt(thermal_energy)*np.linalg.solve(target.reference.lower.T,u))).reshape(source['positions'].shape)
        max_coordinate=max(max_coordinate,float(np.max(abs(position-value['positions']))))
        if target.inside(position)!=bool(value['inside']):raise ValueError('domain changed')
        if bool(value['inside']):
            phi=(float(value['energy'])-float(source['energy']))/thermal_energy-.5*u@u
            gradient=np.linalg.solve(target.reference.lower,-source['basis'].T@value['forces'].ravel())/np.sqrt(thermal_energy)-u
            delta=position[target.free]-common[target.free]
            observation=np.array([float(value['energy']),float((position-common).ravel()@direction),
                float(np.sqrt(np.mean(np.sum(delta*delta,axis=1))))])
            max_phi=max(max_phi,abs(phi-float(value['correction'])))
            max_gradient=max(max_gradient,float(np.max(abs(gradient-value['gradient']))))
            np.testing.assert_allclose(value['observation'],observation,rtol=0,atol=1e-12)
        elif not np.isposinf(value['correction']):raise ValueError('finite outside energy')
    np.testing.assert_allclose(values[0]['energy'],source['energy'],rtol=0,atol=1e-8)
    np.testing.assert_allclose(values[0]['forces'],source['forces'],rtol=0,atol=1e-8)
    current_index=0;u=values[0]['u'].copy();rng=np.random.default_rng(protocol['seeds'][chain])
    retained=[];accepts=[];samples=[];potentials=[];observations=[];errors=[];warm_accept=0;outside_count=0
    cosine,sine=np.cos(protocol['step']),np.sin(protocol['step'])
    for iteration,record in enumerate(records):
        p=rng.standard_normal(len(u));length=int(rng.integers(protocol['steps'][0],protocol['steps'][1]+1))
        if record['iteration']!=iteration or record['steps']!=length or record['start_evaluation']!=current_index:
            raise ValueError('transaction/rng history changed')
        h0=.5*(u@u+p@p)+float(values[current_index]['correction'])
        candidate=u.copy();gradient=values[current_index]['gradient'].copy();outside=False
        index=record['first_evaluation']-1
        for substep in range(length):
            p-=.5*protocol['step']*gradient
            candidate,p=cosine*candidate+sine*p,-sine*candidate+cosine*p
            index+=1;value=values[index]
            if str(value['context'])!='sampling' or int(value['iteration'])!=iteration:
                raise ValueError('evaluation associated with wrong trajectory')
            max_path=max(max_path,float(np.max(abs(candidate-value['u']))))
            if not bool(value['inside']):outside=True;break
            gradient=value['gradient'];p-=.5*protocol['step']*gradient
        if index!=record['last_evaluation']:raise ValueError('path endpoint mismatch')
        if outside:
            accept=False;difference=np.inf;outside_count+=1
            if record['uniform'] is not None or record['delta_h'] is not None:raise ValueError('outside RNG consumption')
        else:
            difference=float(.5*(candidate@candidate+p@p)+float(values[index]['correction'])-h0)
            uniform=float(rng.random());accept=bool(np.log(uniform)<min(0.,-difference))
            max_hamiltonian=max(max_hamiltonian,abs(difference-record['delta_h']))
            if uniform!=record['uniform']:raise ValueError('uniform RNG mismatch')
        if accept!=record['accepted'] or outside!=record['domain_rejected']:raise ValueError('acceptance mismatch')
        if accept:current_index=index;u=candidate.copy()
        if current_index!=record['retained_evaluation']:raise ValueError('rejected-state repeat lost')
        if iteration<protocol['warmup']:warm_accept+=int(accept)
        else:
            retained.append(current_index);accepts.append(accept);samples.append(u.copy());errors.append(difference)
            observations.append(values[current_index]['observation']);potentials.append(values[current_index]['correction'])
    saved=arrays(folder/'chain.npz')
    for key,value in [('samples',samples),('observations',observations),('potentials',potentials),
                      ('accepted_at_saved_step',accepts),('retained_evaluation',retained),('hamiltonian_errors',errors)]:
        np.testing.assert_allclose(saved[key],value,rtol=0,atol=1e-11)
    if rng.bit_generator.state!=metadata['final_rng_state']:raise ValueError('final RNG mismatch')
    if abs(metadata['acceptance_fraction']-np.mean(accepts))>1e-15:raise ValueError('acceptance fraction mismatch')
    if max(max_coordinate,max_phi,max_gradient,max_path,max_hamiltonian)>1e-9:raise ValueError('independent replay tolerance')
    thermals=[];features=[];own_rms=[]
    for index in retained:
        value=values[index];r=value['positions'];f=value['forces']
        terms=canonical_temperature_terms(r,f,target.free,common,source['positions'],temperature_K=300.,
            halfwidth_A=3.,minimum_pair_A=1.,taper_width_A=.4)
        thermals.append([terms['numerator_eV'],terms['denominator'],terms['residual_eV']])
        features.append(value['observation'])
        own_rms.append(float(np.sqrt(np.mean(np.sum((r[target.free]-source['positions'][target.free])**2,axis=1)))))
    thermals=np.asarray(thermals);features=np.asarray(features);blocks=block_statistics(thermals,blocks=8)
    numerator,denominator,residual=thermals.mean(axis=0)
    absolute_errors=np.abs(np.asarray([x for x in errors if np.isfinite(x)]))
    return dict(chain=chain,evaluation_count=len(values),inside_model_calls=sum(bool(v['inside']) for v in values),
        correlated_draws=len(retained),distinct_retained_evaluations=len(set(retained)),acceptance_fraction=float(np.mean(accepts)),
        warmup_acceptance_fraction=warm_accept/protocol['warmup'],domain_rejected_trajectories=outside_count,
        maximum_coordinate_error_A=max_coordinate,maximum_correction_error=max_phi,maximum_score_error=max_gradient,
        maximum_path_error=max_path,maximum_hamiltonian_error=max_hamiltonian,
        mean_total_energy_eV=float(features[:,0].mean()),mean_endpoint_projection_A=float(features[:,1].mean()),
        mean_rms_from_common_center_A=float(features[:,2].mean()),mean_rms_from_own_start_A=float(np.mean(own_rms)),
        mean_identity_numerator_eV=float(numerator),mean_identity_denominator=float(denominator),
        mean_identity_residual_eV=float(residual),
        temperature_ratio_diagnostic_K=float(numerator/denominator/(Boltzmann/electron_volt)) if denominator>0 else None,
        block_residual_standard_error_diagnostic_eV=float(blocks['block_standard_error'][2]),
        first_half_residual_eV=float(blocks['first_half_mean'][2]),second_half_residual_eV=float(blocks['second_half_mean'][2]),
        median_absolute_delta_h=float(np.median(absolute_errors)) if len(absolute_errors) else None,
        retained_evaluation_indices=retained,physical_time_ps=None,equilibrium_certified=False),thermals,features

def main(a):
    if a.output.exists():raise ValueError('fresh audit output required')
    protocol=json.loads((a.pilot/'protocol.json').read_text())
    completion=json.loads((a.pilot/'summary.json').read_text())
    if not completion['complete']:raise ValueError('pilot incomplete')
    for relative,digest in protocol['source_sha256'].items():
        if sha(a.root/relative)!=digest:raise ValueError('input source changed: '+relative)
    if sha(a.root/'results/silicon_wafer_feasibility/sample_force_thermal.py')!=protocol['runner_sha256']:
        raise ValueError('actual runner changed')
    sources=[arrays(a.root/f'results/silicon_initiation_v13/dense_{name}/raw_hessian.npz') for name in ('loading8','return8')]
    common=(sources[0]['positions']+sources[1]['positions'])/2
    results=[];thermals=[];features=[]
    for chain,source in enumerate(sources):
        target=FixedGripTarget(source['positions'],source['basis'],source['hessian'],source['gradient'],300.,common)
        row,thermal,feature=replay(a.root,a.pilot,chain,source,target,protocol,common)
        results.append(row);thermals.append(thermal);features.append(feature)
    if sum(row['inside_model_calls'] for row in results)!=completion['new_model_calls']:
        raise ValueError('actual call count mismatch')
    feature_array=np.asarray(features);rhat=split_rhat(feature_array)
    previous=json.loads((a.root/'results/silicon_crack_thermal_validation/analysis/summary.json').read_text())
    summary=dict(completed_utc=datetime.now(timezone.utc).isoformat(),chains=results,
        feature_names=['total_energy_eV','endpoint_projection_A','rms_from_common_center_A'],
        classical_split_rhat_diagnostic=rhat.tolist(),
        previous_pcn_split_rhat_diagnostic=previous['classical_split_rhat_diagnostic'],
        previous_pcn_acceptance_fractions=[r['acceptance_fraction'] for r in previous['chains']],
        previous_pcn_identity_residuals_eV=[r['mean_identity_residual_eV'] for r in previous['chains']],
        previous_pcn_calls=previous['new_potential_calls'],new_actual_model_calls=completion['new_model_calls'],
        actual_run_elapsed_seconds=completion['elapsed_seconds'],new_calls_in_replay=0,new_MD=0,new_DFT=0,
        comparison_is_not_equal_cost_efficiency_trial=True,short_correlated_chains_are_not_independent_samples=True,
        error_diagnostics_are_not_confidence_intervals=True,equilibrium_certified=False,
        free_energy_estimated=False,kBT_fitted=False,material_approved=False,
        first_crack_states_verified=False,crack_probability_estimated=False,physical_clock=None,production_changed=False,
        audit_source_sha256=sha(Path(__file__)))
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'chain_diagnostics.npz',features=feature_array,thermal_terms=np.asarray(thermals))
    dump(a.output/'summary.json',summary)
    dump(a.output/'input_manifest.json',{str(p.relative_to(a.pilot)).replace('\\','/'):sha(p)
        for p in sorted(a.pilot.rglob('*')) if p.is_file()})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.4),layout='constrained')
    for chain,label in enumerate(['Initial 8% start','Returned 8% start']):
        for axis,values in zip(axes,[feature_array[chain,:,0],feature_array[chain,:,1],np.asarray(thermals)[chain,:,2]]):
            axis.plot(values,label=label)
    for axis,title,ylabel in zip(axes,['Actual energy','Common projection','Boundary-safe thermal identity'],['eV','Angstrom','Residual (eV)']):
        axis.set(title=title,ylabel=ylabel,xlabel='Retained HMC index (NOT time)');axis.grid(alpha=.2)
    axes[2].axhline(0,color='black',lw=.7);axes[0].legend(fontsize=8)
    fig.savefig(a.output/'force_thermal_diagnostics.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--pilot',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);main(parser.parse_args())
