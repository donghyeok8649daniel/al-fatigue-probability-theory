"""Training-only SW shape profiles with the original bulk-energy error guard.

Reuses completed original-SW features. Searches two radial decay shapes at
fixed physical cutoff, then three positive amplitudes. Entire Si(111) families
stay outside each loss and every profile-selection decision. No production use.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from collections import defaultdict
from pathlib import Path
import numpy as np
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames,reference_fields
from .silicon_sw_material_fit import force_features,amplitude_parameters,fit_with_bulk_energy_guard
from .run_silicon_sw_material_fit import (bulk_diagnostic,dump,FORCE_SCALE,ENERGY_SCALE)
from .silicon_bessel_reference import SW111BesselInterface


def design(energy,force,training,frames,offsets,reference_energy,reference_force,anchor):
    counts={kind:sum(frames[i].info['config_type']==kind for i in training)
            for kind in ['dia','surface_001','surface_110']}
    ea=energy[anchor]/len(frames[anchor]);ra=reference_energy[anchor]/len(frames[anchor])
    xf=[];yf=[];xe=[];ye=[];bulk=[]
    for i in training:
        kind=frames[i].info['config_type'];count=counts[kind];n=len(frames[i]);lo,hi=offsets[i:i+2]
        xf.append(force[lo:hi].reshape(-1,3)/(FORCE_SCALE*np.sqrt(3*n*count)))
        yf.append(reference_force[lo:hi].ravel()/(FORCE_SCALE*np.sqrt(3*n*count)))
        xe.append((energy[i]/n-ea)/(ENERGY_SCALE*np.sqrt(count)))
        ye.append((reference_energy[i]/n-ra)/(ENERGY_SCALE*np.sqrt(count)))
        if kind=='dia':bulk.append(len(xe)-1)
    xe,ye=np.array(xe),np.array(ye)
    return np.r_[np.vstack(xf),xe],np.r_[np.concatenate(yf),ye],xe[bulk],ye[bulk]


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);started=time.perf_counter();p,_=source_parameters();rc=p['a']*p['sigma']
    sigma_values=[1.6,2.,p['sigma'],2.4,2.8,3.2];gamma_values=[.6,.9,1.2,1.5,1.8]
    protocol=dict(profile_selection='training objective only, no held-out metric selection',
        sigma_shape_A=sigma_values,gamma_shape=gamma_values,cutoff_A=rc,
        a_equals_cutoff_over_sigma=True,bulk_guard='PW91 dia energy MSE <= original SW on same 489 frames',
        guard_origin='training-only failure of first unconstrained candidate',
        free_amplitudes=3,fixed_shape_profiles=30,global_shape_optimum_certified=False,
        reference_loss_scales=dict(force_eV_A=FORCE_SCALE,energy_eV_atom=ENERGY_SCALE),
        scales_are_uncertainty=False,heldout_used_for_fit=False,
        heldout_previously_audited_for_baseline=True,prior_candidate_diagnostics_are_preserved=True,
        baseline_features_sha256=hashlib.sha256((args.baseline/'features.npz').read_bytes()).hexdigest(),
        new_DFT=0,new_MD=0,new_MACE=0,new_LAMMPS=0,production_enabled=False)
    dump(args.output/'protocol.json',protocol)
    frames,_=read_frames(args.archive,download=False)
    training=json.loads((args.baseline/'fit.json').read_text())['training_frames'];training_set=set(training)
    anchor=json.loads((args.baseline/'fit.json').read_text())['anchor_frame']
    with np.load(args.baseline/'features.npz') as old:
        offsets=old['offsets'];e0=old['energy_basis'];f0=old['force_basis']
        ref_e=old['reference_energy'];ref_f=old['reference_force']
    refs=[reference_fields(a) for a in frames]
    assert np.array_equal(ref_e,np.array([r[0] for r in refs]))
    assert np.array_equal(ref_f,np.vstack([r[1] for r in refs]))
    _,_,b0,t0=design(e0,f0,training,frames,offsets,ref_e,ref_f,anchor)
    bound=float(np.sum((b0@np.ones(3)-t0)**2))
    profiles=[];best=None;evaluations=0
    for sigma in sigma_values:
        for gamma in gamma_values:
            shape=dict(p,sigma=sigma,a=rc/sigma,gamma=gamma)
            cached=sigma==p['sigma'] and gamma==p['gamma']
            if cached:
                energy,force=e0,f0
            else:
                energy=np.zeros_like(e0);force=np.zeros_like(f0)
                for i in training:
                    a=frames[i];lo,hi=offsets[i:i+2]
                    energy[i],force[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,shape)
                    evaluations+=1
            x,y,b,t=design(energy,force,training,frames,offsets,ref_e,ref_f,anchor)
            row=dict(profile=len(profiles),sigma_A=sigma,gamma=gamma,a=rc/sigma,
                     reused_original_features=cached)
            try:
                coefficients,record=fit_with_bulk_energy_guard(x,y,b,t,bound)
                parameters=amplitude_parameters(shape,coefficients)
                row.update(status='converged_feasible',fit=record,parameters=parameters)
                if best is None or record['objective_squared']<best['fit']['objective_squared']:
                    best=dict(row,energy=energy.copy(),force=force.copy(),shape=shape)
            except (ValueError,RuntimeError) as error:
                row.update(status='failed_or_infeasible',reason=str(error))
            profiles.append(row);dump(args.output/'profiles.json',profiles)
            print(json.dumps({k:v for k,v in row.items() if k not in ('parameters','fit')}
                | {'objective_squared':row.get('fit',{}).get('objective_squared')}),flush=True)
    if best is None:raise RuntimeError('no feasible guarded shape profile')
    chosen={k:v for k,v in best.items() if k not in ('energy','force','shape')}
    dump(args.output/'selected_fit.json',chosen)
    e,f=best['energy'],best['force'];coefficients=np.array(best['fit']['amplitudes']);candidate=best['parameters']
    for i,a in enumerate(frames):
        if i not in training_set:
            lo,hi=offsets[i:i+2];e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,best['shape'])
            evaluations+=1
    np.savez_compressed(args.output/'selected_features.npz',offsets=offsets,energy_basis=e,force_basis=f,
                        reference_energy=ref_e,reference_force=ref_f)
    rows=[];ea=e[anchor]/len(frames[anchor]);e0a=e0[anchor]/len(frames[anchor]);ra=ref_e[anchor]/len(frames[anchor])
    for i,a in enumerate(frames):
        lo,hi=offsets[i:i+2];xc=str(a.info.get('xc_functional','UNSPECIFIED'));n=len(a)
        baseline=f0[lo:hi].sum(axis=2);candidate_force=np.einsum('nck,k->nc',f[lo:hi],coefficients)
        row=dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,training=i in training_set,
            atoms=n,baseline_force_RMSE_eV_A=float(np.sqrt(np.mean((baseline-ref_f[lo:hi])**2))),
            candidate_force_RMSE_eV_A=float(np.sqrt(np.mean((candidate_force-ref_f[lo:hi])**2))),
            baseline_relative_energy_error_eV_atom=None,candidate_relative_energy_error_eV_atom=None)
        if xc=='PW91':
            row.update(baseline_relative_energy_error_eV_atom=float((e0[i]/n-e0a).sum()-(ref_e[i]/n-ra)),
                       candidate_relative_energy_error_eV_atom=float((e[i]/n-ea)@coefficients-(ref_e[i]/n-ra)))
        rows.append(row)
    with (args.output/'frame_metrics.csv').open('w',newline='',encoding='utf-8') as output:
        writer=csv.DictWriter(output,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    groups=defaultdict(list)
    for row in rows:groups[(row['config_type'],row['declared_xc'],row['training'])].append(row)
    metrics=[]
    for (kind,xc,trained),group in sorted(groups.items()):
        total=sum(r['atoms'] for r in group)
        record=dict(config_type=kind,declared_xc=xc,training=trained,frames=len(group),atoms=total)
        for model in ['baseline','candidate']:
            record[model+'_force_RMSE_eV_A']=float(np.sqrt(sum(r['atoms']*r[model+'_force_RMSE_eV_A']**2 for r in group)/total))
            err=[r[model+'_relative_energy_error_eV_atom'] for r in group if r[model+'_relative_energy_error_eV_atom'] is not None]
            record[model+'_relative_energy_RMSE_eV_atom']=float(np.sqrt(np.mean(np.square(err)))) if err else None
        metrics.append(record)
    dump(args.output/'group_metrics.json',metrics)
    bulk=bulk_diagnostic(candidate);dump(args.output/'bulk_diagnostic.json',bulk)
    audit=[]
    for kind in ['shuffle','glide']:
        model=SW111BesselInterface(candidate,bulk['lattice_A'],cut_kind=kind)
        for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
            a,b=[model.evaluate(state,method=method) for method in ['bessel','direct']]
            row=dict(cut=kind,state=state,
                bessel=dict(energy=float(a.value),gradient=a.gradient.tolist(),hessian=a.hessian.tolist()),
                direct=dict(energy=float(b.value),gradient=b.gradient.tolist(),hessian=b.hessian.tolist()),
                energy_error_eV=abs(float(a.value-b.value)),gradient_error_eV_A=float(np.max(abs(a.gradient-b.gradient))),
                hessian_error_eV_A2=float(np.max(abs(a.hessian-b.hessian))))
            audit.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('bessel','direct')}),flush=True)
    dump(args.output/'bessel_candidate_audit.json',audit)
    consistency=(max(r['energy_error_eV'] for r in audit)<2e-8 and max(r['gradient_error_eV_A'] for r in audit)<2e-7
                 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6)
    summary=dict(complete=True,profiles=len(profiles),converged_feasible_profiles=sum(r['status']=='converged_feasible' for r in profiles),
        training_frames=len(training),held_out_frames=len(frames)-len(training),selected_profile=best['profile'],
        fit=best['fit'],parameters=candidate,new_feature_evaluations=evaluations,
        cached_training_frames=511,selected_training_features_reused=True,
        candidate_bessel_consistency_passed=consistency,bulk_energy_guard_passed=True,
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        production_enabled=False,new_DFT=0,new_MD=0,new_MACE=0,new_LAMMPS=0,
        elapsed_seconds=time.perf_counter()-started)
    dump(args.output/'summary.json',summary)
    root=Path(__file__).resolve().parents[1]
    sources=['solver_v1/run_silicon_sw_guarded_fit.py','solver_v1/silicon_sw_material_fit.py',
             'solver_v1/test_silicon_sw_material_fit.py','solver_v1/run_silicon_sw_material_fit.py',
             'solver_v1/silicon_bessel_reference.py','solver_v1/silicon_environment_research.py',
             'results/silicon_wafer_feasibility/source_Si.sw']
    dump(args.output/'source_manifest.json',[dict(path=path,sha256=hashlib.sha256((root/path).read_bytes()).hexdigest()) for path in sources])
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',required=True,type=Path)
    parser.add_argument('--baseline',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args())
