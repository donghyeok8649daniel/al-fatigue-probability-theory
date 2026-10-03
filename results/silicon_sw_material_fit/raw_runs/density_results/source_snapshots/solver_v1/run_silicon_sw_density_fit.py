"""Fixed-anchor density-power profiles; preserve perfect-bulk Hessian and EOS.

Only one nonnegative screening exponent is added. Zero controls reuse previous
designs and EOS. All-site density and its derivative enter the same conservative
energy, force, stress, and Fourier-Bessel interface model.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from .silicon_sw_density_screening import (cubic_density_reference,ScreenedDiamondCell,
                                          ScreenedSW111BesselInterface)
from .silicon_sw_material_fit import force_features
from .run_silicon_sw_guarded_fit import design
from .run_silicon_sw_material_fit import dump
from .run_silicon_sw_stress_fit import source_stress
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames,reference_fields


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    sources=['solver_v1/run_silicon_sw_density_fit.py','solver_v1/silicon_sw_density_screening.py',
        'solver_v1/test_silicon_sw_density_screening.py','solver_v1/silicon_sw_material_fit.py',
        'solver_v1/run_silicon_sw_guarded_fit.py','solver_v1/run_silicon_sw_stress_fit.py',
        'solver_v1/silicon_bessel_reference.py','solver_v1/silicon_environment_research.py']
    bindings=[]
    for path in sources:
        raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        bindings.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',bindings)
    anchors=[r for r in json.loads((args.anchors/'profiles.json').read_text()) if r.get('elastic_gate_passed')]
    anchor_protocol=json.loads((args.anchors/'protocol.json').read_text());lattice=anchor_protocol['lattice_anchor_A']
    stress_protocol=json.loads((args.profiles/'protocol.json').read_text());stress_scale=stress_protocol['stress_scale_eV_A3']
    powers=[0.,1.,2.,4.,6.]
    dump(args.output/'protocol.json',dict(known_baseline_failure='elastic-anchored SW worsens Si111 force RMSE 1.04 to 1.54 eV/A',
        new_family='pair SW unchanged; angle_i multiplied by (rho_i/rho_bulk_fixed)**nu',
        rho='sum exp(gamma*sigma/(r-r_cut)); environmental auxiliary, not electron density',
        rho_reference='same fixed shape, perfect cubic lattice at published equilibrium',
        new_free_parameters=1,density_exponents=powers,shape_profiles=len(anchors),profiles=len(anchors)*len(powers),
        angle_bulk_zero_condition='nearest-neighbor tetrahedral cubic branch; no extension claim for longer range',
        all_site_density_before_nonlinearity=True,conservative_density_gradient_included=True,
        exact_bulk_anchors_preserved=True,perfect_cubic_EOS_reused=True,
        zero_exponent_designs_reused=True,bulk_guard='same training dia energy MSE <= original SW',
        selection='training E/F/stress+same calibrated EOS only; entire Si111 and other families excluded',
        stress_scale_eV_A3=stress_scale,reference_uncertainty_claim=False,
        global_shape_optimum_certified=False,material_approved=False,first_initiation_validated=False,
        physical_clock_validated=False,new_DFT=0,new_MD=0,new_LAMMPS=0))
    frames,_=read_frames(args.archive);old=np.load(args.baseline/'features.npz')
    ti=json.loads((args.baseline/'fit.json').read_text());training=ti['training_frames'];anchorframe=ti['anchor_frame']
    offsets=old['offsets'];e0,f0,re,rf=[old[k] for k in ['energy_basis','force_basis','reference_energy','reference_force']]
    dia=[i for i in training if frames[i].info['config_type']=='dia'];diaset=set(dia);ns=len(dia)
    reference=np.stack([source_stress(frames[i]) for i in dia]);rows=[];best=None;diagnostic=None;evaluations=0
    for arow in anchors:
        shape0,_=source_parameters();shape=dict(shape0,sigma=arow['sigma_A'],a=arow['cutoff_A']/arow['sigma_A'],gamma=arow['gamma'])
        alpha=np.array(arow['amplitudes']);rho_ref=cubic_density_reference(shape,lattice)
        for power in powers:
            row=dict(profile=len(rows),anchor_profile=arow['profile'],density_power=power,density_reference=rho_ref,
                     parameters=arow['parameters'],amplitudes=alpha.tolist(),shape=shape)
            if power==0:
                cached=np.load(args.profiles/f"profile_{arow['profile']}_design.npz")
                x,y,b,t=[cached[k] for k in ['design','target','bulk_design','bulk_target']]
                bound=float(cached['maximum_bulk_sq']);stress=cached['stress_basis'];energy=force=None
            else:
                energy=np.zeros_like(e0);force=np.zeros_like(f0);stress=[]
                for i in training:
                    a=frames[i];lo,hi=offsets[i:i+2]
                    energy[i],force[lo:hi],derivative=force_features(a.positions,a.cell.array,a.pbc,shape,
                         strain_derivative=True,density_power=power,density_reference=rho_ref)
                    if i in diaset:stress.append(derivative/a.get_volume())
                    evaluations+=1
                stress=np.stack(stress)
                x,y,b,t=design(energy,force,training,frames,offsets,re,rf,anchorframe)
                x=np.r_[x,stress.reshape(-1,3)/(stress_scale*np.sqrt(9*ns))]
                y=np.r_[y,reference.ravel()/(stress_scale*np.sqrt(9*ns))]
                bound=arow['maximum_bulk_energy_squared']
            np.savez_compressed(args.output/f"profile_{row['profile']}_design.npz",design=x,target=y,
                 bulk_design=b,bulk_target=t,maximum_bulk_sq=bound,stress_basis=stress,reference_stress=reference)
            error=x@alpha-y;bulk=float(np.sum((b@alpha-t)**2));eos=arow['EOS_RMSE_eV_atom']
            loss=float(error@error+(eos/.02)**2)
            row.update(status='completed',training_objective_squared=float(error@error),objective_squared=loss,
                bulk_energy_squared=bulk,maximum_bulk_energy_squared=bound,bulk_guard_passed=bulk<=bound+2e-9,
                EOS_RMSE_eV_atom=eos,bulk_stress_RMSE_eV_A3=float(np.sqrt(np.mean((np.einsum('abck,k->abc',stress,alpha)-reference)**2))),
                reused_zero_design=power==0)
            candidate=dict(row,energy=energy.copy() if energy is not None else None,
                           force=force.copy() if force is not None else None)
            if diagnostic is None or loss<diagnostic['objective_squared']:diagnostic=candidate
            if row['bulk_guard_passed'] and (best is None or loss<best['objective_squared']):best=candidate
            rows.append(row);dump(args.output/'profiles.json',rows)
            print(json.dumps({k:v for k,v in row.items() if k in ['profile','anchor_profile','density_power',
                'objective_squared','bulk_guard_passed','bulk_stress_RMSE_eV_A3']}),flush=True)
    best=best or diagnostic
    chosen={k:v for k,v in best.items() if k not in ['energy','force']};dump(args.output/'selected_fit.json',chosen)
    dump(args.output/'lowest_unrestricted_diagnostic.json',{k:v for k,v in diagnostic.items() if k not in ['energy','force']})
    q=np.zeros(9);bulk=ScreenedDiamondCell(best['parameters'],lattice,density_power=best['density_power'],
                                        density_reference=best['density_reference']).evaluate(q)
    oldanchor=next(r['anchor'] for r in anchors if r['profile']==best['anchor_profile'])
    curvature_error=float(np.max(abs(bulk.hessian-np.array(oldanchor['hessian']))))
    dump(args.output/'bulk_anchor_preservation.json',dict(value=bulk.value,gradient=bulk.gradient.tolist(),
         hessian=bulk.hessian.tolist(),maximum_hessian_difference_eV=curvature_error,
         full_bulk_hessian_preserved=curvature_error<2e-7))
    if best['density_power']==0:
        oldselected=json.loads((args.anchors/'selected_fit.json').read_text())
        if oldselected['profile']!=best['anchor_profile']:
            raise RuntimeError('zero-density winner has no cached full validation; preserve profiles for continuation')
        metrics=json.loads((args.anchors/'frame_metrics.json').read_text())
        audit=json.loads((args.anchors/'bessel_audit.json').read_text());reused_validation=True
    else:
        e,f=best['energy'],best['force'];power=best['density_power'];rho_ref=best['density_reference']
        trainset=set(training)
        for i,a in enumerate(frames):
            if i in trainset:continue
            lo,hi=offsets[i:i+2];e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,best['shape'],
                                              density_power=power,density_reference=rho_ref);evaluations+=1
        np.savez_compressed(args.output/'selected_features.npz',offsets=offsets,energy_basis=e,force_basis=f,
                             reference_energy=re,reference_force=rf)
        alpha=np.array(best['amplitudes']);ea=e[anchorframe]/len(frames[anchorframe]);ra=re[anchorframe]/len(frames[anchorframe])
        metrics=[]
        for i,a in enumerate(frames):
            lo,hi=offsets[i:i+2];n=len(a);xc=str(a.info.get('xc_functional','UNSPECIFIED'))
            diff=np.einsum('nck,k->nc',f[lo:hi],alpha)-rf[lo:hi]
            metrics.append(dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,training=i in trainset,
                atoms=n,force_squared_error=float(np.sum(diff**2)),
                relative_energy_error_eV_atom=float((e[i]/n-ea)@alpha-(re[i]/n-ra)) if xc=='PW91' else None))
        audit=[]
        for cut in ['shuffle','glide']:
            model=ScreenedSW111BesselInterface(best['parameters'],lattice,density_power=power,density_reference=rho_ref,cut_kind=cut)
            for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
                bs,ds=[model.evaluate(state,method=m) for m in ['bessel','direct']]
                audit.append(dict(cut=cut,state=state,energy_error_eV=abs(float(bs.value-ds.value)),
                    gradient_error_eV_A=float(np.max(abs(bs.gradient-ds.gradient))),
                    hessian_error_eV_A2=float(np.max(abs(bs.hessian-ds.hessian))),
                    bessel=dict(energy=bs.value,gradient=bs.gradient.tolist(),hessian=bs.hessian.tolist()),
                    direct=dict(energy=ds.value,gradient=ds.gradient.tolist(),hessian=ds.hessian.tolist())))
                dump(args.output/'bessel_audit.json',audit)
        reused_validation=False
    dump(args.output/'frame_metrics.json',metrics);dump(args.output/'bessel_audit.json',audit)
    groups=defaultdict(list)
    for row in metrics:groups[(row['config_type'],row['declared_xc'],row['training'])].append(row)
    output=[]
    for (kind,xc,trained),group in sorted(groups.items()):
        error=[r['relative_energy_error_eV_atom'] for r in group if r['relative_energy_error_eV_atom'] is not None]
        output.append(dict(config_type=kind,declared_xc=xc,training=trained,frames=len(group),
            force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in group)/(3*sum(r['atoms'] for r in group)))),
            relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(error)))) if error else None))
    dump(args.output/'group_metrics.json',output)
    consistency=max(r['energy_error_eV'] for r in audit)<2e-8 and max(r['gradient_error_eV_A'] for r in audit)<2e-7 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6
    summary=dict(complete=True,profiles=len(rows),selected_profile=best['profile'],density_power=best['density_power'],
        selected_bulk_guard_passed=best['bulk_guard_passed'],new_feature_evaluations=evaluations,
        zero_controls_reused=len(anchors),reused_prior_validation=reused_validation,
        full_bulk_hessian_preserved=curvature_error<2e-7,bessel_consistency_passed=bool(consistency),
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        new_DFT=0,new_MD=0,new_LAMMPS=0,elapsed_seconds=time.perf_counter()-start)
    dump(args.output/'summary.json',summary);print(json.dumps(summary),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['archive','baseline','anchors','profiles','output']:parser.add_argument('--'+name,required=True,type=Path)
    main(parser.parse_args())
