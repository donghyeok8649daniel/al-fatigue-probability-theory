"""Seven predeclared angular profiles on nine known normal-elastic anchors.

Fit selection uses 511 PW91 training frames, their static virial, and the
published perfect cubic EOS only. Si111 and all other families remain outside
selection. Beta=0 designs are reused. Nonzero beta EOS is reevaluated because
compressed second neighbors invalidate the tetrahedral curvature argument.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from .silicon_sw_anharmonic import AngularDiamondCell,AngularSW111BesselInterface
from .silicon_sw_material_fit import force_features
from .run_silicon_sw_guarded_fit import design
from .run_silicon_sw_material_fit import dump
from .run_silicon_sw_stress_fit import source_stress
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    sources=['solver_v1/run_silicon_sw_angular_fit.py','solver_v1/silicon_sw_anharmonic.py',
        'solver_v1/test_silicon_sw_anharmonic.py','solver_v1/silicon_sw_material_fit.py',
        'solver_v1/run_silicon_sw_guarded_fit.py','solver_v1/run_silicon_sw_stress_fit.py',
        'solver_v1/silicon_bessel_reference.py','solver_v1/silicon_environment_research.py']
    bindings=[]
    for path in sources:
        raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        bindings.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',bindings)
    anchors=[r for r in json.loads((args.anchors/'profiles.json').read_text()) if r.get('elastic_gate_passed')]
    protocol=json.loads((args.anchors/'protocol.json').read_text());lattice=protocol['lattice_anchor_A']
    stress_scale=json.loads((args.profiles/'protocol.json').read_text())['stress_scale_eV_A3']
    beta_grid=[-.375,-.25,-.125,0.,.125,.25,.375]
    reference=json.loads(args.reference.read_text());ev=np.array(reference['diamond_E_vs_V'])
    ei=int(np.argmin(ev[:,1]));et=ev[:,1]-ev[ei,1]
    dump(args.output/'protocol.json',dict(
        baseline_failure='normal-elastic SW anchor gives Si111 force RMSE 1.54 vs original 1.04 eV/A; positive density powers worsened training loss',
        family='theta**2*(1-beta*theta)**2, theta=cos(angle)-costheta0',
        beta_grid=beta_grid,additional_parameters=1,shape_profiles=len(anchors),profiles=len(anchors)*len(beta_grid),
        amplitude_anchors_reused=True,bulk_curvature_condition='four tetrahedral neighbors only',
        compressed_EOS_recomputed=True,zero_designs_reused=True,zero_EOS_reused=True,
        tensor_statistics=35,statistics_are_not_fit_parameters=True,J_orders=[0,1,2,3,4],
        self_subtraction=True,all_site_moments_before_squares=True,per_site_bulk_subtraction=True,
        selection='training E/F/stress + known calibration EOS, with original training bulk E guard and previous elastic band',
        training_frames=511,heldout_selection=False,stress_scale_eV_A3=stress_scale,
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        global_shape_optimum_certified=False,new_DFT=0,new_MD=0,new_LAMMPS=0))
    frames,_=read_frames(args.archive);old=np.load(args.baseline/'features.npz')
    fit=json.loads((args.baseline/'fit.json').read_text());training=fit['training_frames'];anchorframe=fit['anchor_frame']
    offsets=old['offsets'];e0,f0,re,rf=[old[k] for k in ['energy_basis','force_basis','reference_energy','reference_force']]
    dia=[i for i in training if frames[i].info['config_type']=='dia'];diaset=set(dia);ns=len(dia)
    rs=np.stack([source_stress(frames[i]) for i in dia]);rows=[];best=None;diagnostic=None;calls=0;eos_calls=0
    for arow in anchors:
        shape0,_=source_parameters();shape=dict(shape0,sigma=arow['sigma_A'],a=arow['cutoff_A']/arow['sigma_A'],gamma=arow['gamma'])
        alpha=np.array(arow['amplitudes'])
        for beta in beta_grid:
            row=dict(profile=len(rows),anchor_profile=arow['profile'],angle_beta=beta,parameters=arow['parameters'],
                     amplitudes=alpha.tolist(),shape=shape)
            if beta==0:
                cached=np.load(args.profiles/f"profile_{arow['profile']}_design.npz")
                x,y,b,t=[cached[k] for k in ['design','target','bulk_design','bulk_target']]
                bound=float(cached['maximum_bulk_sq']);stress=cached['stress_basis'];energy=force=None
                eos=np.array(arow['EOS_energy_eV_atom']);eos_error=np.array(arow['EOS_error_eV_atom'])
            else:
                energy=np.zeros_like(e0);force=np.zeros_like(f0);stress=[]
                for i in training:
                    a=frames[i];lo,hi=offsets[i:i+2]
                    energy[i],force[lo:hi],s=force_features(a.positions,a.cell.array,a.pbc,shape,strain_derivative=True,angle_beta=beta)
                    if i in diaset:stress.append(s/a.get_volume())
                    calls+=1
                stress=np.stack(stress);x,y,b,t=design(energy,force,training,frames,offsets,re,rf,anchorframe)
                x=np.r_[x,stress.reshape(-1,3)/(stress_scale*np.sqrt(9*ns))]
                y=np.r_[y,rs.ravel()/(stress_scale*np.sqrt(9*ns))];bound=arow['maximum_bulk_energy_squared']
                eos=[]
                for volume,_ in ev:
                    eos.append(AngularDiamondCell(row['parameters'],(8*volume)**(1/3),angle_beta=beta).evaluate(np.zeros(9),derivatives=False)/2)
                    eos_calls+=1
                eos=np.array(eos);eos_error=eos-eos[ei]-et
            np.savez_compressed(args.output/f"profile_{row['profile']}_design.npz",design=x,target=y,
                  bulk_design=b,bulk_target=t,maximum_bulk_sq=bound,stress_basis=stress,reference_stress=rs)
            error=x@alpha-y;bulk=float(np.sum((b@alpha-t)**2));eos_rmse=float(np.sqrt(np.mean(eos_error**2)))
            loss=float(error@error+(eos_rmse/.02)**2)
            row.update(status='completed',training_objective_squared=float(error@error),objective_squared=loss,
                bulk_energy_squared=bulk,maximum_bulk_energy_squared=bound,bulk_guard_passed=bulk<=bound+2e-9,
                EOS_RMSE_eV_atom=eos_rmse,EOS_energy_eV_atom=eos.tolist(),EOS_error_eV_atom=eos_error.tolist(),
                bulk_stress_RMSE_eV_A3=float(np.sqrt(np.mean((np.einsum('abck,k->abc',stress,alpha)-rs)**2))),
                reused_zero_design=beta==0)
            candidate=dict(row,energy=energy.copy() if energy is not None else None,force=force.copy() if force is not None else None)
            if diagnostic is None or loss<diagnostic['objective_squared']:diagnostic=candidate
            if row['bulk_guard_passed'] and (best is None or loss<best['objective_squared']):best=candidate
            rows.append(row);dump(args.output/'profiles.json',rows)
            print(json.dumps({k:v for k,v in row.items() if k in ['profile','anchor_profile','angle_beta','objective_squared',
                          'bulk_guard_passed','bulk_stress_RMSE_eV_A3']}),flush=True)
    best=best or diagnostic;dump(args.output/'selected_fit.json',{k:v for k,v in best.items() if k not in ['energy','force']})
    dump(args.output/'lowest_unrestricted_diagnostic.json',{k:v for k,v in diagnostic.items() if k not in ['energy','force']})
    bulk=AngularDiamondCell(best['parameters'],lattice,angle_beta=best['angle_beta']).evaluate(np.zeros(9))
    oldanchor=next(r['anchor'] for r in anchors if r['profile']==best['anchor_profile'])
    herror=float(np.max(abs(bulk.hessian-np.array(oldanchor['hessian']))))
    dump(args.output/'bulk_anchor_preservation.json',dict(value=bulk.value,gradient=bulk.gradient.tolist(),hessian=bulk.hessian.tolist(),
          maximum_hessian_difference_eV=herror,full_bulk_hessian_preserved=herror<2e-7))
    oldselected=json.loads((args.anchors/'selected_fit.json').read_text())
    reuse=best['angle_beta']==0 and best['anchor_profile']==oldselected['profile']
    if reuse:
        metrics=json.loads((args.anchors/'frame_metrics.json').read_text());audit=json.loads((args.anchors/'bessel_audit.json').read_text())
    else:
        e,f=best['energy'],best['force'];beta=best['angle_beta'];trainset=set(training)
        if e is None:e=np.zeros_like(e0);f=np.zeros_like(f0);indices=range(len(frames))
        else:indices=[i for i in range(len(frames)) if i not in trainset]
        for i in indices:
            a=frames[i];lo,hi=offsets[i:i+2];e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,best['shape'],angle_beta=beta);calls+=1
        np.savez_compressed(args.output/'selected_features.npz',offsets=offsets,energy_basis=e,force_basis=f,reference_energy=re,reference_force=rf)
        alpha=np.array(best['amplitudes']);ea=e[anchorframe]/len(frames[anchorframe]);ra=re[anchorframe]/len(frames[anchorframe]);metrics=[]
        for i,a in enumerate(frames):
            lo,hi=offsets[i:i+2];n=len(a);xc=str(a.info.get('xc_functional','UNSPECIFIED'))
            diff=np.einsum('nck,k->nc',f[lo:hi],alpha)-rf[lo:hi]
            metrics.append(dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,training=i in trainset,atoms=n,
                force_squared_error=float(np.sum(diff**2)),relative_energy_error_eV_atom=float((e[i]/n-ea)@alpha-(re[i]/n-ra)) if xc=='PW91' else None))
        audit=[]
        for cut in ['shuffle','glide']:
            model=AngularSW111BesselInterface(best['parameters'],lattice,angle_beta=beta,cut_kind=cut)
            for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
                bs,ds=[model.evaluate(state,method=m) for m in ['bessel','direct']]
                audit.append(dict(cut=cut,state=state,energy_error_eV=abs(float(bs.value-ds.value)),
                    gradient_error_eV_A=float(np.max(abs(bs.gradient-ds.gradient))),hessian_error_eV_A2=float(np.max(abs(bs.hessian-ds.hessian))),
                    bessel=dict(energy=bs.value,gradient=bs.gradient.tolist(),hessian=bs.hessian.tolist()),
                    direct=dict(energy=ds.value,gradient=ds.gradient.tolist(),hessian=ds.hessian.tolist())))
                dump(args.output/'bessel_audit.json',audit)
    dump(args.output/'frame_metrics.json',metrics);dump(args.output/'bessel_audit.json',audit)
    groups=defaultdict(list)
    for row in metrics:groups[(row['config_type'],row['declared_xc'],row['training'])].append(row)
    output=[]
    for (kind,xc,trained),group in sorted(groups.items()):
        ed=[r['relative_energy_error_eV_atom'] for r in group if r['relative_energy_error_eV_atom'] is not None]
        output.append(dict(config_type=kind,declared_xc=xc,training=trained,frames=len(group),
            force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in group)/(3*sum(r['atoms'] for r in group)))),
            relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(ed)))) if ed else None))
    dump(args.output/'group_metrics.json',output)
    consistent=max(r['energy_error_eV'] for r in audit)<2e-8 and max(r['gradient_error_eV_A'] for r in audit)<2e-7 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6
    summary=dict(complete=True,profiles=len(rows),selected_profile=best['profile'],angle_beta=best['angle_beta'],
        selected_bulk_guard_passed=best['bulk_guard_passed'],new_feature_evaluations=calls,new_EOS_energy_evaluations=eos_calls,
        new_bulk_Hessian_evaluations=1,zero_controls_reused=len(anchors),reused_prior_validation=reuse,
        full_bulk_hessian_preserved=herror<2e-7,bessel_consistency_passed=bool(consistent),
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        new_DFT=0,new_MD=0,new_LAMMPS=0,elapsed_seconds=time.perf_counter()-start)
    dump(args.output/'summary.json',summary);print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['archive','baseline','anchors','profiles','reference','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
