"""Joint force/energy/stress/EOS audit after fresh bulk preflight.

Only previously declared stable elastic-band profiles are evaluated. Full
finite-q Hessians are sampled before training to reject nonlocal instability.
No old four-neighbor anchor or EOS is reused across changed support/beta.
"""
import argparse,hashlib,json,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from .silicon_sw_anharmonic import AngularDiamondCell,AngularSW111BesselInterface
from .silicon_sw_phonon import SWBlochMatrix
from .silicon_sw_material_fit import force_features
from .run_silicon_sw_guarded_fit import design
from .run_silicon_sw_stress_fit import source_stress
from .run_silicon_sw_material_fit import dump
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1];manifest=[]
    for name in ['run_silicon_sw_extended_fit.py','silicon_sw_material_fit.py','silicon_sw_anharmonic.py','silicon_sw_phonon.py',
                 'run_silicon_sw_guarded_fit.py','run_silicon_sw_stress_fit.py','silicon_environment_research.py','silicon_bessel_reference.py']:
        path='solver_v1/'+name;raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    inputs=json.loads((args.preflight/'admissible_profiles.json').read_text());reference=json.loads(args.reference.read_text())
    lattice=reference['diamond_a0'];ev=np.array(reference['diamond_E_vs_V']);ei=int(np.argmin(ev[:,1]));et=ev[:,1]-ev[ei,1]
    dump(args.output/'protocol.json',dict(input_profiles=len(inputs),selection='training511 E/F/static stress+known EOS, under original bulk E guard and sampled finite-q stability',
        no_heldout_selection=True,neighborhood_support='actual SW cutoff and angular preference per profile; complete all-site bulk terms',
        preflight_protocol_sha256=hashlib.sha256((args.preflight/'protocol.json').read_bytes()).hexdigest(),
        profile_generation=json.loads((args.preflight/'protocol.json').read_text()),
        new_bulk_anchors_from_preflight=True,sampled_q_cube=9,full_BZ_stability_certified=False,
        numerical_stability_floor='128 machine epsilon times max(1,max Bloch matrix entry); raw eigenvalues saved, no clipping',
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,new_DFT=0,new_MD=0,new_LAMMPS=0))
    frames,_=read_frames(args.archive);old=np.load(args.baseline/'features.npz')
    offsets=old['offsets'];e0,f0,re,rf=[old[k] for k in ['energy_basis','force_basis','reference_energy','reference_force']]
    info=json.loads((args.baseline/'fit.json').read_text());training=info['training_frames'];af=info['anchor_frame'];trainset=set(training)
    dia=[i for i in training if frames[i].info['config_type']=='dia'];diaset=set(dia);rs=np.stack([source_stress(frames[i]) for i in dia]);ns=len(dia)
    _,_,b0,t0=design(e0,f0,training,frames,offsets,re,rf,af);bound=float(np.sum((b0@np.ones(3)-t0)**2));rows=[];best=None;diagnostic=None;calls=0;eos_calls=0
    for prior in inputs:
        row=dict(profile=prior['profile'],shape=prior['shape'],parameters=prior['parameters'],amplitudes=prior['amplitudes'],angle_beta=prior['angle_beta'],anchor=prior['anchor'])
        model=SWBlochMatrix(row['parameters'],lattice,angle_beta=row['angle_beta']);sampled=[];scale=1.
        for q in np.ndindex(9,9,9):
            matrix=model.evaluate(np.array(q)/9);values=np.linalg.eigvalsh(matrix);sampled.append(values.min());scale=max(scale,float(np.max(abs(matrix))))
        floor=128*np.finfo(float).eps*scale;minimum=float(min(sampled));row.update(sampled_minimum_eigenvalue_eV_A2=minimum,numerical_stability_floor_eV_A2=floor)
        if minimum < -floor:
            row.update(status='finite_q_instability',training_evaluated=False)
        else:
            e=np.zeros_like(e0);f=np.zeros_like(f0);stress=[];beta=row['angle_beta'];alpha=np.array(row['amplitudes'])
            for i in training:
                a=frames[i];lo,hi=offsets[i:i+2];e[i],f[lo:hi],s=force_features(a.positions,a.cell.array,a.pbc,row['shape'],strain_derivative=True,angle_beta=beta)
                if i in diaset:stress.append(s/a.get_volume())
                calls+=1
            stress=np.stack(stress);x,y,b,t=design(e,f,training,frames,offsets,re,rf,af)
            x=np.r_[x,stress.reshape(-1,3)/(.05*np.sqrt(9*ns))];y=np.r_[y,rs.ravel()/(.05*np.sqrt(9*ns))]
            np.savez_compressed(args.output/f"profile_{row['profile']}_design.npz",design=x,target=y,bulk_design=b,bulk_target=t,
                maximum_bulk_sq=bound,stress_basis=stress,reference_stress=rs)
            eos=[]
            for volume,_ in ev:
                eos.append(AngularDiamondCell(row['parameters'],(8*volume)**(1/3),angle_beta=beta).evaluate(np.zeros(9),derivatives=False)/2);eos_calls+=1
            eos=np.array(eos);ee=eos-eos[ei]-et;residual=x@alpha-y;bulk=float(np.sum((b@alpha-t)**2));loss=float(residual@residual+np.mean((ee/.02)**2))
            row.update(status='training_completed',training_evaluated=True,training_objective_squared=float(residual@residual),objective_squared=loss,
                bulk_energy_squared=bulk,maximum_bulk_energy_squared=bound,bulk_guard_passed=bool(bulk<=bound+2e-9),
                EOS_energy_eV_atom=eos.tolist(),EOS_error_eV_atom=ee.tolist(),EOS_RMSE_eV_atom=float(np.sqrt(np.mean(ee**2))),
                bulk_stress_RMSE_eV_A3=float(np.sqrt(np.mean((np.einsum('abck,k->abc',stress,alpha)-rs)**2))))
            item=dict(row,energy=e.copy(),force=f.copy())
            if diagnostic is None or loss<diagnostic['objective_squared']:diagnostic=item
            if row['bulk_guard_passed'] and (best is None or loss<best['objective_squared']):best=item
        rows.append(row);dump(args.output/'profiles.json',rows)
        print(json.dumps(dict(profile=row['profile'],status=row['status'],objective=row.get('objective_squared'),guard=row.get('bulk_guard_passed'),qmin=minimum)),flush=True)
    if diagnostic is None:
        dump(args.output/'summary.json',dict(complete=True,input_profiles=len(inputs),all_sampled_q_unstable=True,material_approved=False));return
    eligible=best is not None;best=best or diagnostic;chosen={k:v for k,v in best.items() if k not in ['energy','force']}
    dump(args.output/'selected_fit.json',chosen);dump(args.output/'selection_lock.json',dict(sha256=hashlib.sha256((args.output/'selected_fit.json').read_bytes()).hexdigest(),
        heldout_used=False,selected_eligible=eligible,selection_complete_before_validation=True))
    e,f=best['energy'],best['force'];beta=best['angle_beta'];alpha=np.array(best['amplitudes'])
    for i,a in enumerate(frames):
        if i in trainset:continue
        lo,hi=offsets[i:i+2];e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,best['shape'],angle_beta=beta);calls+=1
    np.savez_compressed(args.output/'selected_features.npz',offsets=offsets,energy_basis=e,force_basis=f,reference_energy=re,reference_force=rf)
    ea=e[af]/len(frames[af]);ra=re[af]/len(frames[af]);metrics=[]
    for i,a in enumerate(frames):
        lo,hi=offsets[i:i+2];xc=str(a.info.get('xc_functional','UNSPECIFIED'));n=len(a);diff=np.einsum('nck,k->nc',f[lo:hi],alpha)-rf[lo:hi]
        metrics.append(dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,training=i in trainset,atoms=n,
            force_squared_error=float(np.sum(diff**2)),relative_energy_error_eV_atom=float((e[i]/n-ea)@alpha-(re[i]/n-ra)) if xc=='PW91' else None))
    dump(args.output/'frame_metrics.json',metrics);groups=defaultdict(list)
    for r in metrics:groups[(r['config_type'],r['declared_xc'],r['training'])].append(r)
    output=[]
    for (kind,xc,trained),g in sorted(groups.items()):
        ed=[r['relative_energy_error_eV_atom'] for r in g if r['relative_energy_error_eV_atom'] is not None]
        output.append(dict(config_type=kind,declared_xc=xc,training=trained,frames=len(g),force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in g)/(3*sum(r['atoms'] for r in g)))),
            relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(ed)))) if ed else None))
    dump(args.output/'group_metrics.json',output);audit=[]
    for cut in ['shuffle','glide']:
        model=AngularSW111BesselInterface(best['parameters'],lattice,angle_beta=beta,cut_kind=cut)
        for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
            bs,ds=[model.evaluate(state,method=m) for m in ['bessel','direct']]
            audit.append(dict(cut=cut,state=state,energy_error_eV=abs(float(bs.value-ds.value)),gradient_error_eV_A=float(np.max(abs(bs.gradient-ds.gradient))),
                hessian_error_eV_A2=float(np.max(abs(bs.hessian-ds.hessian))),bessel=dict(energy=bs.value,gradient=bs.gradient.tolist(),hessian=bs.hessian.tolist()),
                direct=dict(energy=ds.value,gradient=ds.gradient.tolist(),hessian=ds.hessian.tolist())))
            dump(args.output/'bessel_audit.json',audit)
    consistent=max(r['energy_error_eV'] for r in audit)<2e-8 and max(r['gradient_error_eV_A'] for r in audit)<2e-7 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6
    dump(args.output/'summary.json',dict(complete=True,profiles=len(rows),selected_profile=best['profile'],selected_eligible=eligible,
        new_feature_evaluations=calls,new_EOS_energy_evaluations=eos_calls,bessel_consistency_passed=bool(consistent),
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,new_DFT=0,new_MD=0,
        elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['archive','baseline','reference','preflight','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
