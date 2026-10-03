"""Training-only radial shape refinement of the one-parameter angular family.

The previous minimum sigma lies on its search boundary. Extend this radial
shape range while retaining cutoff, exact three bulk anchors, monotone angular
cost, and the original training bulk-energy guard. Beta dependence is quadratic
in observables, so three geometry evaluations determine its quartic loss/guard.
Selection precedes all held-out evaluation and never uses Si111 errors.
"""
import argparse,hashlib,json,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from .silicon_environment_research import DiamondCell
from .silicon_sw_material_fit import force_features
from .silicon_sw_anharmonic import AngularDiamondCell,AngularSW111BesselInterface
from .run_silicon_sw_elastic_anchor import elastic_anchor
from .run_silicon_sw_angular_continuum import real_roots,squared_polynomial
from .run_silicon_sw_guarded_fit import design
from .run_silicon_sw_stress_fit import source_stress
from .run_silicon_sw_material_fit import dump
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    source_names=['run_silicon_sw_angular_shape.py','run_silicon_sw_angular_continuum.py','run_silicon_sw_elastic_anchor.py',
        'run_silicon_sw_guarded_fit.py','run_silicon_sw_stress_fit.py','run_silicon_sw_material_fit.py',
        'silicon_sw_anharmonic.py','silicon_sw_material_fit.py','silicon_environment_research.py','silicon_bessel_reference.py']
    manifest=[]
    for name in source_names:
        path='solver_v1/'+name;raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    reference=json.loads(args.reference.read_text());lattice=reference['diamond_a0'];c11=reference['diamond_c11'];c12=reference['diamond_c12'];c44=reference['diamond_c44']
    ev=np.array(reference['diamond_E_vs_V']);ei=int(np.argmin(ev[:,1]));et=ev[:,1]-ev[ei,1]
    sigmas=[.8,1.2,1.4,1.6,1.8,2.0];gammas=[.9,1.1,1.2,1.3,1.5,1.8];width=.375;controls=[-width,0.,width]
    stress_scale=.05
    dump(args.output/'protocol.json',dict(reason='previous selected training profile sigma=1.6 is on smallest sigma boundary',
        sigmas_A=sigmas,gammas=gammas,fixed_cutoff_A=3.77118,profiles=len(sigmas)*len(gammas),beta_interval=[-width,width],
        beta_search='quadratic conservative observables; quartic loss/guard roots and endpoints',
        extra_angular_parameters=1,bulk_exact_anchors=['a0','C11','C12'],C44_relative_band=.1,
        force_scale_eV_A=.2,energy_scale_eV_atom=.05,stress_scale_eV_A3=stress_scale,EOS_scale_eV_atom=.02,
        bulk_guard='same PW91 dia relative energy MSE <= original SW; no tolerance relaxation',
        heldout_used_for_selection=False,selection_locked_before_validation=True,
        no_global_shape_optimum_claim=True,material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        new_DFT=0,new_MD=0,new_LAMMPS=0))
    dump(args.output/'input_bindings.json',[dict(role=k,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        for k,path in [('archive',args.archive),('baseline_features',args.baseline/'features.npz'),('reference',args.reference),
                       ('prior_profiles',args.prior/'profiles.json')]])
    frames,_=read_frames(args.archive);old=np.load(args.baseline/'features.npz')
    offsets=old['offsets'];e0,f0,re,rf=[old[k] for k in ['energy_basis','force_basis','reference_energy','reference_force']]
    fit=json.loads((args.baseline/'fit.json').read_text());training=fit['training_frames'];af=fit['anchor_frame']
    dia=[i for i in training if frames[i].info['config_type']=='dia'];diaset=set(dia);ns=len(dia);rs=np.stack([source_stress(frames[i]) for i in dia])
    _,_,b0,t0=design(e0,f0,training,frames,offsets,re,rf,af);bound=float(np.sum((b0@np.ones(3)-t0)**2))
    oldrows=json.loads((args.prior/'profiles.json').read_text());p,_=source_parameters();rc=p['sigma']*p['a'];rows=[];best=None
    counts=dict(hessian=0,feature=0,eos=0,identity_controls=0,cached_designs=0)
    for sigma in sigmas:
        for gamma in gammas:
            shape=dict(p,sigma=sigma,a=rc/sigma,gamma=gamma);row=dict(profile=len(rows),sigma_A=sigma,gamma=gamma,cutoff_A=rc)
            try:
                candidate,alpha,anchor=elastic_anchor(shape,lattice,c11,c12,evaluation_counts=counts)
                if abs(anchor['C44_GPa']/c44-1)>.1:raise ValueError('local C44 outside same calibration band')
                designs=[];eos=[]
                for beta in controls:
                    cached=next((r for r in oldrows if r['shape']['sigma']==sigma and r['shape']['gamma']==gamma and r['angle_beta']==beta),None)
                    if cached is not None:
                        d=np.load(args.prior/f"profile_{cached['profile']}_design.npz")
                        x,y,b,t=[d[k] for k in ['design','target','bulk_design','bulk_target']]
                        eose=np.array(cached['EOS_error_eV_atom']);counts['cached_designs']+=1
                        if max(abs(np.array(cached['amplitudes'])-alpha))>2e-8:raise RuntimeError('cached anchor amplitudes differ')
                    else:
                        e=np.zeros_like(e0);f=np.zeros_like(f0);stress=[]
                        for i in training:
                            a=frames[i];lo,hi=offsets[i:i+2]
                            e[i],f[lo:hi],s=force_features(a.positions,a.cell.array,a.pbc,shape,strain_derivative=True,angle_beta=beta)
                            if i in diaset:stress.append(s/a.get_volume())
                            counts['feature']+=1
                        stress=np.stack(stress);x,y,b,t=design(e,f,training,frames,offsets,re,rf,af)
                        x=np.r_[x,stress.reshape(-1,3)/(stress_scale*np.sqrt(9*ns))];y=np.r_[y,rs.ravel()/(stress_scale*np.sqrt(9*ns))]
                        energy=[]
                        for volume,_ in ev:
                            energy.append(AngularDiamondCell(candidate,(8*volume)**(1/3),angle_beta=beta).evaluate(np.zeros(9),derivatives=False)/2)
                            counts['eos']+=1
                        energy=np.array(energy);eose=energy-energy[ei]-et
                    designs.append((x.copy(),y.copy(),b.copy(),t.copy()));eos.append(eose)
                xm,x0,xp=[d[0] for d in designs];bm,bz,bp=[d[2] for d in designs]
                xt=np.stack([x0,(xp-xm)/(2*width),(xp+xm-2*x0)/(2*width**2)])
                bt=np.stack([bz,(bp-bm)/(2*width),(bp+bm-2*bz)/(2*width**2)])
                em,ezz,ep=eos;etensor=np.stack([ezz,(ep-em)/(2*width),(ep+em-2*ezz)/(2*width**2)])
                residual=np.einsum('ijk,k->ij',xt,alpha);residual[0]-=designs[1][1]
                lv=np.concatenate([residual,etensor/(.02*np.sqrt(len(ev)))],axis=1);lp=squared_polynomial(lv)
                br=np.einsum('ijk,k->ij',bt,alpha);br[0]-=designs[1][3];gp=squared_polynomial(br);gc=gp.copy();gc[0]-=bound
                derivative=np.polynomial.polynomial.polyder(lp)
                beta_candidates=[(-width,'endpoint'),(width,'endpoint')]+[(r,'stationary') for r in real_roots(derivative,-width,width)]+[(r,'bulk boundary') for r in real_roots(gc,-width,width)]
                accepted=[]
                for beta,kind in beta_candidates:
                    if not -width<=beta<=width:continue
                    v=lv[0]+beta*lv[1]+beta**2*lv[2];be=br[0]+beta*br[1]+beta**2*br[2]
                    guard=float(be@be)
                    if guard<=bound+2e-9:accepted.append(dict(beta=beta,kind=kind,objective_squared=float(v@v),bulk_energy_squared=guard,
                        loss_derivative=float(np.polynomial.polynomial.polyval(beta,derivative)),bulk_guard_residual=guard-bound))
                if not accepted:
                    row.update(status='no_feasible_beta',parameters=candidate,amplitudes=alpha.tolist(),anchor=anchor,
                               loss_polynomial=lp.tolist(),bulk_polynomial=gp.tolist())
                else:
                    selected=min(accepted,key=lambda r:r['objective_squared'])
                    row.update(status='completed_feasible',parameters=candidate,amplitudes=alpha.tolist(),anchor=anchor,
                        selected=selected,loss_polynomial=lp.tolist(),bulk_polynomial=gp.tolist(),all_feasible_candidates=accepted)
                    if best is None or selected['objective_squared']<best['selected']['objective_squared']:best=dict(row,shape=shape)
                np.savez_compressed(args.output/f"profile_{row['profile']}_polynomial_design.npz",design=xt,bulk_design=bt,
                    target=designs[1][1],bulk_target=designs[1][3],EOS_error_basis=etensor,maximum_bulk_sq=bound,amplitudes=alpha)
                # A fourth beta at fixed direct geometries independently checks
                # observable polynomial interpolation. No held-out frames.
                beta=.137;maxerror=0.
                for i in [training[0],training[len(training)//2],training[-1]]:
                    a=frames[i];actual=force_features(a.positions,a.cell.array,a.pbc,shape,strain_derivative=True,angle_beta=beta);counts['identity_controls']+=1
                    endpoints=[force_features(a.positions,a.cell.array,a.pbc,shape,strain_derivative=True,angle_beta=b) for b in controls];counts['identity_controls']+=3
                    for k in range(3):
                        m,z,pp=[c[k] for c in endpoints];pred=z+beta*(pp-m)/(2*width)+beta**2*(pp+m-2*z)/(2*width**2)
                        maxerror=max(maxerror,float(np.max(abs(pred-actual[k]))))
                if maxerror>3e-9:raise RuntimeError('independent fourth-beta observable identity failed')
                row['fourth_beta_maximum_error']=maxerror
            except (ValueError,RuntimeError) as error:row.update(status='failed_or_inadmissible',reason=str(error))
            rows.append(row);dump(args.output/'profiles.json',rows)
            print(json.dumps(dict(profile=row['profile'],sigma=sigma,gamma=gamma,status=row['status'],selected=row.get('selected'))),flush=True)
    if best is None:raise RuntimeError('no feasible continuous angle/shape candidate')
    dump(args.output/'selected_fit.json',best)
    selection_sha=hashlib.sha256((args.output/'selected_fit.json').read_bytes()).hexdigest()
    dump(args.output/'selection_lock.json',dict(sha256=selection_sha,phase='all training shape selection completed before heldout evaluation',heldout_used=False))
    previous=json.loads((args.prior/'selected_fit.json').read_text());beta=best['selected']['beta'];parameters=best['parameters']
    reused=beta==previous['angle_beta'] and max(abs(parameters[k]-previous['parameters'][k]) for k in parameters)<2e-8
    if reused:
        metrics=json.loads((args.prior/'frame_metrics.json').read_text());audit=json.loads((args.prior/'bessel_audit.json').read_text())
    else:
        e=np.zeros_like(e0);f=np.zeros_like(f0)
        for i,a in enumerate(frames):
            lo,hi=offsets[i:i+2];e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,best['shape'],angle_beta=beta);counts['feature']+=1
        np.savez_compressed(args.output/'selected_features.npz',offsets=offsets,energy_basis=e,force_basis=f,reference_energy=re,reference_force=rf)
        alpha=np.array(best['amplitudes']);ea=e[af]/len(frames[af]);ra=re[af]/len(frames[af]);metrics=[];trainset=set(training)
        for i,a in enumerate(frames):
            lo,hi=offsets[i:i+2];xc=str(a.info.get('xc_functional','UNSPECIFIED'));n=len(a);diff=np.einsum('nck,k->nc',f[lo:hi],alpha)-rf[lo:hi]
            metrics.append(dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,training=i in trainset,atoms=n,
                force_squared_error=float(np.sum(diff**2)),relative_energy_error_eV_atom=float((e[i]/n-ea)@alpha-(re[i]/n-ra)) if xc=='PW91' else None))
        audit=[]
        for cut in ['shuffle','glide']:
            model=AngularSW111BesselInterface(parameters,lattice,angle_beta=beta,cut_kind=cut)
            for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
                bs,ds=[model.evaluate(state,method=m) for m in ['bessel','direct']]
                audit.append(dict(cut=cut,state=state,energy_error_eV=abs(float(bs.value-ds.value)),gradient_error_eV_A=float(np.max(abs(bs.gradient-ds.gradient))),
                    hessian_error_eV_A2=float(np.max(abs(bs.hessian-ds.hessian))),bessel=dict(energy=bs.value,gradient=bs.gradient.tolist(),hessian=bs.hessian.tolist()),
                    direct=dict(energy=ds.value,gradient=ds.gradient.tolist(),hessian=ds.hessian.tolist())))
                dump(args.output/'bessel_audit.json',audit)
    dump(args.output/'frame_metrics.json',metrics);dump(args.output/'bessel_audit.json',audit)
    groups=defaultdict(list)
    for r in metrics:groups[(r['config_type'],r['declared_xc'],r['training'])].append(r)
    output=[]
    for (kind,xc,trained),g in sorted(groups.items()):
        ed=[r['relative_energy_error_eV_atom'] for r in g if r['relative_energy_error_eV_atom'] is not None]
        output.append(dict(config_type=kind,declared_xc=xc,training=trained,frames=len(g),force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in g)/(3*sum(r['atoms'] for r in g)))),
            relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(ed)))) if ed else None))
    dump(args.output/'group_metrics.json',output)
    consistent=max(r['energy_error_eV'] for r in audit)<2e-8 and max(r['gradient_error_eV_A'] for r in audit)<2e-7 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6
    dump(args.output/'summary.json',dict(complete=True,profiles=len(rows),selected_profile=best['profile'],selected_beta=beta,
        counts=counts,reused_prior_validation=reused,selection_sha256=selection_sha,bessel_consistency_passed=bool(consistent),
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,new_DFT=0,new_MD=0,new_LAMMPS=0,
        elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['archive','baseline','reference','prior','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
