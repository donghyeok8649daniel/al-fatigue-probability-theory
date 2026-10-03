"""Joint conservative energy/force/stress SW profiles on the fixed PW91 training set.

Static archived DFT only. Published elastic constants are calibration gates,
not unseen validation. Each completed fit and failed profile remains in output.
"""
from __future__ import annotations
import argparse,json,hashlib,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from .silicon_sw_material_fit import force_features,fit_with_bulk_energy_guard,amplitude_parameters
from .run_silicon_sw_guarded_fit import design
from .run_silicon_sw_material_fit import dump,bulk_diagnostic
from .silicon_bessel_reference import SW111BesselInterface
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames,reference_fields
from results.silicon_wafer_feasibility.run_static_probe import source_parameters


def source_stress(atoms):
    _,_,prefix=reference_fields(atoms)
    value=np.asarray(atoms.info[prefix+'_virial'],float)
    if value.shape!=(9,) or not np.all(np.isfinite(value)):
        raise ValueError('finite source nine-component configurational virial required')
    value=value.reshape(3,3)
    if np.max(abs(value-value.T))>1e-8:
        raise ValueError('source virial symmetry unresolved')
    return -value/atoms.get_volume()


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    sources=['solver_v1/run_silicon_sw_stress_fit.py','solver_v1/silicon_sw_material_fit.py',
             'solver_v1/run_silicon_sw_guarded_fit.py','solver_v1/run_silicon_sw_material_fit.py',
             'solver_v1/silicon_environment_research.py','solver_v1/silicon_bessel_reference.py',
             'results/silicon_wafer_feasibility/source_Si.sw']
    bindings=[]
    for path in sources:
        raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        bindings.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',bindings)
    p,_=source_parameters();rc0=p['a']*p['sigma']
    cutoffs=[rc0,4.2,4.6];sigmas=[1.6,p['sigma'],2.8];gammas=[.9,1.2,1.5]
    elastic_target=np.array([153.3,56.3,72.2]);elastic_band=.1
    protocol=dict(cutoffs_A=cutoffs,sigmas_A=sigmas,gammas=gammas,profiles=27,
        stress_scale_eV_A3=args.stress_scale,scales_are_uncertainty=False,
        stress_loss='dia PW91 only, equal-frame mean Frobenius norm, no slab vacuum normalization',
        source_virial='extxyz: virial=-volume*static stress, no kinetic term',
        virial_reference='https://github.com/libAtoms/extxyz#reading-aseatomsatoms-from-this-format',
        training='same fixed PW91 511 dia/surface001/surface110; all other families excluded',
        held_out_used_for_fit=False,baseline_features_sha256=hashlib.sha256((args.baseline/'features.npz').read_bytes()).hexdigest(),
        bulk_guard='training dia E MSE <= original SW',
        elastic_target_GPa=elastic_target.tolist(),elastic_relative_gate=elastic_band,
        elastic_source='Bartok et al., Phys Rev X 8 041048 (2018), Fig 1, 0 K PW91 DFT',
        elastic_reference_role='known calibration gate; not independent unseen validation',
        selection='lowest joint training loss among profiles passing calibration elastic gate; if none, diagnostic best',
        global_shape_optimum_certified=False,new_DFT=0,new_MD=0,new_LAMMPS=0,material_approved=False)
    dump(args.output/'protocol.json',protocol)
    frames,_=read_frames(args.archive);old=np.load(args.baseline/'features.npz')
    training_info=json.loads((args.baseline/'fit.json').read_text());training=training_info['training_frames']
    anchor=training_info['anchor_frame'];offsets=old['offsets'];e0=old['energy_basis'];f0=old['force_basis']
    re,rf=old['reference_energy'],old['reference_force']
    refs=[reference_fields(a) for a in frames]
    assert np.array_equal(re,np.array([r[0] for r in refs]))
    assert np.array_equal(rf,np.vstack([r[1] for r in refs]))
    _,_,b0,t0=design(e0,f0,training,frames,offsets,re,rf,anchor)
    bound=float(np.sum((b0@np.ones(3)-t0)**2))
    dia=[i for i in training if frames[i].info['config_type']=='dia'];ns=len(dia)
    reference=np.stack([source_stress(frames[i]) for i in dia]);dump(args.output/'stress_source_frames.json',dia)
    dump(args.output/'stress_source_metadata.json',dict(frames=ns,functional='PW91',
        volume_min_A3=min(frames[i].get_volume() for i in dia),
        volume_max_A3=max(frames[i].get_volume() for i in dia),
        virial_min_eV=float(min(np.min(frames[i].info['dft_virial']) for i in dia)),
        virial_max_eV=float(max(np.max(frames[i].info['dft_virial']) for i in dia))))
    profiles=[];best_eligible=None;best_diagnostic=None;evaluations=0
    for cutoff in cutoffs:
        for sigma in sigmas:
            for gamma in gammas:
                index=len(profiles);shape=dict(p,sigma=sigma,a=cutoff/sigma,gamma=gamma)
                energy=np.zeros_like(e0);force=np.zeros_like(f0);stress=[]
                for i in training:
                    a=frames[i];lo,hi=offsets[i:i+2]
                    energy[i],force[lo:hi],derivative=force_features(a.positions,a.cell.array,a.pbc,shape,
                                                                  strain_derivative=True)
                    if i in dia:stress.append(derivative/a.get_volume())
                    evaluations+=1
                stress=np.stack(stress)
                x,y,b,t=design(energy,force,training,frames,offsets,re,rf,anchor)
                x=np.r_[x,stress.reshape(-1,3)/(args.stress_scale*np.sqrt(9*ns))]
                y=np.r_[y,reference.ravel()/(args.stress_scale*np.sqrt(9*ns))]
                # Retain sufficient statistics and exact design to independently
                # recover a failed fit without recomputing atomic geometry.
                np.savez_compressed(args.output/f'profile_{index}_design.npz',design=x,target=y,
                                    bulk_design=b,bulk_target=t,maximum_bulk_sq=bound,
                                    stress_basis=stress,reference_stress=reference)
                row=dict(profile=index,cutoff_A=cutoff,sigma_A=sigma,gamma=gamma)
                try:
                    alpha,fit=fit_with_bulk_energy_guard(x,y,b,t,bound)
                    candidate=amplitude_parameters(shape,alpha)
                    bulk=bulk_diagnostic(candidate)
                    constants=np.array([bulk[k] for k in ['C11_GPa','C12_GPa','C44_GPa']])
                    deviation=constants/elastic_target-1
                    elastic_pass=bool(np.max(abs(deviation))<=elastic_band)
                    row.update(status='completed',fit=fit,parameters=candidate,bulk=bulk,
                        elastic_relative_errors=deviation.tolist(),elastic_calibration_gate_passed=elastic_pass,
                        bulk_stress_RMSE_eV_A3=float(np.sqrt(np.mean((np.einsum('abck,k->abc',stress,alpha)-reference)**2))))
                    saved=dict(row,energy=energy.copy(),force=force.copy(),shape=shape)
                    if best_diagnostic is None or fit['objective_squared']<best_diagnostic['fit']['objective_squared']:
                        best_diagnostic=saved
                    if elastic_pass and (best_eligible is None or fit['objective_squared']<best_eligible['fit']['objective_squared']):
                        best_eligible=saved
                except (ValueError,RuntimeError) as error:
                    row.update(status='failed_or_infeasible',reason=str(error))
                profiles.append(row);dump(args.output/'profiles.json',profiles)
                print(json.dumps({k:v for k,v in row.items() if k not in ['parameters','fit','bulk']}
                      |dict(objective=row.get('fit',{}).get('objective_squared'),
                            elapsed_seconds=time.perf_counter()-start)),flush=True)
    if best_diagnostic is None:raise RuntimeError('no completed joint profile; failed raw designs preserved')
    best=best_eligible or best_diagnostic
    dump(args.output/'selected_fit.json',{k:v for k,v in best.items() if k not in ['energy','force','shape']})
    e,f=best['energy'],best['force'];alpha=np.array(best['fit']['amplitudes']);trainset=set(training)
    for i,a in enumerate(frames):
        if i in trainset:continue
        lo,hi=offsets[i:i+2];e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,best['shape']);evaluations+=1
    np.savez_compressed(args.output/'selected_features.npz',offsets=offsets,energy_basis=e,force_basis=f,
                        reference_energy=re,reference_force=rf)
    rows=[];groups=defaultdict(list);ea=e[anchor]/len(frames[anchor]);ra=re[anchor]/len(frames[anchor])
    for i,a in enumerate(frames):
        lo,hi=offsets[i:i+2];xc=str(a.info.get('xc_functional','UNSPECIFIED'));n=len(a)
        df=np.einsum('nck,k->nc',f[lo:hi],alpha)-rf[lo:hi]
        de=float((e[i]/n-ea)@alpha-(re[i]/n-ra)) if xc=='PW91' else None
        row=dict(frame=i,config_type=str(a.info['config_type']),declared_xc=xc,training=i in trainset,
             atoms=n,force_squared_error=float(np.sum(df*df)),force_RMSE_eV_A=float(np.sqrt(np.mean(df*df))),
             relative_energy_error_eV_atom=de)
        rows.append(row);groups[(row['config_type'],xc,row['training'])].append(row)
    metrics=[]
    for (kind,xc,trained),group in sorted(groups.items()):
        err=[r['relative_energy_error_eV_atom'] for r in group if r['relative_energy_error_eV_atom'] is not None]
        metrics.append(dict(config_type=kind,declared_xc=xc,training=trained,frames=len(group),
            force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in group)/(3*sum(r['atoms'] for r in group)))),
            relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(err)))) if err else None))
    dump(args.output/'frame_metrics.json',rows);dump(args.output/'group_metrics.json',metrics)
    audit=[]
    for cut in ['shuffle','glide']:
        model=SW111BesselInterface(best['parameters'],best['bulk']['lattice_A'],cut_kind=cut)
        for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
            bs,ds=[model.evaluate(state,method=m) for m in ['bessel','direct']]
            row=dict(cut=cut,state=state,energy_error_eV=abs(float(bs.value-ds.value)),
                gradient_error_eV_A=float(np.max(abs(bs.gradient-ds.gradient))),
                hessian_error_eV_A2=float(np.max(abs(bs.hessian-ds.hessian))),
                direct=dict(energy=ds.value,gradient=ds.gradient.tolist(),hessian=ds.hessian.tolist()),
                bessel=dict(energy=bs.value,gradient=bs.gradient.tolist(),hessian=bs.hessian.tolist()))
            audit.append(row);dump(args.output/'bessel_audit.json',audit)
    consistency=bool(max(r['energy_error_eV'] for r in audit)<2e-8 and
                     max(r['gradient_error_eV_A'] for r in audit)<2e-7 and
                     max(r['hessian_error_eV_A2'] for r in audit)<3e-6)
    summary=dict(complete=True,profiles=len(profiles),completed_profiles=sum(r['status']=='completed' for r in profiles),
        elastic_gate_profiles=sum(r.get('elastic_calibration_gate_passed',False) for r in profiles),
        selected_profile=best['profile'],selected_is_calibration_eligible=best_eligible is not None,
        new_feature_evaluations=evaluations,training_frames=len(training),heldout_frames=len(frames)-len(training),
        bessel_consistency_passed=consistency,material_approved=False,first_initiation_validated=False,
        physical_clock_validated=False,new_DFT=0,new_MD=0,new_LAMMPS=0,
        elapsed_seconds=time.perf_counter()-start)
    dump(args.output/'summary.json',summary);print(json.dumps(summary),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['archive','baseline','output']:parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--stress-scale',type=float,default=.05)
    args=parser.parse_args()
    if not np.isfinite(args.stress_scale) or args.stress_scale<=0:parser.error('positive finite stress scale required')
    main(args)
