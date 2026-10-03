"""Anchor SW to published diamond equilibrium and normal elastic benchmarks.

Reuse all completed joint-profile designs; do not reevaluate their training
geometries. Solve three amplitude anchors exactly at each fixed radial shape,
then evaluate the published perfect-diamond EOS. This is static calibration,
not first initiation, finite-T free energy, or physical time validation.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from .silicon_environment_research import DiamondCell
from .silicon_sw_material_fit import amplitude_parameters,force_features
from .silicon_bessel_reference import SW111BesselInterface
from .run_silicon_sw_material_fit import dump,FORCE_SCALE,ENERGY_SCALE
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames,reference_fields


def elastic_anchor(parameters,lattice,c11,c12,*,evaluation_counts=None):
    def evaluate(p):
        if evaluation_counts is not None:evaluation_counts['hessian']+=1
        return DiamondCell(p,lattice).evaluate(np.zeros(9),'direct')
    jets=[evaluate(amplitude_parameters(parameters,a))
          for a in [[1,1,1],[2,1,1],[1,2,1],[1,1,2]]]
    base=jets[0]
    gradients=np.column_stack([j.gradient-base.gradient for j in jets[1:]])
    hessians=np.stack([j.hessian-base.hessian for j in jets[1:]],axis=2)
    volume=lattice**3/4;scale=160.2176634/volume
    matrix=np.array([gradients[0],hessians[0,0],hessians[0,1]])
    target=np.array([0.,c11/scale,c12/scale]);sv=np.linalg.svd(matrix,compute_uv=False)
    if sv[-1]<sv[0]*1e-11:raise ValueError('anchor amplitude rank unresolved')
    alpha=np.linalg.solve(matrix,target)
    if np.any(alpha<=0):raise ValueError('elastic anchors require inadmissible nonpositive SW amplitude')
    p=amplitude_parameters(parameters,alpha)
    exact=evaluate(p)
    if max(abs(exact.gradient))>2e-8:raise ValueError('anchored diamond stationarity unresolved')
    internal=exact.hessian[6:,6:]
    if np.min(np.linalg.eigvalsh(internal))<=0:raise ValueError('anchored internal displacement unstable')
    relaxed=exact.hessian[:6,:6]-exact.hessian[:6,6:]@np.linalg.solve(internal,exact.hessian[6:,:6])
    if np.min(np.linalg.eigvalsh(relaxed))<=0:raise ValueError('anchored cubic branch elastically unstable')
    values=relaxed[[0,0,3],[0,1,3]]*scale
    if max(abs(values[:2]-[c11,c12]))>2e-6:raise ValueError('normal elastic anchors disagree with direct Hessian')
    return p,alpha,dict(matrix=matrix.tolist(),target=target.tolist(),singular_values=sv.tolist(),
        condition_number=float(sv[0]/sv[-1]),gradient=exact.gradient.tolist(),hessian=exact.hessian.tolist(),
        relaxed_hessian=relaxed.tolist(),C11_GPa=float(values[0]),C12_GPa=float(values[1]),
        C44_GPa=float(values[2]),elastic_eigenvalues=np.linalg.eigvalsh(relaxed).tolist())


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    reference=json.loads(args.reference.read_text());lattice=reference['diamond_a0']
    c11,c12,c44=[reference['diamond_'+k] for k in ['c11','c12','c44']]
    ev=np.array(reference['diamond_E_vs_V']);eos_anchor=int(np.argmin(ev[:,1]));eos_target=ev[:,1]-ev[eos_anchor,1]
    eos_scale=.02;elastic_band=.1
    sources=['solver_v1/run_silicon_sw_elastic_anchor.py','solver_v1/silicon_sw_material_fit.py',
             'solver_v1/silicon_environment_research.py','solver_v1/silicon_bessel_reference.py',
             'solver_v1/run_silicon_sw_material_fit.py','results/silicon_wafer_feasibility/source_Si.sw']
    manifest=[]
    for path in sources:
        raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    protocol=dict(scope='three exact amplitudes per fixed shape; reuse prior 27 training designs',
        lattice_anchor_A=lattice,C11_anchor_GPa=c11,C12_anchor_GPa=c12,C44_calibration_gate_GPa=c44,
        elastic_relative_gate=elastic_band,EOS_relative_energy_scale_eV_atom=eos_scale,
        uncertainty_claim=False,EOS_reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        source='https://raw.githubusercontent.com/libAtoms/silicon-testing-framework/fc252cb7d41df7e2bc672d614f3d76a40c9f2ecb/test-results/model-CASTEP_ASE-test-bulk_diamond-properties.json',
        source_DOI='10.1103/PhysRevX.8.041048; 10.5281/zenodo.1250555',
        reference_elasticity='published finite-strain stress fit with internal optimization, not a certified exact zero-strain DFT Hessian',
        current_elasticity='exact local Hessian with static Schur internal relaxation',
        targets_role='calibration, no longer independent validation',EOS_anchor=eos_anchor,
        EOS_same_series_energy_zero=True,thermal_frame_energy_anchor_separate=True,
        train='same 511 PW91 dia/surface001/surface110',heldout_used_for_fit=False,
        selection='lowest training+EOS loss among C44 band and original bulk-energy guard; fallback elastic-only diagnostic',
        new_DFT=0,new_MD=0,new_LAMMPS=0,material_approved=False,global_shape_optimum_certified=False)
    dump(args.output/'protocol.json',protocol)
    oldprofiles=json.loads((args.profiles/'profiles.json').read_text());p,_=source_parameters();rows=[]
    eligible=[];diagnostic=[];counts=dict(hessian=0);eos_calls=0
    for old in oldprofiles:
        i=old['profile'];shape=dict(p,sigma=old['sigma_A'],a=old['cutoff_A']/old['sigma_A'],gamma=old['gamma'])
        row=dict(profile=i,cutoff_A=old['cutoff_A'],sigma_A=old['sigma_A'],gamma=old['gamma'])
        before=counts['hessian']
        try:
            candidate,alpha,anchor=elastic_anchor(shape,lattice,c11,c12,evaluation_counts=counts)
            design=np.load(args.profiles/f'profile_{i}_design.npz');x,y=design['design'],design['target']
            b,t=design['bulk_design'],design['bulk_target'];bound=float(design['maximum_bulk_sq'])
            bulk_error=float(np.sum((b@alpha-t)**2));stress=design['stress_basis'];rs=design['reference_stress']
            eos=[]
            for v,e in ev:
                a=(v*8)**(1/3)
                eos.append(DiamondCell(candidate,a).evaluate(np.zeros(9),'direct',derivatives=False)/2);eos_calls+=1
            eos=np.array(eos);delta=eos-eos[eos_anchor]-eos_target
            residual=x@alpha-y;loss=float(residual@residual+np.mean((delta/eos_scale)**2))
            elastic_pass=abs(anchor['C44_GPa']/c44-1)<=elastic_band
            row.update(status='completed',amplitudes=alpha.tolist(),parameters=candidate,anchor=anchor,
                training_objective_squared=float(residual@residual),EOS_RMSE_eV_atom=float(np.sqrt(np.mean(delta**2))),
                EOS_energy_eV_atom=eos.tolist(),EOS_error_eV_atom=delta.tolist(),objective_squared=loss,
                bulk_energy_squared=bulk_error,maximum_bulk_energy_squared=bound,
                bulk_energy_guard_passed=bulk_error<=bound+2e-9,elastic_gate_passed=bool(elastic_pass),
                bulk_stress_RMSE_eV_A3=float(np.sqrt(np.mean((np.einsum('abck,k->abc',stress,alpha)-rs)**2))))
            if elastic_pass:diagnostic.append(row)
            if elastic_pass and row['bulk_energy_guard_passed']:eligible.append(row)
        except (ValueError,RuntimeError) as error:row.update(status='failed_or_inadmissible',reason=str(error))
        row['new_Hessian_evaluations']=counts['hessian']-before
        rows.append(row);dump(args.output/'profiles.json',rows)
        print(json.dumps({k:v for k,v in row.items() if k in ['profile','status','reason','objective_squared',
             'EOS_RMSE_eV_atom','bulk_energy_guard_passed','elastic_gate_passed']}),flush=True)
    if not diagnostic:
        dump(args.output/'summary.json',dict(complete=True,profiles=len(rows),elastic_gate_profiles=0,
             material_approved=False,elapsed_seconds=time.perf_counter()-start));return
    best=min(eligible or diagnostic,key=lambda r:r['objective_squared']);dump(args.output/'selected_fit.json',best)
    frames,_=read_frames(args.archive);training_info=json.loads((args.baseline/'fit.json').read_text())
    training=training_info['training_frames'];trainset=set(training);anchor_frame=training_info['anchor_frame']
    family_counts=Counter(frames[i].info['config_type'] for i in training)
    cached=np.load(args.profiles/f"profile_{best['profile']}_design.npz")
    x,y=cached['design'],cached['target'];alpha=np.array(best['amplitudes']);metrics=[];offset=0
    efstart=sum(3*len(frames[i]) for i in training)
    for j,i in enumerate(training):
        atoms=frames[i];n=len(atoms);count=family_counts[atoms.info['config_type']];num=3*n
        diff=(x[offset:offset+num]@alpha-y[offset:offset+num])*(FORCE_SCALE*np.sqrt(3*n*count));offset+=num
        ediff=float((x[efstart+j]@alpha-y[efstart+j])*(ENERGY_SCALE*np.sqrt(count)))
        metrics.append(dict(frame=i,config_type=atoms.info['config_type'],declared_xc='PW91',training=True,
             atoms=n,force_squared_error=float(diff@diff),relative_energy_error_eV_atom=ediff))
    priorselected=json.loads((args.profiles/'selected_fit.json').read_text());reused=best['profile']==priorselected['profile']
    shape=dict(p,sigma=best['sigma_A'],a=best['cutoff_A']/best['sigma_A'],gamma=best['gamma'])
    if reused:
        allfeatures=np.load(args.profiles/'selected_features.npz');ea=allfeatures['energy_basis'][anchor_frame]/len(frames[anchor_frame])
    else:ea=force_features(frames[anchor_frame].positions,frames[anchor_frame].cell.array,frames[anchor_frame].pbc,shape)[0]/len(frames[anchor_frame])
    ra=reference_fields(frames[anchor_frame])[0]/len(frames[anchor_frame]);evaluations=0 if reused else 1
    for i,atoms in enumerate(frames):
        if i in trainset:continue
        re,rf,_=reference_fields(atoms);n=len(atoms);xc=str(atoms.info.get('xc_functional','UNSPECIFIED'))
        if reused:
            lo,hi=allfeatures['offsets'][i:i+2];e,f=allfeatures['energy_basis'][i],allfeatures['force_basis'][lo:hi]
        else:e,f=force_features(atoms.positions,atoms.cell.array,atoms.pbc,shape);evaluations+=1
        diff=np.einsum('nck,k->nc',f,alpha)-rf
        metrics.append(dict(frame=i,config_type=str(atoms.info['config_type']),declared_xc=xc,training=False,
             atoms=n,force_squared_error=float(np.sum(diff**2)),
             relative_energy_error_eV_atom=float((e/n-ea)@alpha-(re/n-ra)) if xc=='PW91' else None))
    metrics.sort(key=lambda r:r['frame']);dump(args.output/'frame_metrics.json',metrics)
    groups=defaultdict(list)
    for row in metrics:groups[(row['config_type'],row['declared_xc'],row['training'])].append(row)
    output=[]
    for (kind,xc,trained),group in sorted(groups.items()):
        error=[r['relative_energy_error_eV_atom'] for r in group if r['relative_energy_error_eV_atom'] is not None]
        output.append(dict(config_type=kind,declared_xc=xc,training=trained,frames=len(group),
            force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in group)/(3*sum(r['atoms'] for r in group)))),
            relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(error)))) if error else None))
    dump(args.output/'group_metrics.json',output)
    audit=[]
    for cut in ['shuffle','glide']:
        model=SW111BesselInterface(best['parameters'],lattice,cut_kind=cut)
        for state in [[0.,0.,0.],[.25,.15,-.07],[1.5,.3,.1],[5.,.13,-.07]]:
            bs,ds=[model.evaluate(state,method=m) for m in ['bessel','direct']]
            audit.append(dict(cut=cut,state=state,energy_error_eV=abs(float(bs.value-ds.value)),
                 gradient_error_eV_A=float(np.max(abs(bs.gradient-ds.gradient))),
                 hessian_error_eV_A2=float(np.max(abs(bs.hessian-ds.hessian))),
                 bessel=dict(energy=bs.value,gradient=bs.gradient.tolist(),hessian=bs.hessian.tolist()),
                 direct=dict(energy=ds.value,gradient=ds.gradient.tolist(),hessian=ds.hessian.tolist())))
            dump(args.output/'bessel_audit.json',audit)
    consistency=max(r['energy_error_eV'] for r in audit)<2e-8 and max(r['gradient_error_eV_A'] for r in audit)<2e-7 and max(r['hessian_error_eV_A2'] for r in audit)<3e-6
    dump(args.output/'summary.json',dict(complete=True,profiles=len(rows),elastic_gate_profiles=len(diagnostic),
        elastic_and_bulk_guard_profiles=len(eligible),selected_profile=best['profile'],selected_eligible=bool(eligible),
        new_Hessian_evaluations=counts['hessian'],new_EOS_energy_evaluations=eos_calls,new_force_feature_evaluations=evaluations,
        cached_training_design_profiles=len(oldprofiles),reused_prior_validation_features=reused,
        bessel_consistency_passed=bool(consistency),material_approved=False,first_initiation_validated=False,
        physical_clock_validated=False,new_DFT=0,new_MD=0,new_LAMMPS=0,elapsed_seconds=time.perf_counter()-start))
    print((args.output/'summary.json').read_text(),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['archive','baseline','profiles','reference','output']:parser.add_argument('--'+name,required=True,type=Path)
    main(parser.parse_args())
