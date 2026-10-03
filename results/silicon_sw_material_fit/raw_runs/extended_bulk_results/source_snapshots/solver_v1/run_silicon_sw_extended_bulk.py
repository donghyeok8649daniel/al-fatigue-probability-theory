"""Bulk preflight of longer-support anharmonic SW, with fresh exact anchors.

Second neighbors remove the prior four-neighbor curvature identity. Recompute
the whole centered energy and its full Hessian, including angular bulk terms.
No old anchor transfer, no assumption that changing beta preserves elasticity.
The known DFT elastic benchmarks are calibration; force families are untouched.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .silicon_sw_anharmonic import AngularDiamondCell
from .silicon_sw_material_fit import amplitude_parameters
from .run_silicon_sw_material_fit import dump
from results.silicon_wafer_feasibility.run_static_probe import source_parameters


def angular_anchor(shape,lattice,c11,c12,beta,counts):
    jets=[]
    for a in [[1,1,1],[2,1,1],[1,2,1],[1,1,2]]:
        counts['hessian']+=1;jets.append(AngularDiamondCell(amplitude_parameters(shape,a),lattice,angle_beta=beta).evaluate(np.zeros(9)))
    base=jets[0];grad=np.column_stack([j.gradient-base.gradient for j in jets[1:]])
    h=np.stack([j.hessian-base.hessian for j in jets[1:]],axis=2)
    volume=lattice**3/4;scale=160.2176634/volume;matrix=np.array([grad[0],h[0,0],h[0,1]])
    singular=np.linalg.svd(matrix,compute_uv=False)
    if singular[-1]<singular[0]*1e-11:raise ValueError('bulk anchor amplitudes unresolved')
    alpha=np.linalg.solve(matrix,np.array([0.,c11/scale,c12/scale]))
    parameters=amplitude_parameters(shape,alpha);counts['hessian']+=1
    exact=AngularDiamondCell(parameters,lattice,angle_beta=beta).evaluate(np.zeros(9))
    if max(abs(exact.gradient))>2e-8:raise ValueError('longer-support cubic stationarity unresolved')
    internal=exact.hessian[6:,6:]
    if np.linalg.eigvalsh(internal).min()<=0:raise ValueError('internal relative basis is unstable')
    relaxed=exact.hessian[:6,:6]-exact.hessian[:6,6:]@np.linalg.solve(internal,exact.hessian[6:,:6])
    if np.linalg.eigvalsh(relaxed).min()<=0:raise ValueError('longer-support cubic branch is elastically unstable')
    elastic=relaxed[[0,0,3],[0,1,3]]*scale
    if max(abs(elastic[:2]-[c11,c12]))>2e-6:raise ValueError('bulk normal anchors differ after internal relaxation')
    return parameters,alpha,dict(C11_GPa=float(elastic[0]),C12_GPa=float(elastic[1]),C44_GPa=float(elastic[2]),
        energy_eV=exact.value,gradient=exact.gradient.tolist(),hessian=exact.hessian.tolist(),relaxed_hessian=relaxed.tolist(),
        internal_eigenvalues=np.linalg.eigvalsh(internal).tolist(),elastic_eigenvalues=np.linalg.eigvalsh(relaxed).tolist(),
        anchor_matrix=matrix.tolist(),anchor_singular_values=singular.tolist())


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1];manifest=[]
    for path in ['solver_v1/run_silicon_sw_extended_bulk.py','solver_v1/silicon_sw_anharmonic.py',
                 'solver_v1/silicon_sw_material_fit.py','solver_v1/silicon_environment_research.py']:
        raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    reference=json.loads(args.reference.read_text());lattice=reference['diamond_a0'];c11,c12,c44=[reference['diamond_'+k] for k in ['c11','c12','c44']]
    rc_values=[4.,4.2,4.4,4.6];sigma_values=[1.6,2.1,2.8];gamma_values=[.9,1.2,1.5];beta_values=[0.,.1875,.375]
    dump(args.output/'protocol.json',dict(reason='nearest-neighbor model has only two harmonic stiffnesses; finite-q DFT mismatch persists despite bulk normal anchors',
        cutoffs_A=rc_values,sigmas_A=sigma_values,gammas=gamma_values,betas=beta_values,profiles=108,
        fixed_source_angular_preference=True,all_site_bulk_angular_background_kept=True,fresh_anchor_per_beta=True,
        no_prior_amplitude_transfer=True,DFT_benchmarks_role='known calibration, not independent exact-zero-strain Hessian',
        C44_relative_calibration_band=.1,force_dataset_evaluated=False,heldout_used_for_selection=False,
        material_approved=False,new_DFT=0,new_MD=0,global_shape_optimum_certified=False))
    p,_=source_parameters();rows=[];counts=dict(hessian=0)
    for rc in rc_values:
        for sigma in sigma_values:
            for gamma in gamma_values:
                shape=dict(p,sigma=sigma,a=rc/sigma,gamma=gamma)
                for beta in beta_values:
                    before=counts['hessian'];row=dict(profile=len(rows),cutoff_A=rc,sigma_A=sigma,gamma=gamma,angle_beta=beta,shape=shape)
                    try:
                        parameters,alpha,anchor=angular_anchor(shape,lattice,c11,c12,beta,counts)
                        row.update(status='stable_bulk_anchor',parameters=parameters,amplitudes=alpha.tolist(),anchor=anchor,
                                   elastic_gate_passed=bool(abs(anchor['C44_GPa']/c44-1)<=.1))
                    except (ValueError,RuntimeError) as error:row.update(status='failed_or_inadmissible',reason=str(error))
                    row['new_Hessian_evaluations']=counts['hessian']-before;rows.append(row);dump(args.output/'profiles.json',rows)
                    print(json.dumps(dict(profile=row['profile'],status=row['status'],C44=row.get('anchor',{}).get('C44_GPa'),gate=row.get('elastic_gate_passed'))),flush=True)
    accepted=[r for r in rows if r.get('elastic_gate_passed')]
    dump(args.output/'admissible_profiles.json',accepted)
    dump(args.output/'summary.json',dict(complete=True,profiles=len(rows),stable_bulk_anchors=sum(r['status']=='stable_bulk_anchor' for r in rows),
        calibration_band_passed=len(accepted),new_Hessian_evaluations=counts['hessian'],new_force_geometry_evaluations=0,
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        new_DFT=0,new_MD=0,elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
