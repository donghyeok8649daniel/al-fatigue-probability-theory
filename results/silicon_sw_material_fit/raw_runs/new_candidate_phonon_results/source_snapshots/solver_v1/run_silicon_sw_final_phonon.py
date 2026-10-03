"""Fresh harmonic check of the final locked angular candidate only."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .silicon_sw_phonon import SWBlochMatrix
from .run_silicon_sw_material_fit import dump


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1];manifest=[]
    for name in ['run_silicon_sw_final_phonon.py','silicon_sw_phonon.py','silicon_sw_anharmonic.py','silicon_environment_research.py']:
        path='solver_v1/'+name;raw=(root/path).read_bytes();target=args.output/'source_snapshots'/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    fit=json.loads(args.candidate.read_text());ref=json.loads(args.reference.read_text());old=np.load(args.prior)
    beta=fit.get('angle_beta',fit.get('selected',{}).get('beta'))
    if beta is None:raise ValueError('explicit selected angular beta required')
    model=SWBlochMatrix(fit['parameters'],ref['diamond_a0'],angle_beta=beta)
    paths=old['q_fractional'];expected=old['reference_THz'];values=[model.spectrum(q) for q in paths.reshape(-1,3)]
    eigen=np.array([a for a,b in values]).reshape(3,50,6);bands=np.array([b for a,b in values]).reshape(3,50,6)
    sampled=[];scale=1.
    for q in np.ndindex(9,9,9):
        matrix=model.evaluate(np.array(q)/9);sampled.append(np.linalg.eigvalsh(matrix).min());scale=max(scale,float(np.max(abs(matrix))))
    np.savez_compressed(args.output/'final_candidate.npz',q_fractional=paths,eigenvalues=eigen,frequencies_THz=bands,reference_THz=expected)
    diff=bands-expected
    dump(args.output/'summary.json',dict(complete=True,candidate_profile=fit['profile'],angle_beta=beta,
        new_selected_candidate_harmonic_assembly=True,q_path_points=150,sampled_q_cube=9,
        band_RMSE_THz=float(np.sqrt(np.mean(diff**2))),band_maximum_absolute_error_THz=float(np.max(abs(diff))),
        gamma_optical_THz=bands[0,0,3:].tolist(),reference_gamma_optical_THz=expected[0,0,3:].tolist(),
        maximum_eigenvalue_difference_from_prior_anchored_shape_eV_A2=float(np.max(abs(eigen-old['eigenvalues']))),
        sampled_minimum_eigenvalue_eV_A2=float(min(sampled)),numerical_stability_floor_eV_A2=128*np.finfo(float).eps*scale,
        full_BZ_stability_certified=False,harmonic_frequency_is_not_overdamped_clock=True,
        elapsed_seconds=time.perf_counter()-start,material_approved=False,first_initiation_validated=False,new_DFT=0,new_MD=0))
    dump(args.output/'input_bindings.json',[dict(role=role,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        for role,path in [('selected',args.candidate),('bulk_reference',args.reference),('prior_path_and_reference',args.prior)]])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['candidate','reference','prior','output']:parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
