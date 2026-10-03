"""Known primary phonon bands versus exact harmonic Si research candidates.

Published bands use a 4x4x4 finite displacement (0.03 A) calculation; current
bands are exact local Hessians. These are harmonic validation diagnostics, not
overdamped mobility, free energy, measured lifetime or production Hz.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .silicon_sw_phonon import SWBlochMatrix,SI_MASS_AMU,EV_A2_AMU_TO_THZ
from .run_silicon_sw_material_fit import dump
from results.silicon_wafer_feasibility.run_static_probe import source_parameters


def main(args):
    if args.output.exists():raise ValueError('fresh audit output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    manifest=[]
    for path in ['solver_v1/run_silicon_sw_phonon_audit.py','solver_v1/silicon_sw_phonon.py',
                 'solver_v1/silicon_sw_anharmonic.py','solver_v1/silicon_environment_research.py']:
        raw=(root/path).read_bytes();output=args.output/'source_snapshots'/path;output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    df=json.loads(args.dft.read_text());sw=json.loads(args.sw.read_text());selected=json.loads(args.candidate.read_text())
    parameters=selected['parameters'];beta=selected.get('angle_beta',selected.get('selected',{}).get('beta'))
    if beta is None:raise ValueError('explicit selected angular beta required')
    lattice=df['a0'];p,_=source_parameters();source_p=dict(p,epsilon=2.1675)
    cases=[('original_SW_own_published_lattice',p,sw['a0'],0.,sw),
           ('published_SW_parameters_same_lattice',source_p,sw['a0'],0.,sw),
           ('anchored_pair_angle_control',parameters,lattice,0.,df),
           ('selected_anharmonic',parameters,lattice,beta,df)]
    points=np.array([[0.,0.,0.],[0.,.5,.5],[1.,1.,1.],[.5,.5,.5]])
    paths=np.array([[a+(b-a)*r for r in np.linspace(0.,1.,50)] for a,b in zip(points[:-1],points[1:])]);rows=[];spectra={}
    for name,par,a,b,ref in cases:
        model=SWBlochMatrix(par,a,angle_beta=b);band=[];eigen=[]
        for segment in paths:
            evaluated=[model.spectrum(q) for q in segment];eigen.append([v[0] for v in evaluated]);band.append([v[1] for v in evaluated])
        bands=np.array(band);expected=np.array(ref['phonon_diamond_band_frequencies']);diff=bands-expected;spectra[name]=bands
        translations=np.linalg.eigvalsh(model.evaluate([0.,0.,0.]))[:3]
        sampled=[]
        for q in np.ndindex(9,9,9):sampled.append(np.linalg.eigvalsh(model.evaluate(np.array(q)/9)).min())
        row=dict(case=name,lattice_A=a,beta=b,reference='tagged published SW finite displacement' if ref is sw else 'tagged PW91 DFT finite displacement',
            band_RMSE_THz=float(np.sqrt(np.mean(diff**2))),band_maximum_absolute_error_THz=float(np.max(abs(diff))),
            gamma_optical_THz=bands[0,0,3:].tolist(),reference_gamma_optical_THz=expected[0,0,3:].tolist(),
            X_frequencies_THz=bands[0,-1].tolist(),reference_X_frequencies_THz=expected[0,-1].tolist(),
            L_frequencies_THz=bands[-1,-1].tolist(),reference_L_frequencies_THz=expected[-1,-1].tolist(),
            gamma_translation_eigenvalues_eV_A2=translations.tolist(),sampled_9cube_minimum_eigenvalue_eV_A2=float(min(sampled)),
            analytic_BZ_stability_certified=False)
        rows.append(row);dump(args.output/'metrics.json',rows)
        np.savez_compressed(args.output/(name+'.npz'),q_fractional=paths,eigenvalues=eigen,frequencies_THz=bands,reference_THz=expected)
    invariant=float(np.max(abs(spectra['anchored_pair_angle_control']-spectra['selected_anharmonic'])))
    dump(args.output/'input_bindings.json',[dict(role=role,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for role,path in [('dft',args.dft),('published_sw',args.sw),('selection',args.candidate)]])
    dump(args.output/'summary.json',dict(complete=True,cases=len(cases),harmonic_angular_invariance_maximum_THz=invariant,
        atomic_mass_amu=SI_MASS_AMU,conversion_THz_per_sqrt_eV_A2_amu=EV_A2_AMU_TO_THZ,
        primary_method='tagged phonopy 4x4x4 primitive supercell, displacement .03 A, 50 q per segment; historical phonopy version not pinned',
        current_method='exact per-site Jet Hessian and Hermitian Bloch assembly, signed sqrt for negative curvature',
        reference_frequencies_units='THz under published phonopy default; conversion checked against tagged SW control',
        harmonic_frequencies_are_not_overdamped_clock=True,no_material_approval=True,new_DFT=0,new_MD=0,
        elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['dft','sw','candidate','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
