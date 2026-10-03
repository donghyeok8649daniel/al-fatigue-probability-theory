"""Compare stored relative-basis curvature with the reported Gamma modes.

No new atomic evaluation or candidate selection. The primitive cell contains
two equal-mass atoms; relative displacement has reduced mass m/2. Normal
elastic anchors and a relaxed shear anchor do not fix this optical curvature
independently. DFT-equivalent curvature below is inferred under the stated
same-mass convention, not a new DFT Hessian calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def main(args):
    if args.output.exists():
        raise ValueError('fresh diagnostic output required')
    root = args.package
    read = lambda name: json.loads((root/name).read_text(encoding='utf-8'))
    selection = read('final_selection.json')
    chosen_name = 'raw_runs/'+selection['selected_run']+'/selected_fit.json'
    chosen = read(chosen_name)
    phonon_name = 'raw_runs/final_selected_phonon_recovery_results/summary.json'
    phonon = read(phonon_name)
    ref_name = 'primary_sources/test-results/model-CASTEP_ASE-test-bulk_diamond-properties.json'
    reference = read(ref_name)
    hessian = np.array(chosen['anchor']['hessian'])
    internal = hessian[6:,6:]
    kappa = np.linalg.eigvalsh(internal)
    if min(kappa) <= 0:
        raise RuntimeError('relative basis curvature is unstable')
    mass_amu = 28.0855
    unit = np.sqrt(1.602176634e-19/(1e-20*1.66053906660e-27))/(2*np.pi*1e12)
    frequency = np.sqrt(2*kappa/mass_amu)*unit
    saved = np.sort(phonon['gamma_optical_THz'])
    if np.max(abs(frequency-saved)) > 2e-10:
        raise RuntimeError('stored bulk curvature and Gamma optical modes differ')
    target_frequency = np.sort(phonon['reference_gamma_optical_THz'])
    inferred_target_kappa = (target_frequency/unit)**2*mass_amu/2
    scale = 160.2176634/(reference['diamond_a0']**3/4)
    correction = hessian[:6,6:] @ np.linalg.solve(internal,hessian[6:,:6])
    relaxed = hessian[:6,:6]-correction
    saved_relaxed = np.array(chosen['anchor']['relaxed_hessian'])
    if np.max(abs(relaxed-saved_relaxed)) > 2e-11:
        raise RuntimeError('stored internal relaxation does not reproduce bulk anchors')
    inputs = [chosen_name,phonon_name,ref_name,'diagnose_bulk_optical.py']
    result = dict(complete=True,primitive_cell_atoms=2,relative_coordinate_reduced_mass_amu=mass_amu/2,
        atomic_mass_amu=mass_amu,internal_Hessian_eV_A2=internal.tolist(),
        internal_eigenvalues_eV_A2=kappa.tolist(),reconstructed_gamma_optical_THz=frequency.tolist(),
        recorded_gamma_optical_THz=saved.tolist(),reference_gamma_optical_THz=target_frequency.tolist(),
        inferred_reference_relative_curvature_eV_A2=inferred_target_kappa.tolist(),
        relative_optical_curvature_excess=(kappa/inferred_target_kappa-1).tolist(),
        unrelaxed_shear_curvature_GPa=float(hessian[3,3]*scale),
        internal_relaxation_shear_reduction_GPa=float(correction[3,3]*scale),
        relaxed_C44_GPa=float(relaxed[3,3]*scale),
        optical_mode_maximum_replay_difference_THz=float(np.max(abs(frequency-saved))),
        normal_strain_internal_coupling_maximum_eV_A=float(np.max(abs(hessian[:3,6:]))),
        inputs=[dict(path=n,sha256=hashlib.sha256((root/n).read_bytes()).hexdigest()) for n in inputs],
        new_geometry_evaluations=0,new_DFT=0,new_MD=0,selected_parameters_changed=False,
        actual_DFT_Hessian_inferred_as_direct_calculation=False,material_approved=False,
        physical_clock_calibrated=False,
        interpretation='same-mass stored relative-basis harmonic diagnostic; optical frequencies are not mobility or physical fatigue Hz')
    args.output.write_bytes((json.dumps(result,indent=2,allow_nan=False)+'\n').encode())
    print(json.dumps({k:result[k] for k in ['complete','relaxed_C44_GPa',
        'unrelaxed_shear_curvature_GPa','internal_relaxation_shear_reduction_GPa',
        'optical_mode_maximum_replay_difference_THz','relative_optical_curvature_excess']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
