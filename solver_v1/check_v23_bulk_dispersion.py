"""Static finite-wavevector stiffness, not a vibrational frequency or clock.

Both declared cubic paths are perpendicular to the straight screw line. This
checks whether matching C11/C12/C44 also matches short-wave restoring forces;
it does not attribute a nonlinear core discrepancy to one cause by itself.
"""
import argparse
from pathlib import Path
import time

import numpy as np

from .coordination_screening import CoordinationScreenedBulk
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT, load_current_material
from .run_source_core_reference import load_source_material
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .vector_material_calibration import LENGTH_M


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh static finite-q check required')
    began = time.perf_counter()
    source, _, source_binding = load_source_material()
    operator = SourceEAMBlochHessian(source.source)
    paths = dict(Gamma_L=[(t, t, t) for t in (.005, .05, .15, .3, .5)],
                 Gamma_K=[(t, t, 0.) for t in (.005, .05, .15, .3, .5, .75)])
    model_paths = [('original', None),
        ('v22_wider', ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json'),
        ('v23_radial_joint', ROOT/'results/core_interface_compatibility_v23/radial_full_validation/joint_LS/research_candidate_snapshot.json')]
    save_json(args.out/'definition.json', dict(source_sha256=source_binding['parameter_sha256'],
        paths_fractional_cubic=paths, cubic_basis='corrected +ABC stacked cubic axes',
        candidate_validation_radii_L0=[12., 16.], matrix_units='eV/L0^2',
        source_matrix_conversion='H_eV/angstrom^2 * (L0/angstrom)^2',
        fitting_performed=False, atomic_mass_used=False, physical_Hz=False,
        short_wave_error_alone_proves_core_cause=False, material_accepted=False,
        fit_partition='audit only; overlap with earlier selected bulk projections may exist; not held-out certification'))
    rows, binding_rows = [], []
    for name, path in model_paths:
        model, _, binding = load_current_material() if path is None else load_current_material(path)
        binding_rows.append(dict(model=name, **binding))
        previous = {}
        for radius in (12., 16.):
            bulk = CoordinationScreenedBulk(model.screened_shape, radius=radius, law='power')
            for path_name, points in paths.items():
                for point in points:
                    q = np.asarray(point)
                    target = operator.evaluate(q*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
                    columns, tails = bulk.evaluate(point)
                    actual = np.einsum('c,cij->ij', model.coefficients, columns)
                    bound = float(tails@abs(model.coefficients))
                    k = (path_name, point)
                    rows.append(dict(model=name, path=path_name, radius_L0=radius,
                        qx=point[0], qy=point[1], qz=point[2],
                        matrix_relative_error=float(np.linalg.norm(actual-target)/np.linalg.norm(target)),
                        screw_Hxx_target=target[0, 0], screw_Hxx_prediction=actual[0, 0],
                        screw_relative_error=float(actual[0, 0]/target[0, 0]-1),
                        minimum_H_target=float(np.linalg.eigvalsh(target)[0]),
                        minimum_H_prediction=float(np.linalg.eigvalsh(actual)[0]),
                        candidate_tail_bound=bound, margin=float(np.linalg.eigvalsh(actual)[0]-bound),
                        radius_change=None if k not in previous else float(np.linalg.norm(actual-previous[k], 2)),
                        units='eV/L0^2', physical_Hz=False, held_out_status='not asserted'))
                    previous[k] = actual
        print(name, 'actual source/candidate finite-q check completed', flush=True)
    write_csv(args.out/'static_stiffness_comparison.csv', rows)
    save_json(args.out/'completion.json', dict(completed=True, actual_matrices=len(rows),
        bindings=binding_rows, elapsed_seconds=time.perf_counter()-began,
        candidate_energy_is_infinite_LJ_Bessel=True, finite_neighbor_operator_has_analytic_tail=True,
        physical_time=False, material_accepted=False, nonlinear_core_cause_proved=False))


if __name__ == '__main__':
    main()
