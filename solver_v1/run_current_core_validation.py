"""Independent representation validation, not relaxation or material fitting.

Finite real-atom sums are used ONLY as converging test references. Production
and this research core retain their infinite Poisson/Bessel atomic rows.
"""
from __future__ import annotations

import argparse
import hashlib
import time

import numpy as np

from .current_material_rows import CurrentMaterialRowKernel, CurrentMaterialScrewCore
from .isolated_screw_core import ScrewFarField, row_environment
from .run_current_material_core import load_current_material
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .vector_fcc_validation import direct_row_channels


def direct_core_energy_gradient(core, field, *, images=128):
    """Independent real-atom sum and explicit per-site/bond assembly.

    The unchanged analytic site law has its own interface/FD tests. Only its
    local derivatives and the fixed bulk reference are shared here; no Bessel
    row coefficient or reciprocal derivative is used by this reference.
    """
    full = core.full_field(field)
    reference = [direct_row_channels(core.rows, r, images=images, include_odd_quadrupole=True)
                 for r in core.reference]
    values = []
    derivatives = []
    for source, neighbors in zip(core.energy_ids, core.destination):
        actual = [direct_row_channels(core.rows, r+full[target]-full[source],
                    images=images, include_odd_quadrupole=True)
                  for r, target in zip(core.reference, neighbors)]
        values.append(sum(a['value']-b['value'] for a, b in zip(actual, reference)))
        derivatives.append([a['gradient'] for a in actual])
    energy, weights, _ = core.site_law.evaluate(np.asarray(values), second=False)
    gradient = np.zeros_like(full)
    for source, neighbors, local, w in zip(core.energy_ids, core.destination, derivatives, weights):
        for target, bond in zip(neighbors, local):
            force = w@bond
            gradient[source] -= force
            gradient[target] += force
    return dict(energy=float(energy.sum()), gradient=gradient[core.free_ids],
                all_site_gradient=gradient)


def scaled_error(actual, reference):
    error = abs(actual-reference)
    return dict(maximum_absolute_error=float(np.max(error)),
        maximum_scaled_error=float(np.max(error/np.maximum(1., abs(reference)))))


def main():
    from pathlib import Path
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh validation directory required')
    began = time.perf_counter()
    model, tensor, metadata = load_current_material()
    rows = row_environment(model.base.surface)
    points = np.array([[x, r*np.cos(angle), r*np.sin(angle)]
        for x, r, angle in [(0., .85, 0.), (.17, .72, .3), (-.33, .91, 1.7),
            (.49, 1., 2.9), (.12, 1.4, 4.2), (-.07, 2.8, 5.3), (.38, 4., .7)]])
    kernel = CurrentMaterialRowKernel(rows)
    exact = kernel.evaluate(points, order=2)
    records = []
    for images in (16, 32, 64, 128):
        for index, point in enumerate(points):
            direct = direct_row_channels(rows, point, images=images, include_odd_quadrupole=True)
            for key in ('value', 'gradient', 'hessian'):
                records.append(dict(point=index, x=point[0], y=point[1], z=point[2],
                    real_images_each_side=images, field=key,
                    **scaled_error(direct[key], exact[key][index])))
    write_csv(args.out/'direct_row_convergence.csv', records)
    records = []
    for tolerance in (2e-9, 2e-11, 2e-13):
        result = CurrentMaterialRowKernel(rows, tolerance=tolerance).evaluate(points, order=2)
        for key in ('value', 'gradient', 'hessian'):
            records.append(dict(tolerance=tolerance, field=key,
                modes=result['modes_used'], last_mode_envelope=result['maximum_last_mode_envelope'],
                **scaled_error(result[key], exact[key])))
    write_csv(args.out/'reciprocal_refinement.csv', records)
    far = ScrewFarField(tensor, rows.b, (np.sqrt(3)*rows.b/12, rows.h/2))
    core = CurrentMaterialScrewCore(model, far, free_radius=.85, ring=1)
    field = core.initial+.004*np.sin(np.arange(core.initial.size)*.53).reshape(core.initial.shape)
    value, apply = core.linearize(field)
    records = []
    for images in (32, 64, 128):
        direct = direct_core_energy_gradient(core, field, images=images)
        records.append(dict(real_images_each_side=images,
            energy_error_eV=abs(direct['energy']-value['energy']),
            force_error_eV_L0=float(np.max(abs(direct['gradient']-value['gradient']))),
            all_site_force_sum=float(np.max(abs(direct['all_site_gradient'].sum(axis=0))))))
    write_csv(args.out/'independent_core_assembly.csv', records)
    v = np.cos(np.arange(field.size)*.37).reshape(field.shape); v /= np.linalg.norm(v)
    hessian = apply(v)
    checks = []
    for step in (2e-4, 1e-4, 5e-5, 2.5e-5, 1.25e-5):
        plus, minus = core.evaluate(field+step*v), core.evaluate(field-step*v)
        checks.append(dict(step=step,
            energy_derivative_error=abs((plus['energy']-minus['energy'])/(2*step)-np.sum(value['gradient']*v)),
            hessian_vector_error=float(np.max(abs((plus['gradient']-minus['gradient'])/(2*step)-hessian)))))
    write_csv(args.out/'nonuniform_core_derivatives.csv', checks)
    files = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in args.out.glob('*.csv')}
    save_json(args.out/'validation.json', dict(completed=True, **metadata,
        row_points=len(points), all_22_channels=True, independent_direct_assembly=True,
        maximum_fine_core_energy_error=records[-1]['energy_error_eV'],
        maximum_fine_core_force_error=records[-1]['force_error_eV_L0'],
        actual_relaxations_performed=False, material_validation_passed=False,
        elapsed_seconds=time.perf_counter()-began, result_sha256=files))
    print('independent current-material row/core representation checks completed', flush=True)


if __name__ == '__main__':
    main()
