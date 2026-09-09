"""Re-evaluate saved v9 data with current finite-output guards; no refit."""
import csv
import hashlib
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from .reference_eam_targets import MishinRigidFCCReference
from .run_low_stress_cyclic_diagnostic import build_surface
from .run_vector_registry_audit import OUT, save_json
from .vector_interface_reference import FullRegistryInterface, MishinVectorInterfaceReference, VectorInterfaceEvaluation
from .vector_registry_audit import stationary_state


def main():
    started = time.perf_counter()
    surface, units, meta = build_surface(tolerance=2e-11)
    models = {'analytic_candidate': FullRegistryInterface(surface),
              'Mishin_source_only': MishinVectorInterfaceReference(MishinRigidFCCReference(), units.length_scale_m/1e-10)}
    errors = dict(energy_eV_cell=0., gradient_eV_L0=0., hessian_eV_L0sq=0., grid_J_m2=0.)
    counts = {}
    for filename in ('stationary_points.csv', 'stress_and_unload.csv', 'uniform_interface_folds.csv'):
        with (OUT/filename).open(encoding='utf-8', newline='') as stream:
            rows = list(csv.DictReader(stream))
        counts[filename] = len(rows)
        for row in rows:
            q = [float(row[k]) for k in ('a_L0', 'ux_L0', 'uy_L0')]
            val = models[row['model']].evaluate(q)
            gradient = np.array([float(row['force_'+i+'_eV_L0']) for i in ('a', 'x', 'y')])
            hessian = np.array([[float(row['H_'+i+j+'_eV_L0sq']) for j in ('a', 'x', 'y')]
                                for i in ('a', 'x', 'y')])
            errors['energy_eV_cell'] = max(errors['energy_eV_cell'], abs(val.energy-float(row['energy_eV_cell'])))
            errors['gradient_eV_L0'] = max(errors['gradient_eV_L0'], float(np.max(abs(val.gradient-gradient))))
            errors['hessian_eV_L0sq'] = max(errors['hessian_eV_L0sq'], float(np.max(abs(val.hessian-hessian))))
    with (OUT/'vector_energy_grid.csv').open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    counts['vector_energy_grid.csv'] = len(rows)
    for row in rows:
        val = models[row['model']].evaluate([float(row[k]) for k in ('a_L0', 'ux_L0', 'uy_L0')])
        errors['grid_J_m2'] = max(errors['grid_J_m2'], abs(float(units.energy_to_surface(val.energy))-float(row['energy_J_m2'])))
    for nonfinite in (np.nan, np.inf):
        try:
            VectorInterfaceEvaluation.from_jet(np.full(10, nonfinite))
        except FloatingPointError:
            pass
        else:
            raise AssertionError('nonfinite interface output was not rejected')
    flat = SimpleNamespace(h=1., evaluate=lambda q: VectorInterfaceEvaluation.from_jet(np.zeros(10)))
    if stationary_state(flat, [1., 0., 0.], expected_index=0)['valid']:
        raise AssertionError('a flat zero-force separated state was accepted as a stable intact minimum')
    if max(errors.values()) > 1e-9:
        raise AssertionError(f'saved/current numerical results disagree: {errors}')
    save_json(OUT/'postguard_data_revalidation.json', dict(completed=True,
        wall_seconds=time.perf_counter()-started, state_counts=counts, maximum_errors=errors,
        nonfinite_guard_checked=True, flat_stationary_rejection_checked=True,
        parameter_sha256=meta['parameter_sha256'],
        current_kernel_sha256=hashlib.sha256(Path(__file__).with_name('vector_interface_reference.py').read_bytes()).hexdigest(),
        source_sha256=models['Mishin_source_only'].reference.sha256,
        does_not_prove_absence_of_intermittent_native_or_allocator_errors=True,
        original_single_full_suite_nan_root_cause='not reproduced or identified'))
    print(counts, errors, 'wall seconds', time.perf_counter()-started, flush=True)


if __name__ == '__main__':
    main()
