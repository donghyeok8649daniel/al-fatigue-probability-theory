"""Term budget and exact fixed-shape controls; no new material or PDE.

The old validation points are inspected diagnostics, not fresh blind tests.
Forcing the perfect-interface Haa exactly is an explicit CONTROL showing the
tradeoff, not a claim of zero atomistic uncertainty or a production change.
"""
import argparse
from dataclasses import replace
import hashlib
from pathlib import Path
import time

import numpy as np

from .interface_even_development_targets import even_development_observations
from .isotropic_bulk_validation import IsotropicBulkBasis
from .material_calibration_controls import append_fixed_pair
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .tail_constrained_material import convex_profile, spectral_profile
from .validate_tail_calibration import load_material
from .vector_material_calibration import MaterialObservation
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve the earlier scientific diagnosis')
    began = time.perf_counter()
    source, old, states = source_and_targets()
    observations, _ = even_development_observations(source, old, states)
    budgets, summaries, details, bindings = [], [], [], []
    nonnegative = (0, 1, 2, 4, 5, 8, 9)
    for directory in args.models:
        data, definition, model, cache_type, _ = load_material(directory)
        if not definition.get('quadrupole_saturation_extension') or definition.get('normal_development'):
            raise ValueError('the diagnosis is for the completed v16 family')
        if definition['source_sha256'] != source.reference.sha256:
            raise ValueError('source mismatch')
        bindings.append(dict(model=directory.name, calibration_sha256=hashlib.sha256(
            (directory/'calibration.json').read_bytes()).hexdigest()))
        for label, q in [('perfect', states['perfect']), ('source_saddle', states['saddle']),
                         ('source_fault', states['fault']),
                         ('opening162', (1.62*model.h, 0., 0.)),
                         ('direct039', (model.h, .39, 0.))]:
            candidate, reference = model.evaluate(q), source.evaluate(q)
            for name, jet in candidate.components.items():
                budgets.append(dict(model=directory.name, state=label, component=name,
                    energy_eV_cell=jet[0], force_eV_L0=jet[1], Haa_eV_L0sq=jet[4],
                    Hxx_eV_L0sq=jet[7], source_Haa=reference.hessian[0, 0],
                    candidate_Haa=candidate.hessian[0, 0]))
        shape = data['best']['decays']
        raw = cache_type(observations).matrix(shape)
        matrix, obs = cubic_metric_problem(raw[:, :8], observations)
        extra = raw[:, 8:].copy()
        extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
        matrix = np.column_stack([matrix, extra])
        obs = [replace(o, role='exact') if i < 5 else o for i, o in enumerate(obs)]
        matrix, obs = append_fixed_pair(matrix, obs, definition.get('fixed_pair_control'))
        computed = []
        for stretch in definition['stability_stretches']:
            basis = IsotropicBulkBasis(shape[:3], stretch=stretch,
                                       radius=definition['radius_over_L0'])
            computed.extend(basis.evaluate(q) for q in definition['wavepoints_cubic'])
        columns, tails = map(np.asarray, zip(*computed))
        for control in ('baseline_reprofile', 'perfect_Haa_exact_control',
                        'D1_zero_control', 'omit_spectral_control'):
            M, rows = matrix.copy(), list(obs)
            if control == 'perfect_Haa_exact_control':
                rows = [replace(o, role='exact') if o.name == 'perfect_Haa' else o for o in rows]
            elif control == 'D1_zero_control':
                M = np.vstack([M, np.eye(M.shape[1])[6]])
                rows.append(MaterialObservation('control_D1_zero', 0., 1., 'eV/basis', 'exact'))
            summary = dict(model=directory.name, control=control, fixed_shape=True,
                           new_material_term=False, spectrum_enforced=control != 'omit_spectral_control')
            try:
                fit = (convex_profile(M, rows, nonnegative=nonnegative) if control == 'omit_spectral_control'
                       else spectral_profile(M, rows, columns, tails, nonnegative=nonnegative))
                c = fit['coefficients']
                H = np.einsum('c,qcij->qij', c, columns)
                margin = np.linalg.eigvalsh(H)[:, 0]-tails@abs(c)
                all_residual = (matrix@c-np.array([o.target for o in obs]))/np.array([o.scale for o in obs])
                fit_indices = [i for i, o in enumerate(obs) if o.role == 'fit']
                hold_indices = [i for i, o in enumerate(obs) if o.role == 'heldout']
                summary.update(verified_coefficient_solve=True, loss_on_original_fit_rows=float(
                    all_residual[fit_indices]@all_residual[fit_indices]),
                    inspected_validation_rms=float(np.sqrt(np.mean(all_residual[hold_indices]**2))),
                    minimum_spectral_margin=float(np.min(margin)), D1=c[6],
                    strictly_positive_LJ=fit['strictly_positive_LJ'], kkt=fit['kkt_residual'])
                for i, o in enumerate(obs):
                    if o.name in ('perfect_Haa', 'perfect_Hxx', 'saddle_Haa', 'fault_Haa', 'opening_40h_energy'):
                        summary[o.name] = float(matrix[i]@c)
                details.append(dict(model=directory.name, control=control, result=fit))
            except (ValueError, ArithmeticError) as error:
                summary.update(verified_coefficient_solve=False, rejection=str(error))
            summaries.append(summary)
            save_json(args.out/'progress.json', dict(completed=False, summaries=summaries))
            print(directory.name, control, summary, flush=True)
    keys = list(dict.fromkeys(k for row in summaries for k in row))
    write_csv(args.out/'coefficient_controls.csv', [{k: r.get(k) for k in keys} for r in summaries])
    write_csv(args.out/'normal_curvature_budget.csv', budgets)
    save_json(args.out/'coefficient_profiles.json', details)
    save_json(args.out/'scope.json', dict(completed=True, bindings=bindings,
        source_sha256=source.reference.sha256, actual_shape_optimization=False,
        existing_shape_local_controls_only=True, entire_family_impossibility_proved=False,
        heldout_points_already_inspected=True, new_physical_targets=False,
        production_changed=False, physical_seconds=False, elapsed_seconds=time.perf_counter()-began))
    save_json(args.out/'progress.json', dict(completed=True, summaries=summaries))


if __name__ == '__main__':
    main()
