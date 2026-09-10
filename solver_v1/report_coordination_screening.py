"""Replay completed v19 candidates and their exact-constraint sensitivities.

This reporter does NOT optimize or silently refit a saved parameter vector.
Plots use the actually saved independent curves, not fabricated dynamics.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from .coordination_screening import CoordinationScreenedCache
from .profiled_identifiability import equality_tangent_sensitivity
from .vector_material_calibration import MaterialObservation
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_coordination_shape_refinement import cubic_observations, residual_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--joint-directory', type=Path, required=True)
    parser.add_argument('--validation-directory', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('new report directory required')
    definition = json.loads((args.joint_directory/'definition.json').read_bytes())
    observations = [MaterialObservation(**{**o,
        'state': None if o['state'] is None else tuple(o['state']),
        'jet_weights': None if o['jet_weights'] is None else tuple(o['jet_weights'])})
        for o in definition['observations']]
    observations = [replace(o, role='exact') if i < 5 else o for i, o in enumerate(observations)]
    fit_observations = cubic_observations(observations)
    cache = CoordinationScreenedCache(observations, law='power')

    def matrix(shape):
        raw = cache.matrix(shape)
        M, _ = cubic_metric_problem(raw[:, :8], observations)
        extra = raw[:, 8:].copy(); extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
        return np.column_stack([M, extra])

    summary, all_rows = [], []
    for family in ('old_family', 'power_family'):
        data = json.loads((args.joint_directory/family/'calibration.json').read_bytes())
        if not data['completed'] or data['best'] is None:
            raise ValueError('completed positive-pair fits required for sensitivity replay')
        best = data['best']; c = np.asarray(best['coefficients']); shape = np.asarray(best['shape'])
        M = matrix(shape)
        error = float(np.max(abs(M@c-best['predictions'])))
        if error > 2e-8:
            raise ArithmeticError('saved parameter/observation replay mismatch')
        nshape = 5 if family == 'old_family' else 6
        derivatives = []
        for step in (2e-4, 1e-4):
            dM = []
            for axis in range(nshape):
                if axis < 5:
                    hi, lo = shape.copy(), shape.copy()
                    hi[axis] *= np.exp(step); lo[axis] *= np.exp(-step)
                    delta = (matrix(hi)-matrix(lo))/(2*step)
                elif shape[5] < -1+step:
                    first, second = shape.copy(), shape.copy()
                    first[5] += step; second[5] += 2*step
                    delta = (-3*M+4*matrix(first)-matrix(second))/(2*step)
                elif shape[5] > 1-step:
                    first, second = shape.copy(), shape.copy()
                    first[5] -= step; second[5] -= 2*step
                    delta = (3*M-4*matrix(first)+matrix(second))/(2*step)
                else:
                    hi, lo = shape.copy(), shape.copy()
                    hi[5] += step; lo[5] -= step
                    delta = (matrix(hi)-matrix(lo))/(2*step)
                dM.append(delta)
            derivatives.append(equality_tangent_sensitivity(M, fit_observations, c, dM, exact_rows=list(range(5))))
        last = derivatives[-1]
        last.update(jacobian_step_change=float(np.max(abs(last['jacobian']-derivatives[0]['jacobian']))),
            shape_derivative_steps=[2e-4, 1e-4], shape_coordinates='log five inherited, linear power if present',
            inherited_k_even_bound_active=bool(abs(shape[2]-12.) < 1e-6),
            optimizer_success=data['optimizer']['success'],
            optimizer_message=data['optimizer']['message'], replay_max_error=error,
            fitted_ranges_included=True, statistical_confidence_claimed=False)
        save_json(args.out/(family+'_identifiability.json'), last)
        write_csv(args.out/(family+'_audited_residuals.csv'), residual_table(observations, best))
        residual = np.asarray(best['residuals'])
        row = dict(model=family, loss=best['squared_loss'], profile_count=len(data['profiles']),
            optimizer_success=data['optimizer']['success'], stopping_reason=data['optimizer']['message'],
            source_cohesion_eV=best['predictions'][1], C11_GPa=best['predictions'][2],
            C12_GPa=best['predictions'][3], C44_GPa=best['predictions'][4],
            rank=last['rank'], condition_number=last['condition_number'], screening=shape[5])
        for target_name in ('perfect_Haa', 'perfect_Hxx'):
            i = next(i for i, o in enumerate(observations) if o.name == target_name)
            row[target_name] = best['predictions'][i]
            row[target_name+'_relative_error'] = best['predictions'][i]/observations[i].target-1
        for group in ('v17_', 'v19_new_', 'v19_power_new_'):
            indices = [i for i, o in enumerate(observations) if o.role == 'heldout' and o.name.startswith(group)]
            row[group+'retrospective_rms'] = float(np.sqrt(np.mean(residual[indices]**2)))
        summary.append(row)
        all_rows.extend(dict(model=family, parameter=f'c{i}', value=value) for i, value in enumerate(c))
        all_rows.extend(dict(model=family, parameter=name, value=value) for name, value in zip(
            ('k_scalar', 'k_odd', 'k_even', 'alpha_even', 'k_rank1', 'density_power'), shape))
        print(family, 'exact replay and full equality-tangent sensitivity finished', flush=True)
    write_csv(args.out/'calibration_comparison.csv', summary)
    write_csv(args.out/'parameter_sets.csv', all_rows)
    import csv
    with (args.validation_directory/'interface_curves.csv').open(encoding='utf-8', newline='') as stream:
        curves = list(csv.DictReader(stream))
    with (args.validation_directory/'independent_observations.csv').open(encoding='utf-8', newline='') as stream:
        held = list(csv.DictReader(stream))
    held_summary = []
    for name in dict.fromkeys(row['model'] for row in held):
        for field in ('energy', 'normal_force', 'Haa', 'Hxx'):
            rows = [r for r in held if r['model'] == name and r['field'] == field]
            r = np.array([float(row['normalized_residual']) for row in rows])
            held_summary.append(dict(model=name, field=field, count=len(r), rms=float(np.sqrt(np.mean(r*r))),
                maximum_absolute_normalized_error=float(np.max(abs(r))), within_scale=int(np.sum(abs(r) <= 1))))
    write_csv(args.out/'independent_error_summary.csv', held_summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for name in dict.fromkeys(row['model'] for row in curves):
        rows = [r for r in curves if r['model'] == name and r['path'] == 'opening' and float(r['fraction']) <= 2]
        x = [float(r['fraction']) for r in rows]
        label = name.replace('joint_refinement/', '').replace('source_Al99', 'Al99 source (target only)')
        for ax, field in zip(axes, ('energy_J_m2', 'normal_MPa', 'Haa')):
            ax.plot(x, [float(r[field]) for r in rows], label=label)
            ax.set_xlabel('a / h'); ax.grid(alpha=.2)
    for ax, label in zip(axes, ('Interface energy [J/m²]', 'Ideal coherent traction [MPa]', 'Normal curvature [eV/L0²]')):
        ax.set_ylabel(label)
    axes[0].legend(fontsize=8)
    fig.suptitle('Static calibration comparison; not measured yield or fatigue')
    fig.savefig(args.out/'normal_response_comparison.png', dpi=160)
    plt.close(fig)
    save_json(args.out/'completion.json', dict(completed=True, optimizer_executed_here=False,
        exact_saved_coefficient_replay=True, material_accepted=False, kinetic_calibration=False))


if __name__ == '__main__':
    main()
