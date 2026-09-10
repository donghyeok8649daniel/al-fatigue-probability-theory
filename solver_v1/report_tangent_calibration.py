"""Replay and equality-manifold sensitivity of completed v20 static fits.

No optimization here. Seven source equalities eliminate coefficient directions;
active inequality cones and statistical confidence are not inferred from SVD.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .interface_tangent_calibration import (TangentCalibrationProblem, impose_tangents,
    make_residual_table, shape_selection_status)
from .profiled_identifiability import equality_tangent_sensitivity
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-joint', type=Path, required=True)
    parser.add_argument('--refinement', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--families', nargs='+', choices=('old_family', 'power_family'),
                        default=['old_family', 'power_family'])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--sensitivity-only', action='store_true',
                      help='replay/SVD only; independent response data may still be computing')
    mode.add_argument('--response-only', action='store_true',
                      help='postprocess completed response data without repeating SVD')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh report directory required')
    p = TangentCalibrationProblem(args.parent_joint)
    observations = impose_tangents(p.observations, p.source_tangents)
    exact = [i for i, o in enumerate(observations) if o.role == 'exact']
    lower = np.r_[np.log([.7, 2., 2., 1e-5, 2.]), -1.]
    upper = np.r_[np.log([8., 12., 12., 1e6, 12.]), 1.]
    summary, parameters = [], []
    for family in ([] if args.response_only else args.families):
        data = json.loads((args.refinement/family/'calibration.json').read_bytes())
        if not data['completed'] or data['best'] is None:
            raise ValueError('completed positive-LJ candidate required for report')
        best = data['best']; shape = np.asarray(best['shape']); c = np.asarray(best['coefficients'])
        ndim = 5 if family == 'old_family' else 6
        x = np.r_[np.log(shape[:5]), shape[5]][:ndim]

        def matrix(coords):
            return p.matrix(tuple(np.r_[np.exp(coords[:5]), 0. if ndim == 5 else coords[5]]))

        M = matrix(x)
        replay_error = float(np.max(abs(M@c-best['predictions'])))
        if replay_error > 2e-8:
            raise ArithmeticError('saved vector replay does not agree')
        write_csv(args.out/(family+'_residuals.csv'), make_residual_table(p.observations, observations, best))
        studies = []
        for step in (2e-4, 1e-4):
            derivatives, conventions = [], []
            for j in range(ndim):
                e = np.eye(ndim)[j]*step
                if x[j]-step < lower[j]:
                    derivative = (-3*M+4*matrix(x+e)-matrix(x+2*e))/(2*step)
                    convention = 'forward_second_order_at_bound'
                elif x[j]+step > upper[j]:
                    derivative = (3*M-4*matrix(x-e)+matrix(x-2*e))/(2*step)
                    convention = 'backward_second_order_at_bound'
                else:
                    derivative = (matrix(x+e)-matrix(x-e))/(2*step)
                    convention = 'central_second_order'
                derivatives.append(derivative); conventions.append(convention)
            result = equality_tangent_sensitivity(M, observations, c, derivatives, exact_rows=exact)
            studies.append(result)
        last = studies[-1]
        last.update(step_change_max=float(np.max(abs(last['jacobian']-studies[0]['jacobian']))),
            shape_derivative_steps=[2e-4, 1e-4], derivative_conventions=conventions,
            shape_coordinates='log first five, power linear', replay_error=replay_error,
            physical_kinetics=False, whole_family_identifiability_proved=False)
        save_json(args.out/(family+'_sensitivity.json'), last)
        summary.append(dict(model=family, original_source_loss=best['source_loss_104'],
            common_other_loss=best['other_loss_102'], Haa=best['predictions'][p.tangent_rows['perfect_Haa']],
            Hxx=best['predictions'][p.tangent_rows['perfect_Hxx']],
            rank=last['rank'], dimension=last['jacobian'].shape[1],
            condition_number=last['condition_number'], exact_derivative_residual=last['exact_derivative_residual'],
            optimizer_success=data['optimizer']['success'], optimizer_stop=data['optimizer']['message'],
            optimizer_optimality=data['optimizer'].get('optimality'), **shape_selection_status(data)))
        parameters.extend(dict(model=family, parameter=name, value=value)
            for name, value in zip(('u', 'v', 'A', 'B', 'C', 'D3', 'D1', 'D2', 'D_E', 'K3'), c))
        parameters.extend(dict(model=family, parameter=name, value=value)
            for name, value in zip(('k_scalar', 'k_odd', 'k_even', 'alpha_even', 'k_rank1', 'power'), shape))
        print(family, 'saved-vector replay and seven-equality sensitivity completed', flush=True)
    if not args.response_only:
        write_csv(args.out/'fit_comparison.csv', summary)
        write_csv(args.out/'parameter_sets.csv', parameters)
    if args.sensitivity_only:
        save_json(args.out/'completion.json', dict(completed=True, optimization_executed_here=False,
            sensitivity_replayed=True, response_data_postprocessed=False,
            material_accepted=False, physical_kinetics=False))
        return
    validation_status = json.loads((args.validation/'completion.json').read_bytes())
    if not validation_status['completed']:
        raise ValueError('completed independent response data required')
    with (args.validation/'static_tensor_states.csv').open(encoding='utf-8', newline='') as stream:
        static = list(csv.DictReader(stream))
    reference = {(r['case'], r['fraction']): r for r in static if r['model'] == 'source_Al99_target_only'}
    comparison = []
    primary = dict(pure_interface_normal50MPa='a', pure_interface_shear1_4MPa='x',
                   pure_interface_shear2_4MPa='y')
    for row in static:
        if row['model'] == 'source_Al99_target_only' or abs(float(row['fraction'])) != 1.:
            continue
        source = reference[(row['case'], row['fraction'])]
        for axis in ('a', 'x', 'y'):
            prediction, target = (float(r[f'delta_{axis}_angstrom']) for r in (row, source))
            comparison.append(dict(model=row['model'], case=row['case'], fraction=row['fraction'],
                axis=axis, source_angstrom=target, prediction_angstrom=prediction,
                error_angstrom=prediction-target,
                primary_component_relative_error=(prediction/target-1 if target != 0 and
                    primary.get(row['case']) == axis else None),
                interpretation='all components retained; relative error only on directly forced primary components'))
    write_csv(args.out/'low_stress_comparison.csv', comparison)
    with (args.validation/'interface_curves.csv').open(encoding='utf-8', newline='') as stream:
        curves = list(csv.DictReader(stream))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for name in dict.fromkeys(r['model'] for r in curves):
        rows = [r for r in curves if r['model'] == name and r['path'] == 'opening' and float(r['fraction']) <= 2]
        for ax, field in zip(axes, ('energy_J_m2', 'traction_MPa', 'Haa')):
            ax.plot([float(r['fraction']) for r in rows], [float(r[field]) for r in rows], label=name)
            ax.set_xlabel('a / h'); ax.grid(alpha=.2)
    for ax, label in zip(axes, ('Interface energy [J/m²]', 'Ideal coherent traction [MPa]', 'Haa [eV/L0²]')):
        ax.set_ylabel(label)
    axes[0].legend(fontsize=7)
    fig.suptitle('Exact initial tangents: assess the remaining finite-opening tradeoff')
    fig.savefig(args.out/'opening_tradeoff.png', dpi=160)
    plt.close(fig)
    save_json(args.out/'completion.json', dict(completed=True, optimization_executed_here=False,
        sensitivity_replayed=not args.response_only, response_data_postprocessed=True,
        material_accepted=False, physical_kinetics=False, actual_yield_validated=False))


if __name__ == '__main__':
    main()
