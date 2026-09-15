"""Fit the existing six-shape material family with seven exact anchors eliminated."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from .anchored_shape_calibration import ExactAnchorCoordinates
from .core_interface_compatibility import minimax_compatibility
from .interface_tangent_calibration import NONNEGATIVE
from .profile_shape_sensitivity import difference_stencil
from .run_joint_shape_v31 import decode_shape
from .run_low_frequency_forcing_v29 import sha
from .run_profile_shape_v34 import ShapeProfile
from .run_vector_registry_audit import save_json
from .vector_material_calibration import IDEAL_H


def validation_plan():
    states = [[r*IDEAL_H, 0., 0.] for r in (1.006, 1.061, 1.141, 1.331, 1.681, 1.971)]
    states += [[r*IDEAL_H, x, y] for r, x, y in
               ((1.031, .157, -.053), (1.081, .317, .127),
                (1.181, .417, .237), (.985, .097, .043))]
    return dict(interface_states=states,
        finite_q=[[v, v, 0.] for v in (.2, .5, .7)]+[[v, v, v] for v in (.2, .33, .49)],
        role='declared before this fit; excluded from its loss and selection',
        full_history_blindness_certified=False)


def scan_saturation(parent, resume, out):
    """Check distant existing saturation values with all other shapes fixed."""
    previous = json.loads(Path(resume).read_bytes())
    if not previous['completed'] or not previous['final_profile']['strictly_positive_LJ']:
        raise ValueError('completed positive-LJ previous profile required')
    study = ShapeProfile(parent, out)
    study.upper[[1, 2, 4]] = np.log(24.)
    start = np.asarray(previous['final_shape_coordinates'], float)
    save_json(study.out/'validation_plan.json', validation_plan())
    values = [float(np.exp(start[3])), 1e-5, 1e-2, 1., 1e2, 1e4, 1e5, 1e6]
    save_json(study.out/'saturation_scan_definition.json', dict(resume_sha256=sha(resume),
        saturation_values=values, other_shapes_fixed=True, original_saturation_bounds=[1e-5, 1e6],
        target_scales_changed=False, new_energy_terms=False, production_changed=False))
    records, profiles = [], []
    for alpha in values:
        x = start.copy()
        x[3] = np.log(alpha)
        profile = study.profile(tuple(x))
        if not profiles and abs(profile['minimax_normalized_error']-
                previous['final_profile']['minimax_normalized_error']) > 1e-7:
            raise ArithmeticError('starting coefficient profile did not replay')
        profiles.append(profile)
        records.append(dict(shape_coordinates=x, shape=decode_shape(x, 1., True),
            profile_eta=profile['minimax_normalized_error'], positive_LJ=profile['strictly_positive_LJ'],
            exact_residual=profile['equality_residual'], elapsed_seconds=time.perf_counter()-study.began))
        save_json(study.out/'scan_profiles.json', profiles)
    index = min((i for i, p in enumerate(profiles) if p['strictly_positive_LJ']),
                key=lambda i: profiles[i]['minimax_normalized_error'])
    final = records[index]
    save_json(study.out/'joint_summary.json', dict(completed=True, optimizer_success=None,
        message='deterministic saturation scan; no nonlinear optimizer', iterations=0,
        initial_eta=profiles[0]['minimax_normalized_error'], final_profile=profiles[index],
        final_shape_coordinates=final['shape_coordinates'], final_shape=final['shape'], records=records,
        lower=study.lower, upper=study.upper, elapsed_seconds=time.perf_counter()-study.began,
        global_shape_optimum_certified=False, material_accepted=False, physical_time_calibrated=False,
        production_changed=False))


def run(parent, resume, out, *, maxiter=30, radius=.35, saturation_global=False,
        shape_step=1e-4):
    if maxiter < 1 or not np.isfinite(radius) or not 0 < radius <= .4:
        raise ValueError('positive iteration budget and local log-radius at most .4 required')
    if not np.isfinite(shape_step) or not 0 < shape_step <= 1e-3:
        raise ValueError('positive shape-difference step at most 1e-3 required')
    previous = json.loads(Path(resume).read_bytes())
    if not previous['completed'] or not previous['final_profile']['strictly_positive_LJ']:
        raise ValueError('completed positive-LJ previous profile required')
    study = ShapeProfile(parent, out)
    study.upper[[1, 2, 4]] = np.log(24.)
    start = np.asarray(previous['final_shape_coordinates'], float)
    lower, upper = study.lower.copy(), study.upper.copy()
    lower[:5] = np.maximum(lower[:5], start[:5]-radius)
    upper[:5] = np.minimum(upper[:5], start[:5]+radius)
    if saturation_global:
        # This existing positive shape does not enter the isolated-neighbour
        # decay condition. Retain its ORIGINAL full search bounds.
        lower[3], upper[3] = study.lower[3], study.upper[3]
    # Guarantee the existing isolated-neighbour admissibility throughout this box.
    if 2*np.exp(lower[4])-np.exp(upper[0]) <= 0:
        raise ValueError('local search box does not preserve the isolated-neighbour limit')
    study.lower, study.upper = lower, upper
    study.start = start
    save_json(study.out/'validation_plan.json', validation_plan())
    initial_profile = study.profile(tuple(start))
    if abs(initial_profile['minimax_normalized_error']-
           previous['final_profile']['minimax_normalized_error']) > 1e-7:
        raise ArithmeticError('previous best coefficient profile did not replay')
    c0 = np.asarray(initial_profile['coefficients'])
    D = np.maximum(abs(c0), initial_profile['coefficient_scales'])
    targets, scales = (np.array([getattr(o, key) for o in study.obs]) for key in ('target', 'scale'))
    chart = ExactAnchorCoordinates(study.matrix(tuple(start)), targets, scales,
                                   initial_profile['exact_rows'], D)
    nx = len(start)
    initial = np.r_[start, (c0/D)[chart.free], initial_profile['minimax_normalized_error']]
    save_json(study.out/'anchored_definition.json', dict(resume_sha256=sha(resume),
        start=start, lower=lower, upper=upper, local_log_radius=radius,
        dependent_coefficients=chart.dependent, free_coefficients=chart.free,
        coefficient_scales=D, exact_anchors=len(chart.exact), independent_coefficients=len(chart.free),
        saturation_uses_original_full_bounds=saturation_global,
        maxiter=maxiter, matrix_shape_difference_step=shape_step,
        derivative='implicit differentiation of the exact anchor solve',
        exact_anchor_arithmetic_budget=1e-6, target_scales_changed=False,
        new_energy_terms=False, production_changed=False))

    @lru_cache(maxsize=8)
    def derivative(x):
        return np.array([sum(w*study.matrix(tuple(v)) for v, w in
            difference_stencil(np.array(x), i, lower, upper, shape_step)) for i in range(nx)])

    def evaluate(z, jacobian=False):
        x = tuple(z[:nx])
        M = study.matrix(x)
        dM = derivative(x) if jacobian else np.zeros((nx, *M.shape))
        return chart.evaluate(M, dM, z[nx:-1], z[-1], study.rows, NONNEGATIVE)

    # A bounded direction checks the complete implicit constraint Jacobian using
    # fresh analytic observable matrices, including at the exponent boundary.
    first = evaluate(initial, True)
    direction = np.r_[.2, -.3, -.2, -.15, .1, .2, np.full(len(chart.free), .07), .11]
    audits = []
    for step in (2e-4, 1e-4):
        forward = evaluate(initial+step*direction)['constraints']
        twice = evaluate(initial+2*step*direction)['constraints']
        finite = (-3*first['constraints']+4*forward-twice)/(2*step)
        analytic = first['jacobian']@direction
        error = float(np.linalg.norm(finite-analytic, np.inf))
        audits.append(dict(step=step, maximum_absolute_error=error,
                           relative_error=error/max(1., np.linalg.norm(analytic, np.inf))))
    save_json(study.out/'derivative_audit.json', dict(completed=True, audits=audits,
        initial_exact_residual=first['exact_residual'], initial_pivot_condition=first['pivot_condition']))
    if max(r['relative_error'] for r in audits) > 2e-5:
        raise ArithmeticError('actual implicit calibration Jacobian failed its independent check')

    records = []
    best = dict(shape=start.copy(), profile=initial_profile)

    def checkpoint(z):
        actual = evaluate(z)
        profile = study.profile(tuple(z[:nx]))
        if profile['strictly_positive_LJ'] and profile['minimax_normalized_error'] < best['profile']['minimax_normalized_error']:
            best.update(shape=z[:nx].copy(), profile=profile)
        row = dict(shape_coordinates=z[:nx], joint_eta=float(z[-1]),
            profile_eta=profile['minimax_normalized_error'], exact_residual=actual['exact_residual'],
            inequality_violation=float(max(0., -min(actual['constraints']))),
            pivot_condition=actual['pivot_condition'], positive_LJ=profile['strictly_positive_LJ'],
            elapsed_seconds=time.perf_counter()-study.began)
        records.append(row)
        save_json(study.out/'anchored_checkpoint.json', dict(completed=False, records=records, best=best))
        print('anchored', len(records), 'eta', row['profile_eta'], 'exact', row['exact_residual'], flush=True)

    checkpoint(initial)
    objective_jac = np.r_[np.zeros(len(initial)-1), 1.]
    result = minimize(lambda z:z[-1], initial, jac=lambda z:objective_jac,
        bounds=list(zip(lower, upper))+[(None, None)]*len(chart.free)+[(0., None)],
        constraints=[dict(type='ineq', fun=lambda z:evaluate(z)['constraints'],
                          jac=lambda z:evaluate(z, True)['jacobian'])],
        method='SLSQP', callback=checkpoint, options=dict(maxiter=maxiter, ftol=1e-9))
    checkpoint(result.x)
    # Independently restore and certify coefficients at the best *training* shape.
    final = minimax_compatibility(study.matrix(tuple(best['shape'])), study.obs,
                                   study.rows, nonnegative=NONNEGATIVE)
    final_matrix = study.matrix(tuple(best['shape']))
    final_c = np.asarray(final['coefficients'])
    final_z = np.r_[best['shape'], (final_c/D)[chart.free], final['minimax_normalized_error']]
    final_check = evaluate(final_z, True)
    save_json(study.out/'final_matrix.json', dict(matrix=final_matrix, derivatives=derivative(tuple(best['shape'])),
        target=targets, scales=scales, residuals=(final_matrix@final_c-targets)/scales))
    summary = dict(completed=True, optimizer_success=bool(result.success), message=str(result.message),
        iterations=int(result.nit), initial_eta=initial_profile['minimax_normalized_error'],
        final_profile=final, final_shape_coordinates=best['shape'],
        final_shape=decode_shape(best['shape'], 1., True), records=records,
        maximum_iterate_exact_residual=max(r['exact_residual'] for r in records),
        final_exact_residual=final_check['exact_residual'],
        best_training_point_is_optimizer_endpoint=bool(np.array_equal(best['shape'], result.x[:nx])),
        lower=lower, upper=upper, elapsed_seconds=time.perf_counter()-study.began,
        local_shape_stationarity_certified=False, material_accepted=False,
        physical_time_calibrated=False, production_changed=False)
    save_json(study.out/'joint_summary.json', summary)
    print(json.dumps(dict(initial_eta=summary['initial_eta'], final_eta=final['minimax_normalized_error'],
        optimizer_success=summary['optimizer_success'], message=summary['message'],
        elapsed_seconds=summary['elapsed_seconds']), indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    for name in ('parent', 'resume', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--maxiter', type=int, default=30)
    parser.add_argument('--radius', type=float, default=.35)
    parser.add_argument('--saturation-global', action='store_true')
    parser.add_argument('--scan-saturation', action='store_true')
    parser.add_argument('--shape-step', type=float, default=1e-4)
    args = parser.parse_args()
    if args.scan_saturation:
        scan_saturation(args.parent, args.resume, args.out)
    else:
        run(args.parent, args.resume, args.out, maxiter=args.maxiter, radius=args.radius,
            saturation_global=args.saturation_global, shape_step=args.shape_step)
