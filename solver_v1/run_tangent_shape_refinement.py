"""Actual shape optimization with exact bulk AND source interface tangents.

This is not a replay or a new energy family. The coefficient closure u,v>=0
is a diagnostic search domain; only u,v>0 candidates may be instantiated.
No independent validation row enters the objective or candidate selection.
"""
import argparse
from functools import lru_cache
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares

from .interface_tangent_calibration import TangentCalibrationProblem, tangent_summary
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-joint', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-nfev', type=int, default=30)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh output directory required; preserve earlier fits')
    if args.max_nfev < 1:
        raise ValueError('positive optimizer evaluation budget required')
    started = time.perf_counter()
    problem = TangentCalibrationProblem(args.parent_joint)
    lower = np.r_[np.log([.7, 2., 2., 1e-5, 2.]), -1.]
    upper = np.r_[np.log([8., 12., 12., 1e6, 12.]), 1.]
    save_json(args.out/'definition.json', problem.definition | dict(
        actual_radial_optimization=True, tangent_constraints=problem.source_tangents,
        starts={name: row['shape'] for name, row in problem.parents.items()},
        bounds_transformed=[lower, upper], bounds_inherited_from_v19=True,
        max_nfev=args.max_nfev, transform='log first five, linear power',
        coefficient_zero_LJ_closure_is_material=False,
        optimizer='bounded least_squares; numerical shape derivatives',
        selection='minimum original source loss among verified positive-LJ evaluated profiles',
        global_optimality_claim=False))
    summaries = []
    for family in ('old_family', 'power_family'):
        ndim = 5 if family == 'old_family' else 6
        profiles, failures = [], []

        @lru_cache(maxsize=256)
        def evaluate(coordinates):
            tick = time.perf_counter()
            shape = np.r_[np.exp(coordinates[:5]), 0. if ndim == 5 else coordinates[5]]
            fit = problem.fit(shape, problem.source_tangents)
            if not fit['completed']:
                failures.append(fit)
                save_json(args.out/family/'failed_profiles.json', failures)
                raise ArithmeticError(f"unverified profile: {fit['classification']}")
            fit.pop('residual_table')  # Recompute/export only selected final rows.
            fit['elapsed_seconds'] = time.perf_counter()-tick
            profiles.append(fit)
            save_json(args.out/family/'checkpoint.json', dict(completed=False, profiles=profiles))
            print(family, len(profiles), f"loss={fit['source_loss_104']:.7f}",
                  f"positive_LJ={fit['strictly_positive_LJ']}", shape,
                  f"{fit['elapsed_seconds']:.2f}s", flush=True)
            return fit

        def fun(x):
            fit = evaluate(tuple(x))
            return np.asarray(fit['residuals'])[fit['selected_rows']]

        initial = np.asarray(problem.parents[family]['shape'])
        x0 = np.r_[np.log(initial[:5]), initial[5]][:ndim]
        try:
            run = least_squares(fun, x0, bounds=(lower[:ndim], upper[:ndim]),
                max_nfev=args.max_nfev, jac='2-point', diff_step=1e-4,
                ftol=1e-7, xtol=1e-7, gtol=2e-6)
            optimizer = dict(success=bool(run.success), message=str(run.message),
                nfev=run.nfev, njev=run.njev, last_accepted_coordinates=run.x,
                last_accepted_loss=float(run.fun@run.fun), optimality=run.optimality,
                active_bounds=run.active_mask)
        except (ValueError, ArithmeticError) as error:
            optimizer = dict(success=False, message=str(error), numerical_stop=True)
        if not profiles:
            raise ArithmeticError(f'{family}: no verified profile; inspect failure log')
        minimum = min(profiles, key=lambda r: r['source_loss_104'])
        eligible = [r for r in profiles if r['strictly_positive_LJ']]
        best = min(eligible, key=lambda r: r['source_loss_104']) if eligible else None
        if best is not None:
            verified = problem.fit(best['shape'], problem.source_tangents)
            if not verified['completed']:
                raise ArithmeticError('selected profile failed deterministic replay')
            np.testing.assert_allclose(verified['coefficients'], best['coefficients'], atol=0, rtol=0)
            write_csv(args.out/family/'residuals.csv', verified.pop('residual_table'))
        save_json(args.out/family/'calibration.json', dict(completed=True, family=family,
            best=best, minimum_profile_including_zero_pair=minimum, profiles=profiles,
            optimizer=optimizer, failed_profiles=failures, screening_law='power',
            material_accepted=False, admissible_solution_found=bool(eligible),
            new_heldout_was_used_for_selection=False))
        row = tangent_summary('shape_refit', family, best or dict(completed=False,
            classification='no_positive_LJ_profile_found'), problem)
        row.update(profiles=len(profiles), failures=len(failures),
            minimum_closure_loss=minimum['source_loss_104'],
            optimizer_success=optimizer['success'], termination=optimizer['message'])
        summaries.append(row)
        write_csv(args.out/'summary.csv', summaries)
        evaluate.cache_clear()
    save_json(args.out/'completion.json', dict(completed=True, summary=summaries,
        elapsed_seconds=time.perf_counter()-started, energy_family_changed=False,
        material_accepted=False, physical_kinetics=False,
        whole_family_infeasibility_proved=False))


if __name__ == '__main__':
    main()
