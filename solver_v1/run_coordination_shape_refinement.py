"""Matched actual radial/shape refits of old and screened analytic families.

No frozen-pair CONTROL; all ten energy coefficients are solved subject to the
original exact bulk/spectral constraints. Old data stay development/excluded.
Two deterministic starts differ only by including the extra power shape.
"""
import argparse
from dataclasses import replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares

from .coordination_screening import CoordinationScreenedCache, CoordinationScreenedBulk
from .vector_material_calibration import MaterialObservation
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix
from .tail_constrained_material import spectral_profile
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def cubic_observations(observations):
    """The fit's GPa targets/scales, not the input strain-mode metadata."""
    _, obs = cubic_metric_problem(np.zeros((len(observations), 8)), observations)
    return [replace(o, role='exact') if i < 5 else o for i, o in enumerate(obs)]


def residual_table(observations, best):
    """Audit every exported residual against its displayed target and units."""
    obs = cubic_observations(observations)
    predictions = np.asarray(best['predictions'])
    residuals = (predictions-np.array([o.target for o in obs]))/np.array([o.scale for o in obs])
    np.testing.assert_allclose(residuals, best['residuals'], rtol=1e-10, atol=1e-10)
    return [dict(observable=o.name, role=o.role,
        validation_status='retrospective' if o.role == 'heldout' else 'development',
        target=o.target, prediction=float(p), units=o.units, scale=o.scale,
        normalized_residual=float(r)) for o, p, r in zip(obs, predictions, residuals)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile-directory', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-nfev', type=int, default=20)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh joint-refit directory required')
    began = time.perf_counter()
    definition = json.loads((args.profile_directory/'definition.json').read_bytes())
    if definition.get('screening_law') != 'power':
        raise ValueError('completed power-law comparison required')
    profile = json.loads((args.profile_directory/'free_pair/calibration.json').read_bytes())
    if not profile['completed']:
        raise ValueError('completed initialization required')
    observations = [MaterialObservation(**{**o,
        'state': tuple(o['state']) if o['state'] is not None else None,
        'jet_weights': tuple(o['jet_weights']) if o['jet_weights'] is not None else None})
        for o in definition['observations']]
    initial = np.asarray(profile['shape'])
    # Bounds are inherited numerical windows, not measured material ranges.
    lower = np.log([.7, 2., 2., 1e-5, 2.]); upper = np.log([8., 12., 12., 1e6, 12.])
    bounds = np.array([np.r_[lower, -1.], np.r_[upper, 1.]])
    save_json(args.out/'definition.json', definition | dict(
        joint_refinement=True, actual_radial_optimization=True, fixed_shape_control=False, fixed_pair_control=None,
        previously_inspected_validation=True, optimizer='bounded least_squares; numerical shape derivatives',
        starts=[dict(model='old_family', shape=np.r_[initial[:5], 0.]),
                dict(model='power_family', shape=initial)],
        profile_definition_sha256=hashlib.sha256((args.profile_directory/'definition.json').read_bytes()).hexdigest(),
        max_nfev=args.max_nfev, transformed_bounds=bounds,
        transform='log first five, linear final power; old family fixes power=0'))
    cache = CoordinationScreenedCache(observations, law='power')
    summaries = []
    for family in ('old_family', 'power_family'):
        profiles, failures = [], []
        ndim = 5 if family == 'old_family' else 6

        @lru_cache(maxsize=256)
        def evaluate(coordinates):
            started = time.perf_counter()
            shape = np.r_[np.exp(coordinates[:5]), 0. if ndim == 5 else coordinates[5]]
            try:
                raw = cache.matrix(shape)
                M, obs = cubic_metric_problem(raw[:, :8], observations)
                extra = raw[:, 8:].copy(); extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
                M = np.column_stack([M, extra])
                obs = [replace(o, role='exact') if i < 5 else o for i, o in enumerate(obs)]
                operators, errors = [], []
                for stretch in definition['stability_stretches']:
                    bulk = CoordinationScreenedBulk(shape, radius=definition['radius_over_L0'],
                                                    stretch=stretch, law='power')
                    cols, tails = map(np.asarray, zip(*(bulk.evaluate(q) for q in definition['wavepoints_cubic'])))
                    operators.extend(cols); errors.extend(tails)
                fit = spectral_profile(M, obs, operators, errors, nonnegative=(0, 1, 2, 4, 5, 8, 9))
            except (ValueError, ArithmeticError) as error:
                failures.append(dict(shape=shape, error=str(error)))
                save_json(args.out/family/'failed_profiles.json', failures)
                raise
            row = dict(shape=shape, elapsed_seconds=time.perf_counter()-started, **fit)
            profiles.append(row)
            save_json(args.out/family/'checkpoint.json', dict(completed=False, profiles=profiles))
            print(family, len(profiles), f"loss={fit['squared_loss']:.7f}",
                  f"shape={shape}", f"{row['elapsed_seconds']:.2f}s", flush=True)
            return row

        def fun(x):
            row = evaluate(tuple(x))
            return np.asarray(row['residuals'])[row['selected_rows']]

        x0 = np.r_[np.log(initial[:5]), initial[5]][:ndim]
        optimizer = None
        try:
            run = least_squares(fun, x0, bounds=bounds[:, :ndim], max_nfev=args.max_nfev,
                                jac='2-point', diff_step=1e-4, ftol=1e-7, xtol=1e-7, gtol=2e-6)
            optimizer = dict(success=bool(run.success), message=str(run.message), nfev=run.nfev,
                njev=run.njev, last_accepted_coordinates=run.x, last_accepted_loss=float(run.fun@run.fun),
                optimality=run.optimality)
        except (ValueError, ArithmeticError) as error:
            optimizer = dict(success=False, message=str(error), numerical_stop=True)
        if not profiles:
            raise ArithmeticError(f'{family}: no verified coefficient profile')
        minimum = min(profiles, key=lambda row: row['squared_loss'])
        eligible = [row for row in profiles if row['strictly_positive_LJ']]
        best = min(eligible, key=lambda row: row['squared_loss']) if eligible else None
        save_json(args.out/family/'calibration.json', dict(completed=True, family=family,
            best=best, minimum_profile_including_zero_pair=minimum, profiles=profiles, optimizer=optimizer,
            failed_profiles=failures, screening_law='power', material_accepted=False,
            admissible_solution_found=bool(eligible), new_heldout_was_used_for_selection=False))
        if best is not None:
            write_csv(args.out/family/'residuals.csv', residual_table(observations, best))
        summaries.append(dict(family=family, profiles=len(profiles), failures=len(failures),
            best_positive_LJ_loss=best['squared_loss'] if best else None,
            unconstrained_pair_closure_loss=minimum['squared_loss'], optimizer_success=optimizer['success'],
            termination=optimizer['message']))
        write_csv(args.out/'summary.csv', summaries)
        # Keep the cache of inherited basis columns, but do not keep stale
        # coefficient profiles in the next family's function cache.
        evaluate.cache_clear()
    save_json(args.out/'completion.json', dict(completed=True, summary=summaries,
        elapsed_seconds=time.perf_counter()-began, material_accepted=False,
        kinetics_calibrated=False, whole_family_infeasibility_proved=False))


if __name__ == '__main__':
    main()
