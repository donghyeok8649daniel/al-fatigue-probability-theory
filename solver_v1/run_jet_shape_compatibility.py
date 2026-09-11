"""Deterministic radial-shape minimax audit; no new energy or adopted fit.

Find whether varying existing radial shapes can reduce a certified FIXED-
shape jet incompatibility. Local Powell termination is not global optimality;
coefficient feasibility is not sampled or full bulk stability. Complete
interface/source-core re-evaluation is required after a diagnostic search.
"""
import argparse
from dataclasses import asdict
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from .core_interface_compatibility import ExistingJetSubset, minimax_compatibility
from .current_core_coefficient_basis import NAMES
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .run_current_material_core import ROOT, load_current_material
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


class SearchBudgetReached(RuntimeError):
    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--model', type=Path, default=ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json')
    parser.add_argument('--max-fev', type=int, default=180)
    parser.add_argument('--max-seconds', type=float, default=1200.)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh minimax-shape study required')
    if args.max_fev < 1 or not np.isfinite(args.max_seconds) or args.max_seconds <= 0:
        raise ValueError('positive explicit evaluation and wall budgets required')
    began = time.perf_counter()
    model, _, binding = load_current_material(args.model)
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    original = impose_tangents(problem.observations, problem.source_tangents)
    extra = {'v19_new_4_Haa', 'v19_power_new_2_Haa'}
    selected = [i for i, o in enumerate(original)
                if (o.role == 'fit' and o.units == 'eV/L0^2') or o.name in extra]
    exact = [i for i, o in enumerate(original) if o.role == 'exact']
    subset = ExistingJetSubset(problem, exact+selected)
    observations = [original[i] for i in subset.rows]
    local_selected = [subset.rows.index(i) for i in selected]
    # Previously declared numerical radial windows, NOT measured Al ranges.
    lower = np.log([.7, 2., 2., 1e-5, 2.])
    upper = np.log([8., 12., 12., 1e6, 12.])
    shape0 = np.asarray(model.screened_shape, float)
    save_json(args.out/'definition.json', dict(starting_HEAD='d1ac6695518d3738b8303687d8c6a7b58877203e',
        candidate_sha256=binding['parameter_sha256'], source_sha256=problem.source.reference.sha256,
        original_shape=shape0, subset_rows=subset.rows,
        observations=[asdict(o) for o in observations], selected_rows=local_selected,
        promoted_previously_inspected_constraints=sorted(extra), nonnegative=NONNEGATIVE,
        coefficient_names=NAMES, log_bounds=[lower, upper], screening_exponent_fixed=shape0[5],
        bounds_are_inherited_search_windows_not_material_measurements=True,
        objective='minimum possible worst normalized selected source-jet error at each shape',
        core_forces_in_objective=False, spectral_constraints_included=False,
        optimizer='deterministic bounded Powell', max_fev=args.max_fev, max_seconds=args.max_seconds,
        no_new_energy_term=True, physical_yield_used=False, material_accepted=False))
    # Before optimization require exact replay of full-matrix rows, not a
    # silent bulk-energy/GPa conversion mismatch in the cheaper subset.
    full = problem.matrix(tuple(shape0))
    sub = subset.matrix(tuple(shape0))
    replay = float(np.max(abs(sub-full[subset.rows])))
    if replay > 2e-10:
        raise ArithmeticError(f'subset/full Bessel matrix mismatch: {replay}')
    profiles = []
    def fun(x):
        if time.perf_counter()-began > args.max_seconds:
            raise SearchBudgetReached('declared wall budget reached')
        tick = time.perf_counter()
        shape = np.r_[np.exp(x), shape0[5]]
        try:
            matrix = subset.matrix(tuple(shape))
            fit = minimax_compatibility(matrix, observations, local_selected, nonnegative=NONNEGATIVE)
            if not fit['completed']:
                raise ArithmeticError(f'no verified equality/sign LP: {fit["message"]}')
            profile = dict(fit, shape=shape)
            value = float(fit['minimax_normalized_error'])
        except (ValueError, ArithmeticError) as error:
            profile = dict(completed=False, shape=shape, error=str(error))
            value = np.inf  # an unavailable state, not a clipped physical penalty
        profile['evaluation_seconds'] = time.perf_counter()-tick
        profiles.append(profile)
        save_json(args.out/'checkpoint.json', dict(completed=False, profiles=profiles))
        print(f'jet profile {len(profiles)} eta={value:.7g}, '
              f'positive LJ={profile.get("strictly_positive_LJ")}, '
              f'{profile["evaluation_seconds"]:.2f}s', flush=True)
        return value
    try:
        result = minimize(fun, np.log(shape0[:5]), method='Powell', bounds=list(zip(lower, upper)),
            options=dict(maxfev=args.max_fev, xtol=2e-4, ftol=2e-5))
        status = dict(success=bool(result.success), message=str(result.message), nfev=result.nfev,
                      last_accepted_coordinates=result.x, last_accepted_loss=float(result.fun))
    except (SearchBudgetReached, ArithmeticError, ValueError) as error:
        status = dict(success=False, message=str(error), error_type=type(error).__name__)
    eligible = [p for p in profiles if p['completed'] and p['strictly_positive_LJ']]
    best = min(eligible, key=lambda p: p['minimax_normalized_error']) if eligible else None
    rows = [dict(profile=i, completed=p['completed'], eta=p.get('minimax_normalized_error'),
                 positive_LJ=p.get('strictly_positive_LJ'), elapsed_seconds=p['evaluation_seconds'],
                 error=p.get('error')) for i, p in enumerate(profiles)]
    write_csv(args.out/'profile_summary.csv', rows)
    save_json(args.out/'completion.json', dict(completed=True, optimizer=status, best=best,
        matrix_replay_error=replay, profiles_evaluated=len(profiles),
        positive_profiles=len(eligible), elapsed_seconds=time.perf_counter()-began,
        coefficient_predictions_use_subset_order=True, material_accepted=False,
        whole_family_impossibility_proved=False, physical_yield_validated=False,
        physical_Hz=False, production_changed=False))


if __name__ == '__main__':
    main()
