"""Full-data and frozen-core checks AFTER the v23 selected-jet search.

The local search need not converge; a saved positive trial can be inspected
without calling it an optimized material. Recompute the full Bessel matrix,
sampled bulk operators/tails and excluded source-core forces. All outputs
remain rejected research candidates unless separate material gates pass.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .coordination_screening import CoordinationScreenedInterface
from .core_interface_compatibility import interval_profile, minimax_compatibility
from .frozen_core_force_target import FrozenCoreForceTarget
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .report_current_material_core import restore_case
from .run_current_material_core import ROOT, load_current_material
from .run_source_core_reference import load_source_material, build_source_core
from .tail_constrained_material import spectral_profile
from .vector_material_calibration import MaterialObservation
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh post-search validation directory required')
    began = time.perf_counter()
    data = json.loads((args.search/'completion.json').read_bytes())
    definition = json.loads((args.search/'definition.json').read_bytes())
    if not data['completed'] or data['best'] is None:
        raise ValueError('completed diagnostic search with a positive trial required')
    shape = tuple(data['best']['shape'])
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    if problem.source.reference.sha256 != definition['source_sha256']:
        raise ValueError('source binding changed after search')
    obs = impose_tangents(problem.observations, problem.source_tangents)
    selected = [i for i, o in enumerate(obs) if o.role == 'fit']
    jets = [i for i, o in enumerate(obs) if o.role != 'exact' and o.units == 'eV/L0^2']
    critical = [next(i for i, o in enumerate(obs) if o.name == name)
                for name in ('saddle_Hxx', 'v19_new_4_Haa')]
    source, tensor, binding = load_source_material()
    source_cases = [('R8_development', 'escaped_plus_R8r5'),
                    ('R16_excluded_core_configuration', 'stable_plus_R16r5')]
    frozen, source_hashes = {}, {}
    for label, case in source_cases:
        path = ROOT/'results/current_material_core_v22/source_reference'/case
        core, field, *_ = restore_case(path, source, tensor, binding, core_builder=build_source_core)
        frozen[label] = FrozenCoreForceTarget(core, field, radius=2.)
        source_hashes[label] = dict(case=path.relative_to(ROOT).as_posix(),
            state_sha256=hashlib.sha256((path/'state.csv').read_bytes()).hexdigest())
    save_json(args.out/'definition.json', dict(search=str(args.search.relative_to(ROOT) if args.search.is_absolute() else args.search),
        search_completion_sha256=hashlib.sha256((args.search/'completion.json').read_bytes()).hexdigest(),
        source_sha256=binding['parameter_sha256'], source_cases=source_hashes,
        shape=shape, shape_optimizer=data['optimizer'],
        selected_trial_is_not_claimed_converged=True, coefficients_refitted=True,
        new_blind_target_data=False, previously_inspected_data_explicit=True,
        prior_v22_joint_objective_unchanged=True, physical_yield_used=False,
        no_new_energy_term=True, material_accepted=False))
    M = problem.matrix(shape)
    target = np.array([o.target for o in obs]); scales = np.array([o.scale for o in obs])
    replay = float(np.max(abs(M[definition['subset_rows']]@data['best']['coefficients']
                             -np.asarray(data['best']['predictions']))))
    if replay > 1e-8:
        raise ArithmeticError('full matrix did not reproduce subset trial predictions')
    operators, tails = problem.operators(shape)
    old, _, _ = load_current_material()
    old_M = problem.matrix(tuple(old.screened_shape))
    old_F = frozen['R8_development'].coefficient_matrix(old)
    scale_static = float(np.linalg.norm(((old_M@old.coefficients-target)/scales)[selected]))
    scale_core = float(np.linalg.norm(old_F['design']@old.coefficients-old_F['target']))
    basis_model = CoordinationScreenedInterface(shape, data['best']['coefficients'], law='power')
    force_designs = {label: f.coefficient_matrix(basis_model) for label, f in frozen.items()}
    F = force_designs['R8_development']
    joint_M = np.vstack([M, F['design']])
    joint_obs = [replace(o, scale=o.scale*scale_static) if o.role == 'fit' else o for o in obs]
    joint_obs += [MaterialObservation(f'core_gradient_{i}', float(g), scale_core, 'eV/L0', 'fit')
                  for i, g in enumerate(F['target'])]
    profiles = [('selected_minimax_trial', data['best'])]
    for name, matrix, observations, boxes in (
            ('static_LS', M, obs, []), ('joint_LS', joint_M, joint_obs, []),
            ('joint_critical_boxes', joint_M, joint_obs, critical)):
        try:
            fit = (interval_profile(matrix, observations, boxes, 1., nonnegative=NONNEGATIVE,
                        operators=operators, tails=tails, interval_scales=scales[boxes]) if boxes else
                   spectral_profile(matrix, observations, operators, tails, nonnegative=NONNEGATIVE))
            profiles.append((name, dict(fit, completed=True)))
        except (ArithmeticError, ValueError) as error:
            profiles.append((name, dict(completed=False, error=str(error))))
    summaries, residuals, forces, details = [], [], [], []
    for name, fit in profiles:
        if not fit['completed']:
            details.append(dict(profile=name, **fit))
            continue
        c = np.asarray(fit['coefficients'])
        pred = M@c; r = (pred-target)/scales
        margin = np.linalg.eigvalsh(np.einsum('i,qijk->qjk', c, operators))[:, 0]-tails@abs(c)
        direct_error = None
        if np.all(c[:2] > 0):
            model = CoordinationScreenedInterface(shape, c, law='power')
            errors = []
            for i in critical:
                value = model.evaluate(obs[i].state)
                jet = np.r_[value.energy, value.gradient, value.hessian[0], value.hessian[1, 1:], value.hessian[2, 2]]
                errors.append(abs(float(jet@obs[i].jet_weights)-pred[i]))
            direct_error = max(errors)
            if direct_error > 2e-8:
                raise ArithmeticError('coefficient matrix differs from actual analytic energy jets')
        row = dict(profile=name, static_loss=float(r[selected]@r[selected]),
            maximum_115_jet_error=float(np.max(abs(r[jets]))),
            worst_jet_name=obs[jets[int(np.argmax(abs(r[jets]))) ]].name,
            saddle_Hxx=float(pred[critical[0]]), opening_Haa=float(pred[critical[1]]),
            positive_LJ=bool(np.all(c[:2] > 0)), minimum_robust_bulk_margin=float(np.min(margin)),
            exact_residual=float(np.max(abs(r[[i for i, o in enumerate(obs) if o.role == 'exact']]))),
            independent_jet_replay_error=direct_error, material_accepted=False)
        for label, design in force_designs.items():
            error = design['design']@c-design['target']
            rms = float(np.sqrt(np.mean(error*error)))
            row[f'{label}_force_RMS_eV_L0'] = rms
            forces.append(dict(profile=name, configuration=label, force_RMS_eV_L0=rms,
                force_maximum_eV_L0=float(np.max(abs(error))),
                used_for_coefficient_objective=name.startswith('joint') and label == 'R8_development',
                shape_objective_used_core_forces=False, actual_atomic_relaxation=False))
        for i, (o, y, ri) in enumerate(zip(obs, pred, r)):
            residuals.append(dict(profile=name, observable=o.name, previous_role=o.role,
                units=o.units, target=o.target, scale=o.scale, prediction=float(y),
                normalized_error=float(ri), inspected_before=True,
                source_box_used=name == 'joint_critical_boxes' and i in critical))
        details.append(dict(profile=name, original_profile=fit, full_predictions=pred,
            coefficients=c, shape=shape, full_spectral_margin=margin,
            no_physical_adoption=True))
        if row['positive_LJ'] and row['minimum_robust_bulk_margin'] >= 0:
            # Snapshot stores FULL row order, not the LP's selected-row order.
            snapshot_fit = dict(fit, coefficients=c, shape=shape, predictions=pred,
                strictly_positive_LJ=True, material_accepted=False)
            save_json(args.out/name/'research_candidate_snapshot.json', dict(completed=True,
                best=snapshot_fit, screening_law='power', material_accepted=False,
                shape_optimizer_converged=False, for_static_research_only=True,
                production_changed=False, source_sha256=binding['parameter_sha256']))
        summaries.append(row)
        print(f'{name}: static={row["static_loss"]:.7g}, maxjet={row["maximum_115_jet_error"]:.7g}, '
              f'margin={row["minimum_robust_bulk_margin"]:.6g}', flush=True)
    lp = minimax_compatibility(M, obs, jets, nonnegative=NONNEGATIVE)
    write_csv(args.out/'profile_summary.csv', summaries)
    write_csv(args.out/'all_observable_residuals.csv', residuals)
    write_csv(args.out/'source_force_summary.csv', forces)
    save_json(args.out/'completion.json', dict(completed=True, profiles=details,
        full_115_jet_fixed_shape_minimax=lp, subset_full_prediction_replay_error=replay,
        actual_profiles=len(profiles), failed_profiles=sum(not f['completed'] for _, f in profiles),
        original_static_squared_scale=scale_static**2, original_core_squared_scale=scale_core**2,
        elapsed_seconds=time.perf_counter()-began, material_accepted=False,
        physical_yield_validated=False, physical_Hz=False, whole_family_failure_proved=False))


if __name__ == '__main__':
    main()
