"""Same-energy joint force/interface fit with explicitly audited source-jet boxes.

The two offending inspected curvatures become development constraints. No
target values/scales are changed, no empirical yield is introduced, and every
other residual is saved. A feasible fit is not an adopted material.
"""
import argparse
from dataclasses import replace
import hashlib
from pathlib import Path
import time

import numpy as np

from .core_interface_compatibility import interval_profile
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
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--model', type=Path, default=ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh source-jet-box output required')
    began = time.perf_counter()
    model, _, binding = load_current_material(args.model)
    baseline, _, original_binding = load_current_material()
    source, tensor, source_binding = load_source_material()
    source_path = ROOT/'results/current_material_core_v22/source_reference/escaped_plus_R8r5'
    core, field, *_ = restore_case(source_path, source, tensor, source_binding, core_builder=build_source_core)
    frozen = FrozenCoreForceTarget(core, field, radius=2.)
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    original = impose_tangents(problem.observations, problem.source_tangents)
    lookup = {o.name: i for i, o in enumerate(original)}
    critical = [lookup['saddle_Hxx'], lookup['v19_new_4_Haa']]
    selected = [i for i, o in enumerate(original) if o.role == 'fit']
    base_matrix = problem.matrix(tuple(baseline.screened_shape))
    base_force = frozen.coefficient_matrix(baseline)
    old_residual = (base_matrix@baseline.coefficients-np.array([o.target for o in original]))/np.array([o.scale for o in original])
    scale_static = float(np.linalg.norm(old_residual[selected]))
    scale_core = float(np.linalg.norm(base_force['design']@baseline.coefficients-base_force['target']))
    observations = [replace(o, scale=o.scale*scale_static) if o.role == 'fit' else o for o in original]
    F = frozen.coefficient_matrix(model)
    observations += [MaterialObservation(f'core_gradient_{i}', float(g), scale_core, 'eV/L0', 'fit')
                     for i, g in enumerate(F['target'])]
    M = problem.matrix(tuple(model.screened_shape))
    matrix = np.vstack([M, F['design']])
    operators, tails = problem.operators(tuple(model.screened_shape))
    protocols = [('baseline_replay', []), ('saddle_only', critical[:1]),
                 ('opening_only', critical[1:]), ('critical_pair', critical)]
    save_json(args.out/'definition.json', dict(candidate_sha256=binding['parameter_sha256'],
        original_sha256=original_binding['parameter_sha256'], source_sha256=source_binding['parameter_sha256'],
        shape=model.screened_shape, coefficient_signs=NONNEGATIVE,
        exact_anchors=[o.name for o in original if o.role == 'exact'],
        source_state=source_path.relative_to(ROOT).as_posix(),
        original_static_squared_loss=scale_static**2, original_core_squared_loss=scale_core**2,
        objective='same v22 sum of fractional static/core squared losses',
        protocols=[dict(name=name, constrained=[original[i].name for i in ids]) for name, ids in protocols],
        original_discrepancy_scales_kept=True, inspected_heldout_becomes_development_constraint=True,
        topology_constraint_is_source_target_not_arbitrary_sign_flip=True,
        no_new_energy_term=True, no_nonlinear_optimization=True, material_accepted=False,
        production_changed=False, physical_yield_used=False))
    results, summary, residuals = [], [], []
    for name, rows in protocols:
        tick = time.perf_counter()
        try:
            if rows:
                fit = interval_profile(matrix, observations, rows, 1., nonnegative=NONNEGATIVE,
                    operators=operators, tails=tails, interval_scales=[original[i].scale for i in rows])
            else:
                fit = spectral_profile(matrix, observations, operators, tails, nonnegative=NONNEGATIVE)
            c = fit['coefficients']
            raw = (M@c-np.array([o.target for o in original]))/np.array([o.scale for o in original])
            if name == 'baseline_replay' and np.max(abs(M@c-M@model.coefficients)) > 1e-7:
                raise ArithmeticError('unchanged v22 objective failed mandatory replay')
            row = dict(protocol=name, completed=True, joint_loss=fit['squared_loss'],
                static_loss=float(raw[selected]@raw[selected]),
                force_RMS=float(np.sqrt(np.mean((F['design']@c-F['target'])**2))),
                saddle_Hxx=float(M[critical[0]]@c), opening_Haa=float(M[critical[1]]@c),
                exact_residual=fit['exact_residual'], kkt=fit['kkt_residual'],
                positive_LJ=fit['strictly_positive_LJ'], robust_margin=fit['minimum_robust_margin'],
                error=None, elapsed_seconds=time.perf_counter()-tick, material_accepted=False)
            for i, (o, pred, ri) in enumerate(zip(original, M@c, raw)):
                residuals.append(dict(protocol=name, observable=o.name, previous_role=o.role,
                    target=o.target, scale=o.scale, units=o.units, prediction=pred,
                    normalized_residual=ri, used_in_objective=o.role == 'fit',
                    used_in_constraint=i in rows, independent_validation=False))
            fit.update(shape=model.screened_shape, material_accepted=False)
            if fit['strictly_positive_LJ']:
                save_json(args.out/name/'research_candidate_snapshot.json', dict(completed=True,
                    best=fit, screening_law='power', for_static_research_only=True,
                    source_parameter_sha256=binding['parameter_sha256'],
                    source_state_sha256=hashlib.sha256((source_path/'state.csv').read_bytes()).hexdigest(),
                    material_accepted=False, shape_optimizer_run=False, production_changed=False))
            results.append(dict(protocol=name, completed=True, fit=fit))
        except (ValueError, ArithmeticError) as error:
            row = dict(protocol=name, completed=False, joint_loss=None, static_loss=None,
                force_RMS=None, saddle_Hxx=None, opening_Haa=None, exact_residual=None,
                kkt=None, positive_LJ=None, robust_margin=None, error=str(error),
                elapsed_seconds=time.perf_counter()-tick, material_accepted=False)
            results.append(dict(protocol=name, completed=False, error=str(error)))
        summary.append(row)
        save_json(args.out/'checkpoint.json', dict(completed=False, profiles=results))
        print(f'{name}: {row}', flush=True)
    write_csv(args.out/'profile_summary.csv', summary)
    write_csv(args.out/'all_observable_residuals.csv', residuals)
    save_json(args.out/'completion.json', dict(completed=True, profiles=results,
        actual_profile_attempts=len(summary), elapsed_seconds=time.perf_counter()-began,
        material_accepted=False, production_changed=False, physical_Hz=False))


if __name__ == '__main__':
    main()
