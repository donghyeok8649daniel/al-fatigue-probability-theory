"""Run deterministic material compatibility profiles and independent validation.

One research candidate family, no fatigue/yield/kinetic target, no promotion.
Run stages independently: fit writes actual optimization logs; validate solves
actual new states (not a replay of the old vector registry results).
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize

from .reference_eam_targets import MishinRigidFCCReference
from .run_low_stress_cyclic_diagnostic import FIT, write_csv
from .run_vector_registry_audit import save_json, point_row, stress_probes, relaxed_registry_curve
from .vector_interface_reference import MishinVectorInterfaceReference
from .vector_material_calibration import (
    COEFFICIENTS, IDEAL_H, LENGTH_M, UNITS, VectorCoefficientBasis, matched_observations,
    observation_matrix, fit_coefficients, coefficient_identifiability,
    build_coefficient_surface, exact_constraint_tangent,
)
from .vector_registry_audit import stationary_state, saddle_downhill_endpoints
from .full_fcc_calibration_audit import cubic_constants_gpa
from .static_bulk_stability import StaticBulkHessian


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/fcc111_active_interface/material_strength_v10/joint_fit'
OPENING_RATIOS = (1.02, 1.05, 1.1, 1.2, 1.5, 2., 2.5, 3., 4., 5., 8.)


def source_and_targets():
    reference = MishinRigidFCCReference()
    source = MishinVectorInterfaceReference(reference, LENGTH_M/1e-10)
    observations, states = matched_observations(source)
    return source, observations, states


def execute_fit(out, *, local_evaluations=45):
    started = time.perf_counter(); source, observations, states = source_and_targets()
    raw = FIT.read_bytes(); old = json.loads(raw)['angular_monotone_opening']
    save_json(out/'target_definition.json', dict(observations=[asdict(o) for o in observations],
        source_states=states, source_sha256=source.reference.sha256,
        old_parameter_sha256=hashlib.sha256(raw).hexdigest(),
        temperature_K=0, interface_convention='identical rigid half-crystals; only active gap/registry vary',
        bulk_targets='published rounded zero-K Mishin elastic constants; same nominal 4.05 angstrom lattice',
        normalization='5% bulk curvature,10% interface energy/curvature;.25 GPa force discrepancy, NOT uncertainties',
        reverse_barrier_not_independent=True, experimental_strength_in_loss=False))
    profiles = []; logs = []; matrices = {}
    def evaluate(decays, label):
        key = tuple(map(float, decays))
        if key not in matrices:
            basis = VectorCoefficientBasis(*key)
            mat = observation_matrix(basis, observations)
            inequalities = np.array([basis.jet((r*IDEAL_H, 0., 0.))[1] for r in OPENING_RATIOS])
            matrices[key] = mat, inequalities
        mat, inequalities = matrices[key]
        result = fit_coefficients(mat, observations, opening_inequalities=inequalities)
        row = dict(scalar_decay=key[0], angular_decay=key[1], stage=label, **result)
        profiles.append(row)
        save_json(out/'profile_progress.json', dict(completed=False, profiles=profiles,
            elapsed_seconds=time.perf_counter()-started))
        print(f"profile {len(profiles)} {label}: decays={key}, loss={row['squared_loss']:.6g}, "
              f"feasible={row['admissible']}, {time.perf_counter()-started:.1f}s", flush=True)
        return row
    old_key = old['scalar_decay'], old['angular_decay']
    evaluate(old_key, 'fixed_previous_decays')
    # Full ranges already used by prior studies. All starts declared, no random search.
    scalar = sorted(set([.7, old_key[0], 2.4, 4.5, 8.]))
    angular = sorted(set([2., 3.2, old_key[1], 7.8, 12.]))
    for d in scalar:
        for k in angular:
            evaluate((d, k), 'predeclared_grid')
    eligible = [r for r in profiles if r['admissible'] and r['strictly_positive_LJ_resolved']]
    if not eligible:
        raise RuntimeError('no admissible coefficient profile; no physical result manufactured')
    start = min(eligible, key=lambda r:r['squared_loss'])
    if local_evaluations:
        def objective(x):
            row = evaluate(x, 'local_Powell')
            # Solver success/admissibility is recorded and checked on selection.
            return row['squared_loss']
        result = minimize(objective, [start['scalar_decay'], start['angular_decay']], method='Powell',
            bounds=((.7, 8.), (2., 12.)),
            options=dict(maxfev=local_evaluations, xtol=2e-4, ftol=1e-6))
        logs.append(dict(start=[start['scalar_decay'], start['angular_decay']],
            success=bool(result.success), message=str(result.message), nfev=result.nfev,
            final_objective=float(result.fun)))
    eligible = [r for r in profiles if r['admissible'] and r['strictly_positive_LJ_resolved']]
    best = min(eligible, key=lambda r:r['squared_loss'])
    key = best['scalar_decay'], best['angular_decay']; mat = matrices[key][0]
    old_matrix = matrices[old_key][0]
    target = np.array([o.target for o in observations]); scale = np.array([o.scale for o in observations])
    stages = []
    for label, indices in [('geometry_cohesion_elasticity', set(range(5))),
                           ('add_vector_registry', set(range(12))),
                           ('joint_with_opening', set(range(16)))]:
        # Same optimum radial conventions for transparent coefficient-only stage ablation.
        selected = [i for i in sorted(indices) if observations[i].role != 'heldout']
        from .run_constrained_odd_calibration import constrained_matrix_fit
        keep = np.flatnonzero(np.any(mat[selected] != 0., axis=0))
        lower = np.array([0., 0., 0., -np.inf, 0., 0., -np.inf, -np.inf])
        stage = constrained_matrix_fit(mat[np.ix_(selected,keep)], target[selected], scale[selected],
            coefficient_order=[COEFFICIENTS[i] for i in keep], lower_rest=lower[keep[2:]],
            nonnegative_observation_rows=matrices[key][1][:,keep])
        stages.append(dict(stage=label, selected_rows=selected,
            fixed_zero_unobservable_columns=[COEFFICIENTS[i] for i in range(8) if i not in keep], **stage))
    rows = []
    for name, coeff, matrix in [('old_candidate', np.array(old['coefficients']), old_matrix),
                               ('best_feasible_same_family', best['coefficients'], mat)]:
        prediction = matrix@coeff
        for obs, pred in zip(observations, prediction):
            rows.append(dict(model=name, target=obs.name, reference=obs.target, prediction=pred,
                units=obs.units, absolute_error=pred-obs.target,
                relative_error=(pred-obs.target)/abs(obs.target) if abs(obs.target)>1e-9 else None,
                normalization=obs.scale, normalized_residual=(pred-obs.target)/obs.scale, role=obs.role))
    write_csv(out/'fit_and_heldout_residuals.csv', rows)
    identity = coefficient_identifiability(mat, observations, best['coefficients'])
    # Sensitivity along exact constraints, including both radial decays.
    selected = [i for i, o in enumerate(observations) if o.role == 'fit']
    z = np.asarray(best['coefficients'])[2:]; full_jac = list(identity['jacobian'].T)
    derivative_refinements = []
    for log_step in (4e-4, 2e-4):
        columns = []
        for axis in range(2):
            predictions = []
            for sign in (1, -1):
                shifted = np.array(key); shifted[axis] *= np.exp(sign*log_step)
                shifted_matrix = observation_matrix(VectorCoefficientBasis(*shifted), observations)
                offset, tangent = exact_constraint_tangent(shifted_matrix, target)
                predictions.append(shifted_matrix@(offset+tangent@z))
            columns.append((predictions[0]-predictions[1])[selected]/(2*log_step*scale[selected]))
        derivative_refinements.append(np.array(columns))
    full_jac.extend(derivative_refinements[-1])
    full_jac = np.array(full_jac).T
    _, singular, vt = np.linalg.svd(full_jac, full_matrices=False)
    norms = np.linalg.norm(full_jac, axis=0)
    cosine = full_jac.T@full_jac/(norms[:,None]*norms[None,:])
    identity.update(full_jacobian=full_jac, full_singular_values=singular,
        full_condition_number=float(singular[0]/singular[-1]), full_right_vectors=vt,
        column_cosines=cosine, full_radial_identifiability_checked=True,
        radial_difference_change=float(np.max(abs(derivative_refinements[0]-derivative_refinements[1]))),
        interpretation='scaled local sensitivities, not statistical confidence or a global uniqueness proof')
    save_json(out/'identifiability.json', identity)
    save_json(out/'calibration.json', dict(completed=True, best=best, stages=stages, profiles=profiles,
        local_optimizations=logs, coefficient_order=COEFFICIENTS,
        scalar_grid=scalar, angular_grid=angular, opening_inequalities=OPENING_RATIOS,
        sampled_inequality_not_continuous_certificate=True,
        elapsed_seconds=time.perf_counter()-started, new_energy_terms_added=False,
        production_promoted=False, material_accepted=False))
    save_json(out/'profile_progress.json', dict(completed=True, profiles=profiles,
        elapsed_seconds=time.perf_counter()-started))


def execute_validation(out):
    started = time.perf_counter(); data = json.loads((out/'calibration.json').read_text())
    old = json.loads(FIT.read_bytes())['angular_monotone_opening']; best = data['best']
    source, observations, source_states = source_and_targets()
    static = []; curves = []; extrema = []; bulk_rows = []; barriers = []; folds = []; stress = []
    for name, fit in [('old_candidate', old), ('best_feasible_same_family', best)]:
        print('validate '+name, flush=True)
        model = build_coefficient_surface(fit['scalar_decay'], fit['angular_decay'], fit['coefficients'])
        mat = observation_matrix(VectorCoefficientBasis(fit['scalar_decay'], fit['angular_decay']), observations)
        bulk_rows.append(dict(model=name, **cubic_constants_gpa(mat[:5]@fit['coefficients'])))
        roots = {}
        for label, guess, index in [('perfect', [IDEAL_H,0,0], 0),
                                   ('fault', source_states['fault'], 0),
                                   ('saddle', source_states['saddle'], 1)]:
            r = stationary_state(model, guess, expected_index=index)
            if not r['valid']:
                raise RuntimeError(f'{name}/{label} force/index failure: {r}')
            roots[label] = r
            static.append(point_row(name, label, r['q'], r['evaluation'], UNITS))
        endpoints = saddle_downhill_endpoints(model, roots['saddle'])
        barriers.append(dict(model=name,
            forward_J_m2=float(UNITS.energy_to_surface(roots['saddle']['evaluation'].energy-roots['perfect']['evaluation'].energy)),
            reverse_J_m2=float(UNITS.energy_to_surface(roots['saddle']['evaluation'].energy-roots['fault']['evaluation'].energy)),
            fault_J_m2=float(UNITS.energy_to_surface(roots['fault']['evaluation'].energy-roots['perfect']['evaluation'].energy)),
            endpoint_energy_low=min(r['evaluation'].energy for r in endpoints),
            endpoint_energy_high=max(r['evaluation'].energy for r in endpoints)))
        _, fold = relaxed_registry_curve(model, name, UNITS, 41)
        folds.append(fold)
        stress.extend(stress_probes(model, name, UNITS))
        # Continuous stationary traction candidates from independent W_aa brackets.
        for count in (91, 181):
            grid = np.unique(np.r_[np.linspace(1.00001,2.,count), np.geomspace(2.,12.,count//2)])*IDEAL_H
            values = [model.evaluate([a,0.,0.]) for a in grid]
            for a, v in zip(grid, values):
                if count == 181:
                    curves.append(dict(model=name,a_over_h=a/IDEAL_H,energy_J_m2=float(UNITS.energy_to_surface(v.energy)),
                        normal_traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0])),Haa_eV_L0sq=v.hessian[0,0]))
            for i in range(1,len(grid)):
                if values[i-1].hessian[0,0]*values[i].hessian[0,0] < 0:
                    a = brentq(lambda a:model.evaluate([a,0.,0.]).hessian[0,0],grid[i-1],grid[i],xtol=2e-12)
                    v = model.evaluate([a,0.,0.])
                    extrema.append(dict(model=name,samples=count,a_over_h=a/IDEAL_H,
                        traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0])),
                        Haa_eV_L0sq=v.hessian[0,0]))
        print(f"{name} static completed in {time.perf_counter()-started:.1f}s",flush=True)
    write_csv(out/'bulk_elastic_comparison.csv',bulk_rows)
    write_csv(out/'stationary_states.csv',static)
    write_csv(out/'barrier_comparison.csv',barriers)
    write_csv(out/'ideal_folds_NOT_yield.csv',folds)
    write_csv(out/'static_stress_and_unload.csv',stress)
    write_csv(out/'opening_curve.csv',curves)
    write_csv(out/'opening_traction_extrema.csv',extrema)
    # Independent finite-q validation, never the canonical energy sum.
    model = build_coefficient_surface(best['scalar_decay'],best['angular_decay'],best['coefficients'])
    c = best['coefficients']; finiteq=[]
    for cutoff in (8.,12.):
        audit=StaticBulkHessian(model.face.bulk,cutoff=cutoff,D3=c[5],D1=c[6],D2=c[7],angular_decay=best['angular_decay'])
        for path,end in [('GX',[0.,1.,0.]),('GL',[.5,.5,.5]),('GK',[.75,.75,0.])]:
            for fraction in np.linspace(.02,1.,18):
                q=audit.crystallographic_wavevector(fraction*np.array(end))
                eig=audit.evaluate(q)['eigenvalues']
                finiteq.append(dict(cutoff_L0=cutoff,path=path,fraction=fraction,lambda_min=eig[0],lambda_mid=eig[1],lambda_max=eig[2]))
    write_csv(out/'finite_q_stability.csv',finiteq)
    # Reference-normalization/stationary validation does NOT accept an Al specimen model.
    save_json(out/'validation_status.json', dict(completed=True,elapsed_seconds=time.perf_counter()-started,
        material_accepted=False,finite_source_validated=False,experimental_strength_validated=False,
        physical_seconds=False,physical_Hz=False,production_promoted=False,
        minimum_sampled_finite_q_eigenvalue=min(r['lambda_min'] for r in finiteq),
        most_negative_opening_extremum_MPa=min(r['traction_MPa'] for r in extrema)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['fit','validate'])
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--local-evaluations',type=int,default=45)
    args=parser.parse_args()
    if args.stage=='fit':
        execute_fit(args.output,local_evaluations=args.local_evaluations)
    else:
        execute_validation(args.output)


if __name__=='__main__':
    main()
