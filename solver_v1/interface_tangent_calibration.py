"""Exact initial-interface tangent controls on existing analytic families.

No new energy, force rescaling, target units or production parameters.
See TANGENT_CONSTRAINED_CALIBRATION_V20.md before interpreting a fit.
"""
from dataclasses import asdict, replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from .coordination_screening import CoordinationScreenedCache, CoordinationScreenedBulk
from .vector_material_calibration import MaterialObservation, IDEAL_H, UNITS
from .run_vector_material_calibration import source_and_targets
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix
from .tail_constrained_material import spectral_profile


NONNEGATIVE = (0, 1, 2, 4, 5, 8, 9)
TANGENT_NAMES = ('perfect_Haa', 'perfect_Hxx')


def declared_tangent_validation(source):
    """New state identities declared before v20 coefficient calculations."""
    states = [(r*IDEAL_H, 0., 0.) for r in (1.003, 1.024, 1.077, 1.233, 1.537, 1.821, 2.057)]
    states += [(1.013*IDEAL_H, .083, -.041), (1.047*IDEAL_H, .213, .079),
               (1.103*IDEAL_H, .383, .183), (.993*IDEAL_H, .123, -.021),
               (1.253*IDEAL_H, .443, .163)]
    observations = []
    for number, q in enumerate(states):
        value = source.evaluate(q)
        jet = np.r_[value.energy, value.gradient, value.hessian[0], value.hessian[1, 1:], value.hessian[2, 2]]
        for index, name, unit in ((0, 'energy', 'eV/cell'), (1, 'force', 'eV/L0'),
                                  (4, 'Haa', 'eV/L0^2'), (7, 'Hxx', 'eV/L0^2')):
            target = float(jet[index])
            scale = (float(UNITS.traction_mpa_to_force(250.)) if index == 1 else
                     max(.1*abs(target), .01/float(UNITS.energy_to_surface(1.))) if index == 0 else
                     max(.1*abs(target), .1))
            observations.append(MaterialObservation(f'v20_new_{number}_{name}', target,
                scale, unit, 'heldout', tuple(q), tuple(np.eye(10)[index])))
    return observations


def impose_tangents(observations, imposed):
    """A target CONTROL is explicit; original source observations stay intact."""
    if set(imposed)-set(TANGENT_NAMES):
        raise ValueError('only declared pristine tangent controls are permitted')
    if any(not np.isfinite(v) or v <= 0 for v in imposed.values()):
        raise ValueError('positive finite stable tangent targets required')
    if any(sum(o.name == name for o in observations) != 1 for name in imposed):
        raise ValueError('each tangent observable must occur exactly once')
    return [replace(o, target=float(imposed[o.name]), role='exact') if o.name in imposed else o
            for o in observations]


def equality_sign_feasibility(matrix, observations, nonnegative=NONNEGATIVE):
    """A separate LP for exact rows/signs ONLY, not spectral stability proof."""
    M = np.asarray(matrix, float)
    rows = [i for i, o in enumerate(observations) if o.role == 'exact']
    scales = np.array([observations[i].scale for i in rows])
    E = M[rows]/scales[:, None]
    rhs = np.array([observations[i].target for i in rows])/scales
    norm = np.linalg.norm(M/np.array([o.scale for o in observations])[:, None], axis=0)
    if np.any(norm == 0):
        raise ValueError('zero coefficient column; fix gauge before checking feasibility')
    D = 1/norm
    result = linprog(np.zeros(M.shape[1]), A_eq=E*D, b_eq=rhs,
        bounds=[(0., None) if i in nonnegative else (None, None) for i in range(M.shape[1])],
        method='highs')
    return dict(success=bool(result.success), status=int(result.status), message=str(result.message),
        exact_rank=int(np.linalg.matrix_rank(E*D)), exact_row_count=len(rows),
        null_dimension=M.shape[1]-int(np.linalg.matrix_rank(E*D)),
        maximum_exact_residual=None if not result.success else float(np.max(abs(E@(D*result.x)-rhs))),
        spectral_constraints_included=False, fixed_shape_only=True)


def make_residual_table(reference, imposed_observations, fit):
    """Keep the fit/control residual distinct from the ORIGINAL source error."""
    prediction = np.asarray(fit['predictions'])
    imposed_residual = (prediction-np.array([o.target for o in imposed_observations]))/np.array(
        [o.scale for o in imposed_observations])
    np.testing.assert_allclose(imposed_residual, fit['residuals'], atol=1e-9, rtol=1e-10)
    return [dict(observable=o.name, role=o.role, original_role=source.role,
        validation_status='new_predeclared' if o.name.startswith('v20_new_') else
                          'retrospective' if source.role == 'heldout' else 'development_or_control',
        source_target=source.target, imposed_target=o.target, prediction=float(p),
        units=o.units, scale=o.scale, imposed_normalized_residual=float(r),
        source_normalized_residual=float((p-source.target)/source.scale))
        for source, o, p, r in zip(reference, imposed_observations, prediction, imposed_residual)]


class TangentCalibrationProblem:
    def __init__(self, parent_joint):
        self.parent = Path(parent_joint)
        raw_bytes = (self.parent/'definition.json').read_bytes()
        self.parent_definition = json.loads(raw_bytes)
        self.source, _, _ = source_and_targets()
        if self.parent_definition['source_sha256'] != self.source.reference.sha256:
            raise ValueError('unchanged Al99 source binding required')
        original = [MaterialObservation(**{**o,
            'state': tuple(o['state']) if o['state'] is not None else None,
            'jet_weights': tuple(o['jet_weights']) if o['jet_weights'] is not None else None})
            for o in self.parent_definition['observations']]
        new = declared_tangent_validation(self.source)
        previous_states = {o.state for o in original if o.state is not None}
        # The v19 independent validator's states were also inspected previously.
        from .validate_coordination_screening import independent_states
        previous_states.update(independent_states())
        if previous_states & {o.state for o in new}:
            raise ValueError('new validation overlaps a previously inspected state')
        self.raw_observations = original+new
        _, converted = cubic_metric_problem(np.zeros((len(self.raw_observations), 8)), self.raw_observations)
        self.observations = [replace(o, role='exact') if i < 5 else o for i, o in enumerate(converted)]
        self.tangent_rows = {name: next(i for i, o in enumerate(self.observations) if o.name == name)
                             for name in TANGENT_NAMES}
        self.development = [i for i, o in enumerate(self.observations) if o.role == 'fit']
        self.other_rows = [i for i in self.development if i not in self.tangent_rows.values()]
        if len(self.development) != 104 or len(self.other_rows) != 102:
            raise ValueError('v20 requires the declared unchanged v19 development dataset')
        self.cache = CoordinationScreenedCache(self.raw_observations, law='power')
        self.parents = {}
        for family in ('old_family', 'power_family'):
            parent_bytes = (self.parent/family/'calibration.json').read_bytes()
            data = json.loads(parent_bytes)
            if not data['completed'] or not data['best']['strictly_positive_LJ']:
                raise ValueError('completed positive-LJ parent required')
            self.parents[family] = dict(shape=data['best']['shape'], coefficients=data['best']['coefficients'],
                source_loss=data['best']['squared_loss'], sha256=hashlib.sha256(parent_bytes).hexdigest(),
                optimizer=data['optimizer'])
        self.definition = self.parent_definition | dict(
            experiment='v20 exact initial interface tangent calibration',
            parent_joint_definition_sha256=hashlib.sha256(raw_bytes).hexdigest(),
            parents=self.parents, observations=[asdict(o) for o in self.raw_observations],
            source_tangents=self.source_tangents, screening_law='power',
            original_development_count=104, common_other_count=102,
            new_validation_count=len(new), new_validation_used_for_selection=False,
            energy_family_changed=False, production_changed=False)

    @property
    def source_tangents(self):
        return {name: self.observations[i].target for name, i in self.tangent_rows.items()}

    @lru_cache(maxsize=128)
    def matrix(self, shape):
        raw = self.cache.matrix(np.asarray(shape))
        first, _ = cubic_metric_problem(raw[:, :8], self.raw_observations)
        extra = raw[:, 8:].copy(); extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
        return np.column_stack([first, extra])

    @lru_cache(maxsize=64)
    def operators(self, shape):
        operators, errors = [], []
        for stretch in self.parent_definition['stability_stretches']:
            bulk = CoordinationScreenedBulk(shape, radius=self.parent_definition['radius_over_L0'],
                                            stretch=stretch, law='power')
            columns, tails = map(np.asarray, zip(*(bulk.evaluate(q) for q in self.parent_definition['wavepoints_cubic'])))
            operators.extend(columns); errors.extend(tails)
        return np.array(operators), np.array(errors)

    def fit(self, shape, imposed):
        shape = tuple(map(float, shape))
        M = self.matrix(shape)
        observations = impose_tangents(self.observations, imposed)
        feasibility = equality_sign_feasibility(M, observations)
        if not feasibility['success']:
            return dict(completed=False, classification='equality_sign_LP_failed',
                        feasibility=feasibility, shape=shape, imposed_tangents=imposed)
        operators, tails = self.operators(shape)
        try:
            fit = spectral_profile(M, observations, operators, tails, nonnegative=NONNEGATIVE)
        except (ArithmeticError, ValueError) as error:
            return dict(completed=False, classification='spectral_profile_not_verified',
                error=str(error), feasibility=feasibility, shape=shape, imposed_tangents=imposed)
        original_residual = (fit['predictions']-np.array([o.target for o in self.observations]))/np.array(
            [o.scale for o in self.observations])
        fit.update(completed=True, shape=shape, imposed_tangents=imposed, feasibility=feasibility,
            source_loss_104=float(original_residual[self.development]@original_residual[self.development]),
            other_loss_102=float(original_residual[self.other_rows]@original_residual[self.other_rows]),
            material_accepted=False, energy_family_changed=False)
        fit['residual_table'] = make_residual_table(self.observations, observations, fit)
        return fit


def tangent_summary(case, family, fit, problem):
    row = dict(case=case, family=family, completed=fit['completed'],
        classification=fit.get('classification', 'verified_coefficient_profile'),
        source_loss_104=None, other_loss_102=None, selected_objective=None,
        selected_count=None, Haa=None, Hxx=None, positive_LJ=None,
        exact_residual=None, kkt_residual=None, margin=None, spectral_cuts=None,
        material_accepted=False)
    if not fit['completed']:
        return row
    row.update(source_loss_104=fit['source_loss_104'], other_loss_102=fit['other_loss_102'],
        selected_objective=fit['squared_loss'], selected_count=len(fit['selected_rows']),
        Haa=fit['predictions'][problem.tangent_rows['perfect_Haa']],
        Hxx=fit['predictions'][problem.tangent_rows['perfect_Hxx']],
        positive_LJ=fit['strictly_positive_LJ'], exact_residual=fit['exact_residual'],
        kkt_residual=fit['kkt_residual'], margin=fit['minimum_robust_margin'],
        spectral_cuts=fit['cut_iterations'], material_accepted=False)
    return row


def shape_selection_status(data):
    """A positive trial profile is not necessarily the optimizer endpoint.

    The closed-domain search can terminate with v=0. Retaining an earlier
    positive-LJ trial does not make that trial shape-optimized/converged.
    """
    optimizer = data.get('optimizer', {})
    best = data.get('best')
    coordinates = optimizer.get('last_accepted_coordinates')
    is_endpoint = False
    if best is not None and coordinates is not None:
        x = np.asarray(coordinates, float)
        endpoint = np.r_[np.exp(x[:5]), x[5] if len(x) == 6 else 0.]
        is_endpoint = bool(np.allclose(best['shape'], endpoint, rtol=1e-12, atol=1e-12))
    return dict(positive_profile_found=best is not None,
        selected_profile_is_optimizer_endpoint=is_endpoint,
        positive_material_shape_optimizer_converged=bool(
            best is not None and best['strictly_positive_LJ'] and is_endpoint and optimizer.get('success', False)),
        selection_includes_trial_and_difference_profiles=True,
        material_accepted=False, global_optimality_claimed=False)
