"""Certified fixed-shape compatibility, NOT a new material energy.

Observables are exact coefficient-linear jets of the existing infinite
LJ/Bessel model. LP certificates distinguish a least-squares compromise from
an unattainable target box. All statements are for the supplied fixed shape
and declared coefficient signs; omitted spectral constraints can only shrink
the feasible set. No result proves impossibility of the full radial family.
"""
from functools import lru_cache

import numpy as np
from scipy.optimize import linprog


def minimax_compatibility(matrix, observations, rows, *, nonnegative=(),
                         exact_rows=None):
    """Minimize max_i |(M_i c-y_i)/scale_i| with independently checked LP dual.

    Returns a numerical primal/dual certificate including its residual, not
    merely HiGHS' success flag. Scales are declared discrepancy units, not
    estimated measurement uncertainties. Physical coefficients are not clipped.
    """
    M = np.asarray(matrix, float)
    rows = np.asarray(rows, int)
    if (M.ndim != 2 or M.shape[0] != len(observations) or not M.shape[1]
            or rows.ndim != 1 or not len(rows) or len(set(rows)) != len(rows)
            or np.any(rows < 0) or np.any(rows >= len(M)) or np.any(~np.isfinite(M))):
        raise ValueError('finite matching design and nonempty unique observable indices required')
    p = M.shape[1]
    signs = tuple(nonnegative)
    if len(set(signs)) != len(signs) or any(i < 0 or i >= p for i in signs):
        raise ValueError('unique valid nonnegative coefficient indices required')
    exact = np.asarray([i for i, o in enumerate(observations) if o.role == 'exact']
                       if exact_rows is None else exact_rows, int)
    if (exact.ndim != 1 or not len(exact) or len(set(exact)) != len(exact)
            or np.any(exact < 0) or np.any(exact >= len(M))):
        raise ValueError('explicit nonempty valid equality rows required')
    targets = np.array([o.target for o in observations])
    scales = np.array([o.scale for o in observations])
    if np.any(~np.isfinite(targets)) or np.any(~np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError('finite targets and positive discrepancy scales required')
    norm = np.linalg.norm(M[np.r_[exact, rows]] / scales[np.r_[exact, rows], None], axis=0)
    # A zero column is irrelevant to this LP; keep its original sign and use
    # unit coordinate scale, not a fictitious observable or material prior.
    D = np.ones(p)
    np.divide(1., norm, out=D, where=norm > 0)
    E = M[exact] * D / scales[exact, None]
    e = targets[exact] / scales[exact]
    J = M[rows] * D / scales[rows, None]
    y = targets[rows] / scales[rows]
    Aeq = np.column_stack([E, np.zeros(len(E))])
    A = np.vstack([np.column_stack([J, -np.ones(len(J))]),
                   np.column_stack([-J, -np.ones(len(J))])])
    b = np.r_[y, -y]
    objective = np.r_[np.zeros(p), 1.]
    bounds = [(0., None) if i in signs else (None, None) for i in range(p)] + [(0., None)]
    result = linprog(objective, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=e,
                     bounds=bounds, method='highs',
                     options={'dual_feasibility_tolerance': 1e-9,
                              'primal_feasibility_tolerance': 1e-9})
    common = dict(exact_rows=exact, tested_rows=rows, nonnegative=signs,
                  coefficient_scales=D, fixed_shape_only=True,
                  whole_family_impossibility_proved=False,
                  spectral_constraints_included=False, material_accepted=False)
    if not result.success:
        return dict(**common, completed=False, lp_status=int(result.status),
                    message=str(result.message), certified_lower_bound_available=False)
    x = result.x
    dual_eq = result.eqlin.marginals
    dual_ub = result.ineqlin.marginals
    dual_lb = result.lower.marginals
    equality_error = float(np.max(abs(Aeq @ x - e)))
    violation = max(0., float(np.max(A @ x - b)),
                    max((-x[i] for i in (*signs, p)), default=0.))
    stationarity = float(np.max(abs(objective - Aeq.T @ dual_eq - A.T @ dual_ub - dual_lb)))
    dual_objective = float(e @ dual_eq + b @ dual_ub)  # all finite lower bounds are zero
    duality_gap = float(objective @ x - dual_objective)
    complementarity = max(float(np.max(abs(dual_ub * (A @ x - b)))),
                          float(np.max(abs(dual_lb * x))))
    dual_sign_error = max(0., float(np.max(dual_ub)), float(-np.min(dual_lb)))
    free = [i for i in range(p) if i not in signs]
    free_bound_dual = float(np.max(abs(dual_lb[free]), initial=0.))
    c = D * x[:p]
    actual_error = float(np.max(abs((M[rows] @ c - targets[rows]) / scales[rows])))
    tolerance = 2e-7 * max(1., abs(float(result.fun)))
    if max(equality_error, violation, stationarity, abs(duality_gap),
           complementarity, dual_sign_error, free_bound_dual,
           abs(actual_error - result.fun)) > tolerance:
        raise ArithmeticError('returned minimax LP failed its independent primal/dual certificate')
    return dict(**common, completed=True, lp_status=int(result.status), message=str(result.message),
                coefficients=c, predictions=M @ c, minimax_normalized_error=actual_error,
                dual_lower_bound=dual_objective, duality_gap=duality_gap,
                equality_residual=equality_error, primal_violation=violation,
                stationarity_residual=stationarity, complementarity_residual=complementarity,
                dual_sign_residual=dual_sign_error, certificate_tolerance=tolerance,
                equality_duals=dual_eq, upper_error_duals=dual_ub[:len(rows)],
                lower_error_duals=dual_ub[len(rows):], lower_bound_duals=dual_lb,
                exact_rank=int(np.linalg.matrix_rank(E)),
                target_box_feasible=bool(actual_error <= 1. + tolerance),
                strictly_positive_LJ=bool(p >= 2 and np.all(c[:2] > 0)),
                certified_lower_bound_available=True)


def interval_profile(matrix, observations, rows, widths, *, nonnegative,
                     operators=None, tails=None, interval_scales=None):
    """Existing verified coefficient QP with explicit source-target intervals.

    Homogenize affine bounds by appending c_const=1 as an exact bookkeeping
    coordinate. It is NOT a material parameter or energy term. Existing
    homogeneous spectral constraints get a zero column for this coordinate.
    """
    from .tail_constrained_material import convex_profile, spectral_profile
    from .vector_material_calibration import MaterialObservation

    M = np.asarray(matrix, float)
    rows = np.asarray(rows, int)
    widths = np.broadcast_to(np.asarray(widths, float), rows.shape)
    if (M.ndim != 2 or M.shape[0] != len(observations) or not len(rows)
            or rows.ndim != 1 or len(set(rows)) != len(rows)
            or np.any(rows < 0) or np.any(rows >= len(M))
            or np.any(~np.isfinite(M)) or np.any(~np.isfinite(widths)) or np.any(widths < 0)):
        raise ValueError('finite design, unique target rows and nonnegative interval widths required')
    p = M.shape[1]
    augmented = np.zeros((len(M)+1, p+1))
    augmented[:-1, :-1] = M
    augmented[-1, -1] = 1.
    extra = MaterialObservation('affine_bound_constant_NOT_material', 1., 1., 'bookkeeping', 'exact')
    target = np.array([observations[i].target for i in rows])
    scale = (np.array([observations[i].scale for i in rows]) if interval_scales is None
             else np.asarray(interval_scales, float))
    if scale.shape != rows.shape or np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError('positive matching target-interval scales required')
    lower, upper = target-widths*scale, target+widths*scale
    inequalities = np.vstack([np.column_stack([M[rows], -lower]),
                               np.column_stack([-M[rows], upper])])
    kwargs = dict(nonnegative=tuple(nonnegative)+(p,), inequalities=inequalities)
    if (operators is None) != (tails is None):
        raise ValueError('supply both sampled spectral operators and tails, or neither')
    if operators is None:
        fit = convex_profile(augmented, list(observations)+[extra], **kwargs)
    else:
        op, tail = np.asarray(operators, float), np.asarray(tails, float)
        if op.ndim != 4 or op.shape[1:] != (p, 3, 3) or tail.shape != op.shape[:2]:
            raise ValueError('matching coefficient-linear 3x3 bulk operators required')
        augmented_op = np.pad(op, ((0, 0), (0, 1), (0, 0), (0, 0)))
        augmented_tail = np.pad(tail, ((0, 0), (0, 1)))
        fit = spectral_profile(augmented, list(observations)+[extra],
                               augmented_op, augmented_tail, **kwargs)
    constant = float(fit['coefficients'][-1])
    c = np.asarray(fit['coefficients'][:-1])
    predictions = M @ c
    violation = max(0., float(np.max((lower-predictions[rows])/scale)),
                    float(np.max((predictions[rows]-upper)/scale)))
    if abs(constant-1.) > 1e-9 or violation > 1e-6:
        raise ArithmeticError('affine target-box reconstruction failed')
    fit.update(coefficients=c, predictions=predictions,
               residuals=np.asarray(fit['residuals'][:-1]),
               target_interval_rows=rows, target_interval_widths=widths,
               target_interval_scales=scale,
               target_interval_violation=violation, bookkeeping_constant=constant,
               material_accepted=False, energy_term_added=False)
    return fit


class ExistingJetSubset:
    """Evaluate selected existing jets, with the unchanged cubic-unit map.

    This avoids computing unselected interface states at every radial trial;
    it does not approximate a jet or fit a surrogate. Always retain the first
    five independent raw bulk rows and test against the complete matrix.
    """
    def __init__(self, problem, extra_rows):
        from .coordination_screening import CoordinationScreenedCache
        self.rows = list(range(5)) + sorted(set(extra_rows)-set(range(5)))
        if any(i < 0 or i >= len(problem.observations) for i in self.rows):
            raise ValueError('valid existing-observation indices required')
        self.raw = [problem.raw_observations[i] for i in self.rows]
        self.cache = CoordinationScreenedCache(self.raw, law='power')

    @lru_cache(maxsize=160)
    def matrix(self, shape):
        from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix
        raw = self.cache.matrix(np.asarray(shape, float))
        first, _ = cubic_metric_problem(raw[:, :8], self.raw)
        extra = raw[:, 8:].copy()
        extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
        return np.column_stack([first, extra])
