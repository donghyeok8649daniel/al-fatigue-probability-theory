"""Eliminate exact coefficient-linear anchors before nonlinear shape fitting.

The energy family, target scales and coefficient signs are unchanged. A fixed
QR-selected pivot provides smooth local coordinates; rank loss is rejected,
never regularized into an apparent calibration.
"""
import numpy as np
from scipy.linalg import qr
from scipy.optimize import linprog


class ExactAnchorCoordinates:
    def __init__(self, matrix, targets, scales, exact, coefficient_scales):
        M = np.asarray(matrix, float)
        self.targets = np.asarray(targets, float)
        self.scales = np.asarray(scales, float)
        self.exact = np.asarray(exact, int)
        self.D = np.asarray(coefficient_scales, float)
        if (M.ndim != 2 or self.targets.shape != (len(M),)
                or self.scales.shape != self.targets.shape or self.D.shape != (M.shape[1],)
                or self.exact.ndim != 1 or not 0 < len(self.exact) < M.shape[1]
                or len(set(self.exact)) != len(self.exact)
                or np.any(self.exact < 0) or np.any(self.exact >= len(M))
                or not np.all(np.isfinite(M)) or not np.all(np.isfinite(self.targets))
                or not np.all(np.isfinite(self.scales)) or not np.all(np.isfinite(self.D))
                or np.any(self.D <= 0) or np.any(self.scales <= 0)):
            raise ValueError('finite design, unique anchors and positive numerical scales required')
        self.shape = M.shape
        E = (M*self.D/self.scales[:, None])[self.exact]
        _, _, pivot = qr(E, pivoting=True, mode='economic')
        self.dependent = np.array(pivot[:len(self.exact)])
        self.free = np.array(sorted(set(range(M.shape[1]))-set(self.dependent)))
        self._check_pivot(E[:, self.dependent])

    @staticmethod
    def _check_pivot(pivot):
        condition = float(np.linalg.cond(pivot))
        # A numerical arithmetic budget only, not physical or sampling uncertainty.
        budget = condition*len(pivot)*np.finfo(float).eps
        if not np.isfinite(condition) or budget > 1e-6:
            raise ArithmeticError('exact-anchor pivot exceeds its arithmetic error budget')
        return condition

    def evaluate(self, matrix, derivatives, free_coefficients, eta, tested, nonnegative):
        M, dM = np.asarray(matrix, float), np.asarray(derivatives, float)
        y = np.asarray(free_coefficients, float)
        tested, signs = np.asarray(tested, int), np.asarray(nonnegative, int)
        if (M.shape != self.shape or dM.ndim != 3 or dM.shape[1:] != self.shape
                or y.shape != self.free.shape or not np.all(np.isfinite(M))
                or not np.all(np.isfinite(dM)) or not np.all(np.isfinite(y))
                or not np.isfinite(eta)):
            raise ValueError('matching finite matrix, derivatives and free coefficients required')
        if (tested.ndim != 1 or not len(tested) or len(set(tested)) != len(tested)
                or np.any(tested < 0) or np.any(tested >= len(M))
                or set(tested) & set(self.exact) or signs.ndim != 1
                or len(set(signs)) != len(signs) or np.any(signs < 0) or np.any(signs >= M.shape[1])):
            raise ValueError('valid distinct test rows and coefficient signs required')
        normalized = M*self.D/self.scales[:, None]
        E = normalized[self.exact]
        pivot = E[:, self.dependent]
        condition = self._check_pivot(pivot)
        z = np.zeros(M.shape[1])
        z[self.free] = y
        z[self.dependent] = np.linalg.solve(pivot,
            self.targets[self.exact]/self.scales[self.exact]-E[:, self.free]@y)
        N = np.zeros((M.shape[1], len(self.free)))
        N[self.free] = np.eye(len(self.free))
        N[self.dependent] = -np.linalg.solve(pivot, E[:, self.free])
        dnormalized = dM*self.D/self.scales[None, :, None]
        dz = np.zeros((len(dM), M.shape[1]))
        dz[:, self.dependent] = -np.linalg.solve(pivot,
            (dnormalized[:, self.exact]@z).T).T
        residual = normalized@z-self.targets/self.scales
        jacobian = np.column_stack([(dnormalized@z).T+normalized@dz.T,
                                     normalized@N, np.zeros(len(M))])
        constraints = np.r_[eta-residual[tested], eta+residual[tested], z[signs]]
        J = np.vstack([-jacobian[tested], jacobian[tested],
                       np.column_stack([dz[:, signs].T, N[signs], np.zeros(len(signs))])])
        J[:2*len(tested), -1] = 1.
        return dict(coefficients=self.D*z, scaled_coefficients=z, residuals=residual,
            constraints=constraints, jacobian=J, residual_jacobian=jacobian,
            exact_residual=float(np.max(abs(residual[self.exact]))), pivot_condition=condition)


def tangent_descent(constraints, jacobian, point, lower, upper, active_tolerance):
    """Search the active first-order feasible cone for decreasing epigraph eta.

    Directions have unit infinity norm in the declared numerical coordinates.
    No resolved direction is only a necessary local test with supplied finite
    derivatives and active tolerance, not a global or material certificate.
    """
    c, J, z, lo, hi = (np.asarray(v, float) for v in (constraints, jacobian, point, lower, upper))
    if (c.ndim != 1 or z.ndim != 1 or J.shape != (len(c), len(z))
            or lo.shape != z.shape or hi.shape != z.shape
            or not np.all(np.isfinite(c)) or not np.all(np.isfinite(J))
            or not np.all(np.isfinite(z)) or np.any(np.isnan(lo)) or np.any(np.isnan(hi))
            or np.any(lo > hi) or not np.isfinite(active_tolerance) or active_tolerance <= 0
            or np.any(c < -active_tolerance) or np.any(z < lo-active_tolerance)
            or np.any(z > hi+active_tolerance)):
        raise ValueError('finite feasible constraints, derivatives, bounds and positive active tolerance required')
    active = np.flatnonzero(c <= active_tolerance)
    bounds = [(0. if x-l <= active_tolerance else -1.,
               0. if h-x <= active_tolerance else 1.) for x, l, h in zip(z, lo, hi)]
    objective = np.r_[np.zeros(len(z)-1), 1.]
    result = linprog(objective, A_ub=-J[active], b_ub=np.zeros(len(active)),
                     bounds=bounds, method='highs')
    if not result.success:
        raise ArithmeticError('linearized feasible-cone solve failed')
    return dict(direction=result.x, minimum_linearized_deta=float(result.fun),
        active_constraints=active, active_tolerance=active_tolerance,
        first_order_violation=float(max(0., -np.min(J[active]@result.x, initial=0.))),
        global_shape_optimum_certified=False, material_accepted=False)
