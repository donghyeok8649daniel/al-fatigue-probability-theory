import numpy as np
import pytest
from scipy.optimize import linprog, minimize

from solver_v1.anchored_shape_calibration import ExactAnchorCoordinates, tangent_descent


def matrix(x):
    return np.array([[1.+.1*x, .2, .3], [.1, 1., -.1*x],
                     [.3, 2.+x, 1.], [-.5*x, -.3, 2.]])


def test_implicit_anchor_jacobian_matches_independent_differences():
    targets = np.array([2., 1., .2, -.5])
    scales = np.array([.2, 3., 2., .8])
    D = np.array([.3, 4., 2.])
    chart = ExactAnchorCoordinates(matrix(.2), targets, scales, [0, 1], D)
    dm = (matrix(.3)-matrix(.1))/.2
    z = np.array([.2, .4, 1.])

    def evaluate(v):
        return chart.evaluate(matrix(v[0]), dm[None], v[1:-1], v[-1], [2, 3], [0, 2])

    out = evaluate(z)
    for h in (1e-4, 5e-5):
        numerical = np.column_stack([(evaluate(z+h*e)['constraints']-
                    evaluate(z-h*e)['constraints'])/(2*h) for e in np.eye(3)])
        np.testing.assert_allclose(out['jacobian'], numerical, atol=2e-9, rtol=2e-8)
    for x in (-1., .2, 1.7):
        actual = evaluate(np.array([x, -.2, 3.]))
        assert actual['exact_residual'] < 2e-14
        np.testing.assert_allclose(matrix(x)[:2]@actual['coefficients'], targets[:2], atol=3e-15)


def test_reduced_nonlinear_solver_matches_independent_fixed_shape_lp():
    M = matrix(.2)
    targets = np.array([2., 1., .2, -.5])
    scales = np.array([.2, 3., 2., .8])
    D = np.array([.3, 4., 2.])
    chart = ExactAnchorCoordinates(M, targets, scales, [0, 1], D)
    def evaluate(v):
        return chart.evaluate(M, np.empty((0, *M.shape)), v[:-1], v[-1], [2, 3], [0, 1, 2])
    result = minimize(lambda v:v[-1], [0., 100.], jac=lambda v:np.array([0., 1.]),
        constraints=[dict(type='ineq', fun=lambda v:evaluate(v)['constraints'],
                          jac=lambda v:evaluate(v)['jacobian'])],
        bounds=[(None, None), (0., None)], method='SLSQP', options=dict(ftol=1e-12,maxiter=100))
    J = M[2:]/scales[2:, None]
    rhs = targets[2:]/scales[2:]
    reference = linprog([0., 0., 0., 1.], A_eq=np.column_stack([M[:2], np.zeros(2)]),
        b_eq=targets[:2], A_ub=np.vstack([np.column_stack([J, -np.ones(2)]),
        np.column_stack([-J, -np.ones(2)])]), b_ub=np.r_[rhs, -rhs], bounds=[(0., None)]*4)
    assert result.success and reference.success
    assert result.fun == pytest.approx(reference.fun, abs=1e-9)
    assert evaluate(result.x)['exact_residual'] < 1e-13


def test_singular_anchor_is_rejected_without_regularization():
    M = np.ones((4, 3))
    with pytest.raises(ArithmeticError, match='budget'):
        ExactAnchorCoordinates(M, np.ones(4), np.ones(4), [0, 1], np.ones(3))


def test_pivot_and_physical_coefficients_survive_numerical_rescaling():
    M, t = matrix(.2), np.array([2., 1., .2, -.5])
    coefficients = np.array([1.8, .92, .2])
    t[:2] = M[:2]@coefficients
    for D in (np.ones(3), np.array([.03, 40., .2])):
        chart = ExactAnchorCoordinates(M, t, np.ones(4), [0, 1], D)
        y = (coefficients/D)[chart.free]
        out = chart.evaluate(M, np.empty((0, *M.shape)), y, 10., [2, 3], [])
        np.testing.assert_allclose(out['coefficients'], coefficients, rtol=1e-14, atol=1e-14)


def test_tangent_cone_respects_nonsmooth_cusp_and_active_shape_bounds():
    J = np.array([[-1., 1.], [1., 1.]])
    cusp = tangent_descent([0., 0.], J, [0., 0.], [-2., 0.], [2., np.inf], 1e-8)
    assert cusp['minimum_linearized_deta'] == 0.
    slope = tangent_descent([0., 2.], J, [1., 1.], [-2., 0.], [2., np.inf], 1e-8)
    assert slope['minimum_linearized_deta'] == -1.
    # The same point becomes stationary only under the explicitly imposed x>=1.
    bound = tangent_descent([0., 2.], J, [1., 1.], [1., 0.], [2., np.inf], 1e-8)
    assert bound['minimum_linearized_deta'] == 0.
    assert not bound['global_shape_optimum_certified']
