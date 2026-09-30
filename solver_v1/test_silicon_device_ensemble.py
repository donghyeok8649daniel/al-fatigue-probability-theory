"""Independent matrix and quadrature checks, not material validation."""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.constants import electron_volt
from solver_v1.silicon_conditional_research import KB_EV_K
from solver_v1.silicon_device_ensemble import device_response, FixedGripTarget, EV_A2_TO_N_MM


def test_device_response_matches_direct_coupled_gaussian_inverse():
    h = np.array([[3., .6], [.6, 2.]])
    v = np.array([.2, 1/6])
    for k in (0., .1, 2., 100.):
        r = device_response(h, v, k, 300.)
        direct = np.linalg.inv(h+k*np.outer(v, v))
        assert r['extension_variance_A2'] == pytest.approx(KB_EV_K*300*v@direct@v, rel=2e-14)
        assert r['d_extension_d_actuator'] == pytest.approx(k*v@direct@v, abs=2e-14)
        x = direct@(k*.07*v)
        assert v@x == pytest.approx(.07*r['d_extension_d_actuator'], abs=2e-14)


def test_normalized_grip_coordinate_and_units_are_not_conflated():
    physical_k, norm = 2.3, 6.
    r = device_response([[physical_k/norm**2]], [1/norm], 0., 300.)
    assert r['specimen_stiffness_eV_A2'] == pytest.approx(physical_k)
    assert r['extension_variance_A2'] == pytest.approx(KB_EV_K*300/physical_k)
    assert EV_A2_TO_N_MM*1000 == pytest.approx(electron_volt/1e-20)
    fixed = device_response([[physical_k/norm**2]], [1/norm], np.inf, 300.)
    assert fixed['extension_variance_A2'] == 0.
    assert fixed['d_extension_d_actuator'] == 1.


def test_same_device_formula_matches_independent_nonlinear_energy_derivative():
    # Analytic two-spring minimization obtained by direct energy perturbations.
    ks, km, x = 2.7, 1.3, .15
    energy = lambda d: .5*ks*d*d+.5*km*(x-d)**2
    r = device_response([[ks]], [1.], km, 300.)
    d = r['d_extension_d_actuator']*x
    assert abs((energy(d+1e-5)-energy(d-1e-5))/2e-5) < 1e-10
    for bad in (-1., np.nan):
        with pytest.raises(ValueError): device_response([[ks]], [1.], bad, 300.)
    with pytest.raises(np.linalg.LinAlgError): device_response([[-1.]], [1.], 1., 300.)


def test_nonlinear_gaussian_reference_cancels_from_same_bounded_target():
    positions = np.array([[0., 0., 0.], [3., 0., 0.]])
    b = np.zeros((6, 1)); b[3, 0] = 1.
    energy = lambda x: .7*x*x+.12*x**3+.3*x**4
    t = 1/KB_EV_K
    for curvature, gradient in ((1.4, .2), (3.1, -.4)):
        target = FixedGripTarget(positions, b, [[curvature]], [gradient], t, positions, .7, 1.)
        def mapped(x):
            u = target.reference.whiten([x])
            return np.exp(-.5*u@u-target.correction(u, energy(x), 0.))
        z = quad(mapped, -.7, .7, epsabs=1e-12)[0]
        z_direct = quad(lambda x: np.exp(-energy(x)), -.7, .7, epsabs=1e-12)[0]
        assert z == pytest.approx(z_direct, abs=1e-13)
        inside = target.coordinates(target.reference.whiten([.4]))
        assert target.inside(inside)
        assert np.array_equal(inside[0], positions[0])
        assert not target.inside(target.coordinates(target.reference.whiten([.8])))
        moved = inside.copy(); moved[0, 0] += .01
        assert not target.inside(moved)


def test_common_box_contains_both_starts_without_changing_the_fixed_grips():
    r = np.array([[0., 0., 0.], [3., 0., 0.]])
    b = np.zeros((6, 1)); b[3, 0] = 1.
    center = r.copy(); center[1, 0] += .4
    target = FixedGripTarget(r, b, [[2.]], [0.], 300., center, 1., 1.)
    assert target.inside(r)
    center[0, 0] += .1
    with pytest.raises(ValueError): FixedGripTarget(r, b, [[2.]], [0.], 300., center)
