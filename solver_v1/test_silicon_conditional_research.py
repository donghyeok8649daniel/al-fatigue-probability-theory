"""Independent Gaussian integration, projected dynamics and energy checks."""
import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from solver_v1.silicon_crack_research import AtomicEvaluation
from solver_v1.silicon_conditional_research import (
    AMU_EV_PS2_A2, KB_EV_K, conditional_harmonic, harmonic_free_energy_difference,
    harmonic_memory, harmonic_release, velocity_verlet,
)


def test_gaussian_measure_matches_actual_conditional_integral():
    # Physical q is twice the first Cartesian coordinate: the coarea factor
    # is constant 1/2, not an energy scale or an atom-count factor.
    matrices = [np.array([[4., 1.], [1., 3.]]), np.array([[5., 2.], [2., 7.]])]
    q, temperature = .13, 300.
    reduced = [conditional_harmonic(h, [[0.], [1.]], [[.5], [0.]]) for h in matrices]
    minima = [float(.5*q*q*r.relaxed_curvature[0, 0]) for r in reduced]
    partition = []
    for h in matrices:
        def integrand(z):
            r = np.array([q/2, z])
            return .5*np.exp(-.5*r@h@r/(KB_EV_K*temperature))
        partition.append(quad(integrand, -np.inf, np.inf, epsabs=1e-13)[0])
    exact = -KB_EV_K*temperature*np.log(partition[1]/partition[0])
    candidate = harmonic_free_energy_difference(minima[1]-minima[0],
        reduced[1].logdet_bath, reduced[0].logdet_bath, temperature)
    np.testing.assert_allclose(candidate, exact, rtol=1e-12, atol=1e-14)
    assert harmonic_free_energy_difference(.7, 9, 3, 0.) == .7


def test_two_coordinate_marginal_recovers_one_coordinate_including_metric():
    h = np.array([[5., .4, .8, -.2], [.4, 4., -.3, .6],
                  [.8, -.3, 3., .5], [-.2, .6, .5, 2.]])
    eye = np.eye(4)
    # Two physical relative coordinates each have D.T D=1/2.
    two = conditional_harmonic(h, eye[:, 2:], eye[:, :2]/np.sqrt(2))
    one = conditional_harmonic(h, eye[:, 1:], eye[:, :1]/np.sqrt(2))
    np.testing.assert_allclose(one.logdet_bath,
        two.logdet_bath+np.log(2*two.relaxed_curvature[1, 1]), atol=2e-14)
    k = two.relaxed_curvature
    np.testing.assert_allclose(one.relaxed_curvature, k[0, 0]-k[0, 1]**2/k[1, 1], atol=2e-14)
    angle = .481
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    rotated = conditional_harmonic(h, eye[:, 2:]@rotation, eye[:, :2]/np.sqrt(2))
    np.testing.assert_allclose(rotated.logdet_bath, two.logdet_bath, atol=2e-14)
    np.testing.assert_allclose(rotated.relaxed_curvature, two.relaxed_curvature, atol=2e-14)


@pytest.mark.parametrize('bath', [0., -1.])
def test_unstable_conditional_partition_is_not_regularized(bath):
    with pytest.raises(np.linalg.LinAlgError):
        conditional_harmonic(np.diag([2., bath]), [[0.], [1.]], [[1.], [0.]])


def test_memory_and_full_response_match_independent_coupled_ode():
    h = np.array([[4., 1.], [1., 2.]])
    reduction = conditional_harmonic(h, [[0.], [1.]], [[1.], [0.]])
    mass_amu = 1/AMU_EV_PS2_A2  # unit inertial mass in this synthetic system
    times = np.linspace(0, 2, 101)
    memory = harmonic_memory(reduction, times, mass_amu=mass_amu)
    np.testing.assert_allclose(memory['kernel_eV_A2'][:, 0, 0], .5*np.cos(np.sqrt(2)*times), atol=1e-14)
    np.testing.assert_allclose(memory['q_mass_eV_ps2_A2'][0, 0], 1., atol=2e-15)
    x0 = reduction.relaxed_lift[:, 0]
    solution = solve_ivp(lambda t, y: np.r_[y[2:], -h@y[:2]], [0, 2], np.r_[x0, 0., 0.],
                         method='DOP853', rtol=2e-13, atol=2e-14, dense_output=True)
    prediction = harmonic_release(h, x0, np.array([1., 0.]), times, mass_amu=mass_amu)
    np.testing.assert_allclose(prediction, solution.sol(times)[0], atol=2e-12)
    for t in (.1, .7, 1.3, 2.):
        state = solution.sol(t)
        history = quad(lambda u: .5*np.cos(np.sqrt(2)*(t-u))*solution.sol(u)[2],
                       0, t, epsabs=1e-12)[0]
        # Missing the conditional preparation or the convolution changes this.
        assert abs(-(h@state[:2])[0]+3.5*state[0]+history) < 5e-12


class QuadraticAtoms:
    def evaluate(self, r):
        return AtomicEvaluation(float(np.sum(r*r)), 2*r, np.sum(r*r, axis=1))


def test_verlet_energy_order_fixed_grips_reversal_and_raw_observations():
    positions = np.array([[0., 0., 0.], [.1, 0., 0.]])
    errors = []
    mass = 1/AMU_EV_PS2_A2
    for dt in (.04, .02, .01):
        result = velocity_verlet(QuadraticAtoms(), positions, [True, False],
            dt_ps=dt, steps=round(2/dt), mass_amu=mass, observation=lambda r:r[1, 0])
        np.testing.assert_allclose(result['observations'], .1*np.cos(np.sqrt(2)*result['time_ps']), atol=4e-5)
        assert result['max_fixed_displacement_A'] == 0.
        assert result['production_t0_seconds'] is None
        errors.append(result['max_energy_residual_eV'])
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1
    backward = velocity_verlet(QuadraticAtoms(), result['final_positions_A'], [True, False],
        velocities=-result['final_velocities_A_ps'], dt_ps=.01, steps=200, mass_amu=mass,
        observation=lambda r:r[1, 0])
    np.testing.assert_allclose(backward['final_positions_A'], positions, atol=2e-15)
    np.testing.assert_allclose(backward['final_velocities_A_ps'], 0., atol=2e-15)


def test_nonorthonormal_basis_and_invalid_dynamics_are_rejected():
    with pytest.raises(ValueError):
        conditional_harmonic(np.eye(2), [[0.], [2.]], [[1.], [0.]])
    with pytest.raises(ValueError):
        harmonic_release(np.diag([1., -1.]), [1., 0.], [1., 0.], [0., 1.])
    with pytest.raises(ValueError):
        harmonic_free_energy_difference(1., 2., 3., -10.)
    with pytest.raises(ValueError):
        velocity_verlet(QuadraticAtoms(), [[0., 0., 0.]], [True],
            dt_ps=.1, steps=1, observation=lambda r:r[0])
    with pytest.raises(ValueError):
        velocity_verlet(QuadraticAtoms(), [[0., 0., 0.], [1., 0., 0.]], [True, False],
            velocities=[[.1, 0., 0.], [0., 0., 0.]], dt_ps=.1, steps=1, observation=lambda r:r[0])


def test_retained_coordinate_variance_is_gaussian_conditioning_not_parameter_fit():
    from results.silicon_wafer_feasibility.analyze_conditional_dynamics import explained_force_fraction
    h = np.diag([2., 8.])
    partial = explained_force_fraction(h, [3., 4.], [[1., 0.]])
    np.testing.assert_allclose(partial['explained_force_variance_fraction'], 9/13, atol=2e-15)
    full = explained_force_fraction(h, [3., 4.], np.eye(2))
    np.testing.assert_allclose(full['explained_force_variance_fraction'], 1., atol=2e-15)
    with pytest.raises(np.linalg.LinAlgError):
        explained_force_fraction(h, [3., 4.], [[1., 0.], [1., 0.]])
