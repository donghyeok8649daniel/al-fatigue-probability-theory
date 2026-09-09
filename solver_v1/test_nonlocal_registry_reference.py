"""Analytic test misfit ONLY; production remains the LJ/Bessel energy."""
import numpy as np
import pytest
from solver_v1.nonlocal_registry_reference import PeriodicMisfit, ScrewRegistryFunctional


def sinusoidal_test(m=4., length=128., cells=4096):
    b = 1.
    s = np.arange(64)/64
    potential = PeriodicMisfit.from_samples(.2*(1-np.cos(2*np.pi*s)), period=b)
    return ScrewRegistryFunctional(potential, elastic_factor=m, domain_length=length, cells=cells)


def test_spectral_misfit_derivatives_and_periodicity():
    p = sinusoidal_test().misfit
    s = np.linspace(-1.3,2.1,53)
    np.testing.assert_allclose(p.evaluate(s), .2*(1-np.cos(2*np.pi*s)), atol=1e-14)
    np.testing.assert_allclose(p.evaluate(s,1), .4*np.pi*np.sin(2*np.pi*s), atol=6e-14)
    np.testing.assert_allclose(p.evaluate(s,2), .8*np.pi**2*np.cos(2*np.pi*s), atol=9e-12)
    np.testing.assert_allclose(p.evaluate(s+1), p.evaluate(s), atol=2e-15)


def test_variational_gradient_hessian_and_constant_translation():
    f = sinusoidal_test(cells=256)
    s = .11*np.cos(2*np.pi*f.x/f.length)+.05
    v = np.cos(6*np.pi*f.x/f.length)+.3
    w = np.sin(4*np.pi*f.x/f.length)-.2
    step = 1e-6
    fd = (f.evaluate(s+step*v)["energy"]-f.evaluate(s-step*v)["energy"])/(2*step)
    assert fd == pytest.approx(f.dx*np.dot(f.evaluate(s)["force"],v), rel=3e-9)
    np.testing.assert_allclose((f.evaluate(s+step*v)["force"]-f.evaluate(s-step*v)["force"])/(2*step),
                              f.hessian_action(s,v), rtol=3e-8, atol=1e-8)
    assert np.dot(w,f.hessian_action(s,v)) == pytest.approx(np.dot(v,f.hessian_action(s,w)),abs=1e-12)
    np.testing.assert_allclose(f.elastic_force(np.ones(f.cells)),0,atol=0)
    assert f.evaluate(s+1)["energy"] == pytest.approx(f.evaluate(s)["energy"],abs=1e-12)


def test_nonlocal_linear_response_and_refinement():
    f = sinusoidal_test()
    mode = np.cos(8*np.pi*f.x/f.length)
    amplitude = 1e-6
    solution = f.linear_response(amplitude*mode)
    exact = amplitude/(float(f.misfit.evaluate(0,2))+f.mu/2*8*np.pi/f.length)
    expected = exact*mode-f.misfit.evaluate(0,1)/f.misfit.evaluate(0,2)
    np.testing.assert_allclose(solution,expected,atol=1e-20)
    assert max(abs(f.evaluate(solution,amplitude*mode)["force"])) < 1e-17


def test_periodized_dipole_closed_form_and_conjugate_work():
    f = sinusoidal_test()
    d,w = 20.,.7
    s,sd,sw = f.dipole_trial(d,w)
    assert np.mean(s) == pytest.approx(d/f.length,abs=1e-15)
    assert f.evaluate(s)["elastic"] == pytest.approx(f.dipole_elastic_closed_form(d,w),rel=1e-13)
    for position,derivative in ((0,sd),(1,sw)):
        plus=[d,w]; minus=[d,w]; step=1e-5
        plus[position]+=step; minus[position]-=step
        np.testing.assert_allclose((f.dipole_trial(*plus)[0]-f.dipole_trial(*minus)[0])/(2*step),derivative,atol=5e-11)
    tau=.03
    assert f.evaluate(s,tau)["external_work"] == pytest.approx(-tau*d,abs=1e-14)


def test_nonlinear_intrawell_solution_unloads_without_fabricating_plasticity():
    f=sinusoidal_test(cells=512)
    traction=.15*np.cos(2*np.pi*f.x/f.length)
    solution=f.solve_intrawell_equilibrium(traction)
    assert solution["force_residual"]<1e-11
    assert solution["iterations"]>0
    mirror=f.solve_intrawell_equilibrium(-traction)
    np.testing.assert_allclose(mirror["slip"],-solution["slip"],atol=5e-14)
    unloaded=f.solve_intrawell_equilibrium(0.,initial=solution["slip"])
    assert max(abs(unloaded["slip"]))<1e-12
    with pytest.raises(ValueError,match="nonconvex"):
        f.solve_intrawell_equilibrium(0.,initial=np.full(f.cells,.3))
    # A force-free saddle must not be accepted before the stability check.
    with pytest.raises(ValueError,match="nonconvex"):
        f.solve_intrawell_equilibrium(0.,initial=np.full(f.cells,.5))


def test_width_optimization_is_not_full_equilibrium_or_nucleation_barrier():
    f = sinusoidal_test()
    result = f.optimize_trial_width(25.,width_bounds=(.02,3.))
    assert result["optimizer_success"] and result["width_interior"]
    assert abs(result["width_force"]) < 2e-7
    assert result["balance_traction"] > 0
    # A width-only trial does NOT satisfy the unrestricted Euler equation.
    assert result["full_euler_residual"] > 1e-3


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        PeriodicMisfit.from_samples(np.zeros(13),period=1.)
    with pytest.raises(ValueError):
        sinusoidal_test().dipole_trial(90.,.5)


def test_all_profile_equilibrium_removes_the_trial_shape_force_error():
    f = sinusoidal_test(length=128.,cells=2048)
    trial = f.dipole_trial(20.,.3)[0]
    solved = f.solve_fixed_registry_content(trial,force_tolerance=1e-10)
    assert solved["force_converged"]
    assert solved["projected_force_residual"] < 1e-10
    assert solved["mean_constraint_residual"] < 1e-15
    assert solved["energy"] < f.evaluate(trial)["energy"]
    # The constant constraint reaction is a holding traction, NOT an error or
    # arbitrary fitted Peierls/yield stress. All nonconstant Euler modes vanish.
    assert solved["holding_traction"] > 0
    np.testing.assert_allclose(f.evaluate(solved["slip"])["force"],
                               solved["holding_traction"],atol=1e-10)


def test_registry_content_reaction_is_the_variational_energy_derivative():
    f = sinusoidal_test(length=128.,cells=2048)
    d,h = 20.,.002
    states = [f.solve_fixed_registry_content(f.dipole_trial(v,.3)[0],force_tolerance=2e-10)
              for v in (d-h,d,d+h)]
    assert all(s["force_converged"] for s in states)
    derivative = (states[2]["energy"]-states[0]["energy"])/(2*h)
    assert derivative == pytest.approx(states[1]["holding_traction"]*f.misfit.period,rel=1e-6)


def test_fixed_mean_solution_does_not_turn_optimizer_success_into_a_physical_certificate():
    f = sinusoidal_test(length=64.,cells=1024)
    result = f.solve_fixed_registry_content(f.dipole_trial(8.,.4)[0],force_tolerance=1e-10)
    assert result["force_converged"]
    assert "NOT atomistic core" in result["status"]
    assert "physical_time" not in result and "activation_energy_eV" not in result
