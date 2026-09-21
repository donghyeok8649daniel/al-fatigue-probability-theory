"""Independent integrals and detailed-balance checks for research ensembles."""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtri
from scipy.stats import qmc

from solver_v1.silicon_conditional_research import KB_EV_K
from solver_v1.silicon_thermal_research import (
    GaussianReference, exponential_reweight, pcn_sample, split_rhat, bounded_force_control,
    harmonic_split_proposal, harmonic_hmc_sample,
)


def test_gaussian_whitening_covariance_and_unstable_reference_rejection():
    h = np.array([[3., .7], [.7, 2.]])
    reference = GaussianReference.from_hessian(h, 300., [.2, -.3])
    # Deterministic sigma points have exactly unit mean and covariance.
    u = np.sqrt(2)*np.r_[np.eye(2), -np.eye(2)]
    x = reference.transform(u)
    np.testing.assert_allclose(reference.whiten(x), u, atol=2e-15)
    np.testing.assert_allclose((x-reference.mean).T@(x-reference.mean)/len(x),
                              KB_EV_K*300*np.linalg.inv(h), atol=2e-16)
    with pytest.raises(np.linalg.LinAlgError):
        GaussianReference.from_hessian(np.diag([1., -1.]), 300.)


def test_bounded_nonlinear_partition_matches_independent_quadrature():
    # One-dimensional quartic reference checks normalization, sign, support and
    # weighted moments against actual integration, not a mirrored formula.
    temperature = 1/KB_EV_K
    u = ndtri(qmc.Sobol(1, scramble=True, seed=34).random_base2(17)).ravel()
    residual = .07*u**4+.08*u**3
    inside = abs(u) < 2.
    result = exponential_reweight(residual, temperature, inside=inside, observables=u*u)
    integrand = lambda x: np.exp(-.5*x*x-.07*x**4-.08*x**3)
    z = quad(integrand, -2, 2, epsabs=1e-13)[0]
    exact_f = -np.log(z/np.sqrt(2*np.pi))
    exact_x2 = quad(lambda x:x*x*integrand(x), -2, 2, epsabs=1e-13)[0]/z
    np.testing.assert_allclose(result['correction_eV'], exact_f, atol=3e-5)
    np.testing.assert_allclose(result['mean'], exact_x2, atol=6e-5)


def test_pcn_proposal_preserves_reference_and_metropolis_satisfies_balance():
    beta, x, y = .63, np.array([.3, -.7]), np.array([-.2, .9])
    persistence = np.sqrt(1-beta*beta)
    logproposal = lambda a,b:-np.sum((b-persistence*a)**2)/(2*beta*beta)
    phi = lambda u:.12*np.sum(u**4)
    logpi = lambda u:-.5*u@u-phi(u)
    lhs = logpi(x)+logproposal(x,y)+min(0., phi(x)-phi(y))
    rhs = logpi(y)+logproposal(y,x)+min(0., phi(y)-phi(x))
    np.testing.assert_allclose(lhs, rhs, atol=2e-15)
    result = pcn_sample(lambda u:(0., np.r_[u, u*u]), [0.], proposal_scale=1.,
        draws=20000, warmup=0, rng=np.random.default_rng(512))
    assert result['acceptance_fraction'] == 1.
    assert result['physical_time_ps'] is None
    assert abs(result['observations'][:,0].mean()) < .025
    assert abs(result['observations'][:,1].mean()-1) < .035


def test_pcn_rejections_remain_in_history_and_bounded_target_is_correct():
    result = pcn_sample(lambda u:(0. if abs(u[0]) < .8 else np.inf, u*u),
        [0.], proposal_scale=1., draws=30000, warmup=50, rng=np.random.default_rng(441))
    z = quad(lambda x:np.exp(-x*x/2), -.8, .8)[0]
    moment = quad(lambda x:x*x*np.exp(-x*x/2), -.8, .8)[0]/z
    assert np.all(abs(result['samples']) < .8)
    assert .5 < result['acceptance_fraction'] < .65
    assert abs(result['observations'].mean()-moment) < .008
    assert np.any(np.diff(result['samples'][:,0]) == 0.)
    with pytest.raises(ValueError):
        pcn_sample(lambda u:(np.inf, u), [0.], proposal_scale=.5,
            draws=10, warmup=0, rng=np.random.default_rng(1))


def test_split_rhat_detects_offset_chains_without_calling_them_independent_draws():
    chains = np.random.default_rng(29).standard_normal((4,1000,2))
    assert np.max(split_rhat(chains)) < 1.01
    chains[0,:,0] += 3
    assert split_rhat(chains)[0] > 1.5


def test_force_control_integrates_to_zero_even_with_nonzero_boundary_density():
    radius, temperature = .6, 1/KB_EV_K
    energy = lambda x:.5*2*x*x+.3*x**3+.2*x
    gradient = lambda x:2*x+.9*x*x+.2
    def integrand(x):
        control = bounded_force_control([x],[gradient(x)],[.7],temperature,radius)
        return control*np.exp(-energy(x))
    assert abs(quad(integrand,-radius,radius,epsabs=1e-13)[0]) < 5e-14
    # The usual unbounded-domain v.gradU control is biased on this finite box.
    assert abs(quad(lambda x:.7*gradient(x)*np.exp(-energy(x)),-radius,radius)[0]) > .05


def test_harmonic_split_is_reversible_and_hmc_matches_actual_nonlinear_integral():
    evaluate=lambda u:(.1*np.sum(u**4),.4*u**3,u*u)
    x,p=np.array([.7,-.3]),np.array([.2,.4])
    proposal=harmonic_split_proposal(evaluate,x,p,.2,8)
    reverse=harmonic_split_proposal(evaluate,proposal[0],-proposal[1],.2,8)
    np.testing.assert_allclose(reverse[0],x,atol=3e-15)
    np.testing.assert_allclose(reverse[1],-p,atol=3e-15)
    result=harmonic_hmc_sample(evaluate,[0.],draws=12000,warmup=100,rng=np.random.default_rng(834))
    z=quad(lambda x:np.exp(-.5*x*x-.1*x**4),-np.inf,np.inf)[0]
    exact=quad(lambda x:x*x*np.exp(-.5*x*x-.1*x**4),-np.inf,np.inf)[0]/z
    assert abs(result['observations'].mean()-exact)<.02
    assert result['acceptance_fraction']>.98
    assert result['physical_time_ps'] is None
