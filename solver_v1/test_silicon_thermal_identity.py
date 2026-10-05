import numpy as np
import pytest
from scipy.integrate import quad

from solver_v1.silicon_conditional_research import KB_EV_K
from solver_v1.silicon_thermal_identity import thermal_test_field, canonical_temperature_terms


def terms(x, force, *, center=3., reference=3., radius=1.5, wall=1., width=.4):
    r = np.array([[0.,0.,0.], [x,0.,0.]])
    c = np.array([[0.,0.,0.], [center,0.,0.]])
    ref = np.array([[0.,0.,0.], [reference,0.,0.]])
    f = np.array([[100.,200.,300.], [force,0.,0.]])
    return canonical_temperature_terms(r, f, np.array([False,True]), c, ref,
        temperature_K=1/KB_EV_K, halfwidth_A=radius, minimum_pair_A=wall, taper_width_A=width)


def test_nonlinear_box_identity_by_independent_quadrature_with_nonzero_boundary_density():
    energy = lambda x: .7*x*x+.1*x**3+.2*x**4
    gradient = lambda x: 1.4*x+.3*x*x+.8*x**3
    def integrand(x):
        result = terms(3+x, -gradient(x))
        return (result['numerator_components_eV'][0,0]
                -result['divergence_components'][0,0])*np.exp(-energy(x))
    assert abs(quad(integrand, -1.5, 1.5, epsabs=1e-12)[0]) < 2e-13
    # The unbounded virial formula has a missing box boundary term here.
    naive = quad(lambda x:(x*gradient(x)-1)*np.exp(-energy(x)), -1.5, 1.5)[0]
    assert abs(naive) > .1


def test_pair_wall_requires_the_derivative_of_the_pair_taper():
    energy = lambda x: .3*(x-1.3)**2+.2*(x-1.3)**4
    gradient = lambda x: .6*(x-1.3)+.8*(x-1.3)**3
    def integrand(x):
        result = terms(x, -gradient(x), center=1.4, reference=1.3, radius=.5)
        return (result['numerator_components_eV'][0,0]
                -result['divergence_components'][0,0])*np.exp(-energy(x))
    assert abs(quad(integrand, 1., 1.9, points=[1.4], epsabs=1e-12)[0]) < 2e-13
    # Omitting the pair correction is not cured by having the box factors.
    def missing_pair_derivative(x):
        result = terms(x, -gradient(x), center=1.4, reference=1.3, radius=.5)
        b = 1-((x-1.4)/.5)**2-2*(x-1.3)*(x-1.4)/.5**2
        return (result['numerator_components_eV'][0,0]-result['pair_taper']*b)*np.exp(-energy(x))
    assert abs(quad(missing_pair_derivative, 1., 1.9, points=[1.4])[0]) > .02


def test_full_divergence_matches_independent_coordinate_differences_near_multiple_pair_walls():
    r = np.array([[0.,0.,0.], [1.15,.1,0.], [.25,1.12,.1]])
    mask = np.array([False,True,True]); center = r.copy()
    ref = r.copy(); ref[mask] -= np.array([.08,-.03,.05])
    options = dict(halfwidth_A=.7, minimum_pair_A=1., taper_width_A=.4)
    field = lambda x: thermal_test_field(x, mask, center, ref, **options)
    actual = field(r); finite = 0.; step = 1e-5
    for i in np.flatnonzero(mask):
        for j in range(3):
            plus, minus = r.copy(), r.copy(); plus[i,j] += step; minus[i,j] -= step
            finite += (field(plus)['field_A'][i,j]-field(minus)['field_A'][i,j])/(2*step)
    assert finite == pytest.approx(actual['divergence'], abs=2e-8)
    assert actual['active_pair_tapers'] == 3
    assert np.array_equal(actual['field_A'][~mask], np.zeros((1,3)))


def test_pair_boundary_field_and_divergence_vanish_without_clipping():
    result = terms(1., 4., center=1.4, reference=1.3, radius=.5)
    assert result['pair_taper'] == 0.
    assert result['numerator_eV'] == 0.
    assert result['denominator'] == 0.


def test_fixed_grip_forces_do_not_enter_the_thermal_identity():
    r = np.array([[0.,0.,0.], [3.1,.1,0.]])
    c = np.array([[0.,0.,0.], [3.,0.,0.]])
    mask = np.array([False,True]); f = np.ones_like(r)
    options = dict(temperature_K=300., halfwidth_A=1., minimum_pair_A=1., taper_width_A=.3)
    a = canonical_temperature_terms(r, f, mask, c, c, **options)
    f[0] = 1e10
    b = canonical_temperature_terms(r, f, mask, c, c, **options)
    assert a['residual_eV'] == b['residual_eV']
    assert a['thermal_energy_eV'] == pytest.approx(.025851999786435)


@pytest.mark.parametrize('x,center,radius', [(1.99,3.,1.), (.99,1.4,.5)])
def test_domain_violations_are_rejected(x, center, radius):
    with pytest.raises(ValueError):
        terms(x, 1., center=center, radius=radius)
