import numpy as np
import pytest
from scipy.integrate import quad

from .core_matching import radial_window,log_window_constant,finite_part_samples


@pytest.mark.parametrize('alpha',[.25,.5,.75,.99])
def test_window_constant_against_independent_integral(alpha):
    expected=quad(lambda t:float(radial_window(t,alpha)-1)/t,alpha,1,epsabs=1e-13)[0]
    actual=log_window_constant(alpha)
    assert abs(actual['value']-expected)<actual['estimated_series_tail']+3e-15
    if alpha==.5:
        assert actual['value']==pytest.approx(19/6-5*np.log(2),abs=3e-15)


def test_window_endpoints_and_continuum_reference_cancel():
    np.testing.assert_array_equal(radial_window([0,.5,1,2]),[1,1,0,0])
    # Synthetic continuum log density, NOT an atomistic core-energy target.
    for alpha in (.5,.75):
        for R in (4.,8.,16.):
            integral=quad(lambda r:float(radial_window(r/R,alpha))/r,1.,R,
                          points=[alpha*R],epsabs=2e-13)[0]
            assert abs(integral-np.log(R)-log_window_constant(alpha)['value'])<1e-13


def test_reference_translation_changes_only_declared_log_constant():
    distances=np.array([.5,1.,2.,3.,4.]); energy=np.array([.2,.15,.13,.09,.07])
    options=dict(distances=distances,remainder_energy=energy,radii=[4.],
                 log_coefficient=.3,free_radius=5.)
    first=finite_part_samples(**options,b=1.)
    second=finite_part_samples(**options,b=.8)
    for a,z in zip(first,second):
        assert z['finite_part_at_reference_b_eV']-a['finite_part_at_reference_b_eV']==pytest.approx(.3*np.log(.8),abs=1e-15)
        assert not z['core_radius_fitted'] and not z['convergence_certified']
    with pytest.raises(ValueError):
        finite_part_samples(**{**options,'radii':[6.]},b=1.)
