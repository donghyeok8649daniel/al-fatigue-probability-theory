import numpy as np
import pytest

from .isotropic_bulk_validation import IsotropicBulkBasis
from .quartic_angular_material import QuarticSymmetryTailBulkBasis


def test_reference_limit_preserves_existing_harmonic_columns_and_fixed_gauge():
    decays=(2.8,5.6,10.2)
    old=QuarticSymmetryTailBulkBasis(decays,radius=7.)
    current=IsotropicBulkBasis(decays,stretch=1.,radius=7.)
    assert current.x==pytest.approx(1.,abs=2e-14)
    for key in current.gauge:
        np.testing.assert_array_equal(current.gauge[key],old.gauge[key])
    for q in ((.2,.1,.0),(.5,.5,.5)):
        a,b=current.evaluate(q);ref,tail=old.evaluate(q)
        np.testing.assert_allclose(a,ref,rtol=2e-13,atol=2e-14)
        np.testing.assert_allclose(b,tail,rtol=3e-13,atol=1e-15)
    expanded=IsotropicBulkBasis(decays,stretch=4.065/4.05,radius=7.)
    assert expanded.reference_density==current.reference_density
    for key in current.gauge:
        np.testing.assert_array_equal(expanded.gauge[key],current.gauge[key])
    assert expanded.x<1.
    assert expanded.geometry.b==pytest.approx(4.065/4.05)


def test_strained_analytic_columns_match_independent_direct_site_energies():
    decays=(3.,5.5,9.)
    c=np.array([.33,1.5,10.,16.,.3,4.,8.,.8,.26,1e4,25.])
    obj=IsotropicBulkBasis(decays,stretch=1.007,radius=11.,include_cross=True)
    planes=6;mode=2;q=np.ones(3)*mode/planes
    columns,_=obj.evaluate(q);H=np.einsum('c,cij->ij',c,columns)
    assert obj.x!=1.
    np.testing.assert_allclose(columns[-1],(obj.x-1)*columns[5])
    for v in np.eye(3):
        errors=[];approximations=[]
        for amplitude in (4e-4,2e-4):
            kw=dict(planes=planes,mode=mode,polarization_plane=v,validation_radius=11.)
            zero=obj.direct_sinusoidal_energy(c,amplitude=0.,**kw)
            plus=obj.direct_sinusoidal_energy(c,amplitude=amplitude,**kw)
            minus=obj.direct_sinusoidal_energy(c,amplitude=-amplitude,**kw)
            numeric=2*(plus+minus-2*zero)/amplitude**2
            approximations.append(numeric)
            errors.append(abs(numeric-v@H@v))
        # Equal direct/operator radii isolate derivative refinement; omitted
        # neighbor tails are tested independently below, not mixed into this.
        # A large quartic amplitude gives a visible O(amplitude²) central
        # difference error, especially for the normal polarization. Test that
        # convergence and its fourth-order Richardson extrapolation explicitly.
        assert errors[-1]<.30*errors[0]
        richardson=(4*approximations[1]-approximations[0])/3
        assert abs(richardson-v@H@v)<1e-5


def test_strained_neighbor_tail_controls_radius_change():
    decays=(3.,5.5,9.)
    low=IsotropicBulkBasis(decays,stretch=.99,radius=7.,include_cross=True)
    high=IsotropicBulkBasis(decays,stretch=.99,radius=11.,include_cross=True)
    for q in ((.01,.01,0.),(.5,.5,.5)):
        H,b=low.evaluate(q);Hfine,_=high.evaluate(q)
        assert np.all(np.linalg.norm(Hfine-H,ord=2,axis=(-2,-1))<=b+2e-12)
