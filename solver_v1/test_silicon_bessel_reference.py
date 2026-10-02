"""Direct SW and derivative controls for the Fourier--Bessel research adapter."""
import numpy as np
import pytest
from scipy.optimize import brentq

from results.silicon_wafer_feasibility.run_static_probe import source_parameters, pair, analytic_shuffle_hessian
from .silicon_bessel_reference import SWPlaneBessel, SW111BesselInterface, _point_moments


@pytest.fixture(scope='module')
def reference():
    p, _ = source_parameters()
    bond = brentq(lambda r: float(pair(r,p,derivative=True)),2.2,2.5,xtol=1e-14)
    return p,bond,4*bond/np.sqrt(3)


@pytest.fixture(scope='module')
def plane(reference):
    p,_,lattice = reference
    return SWPlaneBessel(p,lattice,shell_index=96,nodes=512)


@pytest.fixture(scope='module')
def interfaces(reference):
    p,_,lattice = reference
    return {kind: SW111BesselInterface(p,lattice,cut_kind=kind,shell_index=96,nodes=512)
            for kind in ['shuffle','glide']}


@pytest.mark.parametrize('distance,delta', [(2.35,[0.,0.]),(.79,[1.7,.9]),(-1.2,[.31,-.42]),(3.5,[.2,.3])])
def test_plane_all_moments_match_exact_neighbor_sum(plane,distance,delta):
    actual = plane.evaluate(distance,delta)
    expected = sum((_point_moments(v,plane.p) for v in plane.direct_vectors(distance,delta)),np.zeros(15))
    np.testing.assert_allclose(actual.value,expected,rtol=2e-9,atol=3e-9)
    np.testing.assert_allclose(actual.value[6:].reshape(3,3).trace(),0.,atol=2e-11)


def test_phase_periodicity_and_plane_inversion(plane):
    base = plane.evaluate(1.2,[.4,-.3])
    shifted = plane.evaluate(1.2,np.array([.4,-.3])+2*plane.geometry.a1-plane.geometry.a2)
    np.testing.assert_allclose(base.value,shifted.value,atol=2e-11)
    np.testing.assert_allclose(base.gradient,shifted.gradient,atol=5e-10)
    np.testing.assert_allclose(base.hessian,shifted.hessian,atol=2e-9)
    reverse = plane.evaluate(-1.2,[-.4,.3])
    parity = np.ones(15); parity[3:6] = -1
    np.testing.assert_allclose(reverse.value,parity*base.value,atol=2e-11)


def test_plane_derivatives_d_and_lateral_against_value_differences(plane):
    q = np.array([1.2,.31,-.42]); h = 2e-4
    a = plane.evaluate(q[0],q[1:])
    for i in range(3):
        step = np.eye(3)[i]*h
        plus,minus = [plane.evaluate(x[0],x[1:]) for x in [q+step,q-step]]
        np.testing.assert_allclose((plus.value-minus.value)/(2*h),a.gradient[:,i],rtol=5e-6,atol=3e-6)
        np.testing.assert_allclose((plus.gradient-minus.gradient)/(2*h),a.hessian[:,:,i],rtol=2e-5,atol=2e-5)


@pytest.mark.parametrize('kind',['shuffle','glide'])
@pytest.mark.parametrize('q', [[0.,0.,0.],[.21,.15,-.12],[1.1,.32,.21],[5.,.13,-.07]])
def test_interface_all_site_energy_gradient_hessian_against_direct(interfaces,kind,q):
    model = interfaces[kind]
    actual,expected = [model.evaluate(q,method=m) for m in ['bessel','direct']]
    assert abs(actual.value-expected.value) < 2e-8
    np.testing.assert_allclose(actual.gradient,expected.gradient,atol=2e-7)
    np.testing.assert_allclose(actual.hessian,expected.hessian,atol=3e-6)


def test_shuffle_zero_hessian_matches_existing_independent_bond_angle_formula(reference,interfaces):
    p,bond,lattice = reference
    model = interfaces['shuffle']
    actual = model.evaluate([0.,0.,0.])
    expected = analytic_shuffle_hessian(bond,p)
    np.testing.assert_allclose(actual.gradient,0.,atol=2e-9)
    np.testing.assert_allclose(actual.hessian[:2,:2],expected,atol=2e-7)
    np.testing.assert_allclose(actual.hessian[2,2],expected[1,1],atol=2e-7)


def test_self_subtraction_and_sum_before_angular_square(reference):
    from .silicon_bessel_reference import site_energy, PlaneMoments
    p,_,_ = reference
    v = np.array([.2,-.3,2.4])
    value = _point_moments(v,p)
    m = PlaneMoments(value,np.zeros((15,3)),np.zeros((15,3,3)))
    actual = site_energy(m,p)
    assert abs(actual.value-.5*value[0]) < 1e-13
    incomplete = value.copy(); incomplete[2] = 0.
    broken = site_energy(PlaneMoments(incomplete,m.gradient,m.hessian),p)
    assert abs(broken.value-actual.value) > .01


def test_cutoff_d_zero_and_nonordered_planes_are_explicit(reference,plane):
    cutoff = plane.p['a']*plane.p['sigma']
    for d in [cutoff,cutoff+.1,-cutoff]:
        result = plane.evaluate(d,[.4,.2])
        assert np.all(result.value == 0.) and np.all(result.gradient == 0.) and np.all(result.hessian == 0.)
    with pytest.raises(ValueError,match='nonzero'):
        plane.evaluate(0.,[.4,.2])
    p,_,lattice = reference
    model = SW111BesselInterface(p,lattice,cut_kind='glide')
    with pytest.raises(ValueError,match='ordered'):
        model.evaluate([-model.gap,0.,0.])
