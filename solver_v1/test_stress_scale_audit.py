import numpy as np
import pytest

from .interface_static_scenarios import InterfaceUnits, resolved_tensor_tractions, resolved_uniaxial_tractions
from .stress_scale_audit import normal_relaxed_registry


def test_mpa_and_gpa_explicit_roundtrip_and_independent_virtual_work():
    units = InterfaceUnits(2.8e-10,7.1e-20)
    for stress in (4.,50.,2000.):
        force = units.traction_mpa_to_force(stress)
        assert force == pytest.approx(stress*1e6*7.1e-20*2.8e-10/1.602176634e-19)
        assert units.force_to_traction_mpa(force) == pytest.approx(stress)
        assert force == pytest.approx(units.traction_to_force(stress/1000))


def test_tensor_normal_and_two_shears_with_rotation_covariance():
    n,m = np.array([0.,0.,1.]),np.array([1.,0.,0.])
    stress = np.array([[12.,4.,15.],[4.,8.,-7.],[15.,-7.,30.]])
    np.testing.assert_allclose(resolved_tensor_tractions(stress,n,m),[30.,15.,-7.])
    # Proper rotation cycles all three Cartesian axes; input units unchanged.
    R = np.array([[0,1,0],[0,0,1],[1,0,0]])
    np.testing.assert_allclose(resolved_tensor_tractions(R@stress@R.T,R@n,R@m),[30.,15.,-7.])
    e = (n+m)/np.sqrt(2)
    np.testing.assert_allclose(resolved_tensor_tractions(50*np.outer(e,e),n,m)[:2],
        resolved_uniaxial_tractions(50,e,n,m))
    with pytest.raises(ValueError):
        resolved_tensor_tractions(np.triu(stress),n,m)


def test_normal_relaxation_schur_derivative_and_first_spinodal():
    # Analytic potential .5*k*(a-h-c*sin(s))^2 + v*(1-cos(s)).
    h,k,c,v = 1.,10.,.1,2.
    def packed(a,s):
        x=a-h-c*np.sin(s)
        return np.array([.5*k*x*x+v*(1-np.cos(s)), k*x,
            -k*x*c*np.cos(s)+v*np.sin(s), k, -k*c*np.cos(s),
            k*c*c*np.cos(s)**2+k*x*c*np.sin(s)+v*np.cos(s)])
    units=InterfaceUnits(1e-10,1e-20)
    for samples in (32,64):
        rows,peak=normal_relaxed_registry(packed,h=h,period=2*np.pi,units=units,samples=samples)
        assert peak["s_reduced"] == pytest.approx(np.pi/2,abs=1e-10)
        assert peak["shear_mpa"] == pytest.approx(float(units.force_to_traction_mpa(v)))
        assert peak["a_reduced"] == pytest.approx(1.1)
        for row in rows:
            s=row["s_reduced"]
            assert row["schur_curvature"] == pytest.approx(v*np.cos(s),abs=2e-10)
            assert row["relaxed_energy_ev"] == pytest.approx(v*(1-np.cos(s)),abs=2e-12)
