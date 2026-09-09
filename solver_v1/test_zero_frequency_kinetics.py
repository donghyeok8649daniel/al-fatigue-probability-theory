"""Analytic deterministic correlation tests, not fabricated atomistic data."""
import numpy as np
import pytest
from scipy.integrate import quad_vec
from scipy.linalg import expm

from .zero_frequency_kinetics import formal_zero_frequency_response,correlation_integral_audit


def test_noncommuting_overdamped_matrix_integral_recovers_mobility_with_correct_units():
    H=np.array([[3.,.7],[.7,1.]])*10
    M=np.array([[.8,.13],[.13,.2]])*1e-3;kT=4e-21
    C0=kT*np.linalg.inv(H)
    K=np.linalg.solve(M@H,C0)
    result=formal_zero_frequency_response(C0,(K+K.T)/2,kT)
    np.testing.assert_allclose(result['mobility_m2_per_J_second'],M,rtol=1e-13)
    np.testing.assert_allclose(result['friction_J_seconds_per_m2']@M,np.eye(2),atol=1e-14)
    assert not result['production_calibration_available']


def test_inertial_oscillation_can_have_positive_zero_frequency_friction():
    # A nondimensional mathematical control, not mass inserted into our PDE.
    H=np.array([[5.,1.],[1.,2.]]);Gamma=np.array([[.8,.1],[.1,.6]])
    mass=np.diag([1.,1.5]);kT=.4;C0=kT*np.linalg.inv(H)
    operator=np.block([[np.zeros((2,2)),np.eye(2)],[-np.linalg.solve(mass,H),-np.linalg.solve(mass,Gamma)]])
    C=lambda t:expm(operator*t)[:2,:2]@C0
    assert np.linalg.eigvalsh(C(1.5))[0]<0
    integral,error=quad_vec(C,0.,160.,epsabs=1e-13,epsrel=1e-12)
    tail=(-np.linalg.solve(operator,expm(operator*160.)))[:2,:2]@C0
    assert np.linalg.norm(tail,2)<1e-13
    integral+=tail
    assert error<1e-11
    result=formal_zero_frequency_response(C0,(integral+integral.T)/2,kT)
    np.testing.assert_allclose(result['friction_J_seconds_per_m2'],Gamma,rtol=2e-10)


def test_negative_finite_integral_is_not_repaired_or_given_a_mobility():
    with pytest.raises(ValueError,match='do not clip'):
        formal_zero_frequency_response(np.eye(2),np.diag([1.,-1e-12]),1.)
    with pytest.raises(ValueError):formal_zero_frequency_response(np.eye(2),np.eye(2),0.)


def test_undamped_periodic_signal_has_cutoff_dependent_signed_integral():
    t=np.arange(4096.)
    q=np.stack([np.cos(t*2*np.pi/32),np.sin(t*2*np.pi/32)],axis=-1)[:,None]*1e-12
    result=correlation_integral_audit(q,frame_seconds=1e-14,max_lag=126,blocks=4)
    assert np.min(result['integral_eigenvalues_seconds'])<0
    assert np.max(result['integral_eigenvalues_seconds'])>0
    assert not result['production_calibration_available']
    assert np.max(result['half_cutoff_change_seconds'])>2e-14
    with pytest.raises(ValueError,match='two lag windows'):
        correlation_integral_audit(q[:200],frame_seconds=1e-14,max_lag=50,blocks=4)
