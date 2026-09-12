import numpy as np
import pytest
from .collective_forcing import (conjugate_atom_weights,harmonic_response,
                                 susceptibility_from_covariance)


def test_conjugate_force_virtual_work_and_momentum():
    labels=np.array([0,0,1,1,1,2]);w=conjugate_atom_weights(labels)
    direction=np.array([1.,-1.,0.])/np.sqrt(2)
    u=np.arange(18).reshape(6,3)*.002;F=1.7
    q=(u[labels==1].mean(axis=0)-u[labels==0].mean(axis=0))@direction
    assert np.sum(F*w[:,None]*direction*u)==pytest.approx(F*q)
    np.testing.assert_allclose((F*w[:,None]*direction).sum(axis=0),0.,atol=1e-15)
    with pytest.raises(ValueError):conjugate_atom_weights(labels,0,0)


def test_lockin_signed_force_and_phase():
    t=np.arange(10000)*.01;frequency=.2;chi=2.-.3j
    for force in (-.5,1.):
        q=.1+np.real(force*chi*np.exp(2j*np.pi*frequency*t))
        fit=harmonic_response(t,q,frequency,force)
        assert fit['real_A2_eV']==pytest.approx(chi.real)
        assert fit['imag_A2_eV']==pytest.approx(chi.imag)
        assert fit['residual_rms_A']<1e-13


def test_fdt_dho_susceptibility_and_endpoint():
    from .mode_kinetic_calibration import damped_position_correlation
    t=np.arange(0,100,.001);g=.3;wd=2.;omega=.4
    kT=.8;C0=.2;H=kT/C0;mass=H/(g*g+wd*wd)
    c=C0*damped_position_correlation(t,g,wd)
    chi,truncated=susceptibility_from_covariance(t,c,omega/(2*np.pi),kT)
    exact=1/(H-mass*omega**2+2j*g*mass*omega)
    assert abs(chi-exact)<1e-8
    assert abs(chi-truncated)<1e-12
    short,cut=susceptibility_from_covariance(t[:1000],c[:1000],omega/(2*np.pi),kT)
    assert short-cut==pytest.approx(c[999]*np.exp(-1j*omega*t[999])/kT)
