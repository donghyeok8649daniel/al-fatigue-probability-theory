import numpy as np
import pytest
from scipy.integrate import quad
from .mode_kinetic_calibration import (damped_position_correlation, mode_correlations,
    fit_damped_correlation,formal_modal_mobility)


def test_reference_md_sample_protocol_without_optional_engine():
    from .run_reference_thermostat_md import sample_protocol
    assert sample_protocol(.005,.025,25.) == (5,1000)
    assert sample_protocol(.0025,.025,25.) == (10,1000)
    for bad in ((.003,.025,25.),(.005,.025,25.001),(0.,.025,25.)):
        with pytest.raises(ValueError):
            sample_protocol(*bad)


def test_plane_velocity_energy_normalization_and_streaming_removal():
    from .run_reference_thermostat_md import plane_velocity_statistics
    v=np.arange(36,dtype=float).reshape(12,3)/10
    labels=np.repeat(np.arange(3),4)
    mean,raw,internal=plane_velocity_statistics(v,labels,3,mass_amu=26.98)
    factor=.5*26.98*1.66053906660e-27*1e4/1.602176634e-19
    assert raw.mean()==pytest.approx(factor*np.mean(np.sum(v*v,axis=1)))
    np.testing.assert_allclose(raw-internal,factor*np.sum(mean*mean,axis=1),rtol=1e-13)
    shifted=plane_velocity_statistics(v+2.,labels,3,mass_amu=26.98)[2]
    np.testing.assert_allclose(shifted,internal,atol=1e-16)
    np.testing.assert_allclose(plane_velocity_statistics(v,labels,3,mass_amu=26.982)[1],raw*26.982/26.98)


def test_gradient_drag_symbol_and_heldout_prediction():
    from .mode_kinetic_calibration import periodic_difference_symbol, gradient_drag_hypothesis
    N=12;B=np.roll(np.eye(N),-1,axis=0)-np.eye(N)
    w=periodic_difference_symbol(N)
    np.testing.assert_allclose(np.linalg.eigvalsh(B@B.T)[1:],np.sort(w),atol=3e-15)
    np.testing.assert_allclose(B@np.ones(N),0.,atol=1e-15)
    k=np.arange(1,7);fit=gradient_drag_hypothesis(k,.17*w[:6],planes=N,fit_modes=[1,2,3])
    assert fit['coefficient_per_ps']==pytest.approx(.17)
    assert fit['heldout_rmse_per_ps']<1e-15
    assert not fit['production_calibration_available']
    with pytest.raises(ValueError):gradient_drag_hypothesis(k,.17*w[:6],planes=N,fit_modes=k)


def test_experimental_linewidth_provenance_and_frequency_convention():
    import json
    from pathlib import Path
    file=Path(__file__).resolve().parents[1]/'results/modal_kinetic_calibration_v25/experimental_linewidth/reference.json'
    data=json.loads(file.read_bytes())
    assert data['doi']=='10.1103/PhysRevLett.123.235501'
    assert data['external_validation_only'] and not data['fit_target']
    assert data['production_coordinate_mobility'] is None
    assert len(data['records'])==4
    for row in data['records']:
        assert row['damping_g_per_ps']==pytest.approx(np.pi*row['linewidth_THz_from_figure'])
        assert row['envelope_time_ps']*row['damping_g_per_ps']==pytest.approx(1.)
        assert row['experimental_uncertainty_not_reported_by_this_extraction']
    ambient=data['records'][0]
    assert 290<ambient['temperature_K_from_figure']<300
    assert .26<ambient['linewidth_THz_from_figure']<.28


def test_scalar_mixture_mobility_normalization_and_units():
    from .run_modal_mobility_projection import scalar_mixture_mobility
    C=np.arange(1,7)*1e-24;tau=np.ones(6)*2e-12
    weights=np.array([2,2,2,2,2,1])
    r=scalar_mixture_mobility(C,tau,multiplicities=weights,planes=12,temperature_K=300.)
    variance=weights@C/12
    assert r['variance_m2']==pytest.approx(variance,abs=1e-40)
    assert r['mobility_m2_per_J_s']==pytest.approx(variance/(1.380649e-23*300*2e-12))
    assert r['mobility_m2_per_J_s']*r['friction_J_s_per_m2']==pytest.approx(1.)
    # Changing coordinate covariance normalization cannot silently leave M fixed.
    scaled=scalar_mixture_mobility(C/576,tau,multiplicities=weights,planes=12,temperature_K=300.)
    assert scaled['mobility_m2_per_J_s']==pytest.approx(r['mobility_m2_per_J_s']/576)
    with pytest.raises(ValueError):
        scalar_mixture_mobility(C,-tau,multiplicities=weights,planes=12,temperature_K=300.)
    # Energy normalization must preserve BOTH drift and diffusion. It is not
    # permission to put a collective-plane coefficient in a local-cell PDE.
    Np=576;mobility=1.7e10;thermal=1.380649e-23*300
    density=.3;gradient=2e-10;density_gradient=5e10
    flux=-mobility*(density*gradient+thermal*density_gradient)
    normalized=-Np*mobility*(density*gradient/Np+thermal/Np*density_gradient)
    assert normalized==pytest.approx(flux)
    wrong=-Np*mobility*(density*gradient/Np+thermal*density_gradient)
    assert abs(wrong-flux)>abs(flux)


def test_damped_position_identity_and_integral():
    g=.3;w=4.;t=np.arange(0,15,.001)
    c=damped_position_correlation(t,g,w)
    assert c[0]==1
    assert abs((c[2]-2*c[1]+c[0])/.001**2+(g*g+w*w))<.02
    integral,error=quad(lambda x:float(damped_position_correlation(x,g,w)),0,160,epsabs=1e-10,limit=1500)
    assert abs(integral-2*g/(g*g+w*w))<1e-9
    assert error<1e-8


def test_modal_fft_matches_independent_direct_lag_sum():
    t=np.arange(100)[:,None,None];p=np.arange(6)[None,:,None]
    q=1e-12*np.cos(.23*t+2*np.pi*p/6+np.arange(3)[None,None,:])
    q+=.3e-12*np.sin(.37*t+4*np.pi*p/6)
    actual=mode_correlations(q,12)
    z=np.fft.fft(q-q.mean(axis=0,keepdims=True),axis=1,norm='ortho')
    for lag in (0,3,12):
        direct=np.mean(z[lag:]*z[:len(z)-lag].conj(),axis=0).real
        np.testing.assert_allclose(actual[lag],direct,rtol=2e-13,atol=1e-51)
    np.testing.assert_allclose(actual[:,1],actual[:,-1],rtol=1e-13,atol=1e-51)


def test_known_deterministic_dho_fit_and_seconds_units():
    t=np.arange(256)*.025
    c=damped_position_correlation(t,.17,5.)
    fit=fit_damped_correlation(t,c,4.9)
    assert fit['optimizer_success']
    np.testing.assert_allclose([fit['damping_per_ps'],fit['damped_angular_frequency_rad_ps']],[.17,5.],rtol=1e-8)
    assert fit['training_rmse']<1e-9
    assert fit['overdamped_training_rmse']>.2
    M=formal_modal_mobility(2e-24,fit['formal_integral_time_ps']*1e-12,300.)
    np.testing.assert_allclose(M,2e-24*(25+.17**2)/(.34e-12*1.380649e-23*300.))
    assert not fit['production_calibration_available']


def test_invalid_inputs_not_repaired():
    with pytest.raises(ValueError):formal_modal_mobility(1.,-1.,300.)
    with pytest.raises(ValueError):damped_position_correlation([0.,-1.],1.,2.)
    with pytest.raises(ValueError):mode_correlations(np.zeros((20,3,3)),20)
