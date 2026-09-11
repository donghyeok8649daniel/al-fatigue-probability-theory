import numpy as np
import pytest
from .low_frequency_mobility import (multitaper_spectrum,band_integral_proxy,
                                      undamped_leakage_bound)


def test_two_sided_parseval_and_psd():
    n=2048;t=np.arange(n)*.025
    q=np.stack([np.cos(2*np.pi*1.2*t),np.sin(2*np.pi*1.2*t)+.3*np.cos(2*np.pi*.7*t)],axis=1)
    e=multitaper_spectrum(q[:,None,:],.025)
    s=e['spectrum_m2_seconds']
    integrated=(s[0].real+s[-1].real+2*s[1:-1].real.sum(axis=0))/(n*.025)
    np.testing.assert_allclose(integrated,e['taper_weighted_covariance_m2'],atol=1e-14)
    np.testing.assert_allclose(s,s.conj().transpose(0,2,1),atol=1e-14)
    assert np.linalg.eigvalsh(s).min()>-1e-13


def test_phase_leakage_bound_and_no_certification():
    n=4096;dt=.025;w=2*np.pi*1.23;variance=.7
    bound=undamped_leakage_bound(variance,w,n,dt,.05,.15)
    for phase in (0.,.3,1.1,2.5):
        q=np.sqrt(2*variance)*np.cos(w*np.arange(n)*dt+phase)
        e=multitaper_spectrum(q[:,None,None],dt)
        p=band_integral_proxy(e,.05,.15)
        assert p['integral_proxy_m2_seconds'][0,0] <= bound*(1+1e-9)
        assert not p['zero_frequency_limit_certified']
    assert bound>0
    with pytest.raises(ValueError):
        band_integral_proxy(e,.001,.02)


def test_spectral_predictor_heldout_failure_is_not_clipped():
    from .low_frequency_mobility import heldout_spectral_prediction
    training=np.array([[5.,2.],[2.,1.]])
    result=heldout_spectral_prediction(training,training,0,[1])
    assert result['heldout_fraction_explained']==pytest.approx(.8)
    test=np.array([[5.,-2.],[-2.,1.]])
    result=heldout_spectral_prediction(training,test,0,[1])
    assert result['heldout_fraction_explained']<0
    assert not result['causal_identification']
