import numpy as np
import pytest
from .phase_resolution import null_lockin_windows, loss_envelope, harmonic_content


def test_null_windows_use_all_sites_without_claiming_independence():
    t=10+np.arange(1000)*.01
    q=np.column_stack([.4*np.sin(2*np.pi*t),-.2*np.sin(2*np.pi*t)])
    rows=null_lockin_windows(t,q,1.,2.)
    assert len(rows)==10
    np.testing.assert_allclose([x['imag_A'] for x in rows],[-.4,.2]*5,atol=1e-14)
    np.testing.assert_allclose([x['real_A'] for x in rows],0,atol=1e-14)


def test_null_floor_sign_force_and_no_sqrt_two_repair():
    for force in (-2.,2.):
        r=loss_envelope(-.1,[.3,-.4],force)
        assert r['observed_null_envelope_A2_eV']==.2
        assert not r['loss_exceeds_observed_null']
        assert not r['confidence_interval']
        assert not r['production_clock_calibrated']
        assert loss_envelope(-.3,[.3,-.4],force)['loss_exceeds_observed_null']
        assert not loss_envelope(.3,[.3,-.4],force)['loss_exceeds_observed_null']


def test_harmonic_content_does_not_hide_nonlinearity():
    t=np.arange(1234)*.01;w=2*np.pi*.7*t
    q=.2+2*np.cos(w)+.3*np.sin(2*w)+.04*np.cos(3*w)
    r=harmonic_content(t,q,.7)
    np.testing.assert_allclose(r['amplitude_A'],[2,.3,.04],atol=1e-14)
    assert r['residual_rms_A']<1e-14


@pytest.mark.parametrize('duration',[0,-1,np.nan,100])
def test_null_window_invalid_duration(duration):
    with pytest.raises(ValueError):null_lockin_windows(np.arange(100)*.01,np.ones((100,2)),1,duration)
