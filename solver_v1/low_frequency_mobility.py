"""Finite-band mobility diagnostic, not an automatic physical clock.

Two-sided angular-Fourier covariance spectrum sampled at cycle frequencies:
S(0)=2 integral_0^infinity C(t)dt for reversible equilibrium coordinates.
Positive spectral estimates alone do not establish a zero-frequency plateau.
"""
import numpy as np
from scipy.signal.windows import dpss


def multitaper_spectrum(coordinates, dt_seconds, *, nw=3., tapers=2):
    """Return two-sided PSD, without doubling the positive-frequency bins.

    Input (time, equivalent sites, coordinates). Sites are averaged, not
    claimed statistically independent. Remove the temporal mean at each site.
    """
    q = np.asarray(coordinates, dtype=float)
    if (q.ndim != 3 or len(q) < 16 or min(q.shape[1:]) < 1 or not np.all(np.isfinite(q))
            or not np.isfinite(dt_seconds) or dt_seconds <= 0
            or not 0 < nw < len(q)/2 or not 1 <= tapers <= int(2*nw)):
        raise ValueError('finite coordinates, positive sampling and valid tapers required')
    q = q-q.mean(axis=0, keepdims=True)
    windows, concentration = dpss(len(q), nw, Kmax=tapers, sym=False,
                                  return_ratios=True)
    spectrum = np.zeros((len(q)//2+1, q.shape[2], q.shape[2]), complex)
    weighted_covariance = np.zeros((q.shape[2], q.shape[2]))
    for window in windows:
        z = q*window[:, None, None]
        fourier = np.fft.rfft(z, axis=0)
        norm = np.dot(window, window)*q.shape[1]
        spectrum += dt_seconds*np.einsum('fpi,fpj->fij', fourier, fourier.conj())/norm
        weighted_covariance += np.einsum('tpi,tpj->ij', z, z)/norm
    return dict(frequency_hz=np.fft.rfftfreq(len(q), dt_seconds),
                spectrum_m2_seconds=spectrum/tapers,
                taper_weighted_covariance_m2=weighted_covariance/tapers,
                covariance_m2=np.einsum('tpi,tpj->ij', q, q)/(len(q)*q.shape[1]),
                half_bandwidth_hz=nw/(len(q)*dt_seconds),
                taper_concentration=concentration)


def band_integral_proxy(estimate, lower_hz, upper_hz):
    """Measured finite-band S/2 proxy; require band above taper bandwidth.

    This is NOT the demonstrated zero-frequency integral. All bandwidth and
    record-length variations must be retained in the associated study.
    """
    f = estimate['frequency_hz']
    if not estimate['half_bandwidth_hz'] < lower_hz < upper_hz <= f[-1]:
        raise ValueError('band must be resolved above the taper half bandwidth')
    mask = (f >= lower_hz) & (f < upper_hz)
    if mask.sum() < 3:
        raise ValueError('at least three frequency bins required')
    spectral = estimate['spectrum_m2_seconds'][mask].mean(axis=0)
    return dict(integral_proxy_m2_seconds=spectral.real/2,
                imaginary_spectrum_m2_seconds=spectral.imag,
                bins=int(mask.sum()), actual_band_hz=[float(f[mask][0]),float(f[mask][-1])],
                zero_frequency_limit_certified=False)


def undamped_leakage_bound(variance_m2, angular_frequency_rad_s, frames,
                          dt_seconds, lower_hz, upper_hz, *, nw=3., tapers=2):
    """Worst-phase finite-window leakage of ONE undamped sinusoidal mode.

    True infinite-record low-frequency power is zero. Maximize over phase
    analytically (2x2 eigenvalue), not random phase sampling. Variance denotes
    its infinite-record variance. This controls window leakage, not unknown
    anharmonic spectral content or a rigorous total statistical error.
    """
    if (not np.isfinite(variance_m2+angular_frequency_rad_s)
            or variance_m2 < 0 or angular_frequency_rad_s <= 0):
        raise ValueError('nonnegative variance and positive frequency required')
    t = np.arange(frames)*dt_seconds
    basis = np.stack([np.cos(angular_frequency_rad_s*t),
                      np.sin(angular_frequency_rad_s*t)], axis=1)
    result = multitaper_spectrum(basis[:,None,:],dt_seconds,nw=nw,tapers=tapers)
    proxy = band_integral_proxy(result,lower_hz,upper_hz)
    return float(2*variance_m2*np.linalg.eigvalsh(proxy['integral_proxy_m2_seconds'])[-1])


def heldout_spectral_prediction(training, heldout, target, conditioners):
    """Linear spectral predictor fit on one time half, evaluated on another.

    Explained spectral fraction is correlation, NOT causal identification.
    Held-out reduction may be negative; do not clip a failed predictor.
    """
    train=np.asarray(training,complex);test=np.asarray(heldout,complex)
    if (train.ndim!=2 or train.shape[0]!=train.shape[1] or test.shape!=train.shape
            or not np.all(np.isfinite(train+test))
            or not np.allclose(train,train.conj().T,rtol=1e-12,atol=0)
            or not np.allclose(test,test.conj().T,rtol=1e-12,atol=0)):
        raise ValueError('matching Hermitian spectral matrices required')
    indices=np.asarray(conditioners,int)
    if target in indices or len(np.unique(indices))!=len(indices) or len(indices)==0:
        raise ValueError('target and conditioning channels must be distinct')
    ee=train[np.ix_(indices,indices)]
    if np.linalg.eigvalsh(ee).min()<=0 or train[target,target].real<=0 or test[target,target].real<=0:
        raise ValueError('resolved positive channel spectra required')
    weight=np.linalg.solve(ee.T,train[target,indices])
    vector=np.zeros(train.shape[0],complex);vector[target]=1;vector[indices]=-weight
    residual=float((vector@test@vector.conj()).real)
    trained=float((vector@train@vector.conj()).real)
    return dict(weight_real=weight.real,weight_imag=weight.imag,
        training_fraction_explained=float(1-trained/train[target,target].real),
        heldout_fraction_explained=float(1-residual/test[target,target].real),
        heldout_residual_spectrum=residual,causal_identification=False)
