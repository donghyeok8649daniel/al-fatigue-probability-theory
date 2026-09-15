"""Finite-record noise and heating diagnostics; no production time calibration.

Convention: F=F0*cos(omega*t), loss=-Im(chi), q in A, t in ps.
The stationary classical FDT and long-record approximations are conditional:
S_q(omega)=2*kBT*loss/omega (two-sided spectrum),
Var(loss_hat) ~ 2*S_q/(F0**2*T) = 4*kBT*loss/(omega*F0**2*T).
Thus expected input work = 2*kBT*SNR**2 for one sine quadrature.
This is neither a confidence interval nor an MD stationarity certificate.
"""
import numpy as np
from scipy.signal import fftconvolve

KB_EV_K = 1.380649e-23 / 1.602176634e-19


def cycle_windows(time_ps, frequency_per_ps, duration_ps, exclude_ps=0.):
    """Partition the selected record into equal complete-cycle intervals.

    All samples after exclude_ps are used. Neighbouring windows share one
    endpoint, as trapezoidal integrals do. They are not independent replicates.
    Absolute time retains the applied force's phase origin.
    """
    t = np.asarray(time_ps, float)
    values = np.asarray([frequency_per_ps, duration_ps, exclude_ps], float)
    if (t.ndim != 1 or len(t) < 8 or not np.all(np.isfinite(t))
            or not np.all(np.isfinite(values)) or np.any(np.diff(t) <= 0)
            or frequency_per_ps <= 0 or duration_ps <= 0):
        raise ValueError('finite ordered record, positive frequency and duration required')
    dt = t[1]-t[0]
    if not np.allclose(np.diff(t), dt, rtol=1e-8, atol=1e-11):
        raise ValueError('uniform sample spacing required')
    if frequency_per_ps*dt >= .5:
        raise ValueError('forcing frequency must be below sampling Nyquist')
    cycles = duration_ps*frequency_per_ps
    if cycles < 1 or not np.isclose(cycles, round(cycles), atol=1e-8, rtol=0):
        raise ValueError('integer complete cycles required')
    first = int(round((exclude_ps-t[0])/dt))
    width = int(round(duration_ps/dt))
    if (first < 0 or first >= len(t)-1 or width < 7
            or not np.isclose(t[first], exclude_ps, atol=1e-8, rtol=0)
            or not np.isclose(width*dt, duration_ps, atol=1e-8, rtol=0)
            or (len(t)-1-first) % width):
        raise ValueError('complete equal windows with sampled endpoints required')
    return [(i, i+width) for i in range(first, len(t)-1, width)]


def sine_weights(time_ps, frequency_per_ps, force_eV_A):
    """Trapezoidal linear weights for the signed loss estimator, without fit."""
    t = np.asarray(time_ps, float)
    if t.ndim != 1 or len(t) < 8:
        raise ValueError('finite complete-cycle record required')
    if not np.isfinite(force_eV_A) or force_eV_A == 0:
        raise ValueError('nonzero finite signed force required')
    # Validation also checks complete cycles and the sampling frequency.
    cycle_windows(t, frequency_per_ps, t[-1]-t[0], t[0])
    weights = 2*(t[1]-t[0])*np.sin(2*np.pi*frequency_per_ps*t)
    weights[[0, -1]] *= .5
    return weights/(force_eV_A*(t[-1]-t[0]))


def stationary_quadrature_variance(weights, covariance_lags):
    """Exact w^T C w for a supplied stationary covariance C_ij=c[|i-j|].

    FFT contraction is checked against an independent dense matrix in tests.
    Covariance lags must come from a valid stationary model; this routine
    does not infer them from MD or certify their statistical uncertainty.
    A negative contraction is an error, never silently clipped.
    """
    w, c = map(lambda x: np.asarray(x, float), (weights, covariance_lags))
    if (w.ndim != 1 or len(w) < 2 or c.shape != w.shape
            or not np.all(np.isfinite(w)) or not np.all(np.isfinite(c))
            or c[0] < 0):
        raise ValueError('matching finite weights and covariance lags required')
    correlation = fftconvolve(w, w[::-1], mode='full')[len(w)-1:]
    result = float(c[0]*correlation[0]+2*np.dot(c[1:], correlation[1:]))
    if result < 0:
        raise ValueError('negative covariance contraction; do not repair')
    return result


def conditional_precision_budget(force_eV_A, frequency_per_ps, loss_A2_eV,
                                 temperature_K, duration_ps, *, desired_snr=3.):
    """Conditional single-record standard-deviation SNR, not an error envelope.

    This long-time equilibrium prediction assumes resolved spectral density,
    linear response, constant temperature and negligible spectral leakage.
    A finite-band FDT proxy does not prove these assumptions.
    """
    values = np.asarray([force_eV_A, frequency_per_ps, loss_A2_eV,
                         temperature_K, duration_ps, desired_snr], float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError('positive finite conditional planning inputs required')
    theta = KB_EV_K*temperature_K
    omega = 2*np.pi*frequency_per_ps
    power = .5*omega*force_eV_A**2*loss_A2_eV
    variance = 4*theta*loss_A2_eV/(omega*force_eV_A**2*duration_ps)
    required = 4*theta*desired_snr**2/(omega*force_eV_A**2*loss_A2_eV)
    return dict(predicted_loss_std_A2_eV=float(np.sqrt(variance)),
                predicted_snr=float(loss_A2_eV/np.sqrt(variance)),
                expected_power_eV_ps=float(power),
                expected_work_eV=float(power*duration_ps),
                desired_snr=float(desired_snr),
                conditional_duration_ps=float(required),
                conditional_work_at_snr_eV=float(power*required),
                confidence_interval=False, stationarity_certified=False,
                production_clock_calibrated=False)


def thermal_work_summary(time_ps, thermo, cumulative_work_eV=None):
    """Actual temperature trend and U/work bookkeeping on the supplied interval.

    Thermo columns: temperature_K, potential_eV, kinetic_eV, total_eV, pressure.
    Work is external input; U-work is integration error, not dissipation.
    The temperature includes coherent kinetic motion in the original MD.
    """
    t, h = np.asarray(time_ps, float), np.asarray(thermo, float)
    if (t.ndim != 1 or len(t) < 2 or h.shape != (len(t), 5)
            or not np.all(np.isfinite(t)) or not np.all(np.isfinite(h))
            or np.any(np.diff(t) <= 0) or np.any(h[:, 0] <= 0)):
        raise ValueError('matching finite thermo record with positive temperatures required')
    centered = t-t.mean()
    rate = float(np.dot(centered, h[:, 0]-h[:, 0].mean())/np.dot(centered, centered))
    energy = h[:, 3]-h[0, 3]
    result = dict(mean_temperature_K=float(np.mean(h[:, 0])),
                  temperature_linear_trend_K_ns=rate*1000,
                  internal_energy_change_eV=float(energy[-1]),
                  internal_energy_range_eV=float(np.ptp(energy)),
                  work_eV=None, max_work_residual_eV=None,
                  endpoint_work_residual_eV=None)
    if cumulative_work_eV is not None:
        w = np.asarray(cumulative_work_eV, float)
        if w.shape != t.shape or not np.all(np.isfinite(w)):
            raise ValueError('matching finite native-work record required')
        work = w-w[0]
        residual = energy-work
        result.update(work_eV=float(work[-1]),
                      max_work_residual_eV=float(np.max(abs(residual))),
                      endpoint_work_residual_eV=float(residual[-1]))
    return result
