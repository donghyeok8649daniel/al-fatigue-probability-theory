"""Finite-record phase diagnostics. Empirical envelopes are NOT confidence bounds.

Reference-coordinate units: ps, Angstrom, generalized eV/Angstrom.
No production mobility, energy, temperature or clock is inferred here.
"""
import numpy as np
from .collective_forcing import harmonic_response


def null_lockin_windows(time_ps, coordinates_A, frequency_per_ps, duration_ps):
    """All nonoverlapping windows and represented sites; do not count independent.

    Return displacement quadratures using a UNIT virtual force. No physical
    force was applied. Retain absolute timestamps, including noninteger cycles.
    """
    t = np.asarray(time_ps, float)
    q = np.asarray(coordinates_A, float)
    if (t.ndim != 1 or len(t) < 8 or q.ndim != 2 or len(q) != len(t)
            or not np.all(np.isfinite(t)) or not np.all(np.isfinite(q))
            or np.any(np.diff(t) <= 0) or not np.isfinite(duration_ps)
            or duration_ps <= 0):
        raise ValueError('ordered times, sites and positive window required')
    dt = np.median(np.diff(t))
    if not np.allclose(np.diff(t), dt, rtol=1e-7, atol=1e-10):
        raise ValueError('uniform sampling required')
    width = int(round(duration_ps / dt))
    if width < 8 or not np.isclose(width*dt, duration_ps) or width > len(t):
        raise ValueError('resolved complete window required')
    rows = []
    for block in range(len(t)//width):
        sl = slice(block*width, (block+1)*width)
        for site in range(q.shape[1]):
            result = harmonic_response(t[sl], q[sl, site], frequency_per_ps, 1.)
            rows.append(dict(block=block, site=site,
                             real_A=result['real_A2_eV'], imag_A=result['imag_A2_eV']))
    return rows


def loss_envelope(response_imag, null_quadratures_A, force_eV_A):
    """Conservative observed-null envelope, not a probabilistic guarantee.

    For sign-pair averaging, DO NOT divide noise by sqrt(2): the common
    restart is not an independent replicate. This checks loss sign only,
    not agreement with FDT, discretization, weak drive or zero-frequency M.
    """
    values = np.asarray(null_quadratures_A, float)
    if (values.size == 0 or not np.all(np.isfinite(values))
            or not np.isfinite(response_imag+force_eV_A) or force_eV_A == 0):
        raise ValueError('finite observed quadratures and nonzero force required')
    floor = float(np.max(np.abs(values))/abs(force_eV_A))
    return dict(observed_null_envelope_A2_eV=floor,
                loss_exceeds_observed_null=bool(-response_imag > floor),
                loss_to_null_ratio=float(-response_imag/floor) if floor else None,
                confidence_interval=False, production_clock_calibrated=False)


def harmonic_content(time_ps, displacement_A, frequency_per_ps):
    """Joint fundamental/second/third harmonic fit; thermal harmonics retained."""
    t = np.asarray(time_ps, float); q = np.asarray(displacement_A, float)
    if (t.ndim != 1 or q.shape != t.shape or len(t) < 14
            or not np.all(np.isfinite(t+q)) or np.any(np.diff(t) <= 0)
            or not np.isfinite(frequency_per_ps) or frequency_per_ps <= 0):
        raise ValueError('finite ordered data and positive frequency required')
    w = 2*np.pi*frequency_per_ps*t
    X = np.column_stack([np.ones_like(t)] +
                        [fn(k*w) for k in (1, 2, 3) for fn in (np.cos, np.sin)])
    c, _, rank, _ = np.linalg.lstsq(X, q, rcond=None)
    if rank != 7: raise ValueError('unresolved harmonic design')
    return dict(amplitude_A=[float(np.hypot(c[2*k-1], c[2*k])) for k in (1, 2, 3)],
                residual_rms_A=float(np.sqrt(np.mean((q-X@c)**2))))
