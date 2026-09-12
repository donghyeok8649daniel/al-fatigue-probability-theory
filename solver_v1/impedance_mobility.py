"""Research-only scalar conjugate-response mobility; never a production clock.

The inverse of a scalar susceptibility is an effective projected impedance,
NOT an entry of the inverse full response matrix. Input chi uses A^2/eV,
frequency cycles/ps, so drag uses eV ps/A^2. Noise disks are empirical
sensitivity sets, not confidence intervals or rigorous statistical bounds.
"""
import numpy as np

MOBILITY_SI_PER_A2_EV_PS = 1e-20 / (1.602176634e-19 * 1e-12)


def scalar_impedance(chi, frequency_per_ps, *, static_chi=None):
    chi = complex(chi)
    if (not np.isfinite(chi) or abs(chi) == 0
            or not np.isfinite(frequency_per_ps) or frequency_per_ps <= 0):
        raise ValueError('finite nonzero susceptibility and positive frequency required')
    omega = 2*np.pi*frequency_per_ps
    z = 1/chi
    drag = z.imag/omega
    result = dict(storage_impedance_eV_A2=z.real,
        drag_eV_ps_A2=drag,
        mobility_m2_per_J_s=MOBILITY_SI_PER_A2_EV_PS/drag if drag > 0 else None,
        passive_point=drag > 0, production_clock_calibrated=False)
    if static_chi is not None:
        if not np.isfinite(static_chi) or static_chi <= 0:
            raise ValueError('positive static susceptibility required')
        result['dispersion_coefficient_eV_ps2_A2'] = (1/static_chi-z.real)/omega**2
        result['storage_to_static_ratio'] = chi.real/static_chi
    return result


def inverse_disk_bounds(chi, radius, frequency_per_ps):
    """Exact image of |chi-chi_measured|<=radius under z=1/chi.

If radius >= |chi| the set reaches the pole; do not return finite bounds.
If positive drag is feasible but zero is included, mobility has NO finite
upper bound. None denotes unavailable/unbounded, never a silently clipped 0.
"""
    scalar_impedance(chi, frequency_per_ps)
    chi = complex(chi)
    if not np.isfinite(radius) or radius < 0:
        raise ValueError('finite nonnegative radius required')
    flags = dict(confidence_interval=False, production_clock_calibrated=False)
    if radius >= abs(chi):
        return dict(**flags, inverse_set_bounded=False,
            drag_lower_eV_ps_A2=None, drag_upper_eV_ps_A2=None,
            mobility_lower_m2_per_J_s=None, mobility_upper_m2_per_J_s=None,
            positive_drag_resolved=False, passive_drag_feasible=chi.imag-radius < 0)
    denominator = abs(chi)**2-radius**2
    center = chi.conjugate()/denominator
    inverse_radius = radius/denominator
    omega = 2*np.pi*frequency_per_ps
    lower, upper = (center.imag-inverse_radius)/omega, (center.imag+inverse_radius)/omega
    return dict(**flags, inverse_set_bounded=True,
        drag_lower_eV_ps_A2=lower, drag_upper_eV_ps_A2=upper,
        mobility_lower_m2_per_J_s=MOBILITY_SI_PER_A2_EV_PS/upper if upper > 0 else None,
        mobility_upper_m2_per_J_s=MOBILITY_SI_PER_A2_EV_PS/lower if lower > 0 else None,
        positive_drag_resolved=lower > 0, passive_drag_feasible=upper > 0)


def scalar_low_band_drag(covariance_m2, integral_proxy_m2_s, temperature_K):
    """kBT*(S_band/2)/C0^2, under a LOW-FREQUENCY scalar response hypothesis.

This does not replace Re chi(f) with measured finite-frequency storage and
does not certify S_band = S(0). SI drag is J s/m^2.
"""
    values = np.array([covariance_m2, integral_proxy_m2_s, temperature_K], float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError('positive covariance, band integral and temperature required')
    return 1.380649e-23*temperature_K*integral_proxy_m2_s/covariance_m2**2


def fit_constant_drag(drag, sensitivity_scale):
    """Deterministic 1-parameter WLS, scales are documented discrepancies.

No statistical covariance/CI is invented from these sensitivity scales.
Nonpositive unconstrained solution is rejected, not repaired with abs/clip.
"""
    y, scale = np.asarray(drag, float), np.asarray(sensitivity_scale, float)
    if (y.ndim != 1 or len(y) < 2 or scale.shape != y.shape
            or not np.all(np.isfinite(y+scale)) or np.any(scale <= 0)):
        raise ValueError('matching finite vectors and positive scales required')
    reference = float(np.max(np.abs(y)))
    if reference == 0: raise ValueError('nonzero drag scale required')
    weights = (scale.min()/scale)**2
    gamma = float(np.dot(weights, y)/weights.sum())
    residual = (gamma-y)/scale
    return dict(drag=gamma, admissible=gamma > 0,
        normalized_residuals=residual, objective=float(residual@residual),
        scaled_jacobian_singular_value=float(np.linalg.norm(reference/scale)),
        parameter_rank=1, confidence_interval=False, production_clock_calibrated=False)
