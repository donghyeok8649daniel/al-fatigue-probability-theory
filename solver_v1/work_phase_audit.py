"""Independent work/phase audit of reference MD; NOT a production clock.

For F(t)=F0 cos(omega*t), W=integral F dq. Integration by parts gives
W-[Fq]_ends=omega F0 integral q sin(omega*t)dt. On an integer number of
cycles, loss=-Im(chi)=2*(W-[Fq]_ends)/(omega F0**2 T).
The endpoint term is thermal/reversible work, not dissipative phase.
"""
import numpy as np
from scipy.integrate import trapezoid


def work_phase(time_ps, q_angstrom, cumulative_work_eV, force_eV_A,
               frequency_per_ps):
    """Return signed loss estimates without clipping or a confidence claim.

The native MD work and saved-coordinate quadrature are independent numerical
measurements of the SAME response, not independent statistical replicates.
Only uniform samples spanning complete cycles are accepted. Absolute time
must retain the phase origin of the applied force, including after slicing.
"""
    t, q, w = (np.asarray(x, float) for x in
               (time_ps, q_angstrom, cumulative_work_eV))
    if (t.ndim != 1 or len(t) < 8 or q.shape != t.shape or w.shape != t.shape
            or not np.all(np.isfinite(t+q+w)) or np.any(np.diff(t) <= 0)
            or not np.isfinite(force_eV_A+frequency_per_ps)
            or force_eV_A == 0 or frequency_per_ps <= 0):
        raise ValueError('finite matching records and nonzero signed force required')
    dt = np.diff(t)
    if not np.allclose(dt, dt[0], rtol=1e-8, atol=1e-11):
        raise ValueError('uniform sample spacing required')
    if frequency_per_ps*dt[0] >= .5:
        raise ValueError('forcing frequency must be below sampling Nyquist')
    duration = t[-1]-t[0]
    cycles = duration*frequency_per_ps
    if not np.isclose(cycles, round(cycles), rtol=0, atol=1e-8) or cycles < 1:
        raise ValueError('integer complete cycles required')
    omega = 2*np.pi*frequency_per_ps
    force = force_eV_A*np.cos(omega*t)
    endpoint = float(force[-1]*q[-1]-force[0]*q[0])
    work = float(w[-1]-w[0])
    normalization = omega*force_eV_A**2*duration/2
    loss_q = float(2*trapezoid(q*np.sin(omega*t), t)/(force_eV_A*duration))
    loss_w = (work-endpoint)/normalization
    return dict(duration_ps=float(duration), cycles=int(round(cycles)),
        work_eV=work, endpoint_work_eV=endpoint,
        corrected_work_eV=work-endpoint,
        loss_native_A2_eV=float(loss_w), loss_coordinate_A2_eV=loss_q,
        loss_uncorrected_work_A2_eV=work/normalization,
        work_coordinate_discrepancy_eV=work-endpoint-normalization*loss_q,
        loss_discrepancy_A2_eV=float(loss_w-loss_q),
        production_clock_calibrated=False, confidence_interval=False)


def duration_sensitivity(duration_ps, noise_scale, predicted_loss, ratio=3.):
    """Conditional T^-1/2 planning, NOT achieved precision or guaranteed CI.

Uses a positive independent predicted loss, never abs(measured negative loss).
Stationarity, linearity and the scaling assumption need new measurement.
"""
    values = np.array([duration_ps, noise_scale, predicted_loss, ratio], float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError('positive finite planning inputs required')
    return float(duration_ps*(ratio*noise_scale/predicted_loss)**2)
