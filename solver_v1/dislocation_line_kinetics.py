"""Scoped line-translation calibration and pinned-line linear dynamics.

The state is a line position/bow, NOT a/s, and the mobility is not a
PhysicalTimeCalibration. The elastic functional is derived in
LINE_KINETICS_AND_FATIGUE_VALIDATION.md. No yield or fatigue law is fitted.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded


def fit_velocity_response(stress_pa, velocity_m_s, training):
    """Deterministic zero-intercept v=mu_tau*tau fit with disjoint validation.

    Equal velocity residual scales. Point scatter and leave-one-out range
    are empirical diagnostics, NOT independent confidence intervals.
    No Burgers vector, activation length, M_a or M_s is inferred.
    """
    stress = np.asarray(stress_pa, float)
    velocity = np.asarray(velocity_m_s, float)
    mask = np.asarray(training)
    if (stress.ndim != 1 or velocity.shape != stress.shape or mask.shape != stress.shape
            or mask.dtype != bool or np.any(~np.isfinite([stress, velocity]))
            or np.any(stress <= 0) or np.any(velocity < 0)
            or np.count_nonzero(mask) < 3 or np.count_nonzero(~mask) < 2):
        raise ValueError('positive stress, finite velocity and >=3 fit / >=2 held-out points required')

    def fit(which):
        return float(stress[which] @ velocity[which] / (stress[which] @ stress[which]))

    slope = fit(mask)
    if slope <= 0:
        raise ValueError('positive directional mobility not resolved by these data')
    residual = velocity - slope * stress
    jackknife = []
    for index in np.flatnonzero(mask):
        leave = mask.copy()
        leave[index] = False
        jackknife.append(fit(leave))
    return dict(velocity_per_stress_m_per_Pa_s=slope,
        training_rmse_m_s=float(np.sqrt(np.mean(residual[mask]**2))),
        heldout_rmse_m_s=float(np.sqrt(np.mean(residual[~mask]**2))),
        training_points=int(mask.sum()), heldout_points=int((~mask).sum()),
        prediction_m_s=slope*stress, residual_m_s=residual,
        leave_one_out_slopes_m_per_Pa_s=np.array(jackknife),
        all_point_refit_slope_for_sensitivity_only=fit(np.ones(len(mask), bool)),
        source_mean_ratio_for_sensitivity_only=float(velocity.mean()/stress.mean()),
        independent_specimen_validation=False, production_clock_calibrated=False)


@dataclass(frozen=True)
class PinnedLineKinetics:
    """Explicit physical reference inputs; there are deliberately no defaults.

    B_line [Pa s], T_line [J/m], b [m], L [m]. Their material/geometry
    provenance must be reported by the caller. This object only solves the
    stated line model; it does not certify matching material conditions.
    """
    drag_pa_s: float
    stiffness_j_m: float
    burgers_m: float
    span_m: float

    def __post_init__(self):
        values = [self.drag_pa_s, self.stiffness_j_m, self.burgers_m, self.span_m]
        if np.any(~np.isfinite(values)) or min(values) <= 0:
            raise ValueError('explicit positive finite line drag, stiffness, Burgers vector and span required')

    @property
    def slowest_seconds(self):
        return self.drag_pa_s*self.span_m**2/(self.stiffness_j_m*np.pi**2)

    def static_mean_bow_m(self, shear_pa):
        shear = np.asarray(shear_pa, float)
        if np.any(~np.isfinite(shear)):
            raise ValueError('finite resolved shear required')
        return self.burgers_m*shear*self.span_m**2/(12*self.stiffness_j_m)

    def static_maximum_slope(self, shear_pa):
        shear = np.asarray(shear_pa, float)
        if np.any(~np.isfinite(shear)):
            raise ValueError('finite resolved shear required')
        return abs(shear)*self.burgers_m*self.span_m/(2*self.stiffness_j_m)

    def frequency_response(self, frequency_hz, *, tail_tolerance=1e-9):
        """Infinite odd sine series, with explicit conservative absolute tail.

        Physical Hz describes this declared LINE reference, never production
        model frequency. exp(+i omega t) convention gives negative lag phase.
        """
        frequency = np.asarray(frequency_hz, float)
        if (np.any(~np.isfinite(frequency)) or np.any(frequency < 0)
                or not np.isfinite(tail_tolerance) or not 1e-13 <= tail_tolerance < 1):
            raise ValueError('nonnegative finite Hz and 1e-13 <= tail tolerance < 1 required')
        cutoff = max(1, int(np.ceil((32/(np.pi**4*tail_tolerance))**(1/3))))
        if cutoff % 2 == 0:
            cutoff += 1
        modes = np.arange(1, cutoff+1, 2, dtype=float)
        omega_tau = 2*np.pi*frequency*self.slowest_seconds
        response = (96/np.pi**4)*np.sum(
            1/(modes**4*(1+1j*omega_tau[..., None]/modes**2)), axis=-1)
        return dict(transfer=response, modes_used=len(modes), maximum_mode=cutoff,
                    absolute_tail_bound=32/(np.pi**4*cutoff**3))

    def _stiffness_bands(self, segments, diagonal_shift=0., dtype=float):
        if not isinstance(segments, (int, np.integer)) or segments < 4:
            raise ValueError('integer segment count >=4 required')
        dx = self.span_m/segments
        bands = np.zeros((3, segments-1), dtype=dtype)
        edge = self.stiffness_j_m/dx**2
        bands[1] = 2*edge+diagonal_shift
        bands[0, 1:] = -edge
        bands[2, :-1] = -edge
        return dx, bands

    def discrete_harmonic(self, frequency_hz, shear_amplitude_pa, *, segments=64):
        """Independent conservative finite-grid solution and work/dissipation.

        Uses the full complex tridiagonal operator, not the analytic series.
        No small phase is clipped. Units of work and dissipation are joules
        per cycle of this finite line, not per-cell activation energies.
        """
        if (not np.isfinite(frequency_hz) or frequency_hz <= 0
                or not np.isfinite(shear_amplitude_pa) or shear_amplitude_pa == 0):
            raise ValueError('positive Hz and finite nonzero shear amplitude required')
        omega = 2*np.pi*frequency_hz
        dx, bands = self._stiffness_bands(segments, 1j*omega*self.drag_pa_s, complex)
        pressure = shear_amplitude_pa*self.burgers_m
        bow = solve_banded((1, 1), bands, np.full(segments-1, pressure, complex))
        area = dx*np.sum(bow)
        work = -np.pi*pressure*area.imag
        dissipation = np.pi*omega*self.drag_pa_s*dx*np.sum(abs(bow)**2)
        return dict(x_m=np.linspace(0, self.span_m, segments+1),
            bow_hat_m=np.r_[0j, bow, 0j],
            transfer=area/(self.span_m*self.static_mean_bow_m(shear_amplitude_pa)),
            work_per_cycle_J=float(work), dissipation_per_cycle_J=float(dissipation),
            energy_balance_residual_J=float(work-dissipation))

    def integrate_backward_euler(self, times_seconds, shear_pa, *, segments=64):
        """Time-domain independent solve, initially straight, with fixed ends.

        No atomistic timestep or mass is invented. Fully implicit diffusion
        allows independently refined physical times; it is not a new fatigue
        probability path. The user must resolve the transient of interest.
        """
        times = np.asarray(times_seconds, float)
        shear = np.asarray(shear_pa, float)
        if (times.ndim != 1 or times.shape != shear.shape or len(times) < 2
                or np.any(~np.isfinite([times, shear])) or times[0] != 0
                or np.any(np.diff(times) <= 0)):
            raise ValueError('finite increasing seconds starting at zero and matching shear required')
        self._stiffness_bands(segments)  # validate before allocating the state
        bow = np.zeros(segments-1)
        mean = np.zeros(len(times))
        maximum_slope = np.zeros(len(times))
        for index, dt in enumerate(np.diff(times), 1):
            dx, bands = self._stiffness_bands(segments, self.drag_pa_s/dt)
            rhs = self.drag_pa_s/dt*bow+shear[index]*self.burgers_m
            bow = solve_banded((1, 1), bands, rhs)
            mean[index] = dx*np.sum(bow)/self.span_m
            maximum_slope[index] = np.max(abs(np.diff(np.r_[0., bow, 0.])))/dx
        return dict(time_seconds=times.copy(), shear_pa=shear.copy(), mean_bow_m=mean,
            maximum_slope=maximum_slope, final_bow_m=np.r_[0., bow, 0.],
            permanent_slip_generated=False, opening_probability=None,
            physical_time_scope='declared pinned-line reference only',
            production_clock_calibrated=False)
