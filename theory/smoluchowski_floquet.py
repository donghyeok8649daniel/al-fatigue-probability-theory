"""One-cycle survival spectrum of the absorbing Smoluchowski equation.

No fatigue-life law is postulated here.  A periodic load and the already
defined linear absorbing finite-volume evolution determine a positive,
mass-decreasing one-cycle operator K.  Its principal Perron multiplier is the
long-time survival ratio per cycle of that discretized physical model.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from theory.smoluchowski_escape import (
    TransportConfig,
    advance_interval,
    conditional_equilibrium,
    finite_volume_generator,
    transport_grid,
)


@dataclass(frozen=True)
class CycleTrace:
    """Unnormalised evolution over one period from a unit-mass density."""

    time: np.ndarray
    force: np.ndarray
    density: np.ndarray
    survival: np.ndarray
    outflux: np.ndarray
    integrated_outflux: float


@dataclass(frozen=True)
class FloquetResult:
    """Principal one-cycle survival mode and its phase-resolved evolution."""

    period: float
    multiplier: float
    escape_per_cycle: float
    integrated_hazard: float
    mean_hazard_rate: float
    start_density: np.ndarray
    phase_time: np.ndarray
    phase_force: np.ndarray
    phase_density: np.ndarray
    phase_conditional_density: np.ndarray
    phase_survival: np.ndarray
    phase_outflux: np.ndarray
    phase_hazard: np.ndarray
    iterations: int
    residual_l1: float
    multiplier_history: np.ndarray


@dataclass(frozen=True)
class DenseCycleSpectrum:
    """Biorthogonal Perron modes and the leading transient contraction."""

    multiplier: float
    second_eigenvalue_modulus: float
    spectral_ratio: float
    right_density: np.ndarray
    left_survival_weight: np.ndarray
    operator: np.ndarray


def _validate_protocol(time: np.ndarray, force: np.ndarray,
                       c: TransportConfig) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(time, dtype=float)
    f = np.asarray(force, dtype=float)
    if c.boundary != "absorbing":
        raise ValueError("Floquet survival requires an absorbing boundary")
    if (t.ndim != 1 or f.shape != t.shape or t.size < 3
            or np.any(np.diff(t) <= 0)):
        raise ValueError("time and force must be matching, increasing arrays")
    scale = max(1.0, float(np.max(np.abs(f))))
    if abs(float(f[-1] - f[0])) > 1e-12 * scale:
        raise ValueError("one-cycle protocol must have identical endpoint force")
    return t, f


def propagate_cycle(density: np.ndarray, time: np.ndarray, force: np.ndarray,
                    c: TransportConfig, max_dt: float = 0.02,
                    record: bool = False) -> tuple[np.ndarray, float, CycleTrace | None]:
    """Apply the discrete one-cycle evolution operator once.

    ``density`` is not renormalised internally, so linearity and absorbing
    mass balance are retained.  ``integrated_outflux`` is the exact sum of the
    backward-Euler boundary fluxes used by the update.
    """
    t, f = _validate_protocol(time, force, c)
    x, dx = transport_grid(c)
    rho = np.asarray(density, dtype=float).copy()
    if rho.shape != x.shape or np.any(rho < 0) or not np.all(np.isfinite(rho)):
        raise ValueError("density must be finite, nonnegative and match the grid")

    densities = np.empty((t.size, x.size)) if record else None
    survival = np.empty(t.size) if record else None
    outflux = np.zeros(t.size) if record else None
    if record:
        densities[0] = rho
        survival[0] = np.sum(rho) * dx

    escaped = 0.0
    for k in range(1, t.size):
        rho, interval_flux = advance_interval(
            rho, x, dx, float(t[k - 1]), float(t[k]),
            float(f[k - 1]), float(f[k]), c, max_dt=max_dt)
        escaped += interval_flux
        if record:
            densities[k] = rho
            survival[k] = np.sum(rho) * dx
            outflux[k] = interval_flux / (t[k] - t[k - 1])

    trace = None
    if record:
        trace = CycleTrace(t.copy(), f.copy(), densities, survival, outflux, escaped)
    return rho, escaped, trace


def principal_survival_mode(time: np.ndarray, force: np.ndarray,
                            c: TransportConfig,
                            max_dt: float = 0.02,
                            tolerance: float = 1e-11,
                            max_iterations: int = 2000,
                            initial_density: np.ndarray | None = None
                            ) -> FloquetResult:
    """Compute the positive principal mode by mass-normalised power iteration."""
    t, f = _validate_protocol(time, force, c)
    if not (max_dt > 0 and tolerance > 0 and max_iterations >= 2):
        raise ValueError("invalid iteration controls")
    x, dx = transport_grid(c)
    if initial_density is None:
        q = conditional_equilibrium(x, dx, float(f[0]), c)
    else:
        q = np.asarray(initial_density, dtype=float).copy()
        if q.shape != x.shape or np.any(q < 0) or not np.all(np.isfinite(q)):
            raise ValueError("invalid initial density")
        mass = float(np.sum(q) * dx)
        if mass <= 0:
            raise ValueError("initial density has zero mass")
        q /= mass

    history: list[float] = []
    residual = np.inf
    for iteration in range(1, max_iterations + 1):
        end, _, _ = propagate_cycle(q, t, f, c, max_dt=max_dt)
        multiplier = float(np.sum(end) * dx)
        if not (0.0 < multiplier <= 1.0 + 1e-12):
            raise RuntimeError("one-cycle survival multiplier is outside (0,1]")
        q_next = end / multiplier
        residual = float(np.sum(np.abs(q_next - q)) * dx)
        history.append(multiplier)
        q = q_next
        if residual < tolerance:
            break
    else:
        raise RuntimeError("principal survival mode did not converge")

    # Reapply K to the converged start state so the reported multiplier,
    # residual, flux balance and phase trace refer to one identical cycle.
    end, escaped, trace = propagate_cycle(q, t, f, c, max_dt=max_dt, record=True)
    assert trace is not None
    multiplier = float(np.sum(end) * dx)
    residual = float(np.sum(np.abs(end / multiplier - q)) * dx)
    conditional = np.divide(
        trace.density, trace.survival[:, None],
        out=np.zeros_like(trace.density), where=trace.survival[:, None] > 0)
    integrated_hazard = -float(np.log(multiplier))
    period = float(t[-1] - t[0])
    phase_hazard = np.zeros_like(trace.survival)
    phase_hazard[1:] = -np.log(
        trace.survival[1:] / trace.survival[:-1]) / np.diff(t)
    # The conditional mode is periodic. Index zero has no preceding interval
    # in this stored cycle, so use the matching final phase interval.
    phase_hazard[0] = phase_hazard[-1]
    return FloquetResult(
        period=period,
        multiplier=multiplier,
        escape_per_cycle=1.0 - multiplier,
        integrated_hazard=integrated_hazard,
        mean_hazard_rate=integrated_hazard / period,
        start_density=q,
        phase_time=t.copy(),
        phase_force=f.copy(),
        phase_density=trace.density,
        phase_conditional_density=conditional,
        phase_survival=trace.survival,
        phase_outflux=trace.outflux,
        phase_hazard=phase_hazard,
        iterations=iteration,
        residual_l1=residual,
        multiplier_history=np.asarray(history),
    )


def cycle_operator_matrix(time: np.ndarray, force: np.ndarray,
                          c: TransportConfig,
                          max_dt: float = 0.02) -> np.ndarray:
    """Construct K explicitly for independent small-grid spectral checks."""
    _validate_protocol(time, force, c)
    x, _ = transport_grid(c)
    operator = np.empty((x.size, x.size))
    for j in range(x.size):
        basis = np.zeros(x.size)
        basis[j] = 1.0
        operator[:, j], _, _ = propagate_cycle(
            basis, time, force, c, max_dt=max_dt)
    return operator


def dense_cycle_spectrum(time: np.ndarray, force: np.ndarray,
                         c: TransportConfig,
                         max_dt: float = 0.02) -> DenseCycleSpectrum:
    """Return right/left Perron modes of an explicitly assembled cycle map.

    This is intended for verification and moderate grids.  The right mode is
    normalized to unit probability mass.  The left mode ``w`` is normalized
    so ``integral(w*q)=1``.  Consequently a unit-mass initial density ``rho``
    has the asymptotic survival prefactor ``integral(w*rho)``.
    """
    operator = cycle_operator_matrix(time, force, c, max_dt=max_dt)
    x, dx = transport_grid(c)
    eigenvalues = np.linalg.eigvals(operator)
    order = np.argsort(np.abs(eigenvalues))[::-1]
    principal_index = int(order[0])
    multiplier = float(np.real(eigenvalues[principal_index]))
    if abs(float(np.imag(eigenvalues[principal_index]))) > 1e-10:
        raise RuntimeError("principal cycle eigenvalue is unexpectedly complex")
    right = np.ones(x.size)
    right /= np.sum(right) * dx
    for _ in range(10000):
        right_next = operator @ right
        right_next /= np.sum(right_next) * dx
        if np.sum(np.abs(right_next - right)) * dx < 1e-13:
            right = right_next
            break
        right = right_next
    else:
        raise RuntimeError("principal right density did not converge")

    # Direct left eigenvectors of this strongly nonnormal, absorbing matrix
    # can be poorly conditioned near the repulsive boundary.  Positive power
    # iteration on K.T preserves the Perron cone and is numerically stable.
    left = np.ones(x.size)
    left /= np.dot(left, right) * dx
    for _ in range(10000):
        left_next = operator.T @ left
        left_next /= np.dot(left_next, right) * dx
        if np.max(np.abs(left_next - left)) < 1e-13:
            left = left_next
            break
        left = left_next
    else:
        raise RuntimeError("principal left survival weight did not converge")

    second = float(np.abs(eigenvalues[order[1]])) if x.size > 1 else 0.0
    return DenseCycleSpectrum(
        multiplier=multiplier,
        second_eigenvalue_modulus=second,
        spectral_ratio=second / multiplier,
        right_density=right,
        left_survival_weight=left,
        operator=operator,
    )


def asymptotic_survival_prefactor(initial_density: np.ndarray,
                                  spectrum: DenseCycleSpectrum,
                                  c: TransportConfig) -> float:
    """Compute the Perron asymptotic coefficient from an initial density."""
    x, dx = transport_grid(c)
    rho = np.asarray(initial_density, dtype=float)
    if rho.shape != x.shape or np.any(rho < 0):
        raise ValueError("invalid initial density")
    mass = float(np.sum(rho) * dx)
    if mass <= 0:
        raise ValueError("initial density has zero mass")
    return float(np.dot(spectrum.left_survival_weight, rho / mass) * dx)


def frozen_principal_escape_rate(force: float, c: TransportConfig) -> float:
    """Return minus the dominant eigenvalue of the frozen-load subgenerator."""
    if c.boundary != "absorbing":
        raise ValueError("frozen escape rate requires an absorbing boundary")
    eigenvalues = np.linalg.eigvals(finite_volume_generator(float(force), c))
    dominant = eigenvalues[int(np.argmax(np.real(eigenvalues)))]
    if abs(float(np.imag(dominant))) > 1e-9:
        raise RuntimeError("dominant frozen generator eigenvalue is complex")
    rate = -float(np.real(dominant))
    if rate < -1e-11:
        raise RuntimeError("frozen absorbing generator has positive growth")
    return max(rate, 0.0)


def adiabatic_mean_escape_rate(force: np.ndarray,
                               c: TransportConfig) -> float:
    """Phase-average frozen principal rate for a uniformly sampled cycle."""
    values = np.asarray(force, dtype=float)
    if values.ndim != 1 or values.size < 2:
        raise ValueError("force must be a one-dimensional cycle sample")
    rates = np.asarray([frozen_principal_escape_rate(value, c)
                        for value in values])
    return float(np.trapezoid(rates, dx=1.0 / (values.size - 1)))


def direct_cycle_survival_ratios(initial_density: np.ndarray,
                                 time: np.ndarray, force: np.ndarray,
                                 c: TransportConfig, cycles: int,
                                 max_dt: float = 0.02
                                 ) -> tuple[np.ndarray, np.ndarray]:
    """Return boundary survival and successive ratios without renormalisation."""
    if cycles < 1:
        raise ValueError("cycles must be positive")
    x, dx = transport_grid(c)
    rho = np.asarray(initial_density, dtype=float).copy()
    if rho.shape != x.shape or np.any(rho < 0):
        raise ValueError("invalid initial density")
    mass0 = float(np.sum(rho) * dx)
    if mass0 <= 0:
        raise ValueError("initial density has zero mass")
    rho /= mass0
    survival = np.empty(cycles + 1)
    survival[0] = 1.0
    for k in range(1, cycles + 1):
        rho, _, _ = propagate_cycle(rho, time, force, c, max_dt=max_dt)
        survival[k] = np.sum(rho) * dx
    return survival, survival[1:] / survival[:-1]
