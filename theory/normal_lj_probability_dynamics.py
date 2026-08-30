# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 이론 계산에 사용하는 Python 모듈이다.
# - 주요 클래스: SpacingDynamicsParameters, ProbabilityHistory
# - 주요 함수/메서드: validate_spacing_dynamics, _cell_grid, normalize_density
#   metastable_local_equilibrium_density, _chang_cooper_weight, _solve_tridiagonal
#   _implicit_probability_step, _density_observables, solve_spacing_probability_history
#   completed_cycle_hysteresis_areas
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Kinetic probability extension for the one-dimensional normal layer-LJ model.

The physical layer energy is written as

    U(a) = E0 * phi(lambda),       lambda = a / a0,
    E0 = E * A0 * a0,
    f(t) = sigma(t) / E.

Eliminating a fast thermal bath in the overdamped limit gives the reduced
Smoluchowski equation

    partial_tau p
      = partial_lambda[(phi'(lambda)-f)p + chi^(-1) partial_lambda p],

where tau=t/t_r, chi=E0/(k_B T), and t_r is the spacing relaxation time.
The solver below uses a fixed intact computational domain and zero probability
flux at both ends.  It therefore resolves rate-dependent distribution and
energy hysteresis, but deliberately does not turn the upper-tail mass into an
irreversible crack probability.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from theory.normal_lj_chain import (
    critical_stretch,
    normalized_lj_energy,
)


@dataclass(frozen=True)
class SpacingDynamicsParameters:
    """Physical and numerical inputs for the conditional intact distribution."""

    inverse_temperature: float = 2000.0
    relaxation_time_s: float = 0.01
    stretch_min: float = 0.82
    stretch_max: float = 1.25
    grid_cells: int = 240
    substeps_per_interval: int = 2
    max_reduced_substep: float = 0.08
    repulsive_exponent: float = 12.19
    attractive_exponent: float = 6.0


@dataclass(frozen=True)
class ProbabilityHistory:
    """Spacing-distribution history evaluated at the supplied FEM time samples."""

    time_s: np.ndarray
    stress_pa: np.ndarray
    stretch: np.ndarray
    density: np.ndarray
    normalization: np.ndarray
    mean_stretch: np.ndarray
    variance_stretch: np.ndarray
    mean_shifted_lj_energy: np.ndarray
    mean_energy_density_j_m3: np.ndarray
    nonequilibrium_free_energy_density_j_m3: np.ndarray
    cumulative_hysteresis_energy_density_j_m3: np.ndarray
    critical_tail_probability: np.ndarray


def validate_spacing_dynamics(parameters: SpacingDynamicsParameters) -> None:
    p = parameters
    if not (p.repulsive_exponent > p.attractive_exponent > 1.0):
        raise ValueError("require repulsive_exponent > attractive_exponent > 1")
    if not math.isfinite(p.inverse_temperature) or p.inverse_temperature <= 0.0:
        raise ValueError("inverse_temperature must be finite and positive")
    if not math.isfinite(p.relaxation_time_s) or p.relaxation_time_s <= 0.0:
        raise ValueError("relaxation_time_s must be finite and positive")
    if not (0.0 < p.stretch_min < 1.0 < p.stretch_max):
        raise ValueError("stretch domain must contain lambda=1 and remain positive")
    if p.grid_cells < 40:
        raise ValueError("grid_cells must be at least 40")
    if p.substeps_per_interval < 1:
        raise ValueError("substeps_per_interval must be at least 1")
    if p.max_reduced_substep <= 0.0:
        raise ValueError("max_reduced_substep must be positive")


def _cell_grid(parameters: SpacingDynamicsParameters) -> tuple[np.ndarray, float]:
    dx = (parameters.stretch_max - parameters.stretch_min) / parameters.grid_cells
    centers = parameters.stretch_min + (np.arange(parameters.grid_cells) + 0.5) * dx
    return centers, dx


def normalize_density(density: np.ndarray, cell_width: float) -> np.ndarray:
    """Normalize a nonnegative cell-centered density using finite-volume mass."""
    values = np.asarray(density, dtype=float)
    if values.ndim != 1 or values.size == 0 or cell_width <= 0.0:
        raise ValueError("density must be a non-empty 1D array and dx must be positive")
    if np.any(values < 0.0) or not np.all(np.isfinite(values)):
        raise ValueError("density must be finite and nonnegative")
    mass = float(np.sum(values) * cell_width)
    if mass <= 0.0 or not math.isfinite(mass):
        raise ValueError("density has nonpositive or nonfinite mass")
    return values / mass


def metastable_local_equilibrium_density(
    stretch: np.ndarray,
    cell_width: float,
    dimensionless_force: float,
    parameters: SpacingDynamicsParameters,
) -> np.ndarray:
    """Return the normalized Gibbs density on the declared intact domain.

    This is a local/conditional equilibrium initialization.  A tensile LJ
    potential is not globally normalizable on lambda in (0,infinity).
    """
    phi = normalized_lj_energy(
        stretch,
        parameters.repulsive_exponent,
        parameters.attractive_exponent,
    )
    log_weight = -parameters.inverse_temperature * (phi - dimensionless_force * stretch)
    log_weight -= float(np.max(log_weight))
    return normalize_density(np.exp(log_weight), cell_width)


def _chang_cooper_weight(reduced_drift: np.ndarray) -> np.ndarray:
    """Return delta(w)=1/w-1/(exp(w)-1) with a stable small-w expansion."""
    w = np.asarray(reduced_drift, dtype=float)
    delta = np.empty_like(w)
    small = np.abs(w) < 1.0e-5
    ws = w[small]
    delta[small] = 0.5 - ws / 12.0 + ws**3 / 720.0
    positive_large = w > 50.0
    negative_large = w < -50.0
    regular = ~(small | positive_large | negative_large)
    delta[positive_large] = 1.0 / w[positive_large]
    delta[negative_large] = 1.0 + 1.0 / w[negative_large]
    delta[regular] = 1.0 / w[regular] - 1.0 / np.expm1(w[regular])
    return delta


def _solve_tridiagonal(
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    """Thomas solve for one strictly tridiagonal system."""
    n = diagonal.size
    c = np.asarray(upper, dtype=float).copy()
    d = np.asarray(rhs, dtype=float).copy()
    b = np.asarray(diagonal, dtype=float).copy()
    a = np.asarray(lower, dtype=float)
    if n == 0 or d.size != n or a.size != n - 1 or c.size != n - 1:
        raise ValueError("invalid tridiagonal system dimensions")
    for i in range(1, n):
        if abs(b[i - 1]) < 1.0e-300:
            raise RuntimeError("singular probability transport matrix")
        factor = a[i - 1] / b[i - 1]
        b[i] -= factor * c[i - 1]
        d[i] -= factor * d[i - 1]
    if abs(b[-1]) < 1.0e-300:
        raise RuntimeError("singular probability transport matrix")
    x = np.empty(n, dtype=float)
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


def _implicit_probability_step(
    density: np.ndarray,
    stretch: np.ndarray,
    cell_width: float,
    reduced_dt: float,
    dimensionless_force: float,
    parameters: SpacingDynamicsParameters,
) -> np.ndarray:
    """Advance the no-flux Smoluchowski equation by one backward-Euler step."""
    if reduced_dt <= 0.0:
        return density.copy()
    phi = normalized_lj_energy(
        stretch,
        parameters.repulsive_exponent,
        parameters.attractive_exponent,
    )
    # A face-wise discrete energy gradient makes the declared finite-volume
    # Gibbs density an exact zero-current state, not merely an O(dx^2) one.
    drift = dimensionless_force - np.diff(phi) / cell_width
    diffusion = 1.0 / parameters.inverse_temperature
    w = drift * cell_width / diffusion
    delta = _chang_cooper_weight(w)

    # J_(i+1/2) = c_left[i] p_i + c_right[i] p_(i+1).
    c_left = drift * (1.0 - delta) + diffusion / cell_width
    c_right = drift * delta - diffusion / cell_width
    ratio = reduced_dt / cell_width
    n = density.size
    lower = np.empty(n - 1, dtype=float)
    diagonal = np.ones(n, dtype=float)
    upper = np.empty(n - 1, dtype=float)

    diagonal[0] += ratio * c_left[0]
    upper[0] = ratio * c_right[0]
    for i in range(1, n - 1):
        lower[i - 1] = -ratio * c_left[i - 1]
        diagonal[i] += ratio * (c_left[i] - c_right[i - 1])
        upper[i] = ratio * c_right[i]
    lower[-1] = -ratio * c_left[-1]
    diagonal[-1] -= ratio * c_right[-1]

    updated = _solve_tridiagonal(lower, diagonal, upper, density)
    if float(np.min(updated)) < -1.0e-11:
        raise RuntimeError("probability solver produced a materially negative density")
    updated = np.maximum(updated, 0.0)
    return normalize_density(updated, cell_width)


def _density_observables(
    density: np.ndarray,
    stretch: np.ndarray,
    cell_width: float,
    dimensionless_force: float,
    parameters: SpacingDynamicsParameters,
) -> tuple[float, float, float, float, float, float]:
    mass = float(np.sum(density) * cell_width)
    mean = float(np.sum(stretch * density) * cell_width)
    variance = float(np.sum((stretch - mean) ** 2 * density) * cell_width)
    phi = normalized_lj_energy(
        stretch,
        parameters.repulsive_exponent,
        parameters.attractive_exponent,
    )
    psi = phi - float(
        normalized_lj_energy(
            1.0,
            parameters.repulsive_exponent,
            parameters.attractive_exponent,
        )
    )
    mean_energy = float(np.sum(psi * density) * cell_width)
    entropy_integrand = np.zeros_like(density)
    positive = density > 0.0
    entropy_integrand[positive] = density[positive] * np.log(density[positive])
    mean_biased_energy = float(
        np.sum((psi - dimensionless_force * (stretch - 1.0)) * density)
        * cell_width
    )
    free_energy = mean_biased_energy + float(
        np.sum(entropy_integrand) * cell_width / parameters.inverse_temperature
    )
    lam_c = critical_stretch(
        parameters.repulsive_exponent,
        parameters.attractive_exponent,
    )
    tail = float(np.sum(density[stretch >= lam_c]) * cell_width)
    return mass, mean, variance, mean_energy, free_energy, tail


def solve_spacing_probability_history(
    time_s: np.ndarray,
    stress_pa: np.ndarray,
    youngs_modulus_pa: float,
    parameters: SpacingDynamicsParameters = SpacingDynamicsParameters(),
) -> ProbabilityHistory:
    """Solve local p(lambda,t) for one FEM element stress history.

    The hysteresis energy density is the path integral

        w_h(t) = integral_0^t sigma(s) d mean(lambda)(s),

    because mean(lambda)-1 is the reduced normal strain and E0/(A0*a0)=E.
    """
    validate_spacing_dynamics(parameters)
    time = np.asarray(time_s, dtype=float)
    stress = np.asarray(stress_pa, dtype=float)
    if time.ndim != 1 or stress.shape != time.shape or time.size < 2:
        raise ValueError("time and stress must be equal-size 1D arrays with at least 2 values")
    if np.any(np.diff(time) <= 0.0) or not np.all(np.isfinite(time)):
        raise ValueError("time must be finite and strictly increasing")
    if not np.all(np.isfinite(stress)):
        raise ValueError("stress must be finite")
    if not math.isfinite(youngs_modulus_pa) or youngs_modulus_pa <= 0.0:
        raise ValueError("youngs_modulus_pa must be finite and positive")

    stretch, dx = _cell_grid(parameters)
    force = stress / youngs_modulus_pa
    density = metastable_local_equilibrium_density(stretch, dx, float(force[0]), parameters)
    density_history = np.empty((time.size, stretch.size), dtype=float)
    density_history[0] = density

    for k in range(1, time.size):
        physical_dt = float(time[k] - time[k - 1])
        reduced_interval = physical_dt / parameters.relaxation_time_s
        substeps = max(
            parameters.substeps_per_interval,
            int(math.ceil(reduced_interval / parameters.max_reduced_substep)),
        )
        if substeps > 20000:
            raise ValueError("required probability substeps exceed 20000 per FEM interval")
        reduced_dt = reduced_interval / substeps
        for j in range(substeps):
            alpha = (j + 0.5) / substeps
            local_force = (1.0 - alpha) * force[k - 1] + alpha * force[k]
            density = _implicit_probability_step(
                density,
                stretch,
                dx,
                reduced_dt,
                float(local_force),
                parameters,
            )
        density_history[k] = density

    normalization = np.empty(time.size)
    mean = np.empty(time.size)
    variance = np.empty(time.size)
    mean_energy = np.empty(time.size)
    free_energy = np.empty(time.size)
    tail = np.empty(time.size)
    for k in range(time.size):
        (
            normalization[k],
            mean[k],
            variance[k],
            mean_energy[k],
            free_energy[k],
            tail[k],
        ) = _density_observables(
            density_history[k],
            stretch,
            dx,
            float(force[k]),
            parameters,
        )

    cumulative_work = np.zeros(time.size)
    cumulative_work[1:] = np.cumsum(
        0.5 * (stress[1:] + stress[:-1]) * np.diff(mean)
    )
    return ProbabilityHistory(
        time_s=time,
        stress_pa=stress,
        stretch=stretch,
        density=density_history,
        normalization=normalization,
        mean_stretch=mean,
        variance_stretch=variance,
        mean_shifted_lj_energy=mean_energy,
        mean_energy_density_j_m3=youngs_modulus_pa * mean_energy,
        nonequilibrium_free_energy_density_j_m3=youngs_modulus_pa * free_energy,
        cumulative_hysteresis_energy_density_j_m3=cumulative_work,
        critical_tail_probability=tail,
    )


def completed_cycle_hysteresis_areas(
    history: ProbabilityHistory,
    frequency_hz: float,
) -> np.ndarray:
    """Return integral sigma d(mean lambda) for every completed load cycle."""
    if frequency_hz <= 0.0:
        raise ValueError("frequency_hz must be positive")
    period = 1.0 / frequency_hz
    completed = int(math.floor((history.time_s[-1] - history.time_s[0]) / period + 1.0e-10))
    if completed <= 0:
        return np.empty(0, dtype=float)
    boundaries = history.time_s[0] + period * np.arange(completed + 1)
    work_at_boundaries = np.interp(
        boundaries,
        history.time_s,
        history.cumulative_hysteresis_energy_density_j_m3,
    )
    return np.diff(work_at_boundaries)
