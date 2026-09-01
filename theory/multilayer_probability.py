"""Probability observables for the coupled normal-spacing/registry state.

The Smoluchowski PDE is the underlying evolution law.  This module evaluates
the four named governing observables without confusing intrinsic recoverable
energy with cumulative irreversible dissipation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from theory.registry_lattice import plastic_well_index


@dataclass(frozen=True)
class GoverningMetrics:
    """Instantaneous G1, G2, G3-rate and G4 values."""

    mean_spacing: float
    mean_intrinsic_energy: float
    dissipation_rate: float
    normalization_or_survival: float


def probability_currents(
    density: np.ndarray,
    intrinsic_energy: np.ndarray,
    da: float,
    ds: float,
    q_a: float,
    q_s: float,
    mobility_a: float,
    mobility_s: float,
    kbt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the constant-mobility Smoluchowski currents on a regular grid.

    ``intrinsic_energy`` is U0 only.  External loading appears solely through
    Q_a=A0*sigma and Q_s=A0*M*sigma.
    """
    p = np.asarray(density, dtype=float)
    energy = np.asarray(intrinsic_energy, dtype=float)
    if p.shape != energy.shape or p.ndim != 2:
        raise ValueError("density and intrinsic_energy must be matching 2D arrays")
    if min(da, ds, mobility_a, mobility_s, kbt) <= 0.0 or np.any(p < 0.0):
        raise ValueError("grid spacings, mobilities and kBT must be positive")
    dU_da, dU_ds = np.gradient(energy, da, ds, edge_order=2)
    dP_da, dP_ds = np.gradient(p, da, ds, edge_order=2)
    return (
        -mobility_a * (p * (dU_da - q_a) + kbt * dP_da),
        -mobility_s * (p * (dU_ds - q_s) + kbt * dP_ds),
    )


def governing_equations_metrics(
    a: np.ndarray,
    s: np.ndarray,
    density: np.ndarray,
    intrinsic_energy: np.ndarray,
    current_a: np.ndarray,
    current_s: np.ndarray,
    mobility_a: float,
    mobility_s: float,
    reference_energy: float,
) -> GoverningMetrics:
    """Evaluate the repository's official four governing quantities.

    G3 is represented here by its nonnegative rate.  Time integration of that
    rate gives E_hyst; it is not assumed to remain stored in U0.
    """
    a_values, s_values = np.asarray(a), np.asarray(s)
    p = np.asarray(density)
    arrays = (intrinsic_energy, current_a, current_s)
    if a_values.ndim != 1 or s_values.ndim != 1 or p.shape != (a_values.size, s_values.size):
        raise ValueError("a, s and density grid shapes are inconsistent")
    if any(np.asarray(item).shape != p.shape for item in arrays):
        raise ValueError("energy/current shapes must match density")
    if mobility_a <= 0.0 or mobility_s <= 0.0 or np.any(p < 0.0):
        raise ValueError("mobilities must be positive and density nonnegative")
    da, ds = float(a_values[1] - a_values[0]), float(s_values[1] - s_values[0])
    measure = da * ds
    mass = float(np.sum(p) * measure)
    mean_a = float(np.sum(a_values[:, None] * p) * measure)
    mean_u = float(np.sum((np.asarray(intrinsic_energy) - reference_energy) * p) * measure)
    safe = np.maximum(p, 1.0e-300)
    dissipation = float(np.sum(
        np.asarray(current_a) ** 2 / (mobility_a * safe)
        + np.asarray(current_s) ** 2 / (mobility_s * safe)
    ) * measure)
    return GoverningMetrics(mean_a, mean_u, dissipation, mass)


def cumulative_hysteresis(time: np.ndarray, dissipation_rate: np.ndarray) -> np.ndarray:
    """G3: E_hyst(t)=integral dot(D)_irr dt by trapezoidal quadrature."""
    t, rate = np.asarray(time, dtype=float), np.asarray(dissipation_rate, dtype=float)
    if t.ndim != 1 or rate.shape != t.shape or t.size < 2 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time and dissipation arrays are invalid")
    if np.any(rate < -1.0e-13):
        raise ValueError("dissipation rate must be nonnegative")
    result = np.zeros_like(t)
    result[1:] = np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(t))
    return result


def well_populations(
    s: np.ndarray, density: np.ndarray, da: float, ds: float, b: float, s0: float,
) -> tuple[dict[int, float], float]:
    """Return p_z and <z> for unwrapped s=s0+z b+tilde{s}."""
    s_values, p = np.asarray(s), np.asarray(density)
    if p.ndim != 2 or p.shape[1] != s_values.size:
        raise ValueError("density must have s along its second axis")
    z = np.asarray(plastic_well_index(s_values, b, s0))
    populations = {
        int(index): float(np.sum(p[:, z == index]) * da * ds)
        for index in np.unique(z)
    }
    mean_z = float(sum(index * probability for index, probability in populations.items()))
    return populations, mean_z


def plastic_strains(mean_well_index: float, b: float, h_slip: float, schmid: float) -> tuple[float, float]:
    """Map a nonzero residual <z> to shear and tensile plastic strain."""
    if b <= 0.0 or h_slip <= 0.0:
        raise ValueError("b and h_slip must be positive")
    gamma_p = b * mean_well_index / h_slip
    return gamma_p, schmid * gamma_p


def moving_barrier_outflux(
    s: np.ndarray,
    barrier: np.ndarray,
    barrier_velocity: np.ndarray,
    density_on_barrier: np.ndarray,
    current_a_on_barrier: np.ndarray,
    current_s_on_barrier: np.ndarray,
) -> float:
    """Flux relative to a moving graph a=a^dagger(s,t).

    With periodic/no-flux s edges, Reynolds transport gives
    -dS/dt = integral [J_a-J_s*d_s(a^dagger)-P*d_t(a^dagger)] ds.
    At an ideal absorbing boundary P=0, the boundary-motion term vanishes,
    but it is retained for finite-grid/radiation implementations.
    """
    values = [np.asarray(x, dtype=float) for x in (
        s, barrier, barrier_velocity, density_on_barrier,
        current_a_on_barrier, current_s_on_barrier,
    )]
    if any(item.shape != values[0].shape for item in values[1:]):
        raise ValueError("all moving-boundary arrays must have the same shape")
    ds_barrier = np.gradient(values[1], values[0], edge_order=2)
    relative = values[4] - values[5] * ds_barrier - values[3] * values[2]
    return float(np.trapezoid(relative, values[0]))
