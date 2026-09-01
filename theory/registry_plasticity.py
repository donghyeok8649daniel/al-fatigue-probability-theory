# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 이론 계산에 사용하는 Python 모듈이다.
# - 주요 클래스: RegistryTransportConfig, RegistryHistory
# - 주요 함수/메서드: RegistryTransportConfig.validate, registry_grid, metastable_well_distribution
#   _face_coefficients, registry_step, solve_registry
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Finite-volume probability dynamics for the active ideal-registry branch.

The state u=s/b is unwrapped across lattice periods.  Its integer well index
therefore records signed net lattice translations instead of folding every
state back into one period.  A residual shift of the well-index distribution
after unloading is the operational reduced-plasticity criterion.

The model is intentionally modest: it is an ideal one-registry mechanism with
an isothermal bath, not a dislocation-hardening law or a quantitative aluminum
crystal-plasticity model.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from theory.normal_lj_probability_dynamics import _chang_cooper_weight, _solve_tridiagonal
from theory.registry_lattice import (
    RegistryLattice,
    preferred_registry,
    registry_energy,
    registry_energy_derivative,
)


@dataclass(frozen=True)
class RegistryTransportConfig:
    """Dimensionless ideal-registry Smoluchowski parameters.

    Energy is scaled by the generalized-pair coefficient epsilon_c, length by
    b, and time by b^2/(M_s epsilon_c).  Consequently ``generalized_force`` is
    tau*A_repeat*b/epsilon_c and ``inverse_temperature`` is
    epsilon_c/(k_B*T).
    """

    lattice: RegistryLattice = RegistryLattice()
    inverse_temperature: float = 20.0
    u_min: float = -4.0
    u_max: float = 4.0
    cells: int = 480

    def validate(self) -> None:
        self.lattice.validate()
        if self.inverse_temperature <= 0.0:
            raise ValueError("inverse_temperature must be positive")
        if not (self.u_min < -0.5 and self.u_max > 0.5):
            raise ValueError("unwrapped domain must contain at least one full well")
        if self.cells < 80:
            raise ValueError("cells must be at least 80")


@dataclass(frozen=True)
class RegistryHistory:
    time: np.ndarray
    generalized_force: np.ndarray
    registry: np.ndarray
    density: np.ndarray
    mean_registry: np.ndarray
    mean_well_index: np.ndarray
    mean_intrawell_registry: np.ndarray
    variance: np.ndarray
    mean_lattice_energy: np.ndarray
    work: np.ndarray
    entropy_production: np.ndarray
    boundary_probability: np.ndarray
    preferred_registry: float


def registry_grid(c: RegistryTransportConfig) -> tuple[np.ndarray, float]:
    c.validate()
    dx = (c.u_max - c.u_min) / c.cells
    return c.u_min + (np.arange(c.cells) + 0.5) * dx, dx


def metastable_well_distribution(
    u: np.ndarray, dx: float, c: RegistryTransportConfig
) -> tuple[np.ndarray, float]:
    """Conditional Boltzmann density in one registry basin.

    This is an intact/metastable ensemble choice, not a global equilibrium on
    the unwrapped periodic line.  The latter is not normalizable.
    """

    delta0 = preferred_registry(c.lattice)
    # Put the selected minimum in the copy closest to the origin.
    if delta0 > 0.75:
        delta0 -= 1.0
    left_boundary, right_boundary = delta0 - 0.5, delta0 + 0.5
    minimum_energy = float(registry_energy(delta0, c.lattice))
    boundary_energy = min(
        float(registry_energy(left_boundary, c.lattice)),
        float(registry_energy(right_boundary, c.lattice)),
    )
    boundary_slope = max(
        abs(float(registry_energy_derivative(left_boundary, c.lattice))),
        abs(float(registry_energy_derivative(right_boundary, c.lattice))),
    )
    if boundary_energy <= minimum_energy + 1.0e-10 or boundary_slope > 1.0e-8:
        raise ValueError(
            "the selected one-period cell is not a verified metastable basin"
        )
    mask = (u >= left_boundary) & (u < right_boundary)
    energy = np.asarray(registry_energy(u, c.lattice), dtype=float)
    logp = -c.inverse_temperature * energy
    logp -= np.max(logp[mask])
    density = np.where(mask, np.exp(logp), 0.0)
    density /= float(np.sum(density) * dx)
    return density, delta0


def _face_coefficients(
    energy: np.ndarray, dx: float, generalized_force: float, beta: float
) -> tuple[np.ndarray, np.ndarray]:
    drift = generalized_force - np.diff(energy) / dx
    diffusion = 1.0 / beta
    delta = _chang_cooper_weight(drift * dx / diffusion)
    return (
        drift * (1.0 - delta) + diffusion / dx,
        drift * delta - diffusion / dx,
    )


def registry_step(
    density: np.ndarray,
    energy: np.ndarray,
    dx: float,
    dt: float,
    generalized_force: float,
    c: RegistryTransportConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance one backward-Euler reflecting FV step and return face currents."""

    if dt <= 0.0:
        raise ValueError("dt must be positive")
    left, right = _face_coefficients(
        energy, dx, float(generalized_force), c.inverse_temperature
    )
    ratio = dt / dx
    count = density.size
    lower = np.empty(count - 1)
    diagonal = np.ones(count)
    upper = np.empty(count - 1)
    diagonal[0] += ratio * left[0]
    upper[0] = ratio * right[0]
    for i in range(1, count - 1):
        lower[i - 1] = -ratio * left[i - 1]
        diagonal[i] += ratio * (left[i] - right[i - 1])
        upper[i] = ratio * right[i]
    lower[-1] = -ratio * left[-1]
    diagonal[-1] -= ratio * right[-1]
    updated = _solve_tridiagonal(lower, diagonal, upper, density)
    if float(np.min(updated)) < -1.0e-11:
        raise RuntimeError("negative registry density")
    updated = np.maximum(updated, 0.0)
    currents = left * updated[:-1] + right * updated[1:]
    return updated, currents


def solve_registry(
    time: np.ndarray,
    generalized_force: np.ndarray,
    c: RegistryTransportConfig = RegistryTransportConfig(),
    max_dt: float = 0.01,
) -> RegistryHistory:
    """Solve the full unwrapped registry-density evolution.

    Reflecting numerical edges only truncate the unwrapped coordinate.  A run
    is accepted only when ``boundary_probability`` remains negligible under a
    declared tolerance; increasing the domain is a numerical refinement, not
    a physical hardening mechanism.
    """

    c.validate()
    t = np.asarray(time, dtype=float)
    force = np.asarray(generalized_force, dtype=float)
    if (
        t.ndim != 1
        or force.shape != t.shape
        or t.size < 2
        or np.any(np.diff(t) <= 0.0)
        or max_dt <= 0.0
    ):
        raise ValueError("time/force arrays or max_dt are invalid")
    u, dx = registry_grid(c)
    energy = np.asarray(registry_energy(u, c.lattice), dtype=float)
    density, delta0 = metastable_well_distribution(u, dx, c)
    history = np.empty((t.size, u.size), dtype=float)
    history[0] = density
    entropy = np.zeros(t.size)

    for k in range(1, t.size):
        duration = float(t[k] - t[k - 1])
        ratio = duration / max_dt
        substeps = max(1, math.ceil(ratio - 1.0e-12 * max(1.0, abs(ratio))))
        dt = duration / substeps
        accumulated_entropy = 0.0
        for j in range(substeps):
            fraction = (j + 0.5) / substeps
            local_force = force[k - 1] + fraction * (force[k] - force[k - 1])
            density, currents = registry_step(
                density, energy, dx, dt, float(local_force), c
            )
            face_density = np.maximum(
                0.5 * (density[:-1] + density[1:]), 1.0e-300
            )
            accumulated_entropy += dt * float(np.sum(currents**2 / face_density) * dx)
        mass = float(np.sum(density) * dx)
        if abs(mass - 1.0) > 2.0e-10:
            raise RuntimeError("registry probability mass was not conserved")
        history[k] = density
        entropy[k] = accumulated_entropy / duration

    mean = np.sum(history * u[None, :], axis=1) * dx
    variance = np.sum(history * (u[None, :] - mean[:, None]) ** 2, axis=1) * dx
    well_index = np.floor(u - delta0 + 0.5)
    intrawell = u - delta0 - well_index
    mean_z = np.sum(history * well_index[None, :], axis=1) * dx
    mean_intra = np.sum(history * intrawell[None, :], axis=1) * dx
    mean_energy = np.sum(history * energy[None, :], axis=1) * dx
    work = np.zeros(t.size)
    work[1:] = np.cumsum(0.5 * (force[1:] + force[:-1]) * np.diff(mean))
    edge_cells = max(2, c.cells // 100)
    boundary_probability = (
        np.sum(history[:, :edge_cells], axis=1)
        + np.sum(history[:, -edge_cells:], axis=1)
    ) * dx
    return RegistryHistory(
        t,
        force,
        u,
        history,
        mean,
        mean_z,
        mean_intra,
        variance,
        mean_energy,
        work,
        entropy,
        boundary_probability,
        delta0,
    )
