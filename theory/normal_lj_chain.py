# === 한국어 파일 안내 시작 ===
# - 파일 역할: calibrated generalized-LJ layer interaction을 사용하는 보존적 1D normal chain을 적분하고 에너지·spacing·instability 진단값을 계산한다.
# - 주요 클래스: NormalLJParameters, InstabilityEvent, NormalLJResult
# - 주요 함수/메서드: NormalLJResult.energy_balance_relative_error, _validate_exponents, normalized_lj_energy
#   normalized_lj_force, normalized_lj_stiffness, critical_stretch, critical_dimensionless_force
#   stress_to_dimensionless_force, atomic_time_scale, dimensionless_omega_from_frequency
#   physical_frequency_from_dimensionless_omega, simulate_normal_lj_chain
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Normal-opening generalized Lennard-Jones chain.

This module is the main proof-of-principle model for the project's current
normal-deformation direction. It intentionally contains:

- no viscous damping;
- no empirical fatigue damage variable;
- no slip/disregistry coordinate;
- no cycle-evolution law.

The microscopic state is a 1D chain of atom positions x_i. Nearest-neighbor
normal spacings are lambda_i = x_{i+1} - x_i in units of the equilibrium
spacing a0. The left atom is fixed and a prescribed normal cyclic force is
applied to the right atom.

The generalized-LJ energy is normalized so that

    lambda_eq = 1,
    phi'(1) = 0,
    phi''(1) = 1.

Under the earlier 1D stress mapping this makes the dimensionless external
normal force equal to sigma/E.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class NormalLJParameters:
    repulsive_exponent: float = 12.19
    attractive_exponent: float = 6.0
    mean_force: float = 0.0
    force_amplitude: float = 100.0e6 / 69.0e9
    omega: float = 0.02
    ramp_cycles: int = 2


@dataclass(frozen=True)
class InstabilityEvent:
    time: float
    cycle: float
    max_spacing: float


@dataclass
class NormalLJResult:
    time: np.ndarray
    force: np.ndarray
    max_spacing: np.ndarray
    internal_energy: np.ndarray
    external_work: np.ndarray
    cycle_mean_spacing: np.ndarray
    cycle_variance_spacing: np.ndarray
    cycle_max_spacing: np.ndarray
    cycle_min_spacing: np.ndarray
    cycle_snapshots: dict[int, np.ndarray]
    period: float
    first_instability: Optional[InstabilityEvent]

    @property
    def energy_balance_relative_error(self) -> float:
        delta_e = self.internal_energy[-1] - self.internal_energy[0]
        work = self.external_work[-1]
        return abs(delta_e - work) / max(abs(work), 1.0e-14)


def _validate_exponents(m: float, n: float) -> None:
    if not (m > n > 1.0):
        raise ValueError("require repulsive_exponent > attractive_exponent > 1")


def normalized_lj_energy(
    stretch: np.ndarray | float,
    m: float = 12.19,
    n: float = 6.0,
):
    """Normalized generalized-LJ pair energy phi(lambda).

    An additive constant is irrelevant and is intentionally omitted.
    """
    _validate_exponents(m, n)
    lam = np.asarray(stretch)
    return lam ** (-m) / (m * (m - n)) - lam ** (-n) / (n * (m - n))


def normalized_lj_force(
    stretch: np.ndarray | float,
    m: float = 12.19,
    n: float = 6.0,
):
    """Return d phi / d lambda."""
    _validate_exponents(m, n)
    lam = np.asarray(stretch)
    return (lam ** (-n - 1.0) - lam ** (-m - 1.0)) / (m - n)


def normalized_lj_stiffness(
    stretch: np.ndarray | float,
    m: float = 12.19,
    n: float = 6.0,
):
    """Return d^2 phi / d lambda^2."""
    _validate_exponents(m, n)
    lam = np.asarray(stretch)
    return (
        (m + 1.0) * lam ** (-m - 2.0)
        - (n + 1.0) * lam ** (-n - 2.0)
    ) / (m - n)


def critical_stretch(m: float = 12.19, n: float = 6.0) -> float:
    """Stretch where the normalized LJ tangent stiffness first vanishes."""
    _validate_exponents(m, n)
    return ((m + 1.0) / (n + 1.0)) ** (1.0 / (m - n))


def critical_dimensionless_force(m: float = 12.19, n: float = 6.0) -> float:
    lam_c = critical_stretch(m, n)
    return float(normalized_lj_force(lam_c, m, n))


def stress_to_dimensionless_force(stress_pa: float, youngs_modulus_pa: float) -> float:
    """Map normal stress to the dimensionless force used by this model."""
    if youngs_modulus_pa <= 0.0:
        raise ValueError("youngs_modulus_pa must be positive")
    return stress_pa / youngs_modulus_pa


def atomic_time_scale(
    atomic_mass_kg: float,
    a0_m: float,
    youngs_modulus_pa: float,
    reference_area_m2: float,
) -> float:
    """Time scale t0 = sqrt(M a0 / (E A0)) for the normalized coordinate."""
    if min(atomic_mass_kg, a0_m, youngs_modulus_pa, reference_area_m2) <= 0.0:
        raise ValueError("all physical scale inputs must be positive")
    return math.sqrt(
        atomic_mass_kg * a0_m / (youngs_modulus_pa * reference_area_m2)
    )


def dimensionless_omega_from_frequency(frequency_hz: float, t0_s: float) -> float:
    if frequency_hz < 0.0 or t0_s <= 0.0:
        raise ValueError("frequency must be non-negative and t0 positive")
    return 2.0 * math.pi * frequency_hz * t0_s


def physical_frequency_from_dimensionless_omega(
    omega_star: float,
    t0_s: float,
) -> float:
    if omega_star < 0.0 or t0_s <= 0.0:
        raise ValueError("omega_star must be non-negative and t0 positive")
    return omega_star / (2.0 * math.pi * t0_s)


def simulate_normal_lj_chain(
    parameters: NormalLJParameters = NormalLJParameters(),
    *,
    atoms: int = 32,
    dt: float = 0.01,
    cycles: int = 12,
    record_stride: int = 20,
    runaway_spacing: float = 3.0,
) -> NormalLJResult:
    """Integrate the conservative normal chain with velocity Verlet.

    Model:
        V = sum_i phi(lambda_i)
        lambda_i = x_{i+1} - x_i.

    Boundary conditions:
        x_0 is fixed;
        the rightmost atom receives the prescribed normal cyclic force.

    The physical instability diagnostic is the first local crossing of
    phi''(lambda)=0. ``runaway_spacing`` is only a numerical stop used after
    instability to avoid integrating an already separated chain indefinitely.
    """
    p = parameters
    _validate_exponents(p.repulsive_exponent, p.attractive_exponent)
    if atoms < 3:
        raise ValueError("atoms must be at least 3")
    if dt <= 0.0 or cycles <= 0 or record_stride <= 0:
        raise ValueError("dt, cycles, and record_stride must be positive")
    if p.omega <= 0.0:
        raise ValueError("omega must be positive")

    m = p.repulsive_exponent
    n = p.attractive_exponent
    lam_c = critical_stretch(m, n)
    period = 2.0 * math.pi / p.omega
    nsteps = int(round(cycles * period / dt))

    x = np.arange(atoms, dtype=float)
    velocity = np.zeros(atoms, dtype=float)

    def envelope(t: float) -> float:
        if p.ramp_cycles <= 0:
            return 1.0
        ramp_time = p.ramp_cycles * period
        if t >= ramp_time:
            return 1.0
        return 0.5 * (1.0 - math.cos(math.pi * t / ramp_time))

    def external_force(t: float) -> float:
        return envelope(t) * (
            p.mean_force + p.force_amplitude * math.sin(p.omega * t)
        )

    def force_vector(state: np.ndarray, t: float) -> np.ndarray:
        spacing = np.diff(state)
        dphi = normalized_lj_force(spacing, m, n)
        force = np.zeros_like(state)
        force[1:-1] = dphi[1:] - dphi[:-1]
        force[-1] = -dphi[-1] + external_force(t)
        return force

    def internal_energy() -> float:
        kinetic = 0.5 * float(np.dot(velocity[1:], velocity[1:]))
        potential = float(np.sum(normalized_lj_energy(np.diff(x), m, n)))
        return kinetic + potential

    force = force_vector(x, 0.0)
    work = 0.0
    previous_power = 0.0
    next_cycle = 1
    first_instability: Optional[InstabilityEvent] = None

    time_history = [0.0]
    force_history = [0.0]
    max_spacing_history = [1.0]
    energy_history = [internal_energy()]
    work_history = [0.0]

    cycle_mean = []
    cycle_variance = []
    cycle_max = []
    cycle_min = []
    snapshots: dict[int, np.ndarray] = {}

    for step in range(nsteps):
        t = step * dt

        velocity[1:] += 0.5 * dt * force[1:]
        x[1:] += dt * velocity[1:]

        new_force = force_vector(x, t + dt)
        velocity[1:] += 0.5 * dt * new_force[1:]
        force = new_force

        spacing = np.diff(x)
        f_now = external_force(t + dt)
        power = f_now * velocity[-1]
        work += 0.5 * dt * (previous_power + power)
        previous_power = power

        max_spacing = float(np.max(spacing))
        if first_instability is None and max_spacing >= lam_c:
            first_instability = InstabilityEvent(
                time=float(t + dt),
                cycle=float((t + dt) / period),
                max_spacing=max_spacing,
            )

        if t + dt >= next_cycle * period:
            cycle_mean.append(float(np.mean(spacing)))
            cycle_variance.append(float(np.var(spacing)))
            cycle_max.append(max_spacing)
            cycle_min.append(float(np.min(spacing)))
            snapshots[next_cycle] = spacing.copy()
            next_cycle += 1

        if step % record_stride == 0 or step == nsteps - 1:
            time_history.append(float(t + dt))
            force_history.append(float(f_now))
            max_spacing_history.append(max_spacing)
            energy_history.append(internal_energy())
            work_history.append(float(work))

        if max_spacing > runaway_spacing:
            break

    return NormalLJResult(
        time=np.asarray(time_history),
        force=np.asarray(force_history),
        max_spacing=np.asarray(max_spacing_history),
        internal_energy=np.asarray(energy_history),
        external_work=np.asarray(work_history),
        cycle_mean_spacing=np.asarray(cycle_mean),
        cycle_variance_spacing=np.asarray(cycle_variance),
        cycle_max_spacing=np.asarray(cycle_max),
        cycle_min_spacing=np.asarray(cycle_min),
        cycle_snapshots=snapshots,
        period=period,
        first_instability=first_instability,
    )
