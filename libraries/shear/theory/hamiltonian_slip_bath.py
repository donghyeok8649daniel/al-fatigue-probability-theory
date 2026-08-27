"""Hamiltonian nonlinear slip coordinate coupled to a long harmonic lattice bath.

This module is intentionally phenomenology-light:
- no viscous damping term,
- no empirical fatigue damage variable,
- no fitted cycle-evolution equation.

The nonlinear periodic slip potential is a controlled approximation to an atomistic
GSF/gamma-surface and must not be confused with a calibrated Al potential.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SlipBathParameters:
    resolved_mass: float = 1.0
    bath_mass: float = 1.0
    coupling_stiffness: float = 1.0
    bath_stiffness: float = 1.0
    slip_period: float = 1.0
    barrier_height: float = 0.1
    mean_force: float = 0.0
    force_amplitude: float = 0.5
    omega: float = 0.2
    ramp_cycles: int = 2


@dataclass
class SlipBathResult:
    time: np.ndarray
    slip: np.ndarray
    slip_velocity: np.ndarray
    force: np.ndarray
    internal_energy: np.ndarray
    external_work: np.ndarray
    cycle_slip: np.ndarray
    cycle_spacing_variance: np.ndarray
    period: float

    @property
    def final_energy_balance_relative_error(self) -> float:
        delta_e = self.internal_energy[-1] - self.internal_energy[0]
        w = self.external_work[-1]
        return abs(delta_e - w) / max(abs(w), 1.0e-14)


def slip_potential(s: np.ndarray | float, p: SlipBathParameters):
    """One-harmonic periodic approximation to a gamma-surface path."""
    return 0.5 * p.barrier_height * (
        1.0 - np.cos(2.0 * np.pi * np.asarray(s) / p.slip_period)
    )


def slip_potential_gradient(s: float, p: SlipBathParameters) -> float:
    return (
        np.pi * p.barrier_height / p.slip_period
        * np.sin(2.0 * np.pi * s / p.slip_period)
    )


def _envelope(t: float, period: float, ramp_cycles: int) -> float:
    if ramp_cycles <= 0:
        return 1.0
    ramp_time = ramp_cycles * period
    if t >= ramp_time:
        return 1.0
    return 0.5 * (1.0 - np.cos(np.pi * t / ramp_time))


def external_force(t: float, p: SlipBathParameters) -> float:
    period = 2.0 * np.pi / p.omega
    return p.mean_force + p.force_amplitude * _envelope(
        t, period, p.ramp_cycles
    ) * np.sin(p.omega * t)


def simulate_slip_bath(
    parameters: SlipBathParameters = SlipBathParameters(),
    *,
    bath_sites: int = 800,
    dt: float = 0.01,
    cycles: int = 12,
    record_stride: int = 5,
) -> SlipBathResult:
    """Integrate the full conservative equations with velocity Verlet.

    The bath is made long enough that the intended reference runs finish before
    reflected waves can return from the free far boundary.
    """
    p = parameters
    if bath_sites < 3:
        raise ValueError("bath_sites must be at least 3")
    if dt <= 0.0 or cycles <= 0 or record_stride <= 0:
        raise ValueError("dt, cycles, and record_stride must be positive")

    period = 2.0 * np.pi / p.omega
    nsteps = int(round(cycles * period / dt))

    s = 0.0
    vs = 0.0
    u = np.zeros(bath_sites, dtype=float)
    vu = np.zeros(bath_sites, dtype=float)

    def forces(s_: float, u_: np.ndarray, t_: float):
        fs_ = (
            -slip_potential_gradient(s_, p)
            - p.coupling_stiffness * (s_ - u_[0])
            + external_force(t_, p)
        )
        fu_ = np.empty_like(u_)
        fu_[0] = (
            p.coupling_stiffness * (s_ - u_[0])
            + p.bath_stiffness * (u_[1] - u_[0])
        )
        fu_[1:-1] = p.bath_stiffness * (
            u_[2:] - 2.0 * u_[1:-1] + u_[:-2]
        )
        fu_[-1] = p.bath_stiffness * (u_[-2] - u_[-1])
        return fs_, fu_

    def internal_energy():
        return (
            0.5 * p.resolved_mass * vs * vs
            + float(slip_potential(s, p))
            + 0.5 * p.coupling_stiffness * (s - u[0]) ** 2
            + 0.5 * p.bath_mass * np.dot(vu, vu)
            + 0.5 * p.bath_stiffness * np.dot(np.diff(u), np.diff(u))
        )

    fs, fu = forces(s, u, 0.0)

    time = []
    slip = []
    slip_velocity = []
    force = []
    energy = []
    work_history = []
    cycle_slip = []
    cycle_spacing_variance = []

    work = 0.0
    previous_power = external_force(0.0, p) * vs
    next_cycle = 1

    # Store the exact initial state so the energy-balance diagnostic has E(0).
    time.append(0.0)
    slip.append(s)
    slip_velocity.append(vs)
    force.append(external_force(0.0, p))
    energy.append(internal_energy())
    work_history.append(0.0)

    for step in range(nsteps):
        t = step * dt

        vs += 0.5 * dt * fs / p.resolved_mass
        vu += 0.5 * dt * fu / p.bath_mass

        s += dt * vs
        u += dt * vu

        fs, fu = forces(s, u, t + dt)

        vs += 0.5 * dt * fs / p.resolved_mass
        vu += 0.5 * dt * fu / p.bath_mass

        f_now = external_force(t + dt, p)
        power = f_now * vs
        work += 0.5 * dt * (previous_power + power)
        previous_power = power

        if t + dt >= next_cycle * period:
            spacing_like = np.concatenate(([s - u[0]], np.diff(u)))
            cycle_slip.append(s)
            cycle_spacing_variance.append(float(np.var(spacing_like)))
            next_cycle += 1

        if step % record_stride == 0 or step == nsteps - 1:
            time.append(t + dt)
            slip.append(s)
            slip_velocity.append(vs)
            force.append(f_now)
            energy.append(internal_energy())
            work_history.append(work)

    return SlipBathResult(
        time=np.asarray(time),
        slip=np.asarray(slip),
        slip_velocity=np.asarray(slip_velocity),
        force=np.asarray(force),
        internal_energy=np.asarray(energy),
        external_work=np.asarray(work_history),
        cycle_slip=np.asarray(cycle_slip),
        cycle_spacing_variance=np.asarray(cycle_spacing_variance),
        period=period,
    )


def cycle_work(result: SlipBathResult, cycle_index: int) -> float:
    """Return integral F ds over one zero-based cycle index."""
    t0 = cycle_index * result.period
    t1 = (cycle_index + 1) * result.period
    mask = (result.time >= t0) & (result.time < t1)
    if np.count_nonzero(mask) < 3:
        raise ValueError("not enough samples in requested cycle")
    power = result.force[mask] * result.slip_velocity[mask]
    return float(np.trapezoid(power, result.time[mask]))
