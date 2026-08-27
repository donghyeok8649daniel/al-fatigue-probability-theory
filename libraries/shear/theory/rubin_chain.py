from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict

import numpy as np


@dataclass(frozen=True)
class RubinParams:
    """Parameters for a resolved oscillator coupled to a harmonic chain."""

    system_mass: float = 1.0
    bath_mass: float = 1.0
    bath_spring: float = 1.0
    system_spring: float = 1.0


def band_edge(params: RubinParams) -> float:
    """Upper angular-frequency edge of the monatomic chain."""

    return 2.0 * math.sqrt(params.bath_spring / params.bath_mass)


def dynamic_stiffness(
    omega: float,
    params: RubinParams = RubinParams(),
) -> complex:
    """Exact semi-infinite-chain dynamic stiffness.

    Time dependence is exp(i*omega*t) and

        F_hat = Z(omega) * Q_hat.

    No phenomenological damping coefficient is used.
    """

    if omega < 0:
        raise ValueError("omega must be non-negative")

    M = params.system_mass
    m = params.bath_mass
    k = params.bath_spring
    k0 = params.system_spring
    omega_d = band_edge(params)

    if omega == 0.0:
        return complex(k0, 0.0)

    if omega < omega_d:
        q = 2.0 * math.asin(omega / omega_d)
        bath_term = k * (1.0 - np.exp(-1j * q))
        return complex(k0 - M * omega**2, 0.0) + bath_term

    # Above the propagating band, q = pi - i*kappa gives an
    # evanescent lattice field and a purely real dynamic stiffness.
    kappa = 2.0 * math.acosh(omega / omega_d)
    exp_minus_iq = -math.exp(-kappa)
    bath_term = k * (1.0 - exp_minus_iq)
    return complex(k0 - M * omega**2 + bath_term, 0.0)


def analytic_response(
    omega: float,
    force_amplitude: float,
    params: RubinParams = RubinParams(),
) -> Dict[str, float]:
    """Steady harmonic response and analytic hysteresis-loop area."""

    z = dynamic_stiffness(omega, params)
    q_amp = force_amplitude / abs(z)
    phase_lag = math.atan2(z.imag, z.real)
    loop_area = math.pi * z.imag * q_amp**2

    return {
        "omega": omega,
        "band_edge": band_edge(params),
        "z_real": float(z.real),
        "z_imag": float(z.imag),
        "response_amplitude": float(q_amp),
        "phase_lag_rad": float(phase_lag),
        "phase_lag_deg": float(math.degrees(phase_lag)),
        "loop_area": float(loop_area),
    }


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(x)))


def simulate_finite_chain(
    *,
    n_masses: int = 1200,
    omega: float = 0.5,
    force_amplitude: float = 0.1,
    dt: float = 0.02,
    n_periods: int = 60,
    ramp_periods: int = 5,
    params: RubinParams = RubinParams(),
) -> Dict[str, np.ndarray | float]:
    """Integrate the full conservative finite chain with velocity Verlet.

    x[0] is the observed/system coordinate Q with an onsite spring.
    x[1:] are bath masses. There is no phenomenological damping.

    Interpret measured cycles only before the reflected wave from the far end
    returns to x[0].
    """

    if n_masses < 3:
        raise ValueError("n_masses must be at least 3")
    if omega <= 0:
        raise ValueError("omega must be positive")
    if dt <= 0:
        raise ValueError("dt must be positive")

    M = params.system_mass
    m = params.bath_mass
    k = params.bath_spring
    k0 = params.system_spring

    masses = np.full(n_masses, m, dtype=float)
    masses[0] = M

    period = 2.0 * math.pi / omega
    n_steps = int(n_periods * period / dt)

    x = np.zeros(n_masses, dtype=float)
    v = np.zeros(n_masses, dtype=float)

    def envelope(t: float) -> float:
        t_ramp = ramp_periods * period
        if t >= t_ramp:
            return 1.0
        return 0.5 * (1.0 - math.cos(math.pi * t / t_ramp))

    def external_force(t: float) -> float:
        return force_amplitude * envelope(t) * math.sin(omega * t)

    def forces(state: np.ndarray, t: float) -> np.ndarray:
        f = np.empty_like(state)
        f[0] = -k0 * state[0] + k * (state[1] - state[0]) + external_force(t)
        f[1:-1] = k * (state[2:] - 2.0 * state[1:-1] + state[:-2])
        f[-1] = k * (state[-2] - state[-1])
        return f

    time = np.empty(n_steps + 1)
    q = np.empty(n_steps + 1)
    qdot = np.empty(n_steps + 1)
    fext = np.empty(n_steps + 1)
    energy = np.empty(n_steps + 1)
    work = np.empty(n_steps + 1)

    f = forces(x, 0.0)
    time[0] = 0.0
    q[0] = x[0]
    qdot[0] = v[0]
    fext[0] = external_force(0.0)
    energy[0] = 0.0
    work[0] = 0.0

    for step in range(n_steps):
        t = step * dt

        v += 0.5 * dt * f / masses
        x += dt * v

        f_new = forces(x, t + dt)
        v += 0.5 * dt * f_new / masses

        f_new_ext = external_force(t + dt)
        work[step + 1] = work[step] + 0.5 * dt * (
            fext[step] * qdot[step] + f_new_ext * v[0]
        )

        f = f_new
        time[step + 1] = t + dt
        q[step + 1] = x[0]
        qdot[step + 1] = v[0]
        fext[step + 1] = f_new_ext

        kinetic = 0.5 * float(np.dot(masses, v * v))
        onsite = 0.5 * k0 * x[0] ** 2
        springs = 0.5 * k * float(np.dot(np.diff(x), np.diff(x)))
        energy[step + 1] = kinetic + onsite + springs

    return {
        "time": time,
        "q": q,
        "qdot": qdot,
        "force": fext,
        "energy": energy,
        "work": work,
        "period": period,
    }


def cycle_loop_areas(
    result: Dict[str, np.ndarray | float],
    *,
    first_cycle: int,
    last_cycle_exclusive: int,
) -> np.ndarray:
    """Compute integral F dQ for selected cycles."""

    time = np.asarray(result["time"])
    qdot = np.asarray(result["qdot"])
    force = np.asarray(result["force"])
    period = float(result["period"])

    areas = []
    for cycle in range(first_cycle, last_cycle_exclusive):
        mask = (time >= cycle * period) & (time <= (cycle + 1) * period)
        areas.append(_trapz(force[mask] * qdot[mask], time[mask]))

    return np.asarray(areas)


def reference_run() -> Dict[str, float]:
    """Run the repository reference case used in the research note."""

    params = RubinParams()
    analytic = analytic_response(0.5, 0.1, params)
    numeric = simulate_finite_chain(params=params)
    areas = cycle_loop_areas(
        numeric,
        first_cycle=10,
        last_cycle_exclusive=50,
    )

    e_final = float(np.asarray(numeric["energy"])[-1])
    w_final = float(np.asarray(numeric["work"])[-1])
    energy_rel_error = abs(e_final - w_final) / max(abs(w_final), 1e-30)

    area_mean = float(np.mean(areas))

    return {
        **analytic,
        "numeric_loop_area_mean": area_mean,
        "numeric_loop_area_std": float(np.std(areas)),
        "loop_area_relative_error": abs(area_mean - analytic["loop_area"])
        / analytic["loop_area"],
        "final_internal_energy": e_final,
        "final_external_work": w_final,
        "energy_balance_relative_error": energy_rel_error,
    }
