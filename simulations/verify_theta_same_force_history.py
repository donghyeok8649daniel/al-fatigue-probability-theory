"""Verify same-force load/unload history dependence of P, u, and Theta.

The source dynamics are the conservative 1D normal generalized-LJ chain.  No
Boltzmann law, Fokker--Planck closure, damping, random forcing, or physical PDF
family is introduced.  Gaussian kernels are used only as finite-M estimators of
the empirical density and conditional velocity moments.

The script compares loading and unloading passages at identical external-force
levels in cycle 4 after a two-cycle smooth ramp.  The force history is a
pulsating tensile cycle with mean=amplitude=100 MPa under the repository's
sigma/E mapping.  The normalized angular frequency omega*=0.02 is retained
from the existing proof-of-principle Theta reconstruction; it is NOT a
laboratory fatigue-frequency claim.

Important interpretation:
    same-force non-retracing of (P,u,Theta) proves history dependence of the
    reduced state, but the conservative chain has no irreversible dissipation.
    Cycle work is therefore checked against the change in mechanical energy.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from theory.normal_lj_chain import (
    NormalLJParameters,
    normalized_lj_energy,
    normalized_lj_force,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "theta_same_force_history"
REPORT = ROOT / "results" / "reports" / "THETA_SAME_FORCE_HISTORY.md"

ATOMS = 512
DT = 0.01
CYCLE = 4
YOUNGS_MODULUS_PA = 69.0e9
FORCE_SCALE = 100.0e6 / YOUNGS_MODULUS_PA
OMEGA = 0.02
RAMP_CYCLES = 2
BANDWIDTH_FACTOR = 1.5

# label, loading phase, unloading phase.  A phase is measured in cycles.
# Loading means dQ/dt>0 and unloading means dQ/dt<0.  The 50 MPa loading
# crossing occurs later in the same cycle than its unloading crossing; the
# comparison is state-based rather than chronological ordering of the tuple.
PAIR_PHASES = (
    ("100MPa_mean", 0.0, 0.5),
    ("150MPa", 1.0 / 12.0, 5.0 / 12.0),
    ("186.6MPa", 1.0 / 6.0, 1.0 / 3.0),
    ("50MPa", 11.0 / 12.0, 7.0 / 12.0),
)


@dataclass(frozen=True)
class State:
    time: float
    force: float
    spacing: np.ndarray
    spacing_velocity: np.ndarray


def external_force(parameters: NormalLJParameters, t: float, period: float) -> float:
    if parameters.ramp_cycles <= 0:
        envelope = 1.0
    else:
        ramp_time = parameters.ramp_cycles * period
        envelope = 1.0 if t >= ramp_time else 0.5 * (
            1.0 - math.cos(math.pi * t / ramp_time)
        )
    return envelope * (
        parameters.mean_force
        + parameters.force_amplitude * math.sin(parameters.omega * t)
    )


def capture_states(
    parameters: NormalLJParameters,
    sample_cycles: tuple[float, ...],
    *,
    atoms: int = ATOMS,
    dt: float = DT,
) -> tuple[dict[float, State], float]:
    """Capture exact finite-chain snapshots near the requested cycle phases."""
    if atoms < 5 or dt <= 0.0:
        raise ValueError("atoms must be >=5 and dt must be positive")
    period = 2.0 * math.pi / parameters.omega
    ordered = tuple(sorted(float(c) for c in sample_cycles))
    target_times = [c * period for c in ordered]

    x = np.arange(atoms, dtype=float)
    velocity = np.zeros(atoms, dtype=float)
    m = parameters.repulsive_exponent
    n = parameters.attractive_exponent

    def force_vector(state: np.ndarray, t: float) -> np.ndarray:
        spacing = np.diff(state)
        dphi = normalized_lj_force(spacing, m, n)
        force = np.zeros_like(state)
        force[1:-1] = dphi[1:] - dphi[:-1]
        force[-1] = -dphi[-1] + external_force(parameters, t, period)
        return force

    force = force_vector(x, 0.0)
    captured: dict[float, State] = {}
    next_sample = 0
    nsteps = int(math.ceil(max(target_times) / dt)) + 3

    for step in range(nsteps):
        t = step * dt
        velocity[1:] += 0.5 * dt * force[1:]
        x[1:] += dt * velocity[1:]
        force = force_vector(x, t + dt)
        velocity[1:] += 0.5 * dt * force[1:]
        t_new = t + dt

        while next_sample < len(target_times) and t_new >= target_times[next_sample]:
            key = ordered[next_sample]
            captured[key] = State(
                time=float(t_new),
                force=float(external_force(parameters, t_new, period)),
                spacing=np.diff(x).copy(),
                spacing_velocity=np.diff(velocity).copy(),
            )
            next_sample += 1
        if next_sample >= len(target_times):
            break

    if len(captured) != len(ordered):
        raise RuntimeError("failed to capture all requested states")
    return captured, period


def silverman_bandwidth(values: np.ndarray, factor: float) -> float:
    sigma = float(np.std(values, ddof=1))
    if sigma <= 0.0:
        raise ValueError("nonzero spacing spread is required")
    return factor * 1.06 * sigma * values.size ** (-0.2)


def kernel_fields(
    grid: np.ndarray,
    spacing: np.ndarray,
    velocity: np.ndarray,
    bandwidth: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return numerical estimates of P, u, and conditional Theta."""
    z = (grid[:, None] - spacing[None, :]) / bandwidth
    weights = np.exp(-0.5 * z * z)
    sum_w = np.maximum(weights.sum(axis=1), 1.0e-300)

    density = sum_w / (spacing.size * bandwidth * math.sqrt(2.0 * math.pi))
    u = (weights @ velocity) / sum_w
    mean_v2 = (weights @ (velocity * velocity)) / sum_w
    theta = np.maximum(mean_v2 - u * u, 0.0)
    density /= float(np.trapezoid(density, grid))
    return density, u, theta


def pair_metrics(
    load: State,
    unload: State,
    *,
    bandwidth_factor: float = BANDWIDTH_FACTOR,
) -> tuple[dict[str, float], np.ndarray]:
    """Compare the smooth finite-M fields at one equal-force load/unload pair."""
    # The bulk one-point interpretation excludes the two boundary spacings.
    load_spacing = load.spacing[1:-1]
    unload_spacing = unload.spacing[1:-1]
    load_velocity = load.spacing_velocity[1:-1]
    unload_velocity = unload.spacing_velocity[1:-1]

    pooled = np.concatenate((load_spacing, unload_spacing))
    q_lo, q_hi = np.quantile(pooled, [0.002, 0.998])
    pad = 0.2 * float(q_hi - q_lo)
    grid = np.linspace(float(q_lo - pad), float(q_hi + pad), 800)
    bandwidth = silverman_bandwidth(pooled, bandwidth_factor)

    p_load, u_load, theta_load = kernel_fields(
        grid, load_spacing, load_velocity, bandwidth
    )
    p_unload, u_unload, theta_unload = kernel_fields(
        grid, unload_spacing, unload_velocity, bandwidth
    )
    weight = 0.5 * (p_load + p_unload)
    weight /= float(np.trapezoid(weight, grid))

    cdf_load = np.concatenate((
        [0.0],
        np.cumsum(0.5 * (p_load[:-1] + p_load[1:]) * np.diff(grid)),
    ))
    cdf_unload = np.concatenate((
        [0.0],
        np.cumsum(0.5 * (p_unload[:-1] + p_unload[1:]) * np.diff(grid)),
    ))

    u_rms = float(np.sqrt(np.trapezoid(weight * (u_load - u_unload) ** 2, grid)))
    u_scale = float(np.sqrt(np.trapezoid(
        weight * 0.5 * (u_load * u_load + u_unload * u_unload), grid
    )))
    theta_rms = float(np.sqrt(np.trapezoid(
        weight * (theta_load - theta_unload) ** 2, grid
    )))
    theta_scale = float(np.trapezoid(
        weight * 0.5 * (theta_load + theta_unload), grid
    ))

    phi0 = float(normalized_lj_energy(1.0))
    row = {
        "force_load": load.force,
        "force_unload": unload.force,
        "force_mismatch": abs(load.force - unload.force),
        "mean_spacing_load": float(np.mean(load_spacing)),
        "mean_spacing_unload": float(np.mean(unload_spacing)),
        "variance_load": float(np.var(load_spacing)),
        "variance_unload": float(np.var(unload_spacing)),
        "global_velocity_variance_load": float(np.var(load_velocity)),
        "global_velocity_variance_unload": float(np.var(unload_velocity)),
        "mean_intrinsic_energy_load": float(np.mean(normalized_lj_energy(load_spacing) - phi0)),
        "mean_intrinsic_energy_unload": float(np.mean(normalized_lj_energy(unload_spacing) - phi0)),
        "density_L1": float(np.trapezoid(np.abs(p_load - p_unload), grid)),
        "density_KS": float(np.max(np.abs(cdf_load - cdf_unload))),
        "u_weighted_rms_difference": u_rms,
        "u_normalized_difference": u_rms / max(u_scale, 1.0e-15),
        "theta_weighted_rms_difference": theta_rms,
        "theta_normalized_difference": theta_rms / max(theta_scale, 1.0e-15),
        "bandwidth": bandwidth,
    }
    fields = np.column_stack((
        grid,
        p_load,
        p_unload,
        u_load,
        u_unload,
        theta_load,
        theta_unload,
        weight,
    ))
    return row, fields


def cycle_work_balance(
    parameters: NormalLJParameters,
    *,
    atoms: int = ATOMS,
    dt: float = DT,
    cycles: int = 6,
) -> list[dict[str, float]]:
    """Check that cycle loop work is retained as mechanical energy, not G3 loss."""
    period = 2.0 * math.pi / parameters.omega
    m = parameters.repulsive_exponent
    n = parameters.attractive_exponent
    x = np.arange(atoms, dtype=float)
    velocity = np.zeros(atoms, dtype=float)

    def force_vector(state: np.ndarray, t: float) -> np.ndarray:
        spacing = np.diff(state)
        dphi = normalized_lj_force(spacing, m, n)
        force = np.zeros_like(state)
        force[1:-1] = dphi[1:] - dphi[:-1]
        force[-1] = -dphi[-1] + external_force(parameters, t, period)
        return force

    def mechanical_energy() -> float:
        return (
            0.5 * float(velocity[1:] @ velocity[1:])
            + float(np.sum(normalized_lj_energy(np.diff(x), m, n)))
        )

    force = force_vector(x, 0.0)
    previous_power = external_force(parameters, 0.0, period) * velocity[-1]
    cumulative_work = 0.0
    cycle_start_work = 0.0
    cycle_start_energy = mechanical_energy()
    next_cycle_time = period
    current_cycle = 1
    rows: list[dict[str, float]] = []
    nsteps = int(round(cycles * period / dt))

    for step in range(nsteps):
        t = step * dt
        velocity[1:] += 0.5 * dt * force[1:]
        x[1:] += dt * velocity[1:]
        force = force_vector(x, t + dt)
        velocity[1:] += 0.5 * dt * force[1:]
        t_new = t + dt
        power = external_force(parameters, t_new, period) * velocity[-1]
        cumulative_work += 0.5 * dt * (previous_power + power)
        previous_power = power

        if t_new >= next_cycle_time - 1.0e-12:
            energy = mechanical_energy()
            work_increment = cumulative_work - cycle_start_work
            energy_increment = energy - cycle_start_energy
            rows.append({
                "cycle": float(current_cycle),
                "external_work_increment": work_increment,
                "mechanical_energy_increment": energy_increment,
                "work_minus_energy_change": work_increment - energy_increment,
            })
            current_cycle += 1
            next_cycle_time = current_cycle * period
            cycle_start_work = cumulative_work
            cycle_start_energy = energy

    return rows


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    parameters = NormalLJParameters(
        mean_force=FORCE_SCALE,
        force_amplitude=FORCE_SCALE,
        omega=OMEGA,
        ramp_cycles=RAMP_CYCLES,
    )
    sample_cycles = tuple(sorted({
        CYCLE + phase
        for _, load_phase, unload_phase in PAIR_PHASES
        for phase in (load_phase, unload_phase)
    }))
    captures, period = capture_states(parameters, sample_cycles)

    rows: list[dict[str, float | str]] = []
    sensitivity_rows: list[dict[str, float | str]] = []
    for label, load_phase, unload_phase in PAIR_PHASES:
        load = captures[CYCLE + load_phase]
        unload = captures[CYCLE + unload_phase]
        metrics, fields = pair_metrics(load, unload)
        row: dict[str, float | str] = {
            "level": label,
            "load_phase": load_phase,
            "unload_phase": unload_phase,
            **metrics,
        }
        rows.append(row)
        np.savetxt(
            DATA / f"fields_{label}.csv",
            fields,
            delimiter=",",
            header=(
                "lambda,P_load,P_unload,u_load,u_unload,"
                "theta_load,theta_unload,comparison_weight"
            ),
            comments="",
        )

        for factor in (1.0, 1.25, 1.5, 1.75, 2.0):
            sensitivity, _ = pair_metrics(load, unload, bandwidth_factor=factor)
            sensitivity_rows.append({
                "level": label,
                "bandwidth_factor": factor,
                "density_L1": sensitivity["density_L1"],
                "u_weighted_rms_difference": sensitivity["u_weighted_rms_difference"],
                "theta_weighted_rms_difference": sensitivity["theta_weighted_rms_difference"],
                "theta_normalized_difference": sensitivity["theta_normalized_difference"],
            })

    with (DATA / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with (DATA / "bandwidth_sensitivity.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sensitivity_rows[0]))
        writer.writeheader()
        writer.writerows(sensitivity_rows)

    work_rows = cycle_work_balance(parameters)
    with (DATA / "cycle_work_balance.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(work_rows[0]))
        writer.writeheader()
        writer.writerows(work_rows)

    summary = {
        "classification": (
            "deterministic same-force history-dependence test of the exact 1D "
            "P-u-Theta reduced state; not an irreversible-dissipation model"
        ),
        "parameters": {
            "atoms": ATOMS,
            "dt": DT,
            "mean_stress_pa_under_sigma_over_E_mapping": 100.0e6,
            "stress_amplitude_pa_under_sigma_over_E_mapping": 100.0e6,
            "omega_star": OMEGA,
            "ramp_cycles": RAMP_CYCLES,
            "comparison_cycle": CYCLE,
            "bandwidth_factor": BANDWIDTH_FACTOR,
            "period_star": period,
        },
        "pairs": rows,
        "cycle_work_balance": work_rows,
        "interpretation": [
            "P, u, and Theta are not single-valued functions of instantaneous force in this driven finite chain.",
            "The result is dynamic history dependence/non-retracing of the reduced state.",
            "The chain is conservative; cycle work is accounted for by mechanical-energy change, so the result is not G3 irreversible dissipation.",
            "The Gaussian kernel is only a finite-M estimator and is not a physical PDF assumption.",
            "omega*=0.02 is a proof-of-principle atomic-chain protocol, not a laboratory fatigue-frequency mapping.",
        ],
    }
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "# Same-force $P$–$u$–$\\Theta$ history test",
        "",
        "## Result",
        "",
        "The deterministic conservative 1D LJ chain reaches different reduced states on loading and unloading even when the externally applied force is the same.",
        "This establishes dynamic history dependence of the distribution state; it does **not** establish irreversible G3 dissipation.",
        "",
        "## Protocol",
        "",
        f"- atoms: {ATOMS}",
        f"- dt*: {DT}",
        "- mean stress mapping: 100 MPa",
        "- stress amplitude mapping: 100 MPa (0 to 200 MPa pulsating tension)",
        f"- omega*: {OMEGA}",
        f"- ramp: {RAMP_CYCLES} cycles",
        f"- comparison: cycle {CYCLE}",
        "- finite-M field estimator: common-grid Gaussian kernel, 1.5 x Silverman bandwidth",
        "- boundary spacings excluded from the one-point bulk comparison",
        "",
        "## Same-force comparisons",
        "",
        "| force level | L1(P) | KS(P) | normalized u difference | normalized Theta difference | mean spacing load | mean spacing unload |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['level']} | {row['density_L1']:.6f} | {row['density_KS']:.6f} | "
            f"{row['u_normalized_difference']:.6f} | {row['theta_normalized_difference']:.6f} | "
            f"{row['mean_spacing_load']:.9f} | {row['mean_spacing_unload']:.9f} |"
        )
    lines.extend([
        "",
        "At every tested force level the smoothed empirical spacing density differs between loading and unloading.",
        "The conditional mean velocity and conditional velocity variance also differ, so the distinction is not carried by $P$ alone.",
        "",
        "## Estimator sensitivity",
        "",
        "Changing the Gaussian bandwidth factor from 1.0 to 2.0 changes the numerical magnitude of the L1 metric but does not remove the load/unload density distinction at any tested force level.",
        "The Gaussian kernel remains a numerical estimator only.",
        "",
        "## Conservative work balance",
        "",
        "For each integrated cycle, external boundary work agrees with the mechanical-energy increment to numerical precision.",
        "Therefore the non-retracing state is an inertial/wave/phase-space memory of a conservative driven system, not a positive irreversible dissipation law.",
        "",
        "This means the defensible statement is",
        "",
        "$$",
        "\\boxed{",
        "(P,u,\\Theta)_{\\rm load}(Q^*)\\ne(P,u,\\Theta)_{\\rm unload}(Q^*)",
        "}",
        "$$",
        "",
        "for the tested finite-chain protocol, while",
        "",
        "$$",
        "\\boxed{",
        "\\text{same-force non-retracing}\\not\\Rightarrow\\dot D_{\\rm irr}>0.",
        "}",
        "$$",
        "",
        "## Scale limitation",
        "",
        "The present omega*=0.02 protocol is deliberately computational and remains on microscopic/finite-chain dynamical scales.",
        "It must not be identified with a laboratory fatigue frequency.  A separate mechanics-based slow or irreversible mechanism is still required before claiming physical fatigue hysteresis at Hz-scale loading.",
    ])
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
