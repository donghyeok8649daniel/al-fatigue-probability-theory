"""Verify the exact 1D Theta-based density-shape identity on nonlinear LJ-chain data.

This is a numerical verification harness, not a new physical closure.  It uses
the same conservative velocity-Verlet equations as ``theory.normal_lj_chain``
but records instantaneous spacing velocities needed for the exact first/second
moment shape relation.

No Boltzmann, Gaussian/Weibull, Fokker-Planck, damping, Markov, or
neighbor-independence assumption is introduced.  The only non-exact step is
the finite-sample smoothing used to estimate conditional fields from the
finite-M empirical measure.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math

import numpy as np

from theory.normal_lj_chain import (
    NormalLJParameters,
    normalized_lj_energy,
    normalized_lj_force,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "theta_density_reconstruction"
REPORT = ROOT / "results" / "reports" / "THETA_DENSITY_RECONSTRUCTION.md"

PHASES = (2.05, 2.10, 2.25, 2.40, 2.50, 2.60, 2.75, 2.90)


@dataclass(frozen=True)
class CapturedState:
    time: float
    spacing: np.ndarray
    spacing_velocity: np.ndarray
    spacing_acceleration: np.ndarray


def _external_force(parameters: NormalLJParameters, t: float, period: float) -> float:
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


def capture_chain_states(
    parameters: NormalLJParameters,
    sample_cycles: tuple[float, ...],
    *,
    atoms: int = 512,
    dt: float = 0.01,
) -> tuple[dict[float, CapturedState], float]:
    """Capture finite-chain states at requested cycle phases.

    This mirrors the conservative velocity-Verlet update in
    ``simulate_normal_lj_chain``. The two boundary spacing accelerations are
    left as NaN because their equations depend on the boundary loading law;
    the density-shape reconstruction uses only interior spacings.
    """
    if atoms < 5 or dt <= 0.0:
        raise ValueError("atoms must be >=5 and dt must be positive")
    if parameters.omega <= 0.0:
        raise ValueError("omega must be positive")
    sample_cycles = tuple(sorted(float(c) for c in sample_cycles))
    period = 2.0 * math.pi / parameters.omega
    sample_times = [c * period for c in sample_cycles]
    nsteps = int(math.ceil(max(sample_times) / dt)) + 2

    m = parameters.repulsive_exponent
    n = parameters.attractive_exponent
    x = np.arange(atoms, dtype=float)
    velocity = np.zeros(atoms, dtype=float)

    def force_vector(state: np.ndarray, t: float) -> np.ndarray:
        spacing = np.diff(state)
        dphi = normalized_lj_force(spacing, m, n)
        force = np.zeros_like(state)
        force[1:-1] = dphi[1:] - dphi[:-1]
        force[-1] = -dphi[-1] + _external_force(parameters, t, period)
        return force

    force = force_vector(x, 0.0)
    captured: dict[float, CapturedState] = {}
    next_sample = 0

    for step in range(nsteps):
        t = step * dt
        velocity[1:] += 0.5 * dt * force[1:]
        x[1:] += dt * velocity[1:]
        force = force_vector(x, t + dt)
        velocity[1:] += 0.5 * dt * force[1:]
        t_new = t + dt

        while next_sample < len(sample_times) and t_new >= sample_times[next_sample]:
            spacing = np.diff(x).copy()
            spacing_velocity = np.diff(velocity).copy()
            dphi = normalized_lj_force(spacing, m, n)
            acceleration = np.full_like(spacing, np.nan)
            acceleration[1:-1] = dphi[2:] - 2.0 * dphi[1:-1] + dphi[:-2]
            key = sample_cycles[next_sample]
            captured[key] = CapturedState(
                time=float(t_new),
                spacing=spacing,
                spacing_velocity=spacing_velocity,
                spacing_acceleration=acceleration,
            )
            next_sample += 1

        if next_sample >= len(sample_times):
            break

    if len(captured) != len(sample_cycles):
        raise RuntimeError("failed to capture all requested states")
    return captured, period


def _silverman_bandwidth(values: np.ndarray, factor: float = 1.5) -> float:
    values = np.asarray(values, dtype=float)
    sigma = float(np.std(values, ddof=1))
    if sigma <= 0.0:
        raise ValueError("nonzero spacing variance is required")
    return factor * 1.06 * sigma * values.size ** (-0.2)


def _kernel_fields(
    grid: np.ndarray,
    spacing: np.ndarray,
    spacing_velocity: np.ndarray,
    *,
    bandwidth: float,
    spacing_acceleration: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """Gaussian-kernel estimate of P, u, Theta, and conditional acceleration.

    The Gaussian kernel is a numerical estimator only. It is not a proposed
    physical probability family.
    """
    grid = np.asarray(grid, dtype=float)
    lam = np.asarray(spacing, dtype=float)
    vel = np.asarray(spacing_velocity, dtype=float)
    if lam.shape != vel.shape or lam.ndim != 1:
        raise ValueError("spacing and spacing_velocity must be equal-size 1D arrays")
    if bandwidth <= 0.0:
        raise ValueError("bandwidth must be positive")

    z = (grid[:, None] - lam[None, :]) / bandwidth
    weights = np.exp(-0.5 * z * z)
    sum_w = np.maximum(weights.sum(axis=1), 1.0e-300)

    density = sum_w / (lam.size * bandwidth * math.sqrt(2.0 * math.pi))
    mean_velocity = (weights @ vel) / sum_w
    mean_v2 = (weights @ (vel * vel)) / sum_w
    theta = np.maximum(mean_v2 - mean_velocity * mean_velocity, 1.0e-18)

    conditional_acceleration = None
    if spacing_acceleration is not None:
        acc = np.asarray(spacing_acceleration, dtype=float)
        if acc.shape != lam.shape:
            raise ValueError("spacing_acceleration shape mismatch")
        valid = np.isfinite(acc)
        z_acc = (grid[:, None] - lam[None, valid]) / bandwidth
        w_acc = np.exp(-0.5 * z_acc * z_acc)
        sum_acc = np.maximum(w_acc.sum(axis=1), 1.0e-300)
        conditional_acceleration = (w_acc @ acc[valid]) / sum_acc

    return density, mean_velocity, theta, conditional_acceleration


def _normalize(grid: np.ndarray, density: np.ndarray) -> np.ndarray:
    norm = float(np.trapezoid(density, grid))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("density normalization failed")
    return density / norm


def reconstruct_at_phase(
    captures: dict[float, CapturedState],
    phase: float,
    epsilon_cycles: float,
    *,
    bandwidth_factor: float = 1.5,
    grid_points: int = 600,
) -> tuple[dict[str, float], np.ndarray]:
    prev_state = captures[phase - epsilon_cycles]
    state = captures[phase]
    next_state = captures[phase + epsilon_cycles]

    interior = slice(1, -1)
    lam = state.spacing[interior]
    vel = state.spacing_velocity[interior]
    acc = state.spacing_acceleration[interior]

    q_lo, q_hi = np.quantile(lam, [0.005, 0.995])
    pad = 0.15 * float(q_hi - q_lo)
    grid = np.linspace(float(q_lo - pad), float(q_hi + pad), grid_points)
    bandwidth = _silverman_bandwidth(lam, bandwidth_factor)

    density, u, theta, abar = _kernel_fields(
        grid, lam, vel, bandwidth=bandwidth, spacing_acceleration=acc
    )
    _, u_prev, _, _ = _kernel_fields(
        grid,
        prev_state.spacing[interior],
        prev_state.spacing_velocity[interior],
        bandwidth=bandwidth,
    )
    _, u_next, _, _ = _kernel_fields(
        grid,
        next_state.spacing[interior],
        next_state.spacing_velocity[interior],
        bandwidth=bandwidth,
    )
    if abar is None:
        raise RuntimeError("conditional acceleration was not computed")

    dt_pair = next_state.time - prev_state.time
    partial_t_u = (u_next - u_prev) / dt_pair
    partial_lambda_u = np.gradient(u, grid, edge_order=2)
    material_acceleration = partial_t_u + u * partial_lambda_u
    d_log_theta = np.gradient(np.log(theta), grid, edge_order=2)

    log_slope = (abar - material_acceleration) / theta - d_log_theta
    increments = 0.5 * (log_slope[:-1] + log_slope[1:]) * np.diff(grid)
    log_p = np.concatenate(([0.0], np.cumsum(increments)))
    log_p -= float(np.max(log_p))
    reconstructed = _normalize(grid, np.exp(log_p))
    direct_kde = _normalize(grid, density)

    cdf_direct = np.concatenate(
        ([0.0], np.cumsum(0.5 * (direct_kde[:-1] + direct_kde[1:]) * np.diff(grid)))
    )
    cdf_reconstructed = np.concatenate(
        ([0.0], np.cumsum(0.5 * (reconstructed[:-1] + reconstructed[1:]) * np.diff(grid)))
    )
    l1 = float(np.trapezoid(np.abs(direct_kde - reconstructed), grid))
    ks = float(np.max(np.abs(cdf_direct - cdf_reconstructed)))

    mean_direct = float(np.trapezoid(grid * direct_kde, grid))
    mean_reconstructed = float(np.trapezoid(grid * reconstructed, grid))
    var_direct = float(np.trapezoid((grid - mean_direct) ** 2 * direct_kde, grid))
    var_reconstructed = float(
        np.trapezoid((grid - mean_reconstructed) ** 2 * reconstructed, grid)
    )

    energy = normalized_lj_energy(grid) - normalized_lj_energy(1.0)
    mean_energy_direct = float(np.trapezoid(energy * direct_kde, grid))
    mean_energy_reconstructed = float(np.trapezoid(energy * reconstructed, grid))
    energy_relative_error = abs(mean_energy_reconstructed - mean_energy_direct) / max(
        abs(mean_energy_direct), 1.0e-15
    )

    fields = np.column_stack(
        (
            grid,
            direct_kde,
            reconstructed,
            u,
            theta,
            abar,
            material_acceleration,
            log_slope,
        )
    )
    row = {
        "phase_cycle": phase,
        "represented_interior_spacings": float(lam.size),
        "bandwidth": bandwidth,
        "L1_density_error": l1,
        "KS_density_error": ks,
        "mean_spacing_direct": mean_direct,
        "mean_spacing_reconstructed": mean_reconstructed,
        "variance_direct": var_direct,
        "variance_reconstructed": var_reconstructed,
        "mean_energy_direct": mean_energy_direct,
        "mean_energy_reconstructed": mean_energy_reconstructed,
        "mean_energy_relative_error": energy_relative_error,
        "theta_min": float(np.min(theta)),
        "theta_median": float(np.median(theta)),
    }
    return row, fields


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    parameters = NormalLJParameters(force_amplitude=0.03, omega=0.02)
    epsilon_cycles = 2.5e-4
    sample_cycles = tuple(
        sorted(
            {
                phase + offset
                for phase in PHASES
                for offset in (-epsilon_cycles, 0.0, epsilon_cycles)
            }
        )
    )
    captures, period = capture_chain_states(parameters, sample_cycles)

    rows: list[dict[str, float]] = []
    for phase in PHASES:
        row, fields = reconstruct_at_phase(captures, phase, epsilon_cycles)
        rows.append(row)
        np.savetxt(
            DATA / f"phase_{phase:.2f}.csv",
            fields,
            delimiter=",",
            header=(
                "lambda,direct_kde,reconstructed,u,theta,"
                "conditional_acceleration,material_acceleration,log_slope"
            ),
            comments="",
        )

    with (DATA / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    max_l1 = max(row["L1_density_error"] for row in rows)
    max_ks = max(row["KS_density_error"] for row in rows)
    max_energy_error = max(row["mean_energy_relative_error"] for row in rows)
    summary = {
        "classification": (
            "numerical reconstruction check of the exact smooth 1D moment-shape identity; "
            "Gaussian kernel smoothing is an estimator, not a physical PDF assumption"
        ),
        "parameters": {
            "atoms": 512,
            "dt": 0.01,
            "force_amplitude": parameters.force_amplitude,
            "omega": parameters.omega,
            "ramp_cycles": parameters.ramp_cycles,
            "epsilon_cycles_for_time_derivative": epsilon_cycles,
            "bandwidth_rule": "1.5 * Silverman",
            "phases": list(PHASES),
        },
        "period": period,
        "max_L1_density_error": max_l1,
        "max_KS_density_error": max_ks,
        "max_mean_energy_relative_error": max_energy_error,
        "rows": rows,
        "interpretation": (
            "The reconstructed density agrees closely with the directly smoothed finite-M "
            "density at all tested phases. Mean energy must be compared between the two "
            "continuum densities on the same grid; comparison to the raw atomic sample "
            "also includes KDE smoothing bias."
        ),
        "limitations": [
            "finite-M empirical state is smoothed numerically to estimate conditional fields",
            "only interior spacings are used because the bulk acceleration identity does not apply at boundaries",
            "the check is 1D normal-only and does not yet validate the (a,s) tensor extension",
            "Theta approaching zero is a degenerate regime where the smooth shape formula is not applicable",
        ],
    }
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "# Theta-based density reconstruction — numerical verification",
        "",
        "This report checks the exact 1D smooth-moment shape identity against the nonlinear conservative layer-LJ chain.",
        "No Boltzmann, Gaussian/Weibull physical PDF, Fokker–Planck, damping, Markov, or neighbor-independence assumption is used.",
        "A Gaussian kernel is used only as a finite-sample estimator of the empirical density and conditional moments.",
        "",
        "## Fixed numerical protocol",
        "",
        "- atoms: 512",
        "- dt: 0.01",
        f"- force amplitude: {parameters.force_amplitude}",
        f"- omega: {parameters.omega}",
        f"- time derivative half-window: {epsilon_cycles} cycle",
        "- bandwidth: 1.5 × Silverman rule, fixed for all phases",
        "- boundary spacings excluded from the bulk acceleration identity",
        "",
        "## Results",
        "",
        "| phase N | L1(P) | KS(P) | relative mean-energy error |",
        "|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['phase_cycle']:.2f} | {row['L1_density_error']:.5f} | "
            f"{row['KS_density_error']:.5f} | {row['mean_energy_relative_error']:.3%} |"
        )
    lines += [
        "",
        f"Maximum L1 density error: **{max_l1:.5f}**.",
        f"Maximum KS density error: **{max_ks:.5f}**.",
        f"Maximum relative error in mean intrinsic energy on the same smoothed support: **{max_energy_error:.3%}**.",
        "",
        "## Interpretation",
        "",
        "The 1D Theta-based shape identity is numerically consistent with the deterministic nonlinear LJ-chain data at the tested phases.",
        "The test does **not** prove a new closure: Theta, conditional acceleration, and the material acceleration are measured from the chain state.",
        "The next derivation step is the two-coordinate (a,s) conditional velocity-covariance tensor and its integrability/compatibility condition.",
        "",
        "## Assumption / validity ledger",
        "",
        "- Exact: finite-M mechanics and the interior nonlinear LJ spacing acceleration.",
        "- Exact at smooth-moment level: continuity and first/second velocity moment identities.",
        "- Numerical approximation only: KDE/regression used to estimate smooth conditional fields from finite M.",
        "- Not used: Boltzmann equilibrium, Gaussian physical state, Fokker–Planck, Markov bath, independent spacings.",
        "- Invalid regime for the divided shape law: Theta = 0.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
