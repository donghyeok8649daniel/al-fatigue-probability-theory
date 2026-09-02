"""Run the deterministic reduced (a_i,s_i) spatial chain from theory-core.

This script is intentionally self-contained because the pruned numerical-fem
branch does not own the theory package. It mirrors the Milestone-21 reduced
model

    V_M = sum_i U0(a_i,s_i),  a_i = x_{i+1}-x_i,

using the current multiplicity-free row/layer generalized-LJ energy

    U0(a,s)=sum_{k>=1}sum_{p in Z} v(sqrt((k a)^2+(p b+s)^2)).

The (k,p) sums are numerically truncated; no physical interaction cutoff,
random forcing, Gaussian initial spread, damping, Boltzmann law, or
Fokker--Planck closure is introduced.

Two diagnostics are run:

1. normal_only: q_s(t)=0. This tests whether boundary-driven axial diversity
   alone can move a system started exactly at the symmetric registry well.
2. registry_driven: a declared uniform oscillatory registry force displaces the
   state from registry symmetry so the existing mixed U0(a,s) coupling can
   transfer axial spatial diversity into s. This second case is a normalized
   mechanism diagnostic, not a calibrated Al loading law.

The reported Theta entries are GLOBAL finite-cell velocity covariances. They
must not be confused with the conditional field Theta(a,s,t) used in the
exact density-shape identity.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

M_EXP = 12.0
N_EXP = 6.0
B = 1.0
EPSILON = 1.0
SIGMA_LJ = 1.0
KMAX = 12
PMAX = 30
CELLS = 48
DT = 0.004
OMEGA = 0.35
NORMAL_MEAN = 0.2
NORMAL_AMPLITUDE = 0.2
REGISTRY_AMPLITUDE = 0.8
RAMP_CYCLES = 1

C_MN = M_EXP / (M_EXP - N_EXP) * (M_EXP / N_EXP) ** (N_EXP / (M_EXP - N_EXP))
_K = np.arange(1, KMAX + 1, dtype=float)[:, None, None]
_P = np.arange(-PMAX, PMAX + 1, dtype=float)[None, None, :]


def lattice_fields(a: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return truncated direct (U0, dU/da, dU/ds) for equal-size cell arrays."""
    aa = np.asarray(a, dtype=float)
    ss = np.asarray(s, dtype=float)
    if aa.ndim != 1 or ss.shape != aa.shape or np.any(aa <= 0.0):
        raise ValueError("a and s must be equal-size 1D arrays and a>0")

    x = _K * aa[None, :, None]
    y = _P * B + ss[None, :, None]
    r2 = x * x + y * y
    if np.any(r2 <= 0.0):
        raise ValueError("singular pair distance")

    rep = (SIGMA_LJ * SIGMA_LJ / r2) ** (0.5 * M_EXP)
    att = (SIGMA_LJ * SIGMA_LJ / r2) ** (0.5 * N_EXP)
    energy = C_MN * EPSILON * np.sum(rep - att, axis=(0, 2))

    radial = C_MN * EPSILON * (
        -M_EXP * SIGMA_LJ**M_EXP * r2 ** (-0.5 * (M_EXP + 2.0))
        + N_EXP * SIGMA_LJ**N_EXP * r2 ** (-0.5 * (N_EXP + 2.0))
    )
    d_u_da = np.sum(radial * (x * _K), axis=(0, 2))
    d_u_ds = np.sum(radial * y, axis=(0, 2))
    return energy, d_u_da, d_u_ds


def scalar_u(a: float, s: float) -> float:
    return float(lattice_fields(np.array([a]), np.array([s]))[0][0])


def golden_minimum(function, lower: float, upper: float, tolerance: float = 1.0e-11) -> float:
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = upper - ratio * (upper - lower)
    x2 = lower + ratio * (upper - lower)
    f1, f2 = function(x1), function(x2)
    for _ in range(250):
        if abs(upper - lower) <= tolerance:
            break
        if f1 > f2:
            lower, x1, f1 = x1, x2, f2
            x2 = lower + ratio * (upper - lower)
            f2 = function(x2)
        else:
            upper, x2, f2 = x2, x1, f1
            x1 = upper - ratio * (upper - lower)
            f1 = function(x1)
    return 0.5 * (lower + upper)


def reference_state() -> tuple[float, float, float]:
    # In the present p+s phase convention, s/b=1/2 is the stable symmetric
    # registry well. One may equivalently shift the registry origin by b/2.
    s0 = 0.5 * B
    a0 = golden_minimum(lambda a: scalar_u(a, s0), 0.85, 1.15)
    return a0, s0, scalar_u(a0, s0)


def run_case(name: str, cycles: int, registry_amplitude: float) -> dict:
    a0, s0, u_ref = reference_state()
    period = 2.0 * math.pi / OMEGA
    steps = int(round(cycles * period / DT))

    x = np.arange(CELLS + 1, dtype=float) * a0
    vx = np.zeros(CELLS + 1)
    s = np.full(CELLS, s0)
    vs = np.zeros(CELLS)

    def envelope(t: float) -> float:
        ramp_time = RAMP_CYCLES * period
        if t >= ramp_time:
            return 1.0
        return 0.5 * (1.0 - math.cos(math.pi * t / ramp_time))

    def q_a(t: float) -> float:
        return envelope(t) * (NORMAL_MEAN + NORMAL_AMPLITUDE * math.sin(OMEGA * t))

    def q_s(t: float) -> float:
        return envelope(t) * registry_amplitude * math.sin(OMEGA * t)

    def accelerations(t: float):
        spacing = np.diff(x)
        energy, grad_a, grad_s = lattice_fields(spacing, s)
        ax = np.zeros(CELLS + 1)
        ax[1:CELLS] = grad_a[1:] - grad_a[:-1]
        ax[CELLS] = -grad_a[-1] + q_a(t)
        ass = q_s(t) - grad_s
        return ax, ass, energy

    def row(t: float, energy: np.ndarray, work: float) -> tuple[dict, tuple[np.ndarray, ...]]:
        a = np.diff(x)
        va = np.diff(vx)
        ca = a - np.mean(a)
        cs = s - np.mean(s)
        cva = va - np.mean(va)
        cvs = vs - np.mean(vs)
        mechanical_energy = (
            0.5 * float(vx[1:] @ vx[1:])
            + 0.5 * float(vs @ vs)
            + float(np.sum(energy))
        )
        values = {
            "t": t,
            "cycle": t / period,
            "Q_a": q_a(t),
            "Q_s": q_s(t),
            "mean_a": float(np.mean(a)),
            "var_a": float(np.mean(ca * ca)),
            "mean_s": float(np.mean(s)),
            "var_s": float(np.mean(cs * cs)),
            "cov_as": float(np.mean(ca * cs)),
            "mean_delta_U": float(np.mean(energy - u_ref)),
            "theta_aa_global": float(np.mean(cva * cva)),
            "theta_as_global": float(np.mean(cva * cvs)),
            "theta_ss_global": float(np.mean(cvs * cvs)),
            "a_min": float(np.min(a)),
            "a_max": float(np.max(a)),
            "s_min": float(np.min(s)),
            "s_max": float(np.max(s)),
            "mechanical_energy": mechanical_energy,
            "external_work": work,
        }
        return values, (a.copy(), s.copy(), va.copy(), vs.copy())

    ax, ass, energy = accelerations(0.0)
    initial_energy = 0.5 * float(vx[1:] @ vx[1:]) + 0.5 * float(vs @ vs) + float(np.sum(energy))
    work = 0.0
    previous_power = 0.0
    sample_stride = max(1, int(round(period / DT / 100.0)))

    history: list[dict] = []
    first, state = row(0.0, energy, work)
    history.append(first)
    peak_var_s = (first["var_s"], first["t"], state)

    for step in range(steps):
        t = step * DT
        vx[1:] += 0.5 * DT * ax[1:]
        vs += 0.5 * DT * ass
        x[1:] += DT * vx[1:]
        s += DT * vs

        ax_new, ass_new, energy = accelerations(t + DT)
        vx[1:] += 0.5 * DT * ax_new[1:]
        vs += 0.5 * DT * ass_new
        ax, ass = ax_new, ass_new

        power = q_a(t + DT) * vx[-1] + q_s(t + DT) * float(np.sum(vs))
        work += 0.5 * DT * (previous_power + power)
        previous_power = power

        if (step + 1) % sample_stride == 0 or step == steps - 1:
            values, state = row(t + DT, energy, work)
            history.append(values)
            if values["var_s"] > peak_var_s[0]:
                peak_var_s = (values["var_s"], values["t"], state)

    final_energy = history[-1]["mechanical_energy"]
    energy_residual = abs((final_energy - initial_energy) - work) / max(abs(work), 1.0e-14)

    def peak(key: str, absolute: bool = False):
        candidate = max(history, key=lambda r: abs(r[key]) if absolute else r[key])
        return (abs(candidate[key]) if absolute else candidate[key], candidate["cycle"])

    peak_a = peak("var_a")
    peak_s = peak("var_s")
    peak_taa = peak("theta_aa_global")
    peak_tas = peak("theta_as_global", absolute=True)
    peak_tss = peak("theta_ss_global")
    peak_u = peak("mean_delta_U")

    a_peak, s_peak, va_peak, vs_peak = peak_var_s[2]
    return {
        "summary": {
            "case": name,
            "cells": CELLS,
            "cycles": cycles,
            "dt": DT,
            "omega": OMEGA,
            "normal_mean_force": NORMAL_MEAN,
            "normal_force_amplitude": NORMAL_AMPLITUDE,
            "registry_force_amplitude": registry_amplitude,
            "equilibrium_a0": a0,
            "equilibrium_s0": s0,
            "kmax": KMAX,
            "pmax": PMAX,
            "peak_var_a": peak_a[0],
            "peak_var_a_cycle": peak_a[1],
            "peak_var_s": peak_s[0],
            "peak_var_s_cycle": peak_s[1],
            "peak_theta_aa_global": peak_taa[0],
            "peak_theta_aa_cycle": peak_taa[1],
            "peak_abs_theta_as_global": peak_tas[0],
            "peak_abs_theta_as_cycle": peak_tas[1],
            "peak_theta_ss_global": peak_tss[0],
            "peak_theta_ss_cycle": peak_tss[1],
            "peak_mean_delta_U": peak_u[0],
            "peak_mean_delta_U_cycle": peak_u[1],
            "a_min_over_run": min(r["a_min"] for r in history),
            "a_max_over_run": max(r["a_max"] for r in history),
            "s_min_over_run": min(r["s_min"] for r in history),
            "s_max_over_run": max(r["s_max"] for r in history),
            "final_mean_delta_U": history[-1]["mean_delta_U"],
            "work_energy_relative_residual": energy_residual,
        },
        "history": history,
        "peak_var_s_state": [
            {"cell": i, "a": float(a_peak[i]), "s": float(s_peak[i]),
             "a_dot": float(va_peak[i]), "s_dot": float(vs_peak[i])}
            for i in range(CELLS)
        ],
    }


def write_outputs(results: list[dict], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = [result["summary"] for result in results]
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)
    with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    # Keep one raw cloud for the main coupled 5-cycle diagnostic.
    coupled = next(r for r in results if r["summary"]["case"] == "registry_driven_5cycle")
    cloud = coupled["peak_var_s_state"]
    with (output_dir / "peak_state_registry_driven.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(cloud[0]))
        writer.writeheader()
        writer.writerows(cloud)


if __name__ == "__main__":
    cases = [
        run_case("normal_only_5cycle", cycles=5, registry_amplitude=0.0),
        run_case("registry_driven_5cycle", cycles=5, registry_amplitude=REGISTRY_AMPLITUDE),
        run_case("registry_driven_10cycle", cycles=10, registry_amplitude=REGISTRY_AMPLITUDE),
    ]
    target = Path("results/data/reduced_as_spatial_chain")
    write_outputs(cases, target)
    for result in cases:
        print(json.dumps(result["summary"], indent=2))
