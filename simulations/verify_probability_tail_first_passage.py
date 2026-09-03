"""Verify snapshot-tail flux balance and cumulative first passage in the 1D normal LJ chain.

This numerical audit deliberately distinguishes three objects:

1. the instantaneous/snapshot tail mass above a chosen spacing threshold;
2. cumulative first-passage mass, which remembers whether a spacing has ever crossed;
3. irreversible dissipation, which is not introduced here.

The script uses the same conservative generalized-LJ force functions as the
active theory branch after branch integration.
"""

from __future__ import annotations

from pathlib import Path
import csv
import json
import math

import numpy as np

from theory.normal_lj_chain import normalized_lj_energy, normalized_lj_force, critical_stretch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "data" / "probability_tail_first_passage"

ATOMS = 512
DT = 0.01
CYCLES = 6
OMEGA = 0.02
RAMP_CYCLES = 2
YOUNGS_MODULUS = 69.0e9
MEAN_STRESS = 100.0e6
STRESS_AMPLITUDE = 100.0e6
MEAN_FORCE = MEAN_STRESS / YOUNGS_MODULUS
FORCE_AMPLITUDE = STRESS_AMPLITUDE / YOUNGS_MODULUS
THRESHOLDS = (1.002, 1.003, 1.004, 1.005, 1.006)
AUDIT_THRESHOLD = 1.004


def run() -> tuple[list[dict[str, float | int]], dict]:
    period = 2.0 * math.pi / OMEGA
    nsteps = int(round(CYCLES * period / DT))

    x = np.arange(ATOMS, dtype=float)
    velocity = np.zeros(ATOMS, dtype=float)

    def envelope(t: float) -> float:
        ramp_time = RAMP_CYCLES * period
        if t >= ramp_time:
            return 1.0
        return 0.5 * (1.0 - math.cos(math.pi * t / ramp_time))

    def external_force(t: float) -> float:
        return envelope(t) * (
            MEAN_FORCE + FORCE_AMPLITUDE * math.sin(OMEGA * t)
        )

    def force_vector(state: np.ndarray, t: float) -> np.ndarray:
        spacing = np.diff(state)
        dphi = normalized_lj_force(spacing)
        force = np.zeros_like(state)
        force[1:-1] = dphi[1:] - dphi[:-1]
        force[-1] = -dphi[-1] + external_force(t)
        return force

    def mechanical_energy() -> float:
        spacing = np.diff(x)
        kinetic = 0.5 * float(np.dot(velocity[1:], velocity[1:]))
        potential = float(np.sum(normalized_lj_energy(spacing)))
        return kinetic + potential

    force = force_vector(x, 0.0)
    initial_energy = mechanical_energy()
    work = 0.0
    previous_power = 0.0

    spacing = np.diff(x)
    maximum_ever = spacing.copy()
    previous_audit_state = spacing >= AUDIT_THRESHOLD
    upward_crossings = 0
    downward_crossings = 0

    next_cycle = 1
    rows: list[dict[str, float | int]] = []

    for step in range(nsteps):
        t = step * DT

        velocity[1:] += 0.5 * DT * force[1:]
        x[1:] += DT * velocity[1:]

        new_force = force_vector(x, t + DT)
        velocity[1:] += 0.5 * DT * new_force[1:]
        force = new_force

        spacing = np.diff(x)
        maximum_ever = np.maximum(maximum_ever, spacing)

        current_audit_state = spacing >= AUDIT_THRESHOLD
        upward_crossings += int(np.count_nonzero(~previous_audit_state & current_audit_state))
        downward_crossings += int(np.count_nonzero(previous_audit_state & ~current_audit_state))
        previous_audit_state = current_audit_state

        q_now = external_force(t + DT)
        power = q_now * velocity[-1]
        work += 0.5 * DT * (previous_power + power)
        previous_power = power

        if t + DT >= next_cycle * period:
            energy_change = mechanical_energy() - initial_energy
            row: dict[str, float | int] = {
                "cycle": next_cycle,
                "mean": float(np.mean(spacing)),
                "variance": float(np.var(spacing)),
                "q99": float(np.quantile(spacing, 0.99)),
                "max": float(np.max(spacing)),
                "energy_change": float(energy_change),
                "external_work": float(work),
                "work_minus_energy": float(work - energy_change),
                "audit_up": upward_crossings,
                "audit_down": downward_crossings,
                "audit_net": upward_crossings - downward_crossings,
                "audit_snapshot_count": int(np.count_nonzero(current_audit_state)),
            }
            for threshold in THRESHOLDS:
                row[f"snapshot_ge_{threshold:.3f}"] = float(np.mean(spacing >= threshold))
                row[f"firstpass_ge_{threshold:.3f}"] = float(
                    np.mean(maximum_ever >= threshold)
                )
            rows.append(row)
            upward_crossings = 0
            downward_crossings = 0
            next_cycle += 1

    snapshot_004 = [float(row["snapshot_ge_1.004"]) for row in rows]
    work_errors = [abs(float(row["work_minus_energy"])) for row in rows]

    crossing_balance = True
    previous_count = 0
    for row in rows:
        count = int(row["audit_snapshot_count"])
        crossing_balance &= count - previous_count == int(row["audit_net"])
        previous_count = count

    all_firstpass_monotone = True
    for threshold in THRESHOLDS:
        key = f"firstpass_ge_{threshold:.3f}"
        values = [float(row[key]) for row in rows]
        all_firstpass_monotone &= all(
            b + 1.0e-15 >= a for a, b in zip(values[:-1], values[1:])
        )

    summary = {
        "classification": "deterministic conservative 1D normal-chain tail/first-passage audit",
        "parameters": {
            "atoms": ATOMS,
            "dt": DT,
            "cycles": CYCLES,
            "mean_stress_pa_under_sigma_over_E_mapping": MEAN_STRESS,
            "stress_amplitude_pa_under_sigma_over_E_mapping": STRESS_AMPLITUDE,
            "omega_star": OMEGA,
            "ramp_cycles": RAMP_CYCLES,
            "thresholds": list(THRESHOLDS),
            "crossing_balance_threshold": AUDIT_THRESHOLD,
        },
        "checks": {
            "max_abs_work_minus_energy": max(work_errors),
            "crossing_balance_exact_at_each_cycle": bool(crossing_balance),
            "snapshot_tail_monotone": bool(
                all(b + 1.0e-15 >= a for a, b in zip(snapshot_004[:-1], snapshot_004[1:]))
            ),
            "first_passage_tail_monotone_for_each_threshold": bool(all_firstpass_monotone),
            "critical_stretch": critical_stretch(),
            "critical_first_passage_fraction_through_cycle_6": float(
                np.mean(maximum_ever >= critical_stretch())
            ),
        },
        "interpretation": [
            "Snapshot tail mass is not a cumulative damage variable: it can grow and later decrease because cells can cross back below the threshold.",
            "Cumulative first-passage mass is nondecreasing by construction and remains distinct from irreversible thermodynamic dissipation.",
            "At the audit threshold the cycle-end change in the number of cells above threshold equals upward crossings minus downward crossings.",
            "External work matches the change in conservative mechanical energy to numerical integration accuracy.",
            "The omega*=0.02 protocol is a proof-of-principle atomic-chain calculation and is not a direct mapping to laboratory fatigue frequency.",
        ],
    }

    return rows, summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, summary = run()

    csv_path = OUT / "tail_audit.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
