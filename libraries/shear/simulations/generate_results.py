"""Generate all reference figures and machine-readable result files.

Run from repository root:
    python -m simulations.generate_results

The models are proof-of-principle reduced mechanics models, not calibrated Al fatigue-life predictors.
"""
from __future__ import annotations

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from theory.rubin_chain import (
    RubinParams,
    analytic_response,
    cycle_loop_areas,
    simulate_finite_chain,
)
from theory.hamiltonian_slip_bath import (
    SlipBathParameters,
    simulate_slip_bath,
)

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures"
DATA = ROOT / "results" / "data"


def _save_figure(stem: str) -> None:
    plt.savefig(FIG / f"{stem}.png", dpi=180)
    plt.savefig(FIG / f"{stem}.svg")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    params = RubinParams()
    analytic = analytic_response(0.5, 0.1, params)
    rubin = simulate_finite_chain(params=params)
    areas = cycle_loop_areas(rubin, first_cycle=10, last_cycle_exclusive=50)

    e_final = float(np.asarray(rubin["energy"])[-1])
    w_final = float(np.asarray(rubin["work"])[-1])
    rubin_summary = {
        **analytic,
        "numeric_loop_area_mean": float(np.mean(areas)),
        "numeric_loop_area_std": float(np.std(areas)),
        "loop_area_relative_error": float(
            abs(np.mean(areas) - analytic["loop_area"]) / analytic["loop_area"]
        ),
        "final_internal_energy": e_final,
        "final_external_work": w_final,
        "energy_balance_relative_error": float(
            abs(e_final - w_final) / max(abs(w_final), 1.0e-30)
        ),
    }

    period = float(rubin["period"])
    time = np.asarray(rubin["time"])
    q = np.asarray(rubin["q"])
    force = np.asarray(rubin["force"])

    cycle = 30
    mask = (time >= cycle * period) & (time < (cycle + 1) * period)
    plt.figure(figsize=(7, 5))
    plt.plot(q[mask], force[mask])
    plt.xlabel("Resolved coordinate Q")
    plt.ylabel("External force F")
    plt.title("Rubin-chain reduced hysteresis (cycle 30)")
    plt.grid(True)
    plt.tight_layout()
    _save_figure("rubin_hysteresis_cycle30")
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(np.arange(10, 50), areas, marker="o", markersize=3)
    plt.axhline(
        analytic["loop_area"],
        linestyle="--",
        label="Analytic semi-infinite result",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Loop area ∮F dQ")
    plt.title("Rubin-chain loop-area convergence")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    _save_figure("rubin_loop_area_convergence")
    plt.close()

    (DATA / "rubin_reference.json").write_text(
        json.dumps(rubin_summary, indent=2), encoding="utf-8"
    )
    with (DATA / "rubin_loop_areas.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cycle", "loop_area"])
        for cycle, area in zip(range(10, 50), areas):
            writer.writerow([cycle, float(area)])

    slip_results = {}
    for amplitude in (0.34, 0.40, 0.50):
        p = SlipBathParameters(force_amplitude=amplitude)
        slip_results[amplitude] = simulate_slip_bath(p)

    plt.figure(figsize=(8, 5))
    for amplitude, result in slip_results.items():
        n = np.arange(1, len(result.cycle_slip) + 1)
        plt.plot(n, result.cycle_slip, marker="o", label=f"Fa={amplitude:.2f}")
    plt.xlabel("Cycle number N")
    plt.ylabel("Resolved slip state s(NT)")
    plt.title("Cycle-to-cycle structural evolution")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    _save_figure("slip_cycle_accumulation")
    plt.close()

    running = slip_results[0.50]
    cycle = 8
    mask = (
        (running.time >= cycle * running.period)
        & (running.time < (cycle + 1) * running.period)
    )
    plt.figure(figsize=(7, 5))
    plt.plot(running.slip[mask], running.force[mask])
    plt.xlabel("Slip coordinate s")
    plt.ylabel("External generalized force F")
    plt.title("Nonlinear slip-bath hysteresis (Fa=0.50, cycle 8)")
    plt.grid(True)
    plt.tight_layout()
    _save_figure("slip_hysteresis_running")
    plt.close()

    plt.figure(figsize=(8, 5))
    for amplitude, result in slip_results.items():
        n = np.arange(1, len(result.cycle_spacing_variance) + 1)
        plt.plot(
            n,
            result.cycle_spacing_variance,
            marker="o",
            label=f"Fa={amplitude:.2f}",
        )
    plt.xlabel("Cycle number N")
    plt.ylabel("Variance of spacing-like relative displacement")
    plt.title("Redistribution into unresolved lattice modes")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    _save_figure("spacing_variance_by_cycle")
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(
        running.time,
        running.internal_energy - running.internal_energy[0],
        label="Δ internal energy",
    )
    plt.plot(
        running.time,
        running.external_work,
        label="Integrated external work",
    )
    plt.xlabel("Time")
    plt.ylabel("Energy / work (nondimensional)")
    plt.title("Energy-balance check, Fa=0.50")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    _save_figure("slip_energy_balance")
    plt.close()

    slip_summary = {}
    with (DATA / "slip_cycle_history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "force_amplitude",
                "cycle",
                "cycle_slip",
                "spacing_like_variance",
            ]
        )
        for amplitude, result in slip_results.items():
            for cycle, (s, variance) in enumerate(
                zip(result.cycle_slip, result.cycle_spacing_variance),
                start=1,
            ):
                writer.writerow([amplitude, cycle, float(s), float(variance)])

            slip_summary[str(amplitude)] = {
                "last_six_cycle_states": [
                    float(x) for x in result.cycle_slip[-6:]
                ],
                "final_energy_balance_relative_error":
                    result.final_energy_balance_relative_error,
                "cycle_2_variance": float(result.cycle_spacing_variance[1]),
                "final_cycle_variance": float(
                    result.cycle_spacing_variance[-1]
                ),
            }

    (DATA / "slip_reference.json").write_text(
        json.dumps(slip_summary, indent=2), encoding="utf-8"
    )

    manifest = {
        "model_status": "proof-of-principle, not calibrated Al fatigue life",
        "rubin": rubin_summary,
        "slip": slip_summary,
        "figures": sorted(
            p.name for p in FIG.glob("*.*") if p.suffix in {".png", ".svg"}
        ),
    }
    (DATA / "result_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
