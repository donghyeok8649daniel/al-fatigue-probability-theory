"""Direct 1D layer-LJ mechanics-versus-closure falsification test.

The deterministic chain is sampled at phase-locked physical times. For each
snapshot, only its measured mean stretch and mean configurational energy are
passed to the large-M distribution closure. The resulting density is then
compared with the actual finite spacing sample.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from theory.normal_lj_chain import NormalLJParameters, critical_stretch, simulate_normal_lj_chain
from theory.normal_lj_closure_validation import compare_snapshot_to_closure
from theory.normal_lj_distribution import closure_density, solve_distribution_closure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIG = ROOT / "results" / "figures"

CASES = (
    {
        "name": "subcritical_slow_t10T",
        "parameters": NormalLJParameters(force_amplitude=0.03, omega=0.01),
        "integration_cycles": 12,
        "sample_index": 10,
    },
    {
        "name": "subcritical_dynamic_t2T",
        "parameters": NormalLJParameters(force_amplitude=0.03, omega=0.02),
        "integration_cycles": 3,
        "sample_index": 2,
    },
)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows = []
    lam_c = critical_stretch()

    for case in CASES:
        result = simulate_normal_lj_chain(case["parameters"], cycles=case["integration_cycles"])
        values = result.cycle_snapshots[case["sample_index"]]
        sample_time = case["sample_index"] * result.period
        comparison = compare_snapshot_to_closure(
            values,
            closure_quadrature_order=640,
            cdf_quadrature_order=128,
        )
        row = {"case": case["name"], "time_star": sample_time, **asdict(comparison)}
        rows.append(row)

        solution = solve_distribution_closure(
            comparison.empirical_mean_stretch,
            comparison.empirical_mean_energy,
            quadrature_order=640,
        )
        lo = max(0.90, float(np.min(values)) - 0.02)
        hi = min(1.20, max(float(np.max(values)) + 0.06, 1.04))
        grid = np.linspace(lo, hi, 1800)
        plt.figure(figsize=(8, 5))
        plt.hist(values, bins=12, density=True, histtype="step", linewidth=2, label="deterministic 1D layer-LJ")
        plt.plot(
            grid,
            closure_density(
                grid,
                solution.moments.alpha,
                solution.moments.beta,
                quadrature_order=640,
            ),
            label="same mean + energy closure",
        )
        plt.axvline(lam_c, linestyle="--", label="lambda_c")
        plt.xlabel("Layer stretch lambda")
        plt.ylabel("Density")
        plt.title(case["name"])
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIG / f"normal_lj_closure_{case['name']}.svg")
        plt.close()

    with (DATA / "normal_lj_closure_falsification.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "classification": "direct numerical falsification test of a controlled closure",
        "lambda_c": lam_c,
        "comparison_rule": "The closure receives only the deterministic snapshot's measured mean stretch and mean shifted configurational energy.",
        "rows": rows,
    }
    (DATA / "normal_lj_closure_falsification.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
