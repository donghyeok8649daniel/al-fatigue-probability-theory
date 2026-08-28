"""Generate the continuous-time 1D normal-LJ energy-feasibility reference result.

The lower compression bounds used for the plotted parameter sweep are
mathematical illustrations only. They are not calibrated Al material values.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from theory.normal_lj_chain import critical_stretch
from theory.normal_lj_energy_feasibility import (
    no_compression_bound_counterexample_energy,
    safe_energy_interval,
    shifted_lj_energy,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIG = ROOT / "results" / "figures"


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    lam_c = critical_stretch()
    illustrative_bounds = (0.90, 0.95, 0.98, 0.99)

    summary = {
        "classification": "exact convex-measure bound under stated support constraint",
        "active_dimension": "1D",
        "state_variable": "P(a,t)",
        "lambda_c": lam_c,
        "shifted_energy_at_lambda_c": float(shifted_lj_energy(lam_c)),
        "illustrative_lower_bounds_are_material_inputs": False,
        "illustrative_safe_ceiling_at_mean_1": {
            str(lower): safe_energy_interval(1.0, lower).maximum_energy
            for lower in illustrative_bounds
        },
        "no_lower_bound_counterexample": {
            "epsilon_0.8": no_compression_bound_counterexample_energy(0.8, 1.0),
            "epsilon_0.5": no_compression_bound_counterexample_energy(0.5, 1.0),
            "epsilon_0.3": no_compression_bound_counterexample_energy(0.3, 1.0),
        },
    }

    (DATA / "normal_lj_energy_feasibility_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    means = np.linspace(0.995, 1.07, 300)
    plt.figure(figsize=(8, 5))
    for lower in illustrative_bounds:
        values = [safe_energy_interval(mu, lower).maximum_energy for mu in means]
        plt.plot(means, values, label=fr"$\lambda_L={lower:.2f}$")

    plt.plot(
        means,
        shifted_lj_energy(means),
        linestyle="--",
        label=r"minimum $\psi(\mu)$",
    )
    plt.axvline(lam_c, linestyle=":", label=fr"$\lambda_c={lam_c:.4f}$")
    plt.xlabel(r"Mean stretch $\mu(t)$")
    plt.ylabel("Shifted LJ energy per bond (dimensionless)")
    plt.title("Exact crack-free energy interval for 1D LJ")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "normal_lj_safe_energy_ceiling.svg")
    plt.close()

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
