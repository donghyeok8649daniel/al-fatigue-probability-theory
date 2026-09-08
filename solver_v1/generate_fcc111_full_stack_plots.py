"""Generate compact static plots from the FCC(111) reference audit."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .energy_model_registry import (
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
)
from .fcc111_full_energy import FullFCC111StackEnergy


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "results" / "fcc111_full_stack"
    output.mkdir(parents=True, exist_ok=True)
    models = {
        "Full FCC LJ (unrefitted)": FullFCC111StackEnergy.from_reduced_model(
            build_energy_model(TWO_ROW_LJ_REFERENCE)
        ),
        "Full FCC hypothetical hybrid (unrefitted)":
            FullFCC111StackEnergy.from_reduced_model(
                build_energy_model(ANALYTIC_LJ_EAM_HYPOTHETICAL)
            ),
    }

    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    for label, model in models.items():
        normal = np.linspace(0.82 * model.a0, 1.8 * model.a0, 100)
        normal_energy = np.array([model.energy(float(a), 0.0) for a in normal])
        axes[0].plot(normal, normal_energy - np.min(normal_energy), label=label)
        registry = np.linspace(0.0, model.p.b, 101)
        registry_energy = np.array([model.energy(model.a0, float(s)) for s in registry])
        axes[1].plot(registry / model.p.b, registry_energy - np.min(registry_energy), label=label)
    axes[0].set(xlabel="a / b", ylabel="W - min(W)", title="Normal energy slice")
    axes[1].set(xlabel="s / b", ylabel="W - min(W)", title="Registry energy slice")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(output / "normal_and_registry_slices.png", dpi=180)
    plt.close(figure)

    with (output / "barrier_comparison.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    for label, model_id in (
        ("Full FCC LJ", TWO_ROW_LJ_REFERENCE),
        ("Full FCC hypothetical hybrid", ANALYTIC_LJ_EAM_HYPOTHETICAL),
    ):
        selected = [
            row for row in rows
            if row["model"] == model_id and row["geometry"].startswith("full_")
        ]
        x = np.array([float(row["force"]) for row in selected])
        opening = np.array([float(row["opening_barrier"]) for row in selected])
        configuration = np.array([
            np.nan if not row["configurational_barrier"]
            else float(row["configurational_barrier"])
            for row in selected
        ])
        ratio = np.array([
            np.nan if not row["barrier_ratio"] else float(row["barrier_ratio"])
            for row in selected
        ])
        axes[0].plot(x, opening, "o--", label=f"{label}: opening")
        axes[0].plot(x, configuration, "s-", label=f"{label}: registry")
        axes[1].plot(x, ratio, "o-", label=label)
    axes[0].set(xlabel="dimensionless generalized force f*", ylabel="barrier",
                title="Static barriers")
    axes[1].set(xlabel="dimensionless generalized force f*", ylabel="registry / opening",
                title="Static barrier ratio")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(output / "barrier_hierarchy.png", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
