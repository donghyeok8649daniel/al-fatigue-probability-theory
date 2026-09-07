"""Reproduce the compact calibration plots from committed CSV results."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def generate(output: Path | None = None) -> tuple[Path, Path]:
    root = Path(__file__).resolve().parents[1]
    output = output or root / "results" / "aluminum_calibration"

    with (output / "fit_residuals.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    fitted = [row for row in rows if row["used_for_fit"] == "yes"]
    labels = [row["target"] for row in fitted[:5]]
    square = [float(row["normalized_residual"]) for row in fitted[:5]]
    linear = [float(row["normalized_residual"]) for row in fitted[5:10]]
    x = np.arange(len(labels))
    fig, axis = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    axis.bar(x - 0.19, square, 0.38, label="Square-root")
    axis.bar(x + 0.19, linear, 0.38, label="Linear extension")
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xticks(x, labels, rotation=22, ha="right")
    axis.set_ylabel("Normalized residual")
    axis.set_title("Pure-Al reduced-target fit residuals")
    axis.legend()
    residual_path = output / "fit_residuals.png"
    fig.savefig(residual_path, dpi=160)
    plt.close(fig)

    with (output / "barrier_hierarchy.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        barrier_rows = [
            row for row in csv.DictReader(stream)
            if row["model"] == "linear_best_feasible"
            and row["configurational_barrier_eV"]
        ]
    force = np.array([float(row["force"]) for row in barrier_rows])
    slip = np.array(
        [float(row["configurational_barrier_eV"]) for row in barrier_rows]
    )
    opening = np.array(
        [float(row["opening_barrier_at_s_saddle_eV"]) for row in barrier_rows]
    )
    fig, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    axis.plot(force, slip, "o-", label="Configurational barrier")
    axis.plot(force, opening, "s-", label="Opening barrier at s saddle")
    axis.set_xlabel("Generalized force f*")
    axis.set_ylabel("Barrier (eV / reduced cell)")
    axis.set_title("Best-feasible hybrid, hypothetical chi=0.5")
    axis.legend()
    barrier_path = output / "barrier_hierarchy.png"
    fig.savefig(barrier_path, dpi=160)
    plt.close(fig)
    return residual_path, barrier_path


if __name__ == "__main__":
    for path in generate():
        print(path)
