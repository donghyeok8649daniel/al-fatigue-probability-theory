# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: _element_histories, write_probability_element_history, _last_cycle_indices
#   plot_probability_hysteresis_summary, plot_mesh_probability_fields, run_demo, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Run the tensile FEM -> local spacing probability demonstration."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import numpy as np

from simulations.fem_tension_app import TensionRunConfig, run_fem_solver
from simulations.fem_tension_ui import load_fem_history, save_preview_images
from theory.normal_lj_probability_dynamics import (
    ProbabilityHistory,
    SpacingDynamicsParameters,
    completed_cycle_hysteresis_areas,
    solve_spacing_probability_history,
)


def _element_histories(
    elements: np.ndarray,
    youngs_modulus_pa: float,
    parameters: SpacingDynamicsParameters,
) -> dict[int, ProbabilityHistory]:
    histories: dict[int, ProbabilityHistory] = {}
    cache: dict[bytes, ProbabilityHistory] = {}
    for element in np.unique(elements["element"]).astype(int):
        rows = elements[elements["element"] == element]
        rows = rows[np.argsort(rows["step"])]
        stress = np.asarray(rows["stress_pa"], dtype=float)
        key = np.round(stress, decimals=6).tobytes()
        if key not in cache:
            cache[key] = solve_spacing_probability_history(
                np.asarray(rows["time_s"], dtype=float),
                stress,
                youngs_modulus_pa,
                parameters,
            )
        histories[element] = cache[key]
    return histories


def write_probability_element_history(
    path: Path,
    elements: np.ndarray,
    histories: dict[int, ProbabilityHistory],
    equilibrium_spacing_m: float,
) -> None:
    """Write the FEM-element/probability interface as one tidy CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    x_by_element = {
        int(element): float(np.mean(elements["x_mid_m"][elements["element"] == element]))
        for element in np.unique(elements["element"]).astype(int)
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_s",
                "step",
                "element",
                "x_mid_m",
                "stress_pa",
                "normalization",
                "mean_stretch",
                "mean_spacing_m",
                "variance_stretch",
                "mean_lj_energy_density_j_m3",
                "nonequilibrium_free_energy_density_j_m3",
                "cumulative_hysteresis_energy_density_j_m3",
                "critical_tail_probability",
            ]
        )
        for element, history in sorted(histories.items()):
            for step in range(history.time_s.size):
                writer.writerow(
                    [
                        f"{history.time_s[step]:.17g}",
                        step,
                        element,
                        f"{x_by_element[element]:.17g}",
                        f"{history.stress_pa[step]:.17g}",
                        f"{history.normalization[step]:.17g}",
                        f"{history.mean_stretch[step]:.17g}",
                        f"{equilibrium_spacing_m * history.mean_stretch[step]:.17g}",
                        f"{history.variance_stretch[step]:.17g}",
                        f"{history.mean_energy_density_j_m3[step]:.17g}",
                        f"{history.nonequilibrium_free_energy_density_j_m3[step]:.17g}",
                        f"{history.cumulative_hysteresis_energy_density_j_m3[step]:.17g}",
                        f"{history.critical_tail_probability[step]:.17g}",
                    ]
                )


def _last_cycle_indices(history: ProbabilityHistory, frequency_hz: float) -> np.ndarray:
    period = 1.0 / frequency_hz
    start = history.time_s[-1] - period - 1.0e-12
    return np.flatnonzero(history.time_s >= start)


def plot_probability_hysteresis_summary(
    history: ProbabilityHistory,
    frequency_hz: float,
    equilibrium_spacing_m: float,
    output_path: Path,
) -> None:
    """Plot the four active equations as one numerical diagnostic."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cycle = _last_cycle_indices(history, frequency_hz)
    n = history.time_s.size
    steps_per_cycle = cycle.size - 1
    cycle_start = cycle[0]
    phase_indices = [
        cycle_start,
        min(cycle_start + steps_per_cycle // 4, n - 1),
        min(cycle_start + steps_per_cycle // 2, n - 1),
        min(cycle_start + 3 * steps_per_cycle // 4, n - 1),
    ]
    labels = ("mean/loading", "peak tension", "mean/unloading", "minimum tension")

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0))
    ax = axes[0, 0]
    ax.plot(history.time_s, history.stress_pa / 1.0e6, color="tab:red", label="stress")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("stress [MPa]", color="tab:red")
    ax.tick_params(axis="y", labelcolor="tab:red")
    ax.grid(True, alpha=0.3)
    twin = ax.twinx()
    twin.plot(
        history.time_s,
        equilibrium_spacing_m * history.mean_stretch * 1.0e10,
        color="tab:blue",
        label="mean spacing",
    )
    twin.set_ylabel("mean spacing [angstrom]", color="tab:blue")
    twin.tick_params(axis="y", labelcolor="tab:blue")
    ax.set_title("Stress input and mean-spacing response")

    ax = axes[0, 1]
    for index, label in zip(phase_indices, labels):
        ax.plot(history.stretch, history.density[index], label=label)
    ax.set_xlabel("normalized spacing lambda")
    ax.set_ylabel("p(lambda,t)")
    ax.set_title("Loading/unloading distributions in the last cycle")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    mean_strain = history.mean_stretch[cycle] - 1.0
    ax.plot(mean_strain, history.stress_pa[cycle] / 1.0e6, color="tab:purple")
    ax.fill(mean_strain, history.stress_pa[cycle] / 1.0e6, alpha=0.15, color="tab:purple")
    ax.set_xlabel("mean spacing strain = mean(lambda)-1")
    ax.set_ylabel("stress [MPa]")
    area = completed_cycle_hysteresis_areas(history, frequency_hz)[-1]
    ax.set_title(f"Energy hysteresis — last-cycle area = {area:.4g} J/m^3")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(
        history.time_s,
        history.mean_energy_density_j_m3,
        label="mean LJ energy density",
        color="tab:green",
    )
    ax.plot(
        history.time_s,
        history.cumulative_hysteresis_energy_density_j_m3,
        label="cumulative path work",
        color="tab:orange",
    )
    ax.set_xlabel("time [s]")
    ax.set_ylabel("energy density [J/m^3]")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(
        history.time_s,
        history.critical_tail_probability,
        linestyle="--",
        color="tab:red",
        label="tail above lambda_c",
    )
    ax2.set_ylabel("critical-tail probability", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    lines, line_labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, line_labels + labels2, fontsize=8, loc="upper left")
    ax.set_title("Energy state and instability-tail diagnostic")

    fig.suptitle(
        "1D layer-LJ probability dynamics — candidate kinetic extension, not calibrated Al life prediction",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_mesh_probability_fields(
    nodes: np.ndarray,
    elements: np.ndarray,
    histories: dict[int, ProbabilityHistory],
    output_path: Path,
) -> None:
    """Render the actual 1D FEM partition and probability-derived element fields."""
    first_history = histories[min(histories)]
    peak_step = int(np.argmax(first_history.stress_pa))
    node_rows = nodes[nodes["step"] == peak_step]
    element_rows = elements[elements["step"] == peak_step]
    order = np.argsort(element_rows["element"])
    element_rows = element_rows[order]
    element_ids = np.asarray(element_rows["element"], dtype=int)
    mean_stretch = np.asarray([histories[int(e)].mean_stretch[peak_step] for e in element_ids])
    tail = np.asarray([histories[int(e)].critical_tail_probability[peak_step] for e in element_ids])
    stress_mpa = np.asarray(element_rows["stress_pa"], dtype=float) / 1.0e6
    x_nodes = np.asarray(node_rows["x_m"], dtype=float)
    x_mid = np.asarray(element_rows["x_mid_m"], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(11.0, 6.8), sharex=True)
    ax = axes[0]
    tail_span = float(np.ptp(tail))
    tail_scale = max(abs(float(np.mean(tail))), 1.0e-12)
    if tail_span <= 1.0e-9 * tail_scale:
        tail = np.full_like(tail, float(np.mean(tail)))
    color_half_span = max(0.05 * abs(float(np.mean(tail))), 1.0e-12)
    norm = Normalize(
        vmin=float(np.min(tail)) - color_half_span,
        vmax=float(np.max(tail)) + color_half_span,
    )
    cmap = plt.get_cmap("magma")
    height = 1.0
    for i, element in enumerate(element_ids):
        width = x_nodes[i + 1] - x_nodes[i]
        ax.add_patch(
            Rectangle(
                (x_nodes[i], -height / 2.0),
                width,
                height,
                facecolor=cmap(norm(tail[i])),
                edgecolor="black",
                linewidth=0.65,
            )
        )
        if len(element_ids) <= 30:
            ax.text(x_mid[i], 0.0, str(element), ha="center", va="center", fontsize=6)
    ax.scatter(x_nodes, np.full_like(x_nodes, -height / 2.0), s=9, color="black", zorder=3)
    ax.set_ylim(-0.8, 0.8)
    ax.set_yticks([])
    ax.set_title(f"Actual 1D FEM mesh at peak tension — {len(element_ids)} elements")

    axes[1].plot(x_mid, stress_mpa, marker="o", markersize=3)
    axes[1].set_ylabel("axial stress [MPa]")
    axes[1].grid(True, alpha=0.3)
    if float(np.ptp(stress_mpa)) <= 1.0e-9 * max(abs(float(np.mean(stress_mpa))), 1.0):
        stress_margin = max(0.01 * abs(float(np.mean(stress_mpa))), 1.0)
        axes[1].set_ylim(
            float(np.mean(stress_mpa)) - stress_margin,
            float(np.mean(stress_mpa)) + stress_margin,
        )
        axes[1].text(
            0.5,
            0.12,
            "uniform-area bar: axial stress is spatially uniform",
            transform=axes[1].transAxes,
            ha="center",
            fontsize=8,
        )
    mean_line = axes[2].plot(
        x_mid,
        mean_stretch,
        marker="o",
        markersize=3,
        color="tab:blue",
        label="mean stretch",
    )
    axes[2].set_ylabel("mean stretch", color="tab:blue")
    axes[2].tick_params(axis="y", labelcolor="tab:blue")
    tail_ax = axes[2].twinx()
    tail_line = tail_ax.plot(
        x_mid,
        tail,
        marker="s",
        markersize=3,
        color="tab:orange",
        label="critical tail",
    )
    tail_ax.set_ylabel("critical-tail probability", color="tab:orange")
    tail_ax.tick_params(axis="y", labelcolor="tab:orange")
    axes[2].set_xlabel("tensile coordinate x [m]")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(mean_line + tail_line, ["mean stretch", "critical tail"], fontsize=8)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_demo(
    data_dir: Path,
    figure_dir: Path,
    elements_count: int = 32,
) -> dict[str, float | int | str | list[float]]:
    """Run a reproducible tension-tension demonstration and save all outputs."""
    fem_dir = data_dir / "fem"
    config = TensionRunConfig(
        length_mm=50.0,
        width_mm=10.0,
        thickness_mm=1.0,
        young_gpa=69.0,
        elements=elements_count,
        stress_mean_mpa=120.0,
        stress_amplitude_mpa=80.0,
        frequency_hz=5.0,
        cycles=3,
        steps_per_cycle=80,
        deformation_scale=25.0,
    )
    probability_parameters = SpacingDynamicsParameters(
        inverse_temperature=2000.0,
        relaxation_time_s=0.03,
        grid_cells=220,
        substeps_per_interval=2,
    )
    equilibrium_spacing_m = 2.86e-10

    run_fem_solver(config, fem_dir)
    nodes, elements = load_fem_history(fem_dir)
    histories = _element_histories(elements, config.young_pa, probability_parameters)
    write_probability_element_history(
        data_dir / "probability_elements.csv",
        elements,
        histories,
        equilibrium_spacing_m,
    )

    representative = histories[int(np.unique(elements["element"])[len(histories) // 2])]
    areas = completed_cycle_hysteresis_areas(representative, config.frequency_hz)
    plot_probability_hysteresis_summary(
        representative,
        config.frequency_hz,
        equilibrium_spacing_m,
        figure_dir / "probability_hysteresis_summary.png",
    )
    plot_mesh_probability_fields(
        nodes,
        elements,
        histories,
        figure_dir / "fem_mesh_probability_fields.png",
    )
    save_preview_images(
        fem_dir,
        figure_dir,
        half_width_m=config.width_m / 2.0,
        half_thickness_m=config.thickness_m / 2.0,
        deformation_scale=config.deformation_scale,
        field="stress",
    )

    summary: dict[str, float | int | str | list[float]] = {
        "status": "candidate kinetic demonstration; chi and relaxation time are not calibrated aluminum values",
        "elements": elements_count,
        "probability_grid_cells": probability_parameters.grid_cells,
        "inverse_temperature_chi": probability_parameters.inverse_temperature,
        "relaxation_time_s": probability_parameters.relaxation_time_s,
        "max_normalization_error": float(np.max(np.abs(representative.normalization - 1.0))),
        "last_cycle_hysteresis_energy_density_j_m3": float(areas[-1]),
        "cycle_hysteresis_energy_density_j_m3": [float(value) for value in areas],
        "peak_mean_spacing_m": float(equilibrium_spacing_m * np.max(representative.mean_stretch)),
        "peak_critical_tail_probability": float(np.max(representative.critical_tail_probability)),
    }
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run FEM -> layer-LJ probability coupling demo")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("results/data/fem_probability_demo"),
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("results/figures/fem_probability_demo"),
    )
    parser.add_argument("--elements", type=int, default=32)
    args = parser.parse_args()
    summary = run_demo(args.data_dir, args.figure_dir, args.elements)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
