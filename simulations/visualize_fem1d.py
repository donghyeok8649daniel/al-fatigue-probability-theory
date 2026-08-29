# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: load_numeric_csv, select_snapshot_step, _save_line, plot_fem1d_results, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Visualizer for the standalone 1D bar FEM scaffold."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_numeric_csv(path: Path) -> np.ndarray:
    """Load one headered numeric CSV as a one-dimensional structured array."""
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    return np.atleast_1d(data)


def select_snapshot_step(elements: np.ndarray, mode: str = "peak-tension") -> int:
    """Choose a representative time step from element history."""
    steps = np.unique(elements["step"]).astype(int)
    mean_stress = []
    for step in steps:
        values = elements["stress_pa"][elements["step"] == step]
        mean_stress.append(float(np.mean(values)))
    mean_stress = np.asarray(mean_stress)

    if mode == "peak-tension":
        return int(steps[int(np.argmax(mean_stress))])
    if mode == "peak-absolute":
        return int(steps[int(np.argmax(np.abs(mean_stress)))])
    if mode == "final":
        return int(steps[-1])
    raise ValueError(f"unknown snapshot mode: {mode}")


def _save_line(x, y, xlabel: str, ylabel: str, title: str, path: Path) -> None:
    plt.figure(figsize=(7.2, 4.6))
    plt.plot(x, y, marker="o")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_fem1d_results(
    input_dir: Path,
    output_dir: Path,
    snapshot_mode: str = "peak-tension",
) -> dict[str, float | int]:
    """Create basic 1D FEM plots and return a compact numerical summary."""
    nodes = load_numeric_csv(input_dir / "nodes.csv")
    elements = load_numeric_csv(input_dir / "elements.csv")
    output_dir.mkdir(parents=True, exist_ok=True)

    step = select_snapshot_step(elements, snapshot_mode)
    node_snapshot = nodes[nodes["step"] == step]
    element_snapshot = elements[elements["step"] == step]

    _save_line(
        node_snapshot["x_m"],
        node_snapshot["displacement_m"],
        "x [m]",
        "displacement u [m]",
        f"1D FEM displacement — step {step}",
        output_dir / "displacement_snapshot.png",
    )
    _save_line(
        element_snapshot["x_mid_m"],
        element_snapshot["strain"],
        "x [m]",
        "axial strain",
        f"1D FEM strain — step {step}",
        output_dir / "strain_snapshot.png",
    )
    _save_line(
        element_snapshot["x_mid_m"],
        element_snapshot["stress_pa"] / 1.0e6,
        "x [m]",
        "axial stress [MPa]",
        f"1D FEM stress — step {step}",
        output_dir / "stress_snapshot.png",
    )

    element_ids = np.unique(elements["element"]).astype(int)
    selected_element = int(element_ids[len(element_ids) // 2])
    history = elements[elements["element"] == selected_element]
    _save_line(
        history["time_s"],
        history["stress_pa"] / 1.0e6,
        "time [s]",
        "axial stress [MPa]",
        f"Stress history — element {selected_element}",
        output_dir / "stress_history.png",
    )

    return {
        "snapshot_step": step,
        "snapshot_time_s": float(element_snapshot["time_s"][0]),
        "max_snapshot_stress_pa": float(np.max(element_snapshot["stress_pa"])),
        "min_snapshot_stress_pa": float(np.min(element_snapshot["stress_pa"])),
        "tip_displacement_m": float(node_snapshot["displacement_m"][-1]),
        "selected_history_element": selected_element,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize the standalone 1D FEM CSV output.")
    parser.add_argument("--input-dir", type=Path, default=Path("fem1d_output"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/fem1d_demo"))
    parser.add_argument(
        "--snapshot",
        choices=("peak-tension", "peak-absolute", "final"),
        default="peak-tension",
    )
    args = parser.parse_args()

    summary = plot_fem1d_results(args.input_dir, args.output_dir, args.snapshot)
    for key, value in summary.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
