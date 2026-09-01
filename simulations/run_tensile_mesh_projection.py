# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D FEM/확률 이력의 normal-only scalar를 실제 2D/3D 또는 CAD mesh cell에 매핑해 저장한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: load_element_snapshot, build_or_load_mesh, write_projection_csv, run_projection
#   _parse_axis, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Map a one-dimensional axial result to an actual 2D/3D geometry mesh.

The generated field is explicitly a normal-only projection.  It is suitable
for mesh inspection and for carrying ``sigma_nn`` or a probability-derived
scalar to each cell, but it is not a multidimensional elasticity solve.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from simulations.fem_geometry_mesh import (
    GeometryMesh,
    load_geometry_mesh,
    map_axial_element_field,
    save_mesh_npz,
    structured_box_mesh,
    structured_rectangle_mesh,
)
from simulations.fem_mesh_ui import save_geometry_mesh_preview
from simulations.visualize_fem1d import load_numeric_csv


def load_element_snapshot(
    history_csv: Path,
    field: str,
    step: str | int,
) -> tuple[np.ndarray, np.ndarray, int, float]:
    """Load one element-centered scalar snapshot from FEM/probability CSV."""
    history = load_numeric_csv(Path(history_csv))
    required = {"step", "time_s", "element", "x_mid_m", field}
    missing = sorted(required.difference(history.dtype.names or ()))
    if missing:
        raise ValueError(f"history CSV is missing columns: {', '.join(missing)}")
    steps = np.unique(history["step"]).astype(int)
    if isinstance(step, str):
        if step != "peak-tension":
            raise ValueError("step string must be 'peak-tension'")
        if "stress_pa" not in (history.dtype.names or ()):
            raise ValueError("peak-tension selection requires stress_pa")
        means = np.array(
            [np.mean(history["stress_pa"][history["step"] == candidate]) for candidate in steps]
        )
        selected = int(steps[int(np.argmax(means))])
    else:
        selected = int(step)
        if selected not in steps:
            raise ValueError(f"step {selected} is absent from history")
    rows = history[history["step"] == selected]
    rows = rows[np.argsort(rows["element"])]
    return (
        np.asarray(rows["x_mid_m"], dtype=float),
        np.asarray(rows[field], dtype=float),
        selected,
        float(rows["time_s"][0]),
    )


def build_or_load_mesh(
    *,
    dimension: int,
    length_m: float,
    width_m: float,
    thickness_m: float,
    nx: int,
    ny: int,
    nz: int,
    geometry_path: Path | None,
    coordinate_scale_to_m: float,
    characteristic_length_m: float | None,
) -> GeometryMesh:
    """Create a structured mesh or load/mesh the supplied geometry file."""
    if geometry_path is not None:
        return load_geometry_mesh(
            geometry_path,
            coordinate_scale_to_m=coordinate_scale_to_m,
            target_dimension=dimension,
            characteristic_length_m=characteristic_length_m,
        )
    if dimension == 2:
        return structured_rectangle_mesh(length_m, width_m, nx, ny)
    if dimension == 3:
        return structured_box_mesh(length_m, width_m, thickness_m, nx, ny, nz)
    raise ValueError("dimension must be 2 or 3")


def write_projection_csv(
    path: Path,
    centers_m: np.ndarray,
    axial_coordinate: np.ndarray,
    values: np.ndarray,
    field: str,
) -> None:
    """Write cell centers and projected scalar without inventing tensor fields."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    centers = np.asarray(centers_m, dtype=float)
    padded = (
        np.column_stack([centers, np.zeros(centers.shape[0])])
        if centers.shape[1] == 2
        else centers
    )
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["cell", "x_m", "y_m", "z_m", "normalized_axial_coordinate", field])
        for cell, (center, xi, value) in enumerate(zip(padded, axial_coordinate, values)):
            writer.writerow(
                [
                    cell,
                    f"{center[0]:.17g}",
                    f"{center[1]:.17g}",
                    f"{center[2]:.17g}",
                    f"{xi:.17g}",
                    f"{value:.17g}",
                ]
            )


def run_projection(
    *,
    history_csv: Path,
    field: str,
    step: str | int,
    output_dir: Path,
    dimension: int,
    length_m: float,
    width_m: float,
    thickness_m: float,
    nx: int,
    ny: int,
    nz: int,
    geometry_path: Path | None = None,
    coordinate_scale_to_m: float = 1.0,
    characteristic_length_m: float | None = None,
    tensile_axis: np.ndarray | None = None,
    preview_path: Path | None = None,
) -> dict[str, float | int | str]:
    """Execute and persist one explicitly labeled axial mesh projection."""
    x_mid, source_values, selected_step, time_s = load_element_snapshot(history_csv, field, step)
    mesh = build_or_load_mesh(
        dimension=dimension,
        length_m=length_m,
        width_m=width_m,
        thickness_m=thickness_m,
        nx=nx,
        ny=ny,
        nz=nz,
        geometry_path=geometry_path,
        coordinate_scale_to_m=coordinate_scale_to_m,
        characteristic_length_m=characteristic_length_m,
    )
    if tensile_axis is None:
        tensile_axis = np.eye(mesh.embedding_dimension, dtype=float)[0]
    projection = map_axial_element_field(
        mesh,
        x_mid,
        source_values,
        tensile_axis=tensile_axis,
    )
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    save_mesh_npz(output / "geometry_mesh.npz", mesh)
    write_projection_csv(
        output / "projected_cell_field.csv",
        projection.cell_centers_m,
        projection.normalized_axial_coordinate,
        projection.values,
        field,
    )
    resolved_preview_path = output / "mesh_field_preview.png" if preview_path is None else Path(preview_path)
    preview = save_geometry_mesh_preview(
        resolved_preview_path,
        mesh,
        projection.values,
        field_label=field,
        title=(
            f"{dimension}D mesh — {field.replace('_', ' ')}\n"
            "normal-only axial projection; not multidimensional elasticity"
        ),
    )
    summary: dict[str, float | int | str] = {
        **preview,
        "history_csv": str(history_csv),
        "field": field,
        "step": selected_step,
        "time_s": time_s,
        "tensile_axis": ",".join(f"{value:.17g}" for value in projection.tensile_axis),
        "mapping": projection.source_role,
        "preview_path": str(resolved_preview_path),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return summary


def _parse_axis(text: str) -> np.ndarray:
    values = np.asarray([float(value.strip()) for value in text.split(",")], dtype=float)
    if values.shape not in ((2,), (3,)):
        raise argparse.ArgumentTypeError("axis must contain 2 or 3 comma-separated values")
    return values


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Project a 1D normal-only FEM/probability field onto an actual geometry mesh"
    )
    parser.add_argument("--history-csv", type=Path, required=True)
    parser.add_argument("--field", default="stress_pa")
    parser.add_argument("--step", default="peak-tension")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dimension", type=int, choices=(2, 3), default=3)
    parser.add_argument("--geometry", type=Path, default=None)
    parser.add_argument("--coordinate-scale-to-m", type=float, default=1.0)
    parser.add_argument("--characteristic-length-m", type=float, default=None)
    parser.add_argument("--length-m", type=float, default=0.05)
    parser.add_argument("--width-m", type=float, default=0.01)
    parser.add_argument("--thickness-m", type=float, default=0.001)
    parser.add_argument("--nx", type=int, default=32)
    parser.add_argument("--ny", type=int, default=6)
    parser.add_argument("--nz", type=int, default=2)
    parser.add_argument("--axis", default=None, help="comma-separated tensile-axis components")
    parser.add_argument("--preview-path", type=Path, default=None)
    args = parser.parse_args(argv)
    try:
        selected_step: str | int = int(args.step)
    except ValueError:
        selected_step = args.step
    axis = None if args.axis is None else _parse_axis(args.axis)
    summary = run_projection(
        history_csv=args.history_csv,
        field=args.field,
        step=selected_step,
        output_dir=args.output_dir,
        dimension=args.dimension,
        length_m=args.length_m,
        width_m=args.width_m,
        thickness_m=args.thickness_m,
        nx=args.nx,
        ny=args.ny,
        nz=args.nz,
        geometry_path=args.geometry,
        coordinate_scale_to_m=args.coordinate_scale_to_m,
        characteristic_length_m=args.characteristic_length_m,
        tensile_axis=axis,
        preview_path=args.preview_path,
    )
    for key, value in summary.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
