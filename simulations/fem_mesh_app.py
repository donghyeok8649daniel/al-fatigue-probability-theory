# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D/2D/3D FEM geometry mesh 생성·CAD 입력·가시성 조절을 제공하는 경량 Matplotlib UI다.
# - 주요 클래스: MeshAppConfig, FEMMeshApp
# - 주요 함수/메서드: MeshAppConfig.length_m, MeshAppConfig.width_m, MeshAppConfig.thickness_m
#   MeshAppConfig.mesh_size_m, validate_mesh_app_config, create_geometry_mesh, axial_visibility_field
#   clipped_mesh_view, FEMMeshApp.__init__, FEMMeshApp._create_controls, FEMMeshApp._read_config
#   FEMMeshApp._set_status, FEMMeshApp._on_dimension, FEMMeshApp._on_visibility, FEMMeshApp._browse
#   FEMMeshApp._generate, FEMMeshApp._clear_axes, FEMMeshApp._redraw, FEMMeshApp._save, FEMMeshApp.show
#   run_headless_smoke, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Lightweight 1D/2D/3D FEM geometry-meshing application.

The application depends only on NumPy and Matplotlib in its default mode.
Optional Gmsh/meshio backends are loaded only when a corresponding CAD or mesh
file is requested.  Cell colors show normalized axial position so connectivity
and clipping are visible before a mechanics result is attached.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider, TextBox
import numpy as np

from simulations.fem_geometry_mesh import (
    GeometryMesh,
    cell_centers,
    load_geometry_mesh,
    save_mesh_npz,
    structured_box_mesh,
    structured_line_mesh,
    structured_rectangle_mesh,
    subset_cells,
)
from simulations.fem_mesh_ui import (
    plot_geometry_mesh_1d,
    plot_geometry_mesh_2d,
    plot_geometry_mesh_3d,
    save_geometry_mesh_preview,
)


@dataclass(frozen=True)
class MeshAppConfig:
    """Geometry and discretization inputs independent of fatigue parameters."""

    dimension: int = 3
    length_mm: float = 50.0
    width_mm: float = 10.0
    thickness_mm: float = 1.0
    nx: int = 24
    ny: int = 4
    nz: int = 2
    geometry_path: str = ""
    coordinate_scale_to_m: float = 1.0
    mesh_size_mm: float = 1.0

    @property
    def length_m(self) -> float:
        return self.length_mm * 1.0e-3

    @property
    def width_m(self) -> float:
        return self.width_mm * 1.0e-3

    @property
    def thickness_m(self) -> float:
        return self.thickness_mm * 1.0e-3

    @property
    def mesh_size_m(self) -> float:
        return self.mesh_size_mm * 1.0e-3


def validate_mesh_app_config(config: MeshAppConfig) -> None:
    if config.dimension not in (1, 2, 3):
        raise ValueError("dimension must be 1, 2, or 3")
    for name, value in {
        "length_mm": config.length_mm,
        "width_mm": config.width_mm,
        "thickness_mm": config.thickness_mm,
        "coordinate_scale_to_m": config.coordinate_scale_to_m,
        "mesh_size_mm": config.mesh_size_mm,
    }.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if config.nx < 1 or config.ny < 1 or config.nz < 1:
        raise ValueError("nx, ny, and nz must be at least 1")


def create_geometry_mesh(config: MeshAppConfig) -> GeometryMesh:
    """Generate or import the requested topological mesh dimension."""
    validate_mesh_app_config(config)
    if config.geometry_path.strip():
        return load_geometry_mesh(
            Path(config.geometry_path.strip()),
            coordinate_scale_to_m=config.coordinate_scale_to_m,
            target_dimension=config.dimension,
            characteristic_length_m=config.mesh_size_m,
        )
    if config.dimension == 1:
        return structured_line_mesh(config.length_m, config.nx)
    if config.dimension == 2:
        return structured_rectangle_mesh(config.length_m, config.width_m, config.nx, config.ny)
    return structured_box_mesh(
        config.length_m,
        config.width_m,
        config.thickness_m,
        config.nx,
        config.ny,
        config.nz,
    )


def axial_visibility_field(mesh: GeometryMesh) -> np.ndarray:
    """Return normalized x-position solely for mesh-visibility coloring."""
    centers = cell_centers(mesh)
    x = centers[:, 0]
    span = float(np.ptp(x))
    if span <= 0.0:
        if mesh.cell_count == 1:
            return np.zeros(1, dtype=float)
        return np.linspace(0.0, 1.0, mesh.cell_count)
    return (x - float(np.min(x))) / span


def clipped_mesh_view(
    mesh: GeometryMesh,
    values: np.ndarray,
    fraction: float,
) -> tuple[GeometryMesh, np.ndarray]:
    """Clip cells by axial centroid to expose interior volume connectivity."""
    scalar = np.asarray(values, dtype=float)
    if scalar.shape != (mesh.cell_count,):
        raise ValueError("values must match mesh cell count")
    fraction = float(np.clip(fraction, 0.0, 1.0))
    if fraction >= 1.0 - 1.0e-12:
        return mesh, scalar
    centers = cell_centers(mesh)
    x = centers[:, 0]
    threshold = float(np.min(x) + fraction * (np.max(x) - np.min(x)))
    keep = x <= threshold + 1.0e-15
    if not np.any(keep):
        keep[int(np.argmin(x))] = True
    clipped, original_indices = subset_cells(mesh, keep)
    return clipped, scalar[original_indices]


class FEMMeshApp:
    """Matplotlib UI with lightweight mesh generation and visibility controls."""

    _INPUTS = (
        ("length_mm", "Length [mm]", "50"),
        ("width_mm", "Width [mm]", "10"),
        ("thickness_mm", "Thickness [mm]", "1"),
        ("nx", "Nx", "24"),
        ("ny", "Ny", "4"),
        ("nz", "Nz", "2"),
        ("geometry_path", "CAD/mesh path", ""),
        ("coordinate_scale_to_m", "CAD scale -> m", "1"),
        ("mesh_size_mm", "CAD mesh size [mm]", "1"),
    )

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = Path(output_dir)
        self.dimension = 3
        self.mesh: GeometryMesh | None = None
        self.values: np.ndarray | None = None
        self.main_ax = None
        self.colorbar = None
        self.show_nodes = False
        self.show_edges = True

        self.fig = plt.figure(figsize=(13.2, 7.6))
        self.fig.canvas.manager.set_window_title("Lightweight FEM Mesh — 1D / 2D / 3D")
        self.textboxes: dict[str, TextBox] = {}
        self._create_controls()
        self.status = self.fig.text(
            0.255,
            0.975,
            "Ready — geometry mesh only; active fatigue theory remains 1D normal-only",
            ha="left",
            va="top",
            fontsize=8.5,
        )
        self._generate(None)

    def _create_controls(self) -> None:
        self.fig.text(0.025, 0.955, "Mesh inputs", weight="bold", va="top")
        dim_ax = self.fig.add_axes([0.025, 0.835, 0.18, 0.095])
        self.dimension_radio = RadioButtons(dim_ax, ("1D", "2D", "3D"), active=2)
        self.dimension_radio.on_clicked(self._on_dimension)
        top = 0.785
        spacing = 0.056
        for index, (key, label, initial) in enumerate(self._INPUTS):
            y = top - index * spacing
            self.fig.text(0.025, y + 0.031, label, fontsize=8, va="center")
            ax = self.fig.add_axes([0.025, y, 0.18, 0.029])
            self.textboxes[key] = TextBox(ax, "", initial=initial)

        generate_ax = self.fig.add_axes([0.025, 0.245, 0.084, 0.042])
        self.generate_button = Button(generate_ax, "Generate")
        self.generate_button.on_clicked(self._generate)
        save_ax = self.fig.add_axes([0.12, 0.245, 0.084, 0.042])
        self.save_button = Button(save_ax, "Save")
        self.save_button.on_clicked(self._save)
        browse_ax = self.fig.add_axes([0.025, 0.192, 0.179, 0.038])
        self.browse_button = Button(browse_ax, "Browse CAD/mesh")
        self.browse_button.on_clicked(self._browse)

        visibility_ax = self.fig.add_axes([0.025, 0.088, 0.18, 0.08])
        self.visibility_checks = CheckButtons(
            visibility_ax,
            ("show nodes", "show edges"),
            (self.show_nodes, self.show_edges),
        )
        self.visibility_checks.on_clicked(self._on_visibility)

        opacity_ax = self.fig.add_axes([0.33, 0.055, 0.25, 0.025])
        self.opacity_slider = Slider(opacity_ax, "opacity", 0.15, 1.0, valinit=0.9)
        self.opacity_slider.on_changed(self._redraw)
        clip_ax = self.fig.add_axes([0.68, 0.055, 0.25, 0.025])
        self.clip_slider = Slider(clip_ax, "axial clip", 0.04, 1.0, valinit=1.0)
        self.clip_slider.on_changed(self._redraw)
        self.fig.text(
            0.025,
            0.025,
            "Core: NumPy + Matplotlib | STEP/IGES: optional Gmsh",
            fontsize=8,
        )

    def _read_config(self) -> MeshAppConfig:
        def floating(key: str) -> float:
            return float(self.textboxes[key].text.strip())

        def integer(key: str) -> int:
            value = floating(key)
            if not value.is_integer():
                raise ValueError(f"{key} must be an integer")
            return int(value)

        config = MeshAppConfig(
            dimension=self.dimension,
            length_mm=floating("length_mm"),
            width_mm=floating("width_mm"),
            thickness_mm=floating("thickness_mm"),
            nx=integer("nx"),
            ny=integer("ny"),
            nz=integer("nz"),
            geometry_path=self.textboxes["geometry_path"].text.strip(),
            coordinate_scale_to_m=floating("coordinate_scale_to_m"),
            mesh_size_mm=floating("mesh_size_mm"),
        )
        validate_mesh_app_config(config)
        return config

    def _set_status(self, message: str) -> None:
        self.status.set_text(message)
        self.fig.canvas.draw_idle()

    def _on_dimension(self, label: str) -> None:
        self.dimension = int(label[0])
        self._generate(None)

    def _on_visibility(self, _label: str) -> None:
        self.show_nodes, self.show_edges = self.visibility_checks.get_status()
        self._redraw(None)

    def _browse(self, _event) -> None:
        try:
            from tkinter import Tk, filedialog

            root = Tk()
            root.withdraw()
            selected = filedialog.askopenfilename(
                title="Select CAD or mesh",
                filetypes=[
                    ("CAD/mesh", "*.stl *.obj *.step *.stp *.iges *.igs *.brep *.msh *.vtk *.vtu *.xdmf *.npz"),
                    ("All files", "*.*"),
                ],
            )
            root.destroy()
            if selected:
                self.textboxes["geometry_path"].set_val(selected)
        except Exception as exc:
            self._set_status(f"Browse unavailable: enter the path manually ({exc})")

    def _generate(self, _event) -> None:
        try:
            config = self._read_config()
            self.mesh = create_geometry_mesh(config)
            self.values = axial_visibility_field(self.mesh)
            self.clip_slider.set_val(1.0)
            self._set_status(
                f"Loaded {self.mesh.topological_dimension}D mesh: "
                f"{self.mesh.points_m.shape[0]} nodes, {self.mesh.cell_count} cells — {self.mesh.source}"
            )
            self._redraw(None)
        except Exception as exc:
            self._set_status(f"ERROR: {exc}")

    def _clear_axes(self) -> None:
        if self.colorbar is not None:
            try:
                self.colorbar.remove()
            except Exception:
                pass
            self.colorbar = None
        if self.main_ax is not None:
            try:
                self.fig.delaxes(self.main_ax)
            except Exception:
                pass
            self.main_ax = None

    def _redraw(self, _event) -> None:
        if self.mesh is None or self.values is None:
            return
        visible_mesh, visible_values = clipped_mesh_view(
            self.mesh,
            self.values,
            float(self.clip_slider.val),
        )
        self._clear_axes()
        spatial_3d = visible_mesh.topological_dimension == 3 or visible_mesh.embedding_dimension == 3
        self.main_ax = self.fig.add_axes(
            [0.25, 0.14, 0.72, 0.70],
            projection="3d" if spatial_3d else None,
        )
        if visible_mesh.topological_dimension == 1:
            artist = plot_geometry_mesh_1d(
                self.main_ax,
                visible_mesh,
                visible_values,
                show_nodes=self.show_nodes,
            )
        elif visible_mesh.topological_dimension == 2 and visible_mesh.embedding_dimension == 2:
            artist = plot_geometry_mesh_2d(
                self.main_ax,
                visible_mesh,
                visible_values,
                show_edges=self.show_edges,
                show_nodes=self.show_nodes,
                alpha=float(self.opacity_slider.val),
            )
        else:
            artist = plot_geometry_mesh_3d(
                self.main_ax,
                visible_mesh,
                visible_values,
                show_edges=self.show_edges,
                show_nodes=self.show_nodes,
                alpha=float(self.opacity_slider.val),
            )
        self.main_ax.set_title(
            f"{visible_mesh.topological_dimension}D FEM geometry mesh — "
            f"showing {visible_mesh.cell_count}/{self.mesh.cell_count} cells\n"
            "color = normalized axial position; not a mechanics result",
            fontsize=11,
        )
        self.colorbar = self.fig.colorbar(artist, ax=self.main_ax, shrink=0.75, pad=0.08)
        self.colorbar.set_label("normalized axial cell coordinate")
        self.fig.canvas.draw_idle()

    def _save(self, _event) -> None:
        if self.mesh is None or self.values is None:
            self._set_status("Generate or load a mesh before saving.")
            return
        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            save_mesh_npz(self.output_dir / "geometry_mesh.npz", self.mesh)
            summary = save_geometry_mesh_preview(
                self.output_dir / "mesh_preview.png",
                self.mesh,
                self.values,
                field_label="normalized axial cell coordinate",
                title=(
                    f"{self.mesh.topological_dimension}D FEM geometry mesh\n"
                    "visibility field only; not a mechanics result"
                ),
            )
            (self.output_dir / "summary.json").write_text(
                json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            self._set_status(f"Saved mesh bundle and preview to {self.output_dir}")
        except Exception as exc:
            self._set_status(f"ERROR saving mesh: {exc}")

    def show(self) -> None:
        plt.show()


def run_headless_smoke(
    output_dir: Path,
    preview_dir: Path | None = None,
) -> dict[str, int]:
    """Generate and render all supported dimensions without opening a window."""
    output = Path(output_dir)
    previews = output if preview_dir is None else Path(preview_dir)
    cell_counts: dict[str, int] = {}
    configs = {
        1: MeshAppConfig(dimension=1, nx=8),
        2: MeshAppConfig(dimension=2, nx=8, ny=3),
        3: MeshAppConfig(dimension=3, nx=8, ny=3, nz=2),
    }
    for dimension, config in configs.items():
        mesh = create_geometry_mesh(config)
        values = axial_visibility_field(mesh)
        directory = output / f"mesh_{dimension}d"
        save_mesh_npz(directory / "geometry_mesh.npz", mesh)
        save_geometry_mesh_preview(
            directory / "mesh_preview.png"
            if preview_dir is None
            else previews / f"mesh_{dimension}d.png",
            mesh,
            values,
            field_label="normalized axial cell coordinate",
            title=f"{dimension}D lightweight mesh smoke test",
        )
        cell_counts[f"cells_{dimension}d"] = mesh.cell_count
    app = FEMMeshApp(output / "ui_saved_mesh")
    if not app.show_nodes:
        app.visibility_checks.set_active(0)
    app.clip_slider.set_val(0.55)
    previews.mkdir(parents=True, exist_ok=True)
    app.fig.savefig(previews / "lightweight_mesh_ui.png", dpi=160)
    plt.close(app.fig)
    return cell_counts


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Lightweight 1D/2D/3D FEM geometry mesh UI")
    parser.add_argument("--output-dir", type=Path, default=Path("results/data/fem_mesh_ui"))
    parser.add_argument("--preview-dir", type=Path, default=None)
    parser.add_argument("--headless-smoke", action="store_true")
    args = parser.parse_args(argv)
    if args.headless_smoke:
        for key, value in run_headless_smoke(args.output_dir, args.preview_dir).items():
            print(f"{key}={value}")
        return
    FEMMeshApp(args.output_dir).show()


if __name__ == "__main__":
    main()
