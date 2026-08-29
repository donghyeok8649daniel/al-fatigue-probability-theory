# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: TensionViewer
# - 주요 함수/메서드: load_fem_history, axial_snapshot, extruded_scalar_2d, _deformed_x, plot_tension_2d
#   _element_prism_faces, plot_tension_3d, _field_range, save_preview_images, TensionViewer.__init__
#   TensionViewer._on_slider, TensionViewer._on_view, TensionViewer._on_field, TensionViewer.redraw
#   TensionViewer.show, launch_ui, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""2D/3D presentation layer for the strictly one-dimensional tensile FEM.

The geometry is extruded only for visualization.  No transverse degree of
freedom, shear stress, von-Mises stress, Poisson contraction, or multiaxial
failure criterion is introduced here.  Each displayed cross-section simply
inherits the scalar axial field from the underlying 1D element.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.widgets import RadioButtons, Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from simulations.visualize_fem1d import load_numeric_csv, select_snapshot_step


def load_fem_history(input_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load node and element histories emitted by the C 1D FEM solver."""
    return (
        load_numeric_csv(input_dir / "nodes.csv"),
        load_numeric_csv(input_dir / "elements.csv"),
    )


def axial_snapshot(
    nodes: np.ndarray,
    elements: np.ndarray,
    step: int,
    field: str = "stress",
) -> dict[str, np.ndarray | float | int | str]:
    """Return one snapshot containing only axial 1D mechanical quantities."""
    node_rows = nodes[nodes["step"] == step]
    elem_rows = elements[elements["step"] == step]
    if node_rows.size == 0 or elem_rows.size == 0:
        raise ValueError(f"step {step} is absent from FEM history")

    if field == "stress":
        scalar = np.asarray(elem_rows["stress_pa"], dtype=float) / 1.0e6
        label = "axial stress [MPa]"
    elif field == "strain":
        scalar = np.asarray(elem_rows["strain"], dtype=float)
        label = "axial strain"
    else:
        raise ValueError("field must be 'stress' or 'strain'")

    return {
        "step": int(step),
        "time_s": float(node_rows["time_s"][0]),
        "x_nodes_m": np.asarray(node_rows["x_m"], dtype=float),
        "u_nodes_m": np.asarray(node_rows["displacement_m"], dtype=float),
        "x_mid_m": np.asarray(elem_rows["x_mid_m"], dtype=float),
        "scalar": scalar,
        "field": field,
        "label": label,
    }


def extruded_scalar_2d(
    x_nodes_m: np.ndarray,
    element_scalar: np.ndarray,
    half_width_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extrude a 1D element scalar uniformly across a visual 2D width."""
    x_nodes = np.asarray(x_nodes_m, dtype=float)
    values = np.asarray(element_scalar, dtype=float)
    if x_nodes.ndim != 1 or values.ndim != 1 or x_nodes.size != values.size + 1:
        raise ValueError("x_nodes must contain exactly one more entry than element_scalar")
    if half_width_m <= 0.0:
        raise ValueError("half_width_m must be positive")
    x_grid = np.vstack([x_nodes, x_nodes])
    y_grid = np.vstack([
        np.full_like(x_nodes, -half_width_m),
        np.full_like(x_nodes, half_width_m),
    ])
    field_grid = values[np.newaxis, :]
    return x_grid, y_grid, field_grid


def _deformed_x(snapshot: dict, deformation_scale: float) -> np.ndarray:
    return np.asarray(snapshot["x_nodes_m"]) + deformation_scale * np.asarray(snapshot["u_nodes_m"])


def plot_tension_2d(
    ax,
    snapshot: dict,
    half_width_m: float,
    deformation_scale: float,
    norm: Normalize,
    cmap,
):
    """Render a 2D extruded tensile bar without introducing 2D mechanics."""
    x_def = _deformed_x(snapshot, deformation_scale)
    x_grid, y_grid, field_grid = extruded_scalar_2d(
        x_def,
        np.asarray(snapshot["scalar"]),
        half_width_m,
    )
    mesh = ax.pcolormesh(x_grid, y_grid, field_grid, shading="flat", cmap=cmap, norm=norm)
    ax.plot(x_def, np.full_like(x_def, -half_width_m), linewidth=0.8)
    ax.plot(x_def, np.full_like(x_def, half_width_m), linewidth=0.8)
    ax.set_xlabel("tensile axis x [m]")
    ax.set_ylabel("visual width y [m]")
    ax.set_title(f"2D tensile-only view — t={snapshot['time_s']:.6g} s")
    ax.set_aspect("auto")
    return mesh


def _element_prism_faces(x0: float, x1: float, hw: float, ht: float):
    a = (x0, -hw, -ht)
    b = (x1, -hw, -ht)
    c = (x1, hw, -ht)
    d = (x0, hw, -ht)
    e = (x0, -hw, ht)
    f = (x1, -hw, ht)
    g = (x1, hw, ht)
    h = (x0, hw, ht)
    return [
        [a, b, c, d],
        [e, f, g, h],
        [a, b, f, e],
        [d, c, g, h],
        [a, d, h, e],
        [b, c, g, f],
    ]


def plot_tension_3d(
    ax,
    snapshot: dict,
    half_width_m: float,
    half_thickness_m: float,
    deformation_scale: float,
    norm: Normalize,
    cmap,
):
    """Render a rectangular 3D extrusion colored only by the axial scalar."""
    if half_width_m <= 0.0 or half_thickness_m <= 0.0:
        raise ValueError("visual half-width and half-thickness must be positive")
    x_def = _deformed_x(snapshot, deformation_scale)
    values = np.asarray(snapshot["scalar"], dtype=float)

    all_faces = []
    all_colors = []
    for i, value in enumerate(values):
        faces = _element_prism_faces(
            float(x_def[i]),
            float(x_def[i + 1]),
            half_width_m,
            half_thickness_m,
        )
        color = cmap(norm(float(value)))
        all_faces.extend(faces)
        all_colors.extend([color] * len(faces))

    collection = Poly3DCollection(
        all_faces,
        facecolors=all_colors,
        linewidths=0.25,
        edgecolors="k",
        alpha=0.95,
    )
    ax.add_collection3d(collection)
    span = max(float(x_def[-1] - x_def[0]), 1.0e-12)
    ax.set_xlim(float(x_def[0] - 0.02 * span), float(x_def[-1] + 0.02 * span))
    ax.set_ylim(-1.3 * half_width_m, 1.3 * half_width_m)
    ax.set_zlim(-1.3 * half_thickness_m, 1.3 * half_thickness_m)
    ax.set_xlabel("tensile axis x [m]")
    ax.set_ylabel("visual y [m]")
    ax.set_zlabel("visual z [m]")
    ax.set_title(f"3D tensile-only view — t={snapshot['time_s']:.6g} s")
    return cm.ScalarMappable(norm=norm, cmap=cmap)


def _field_range(elements: np.ndarray, field: str) -> tuple[float, float]:
    if field == "stress":
        values = np.asarray(elements["stress_pa"], dtype=float) / 1.0e6
    elif field == "strain":
        values = np.asarray(elements["strain"], dtype=float)
    else:
        raise ValueError("field must be 'stress' or 'strain'")
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if np.isclose(vmin, vmax):
        delta = max(abs(vmin), 1.0) * 1.0e-9
        vmin -= delta
        vmax += delta
    return vmin, vmax


def save_preview_images(
    input_dir: Path,
    output_dir: Path,
    half_width_m: float = 0.0025,
    half_thickness_m: float = 0.001,
    deformation_scale: float = 1.0,
    field: str = "stress",
) -> dict[str, float | int]:
    """Generate deterministic 2D and 3D peak-tension preview images."""
    nodes, elements = load_fem_history(input_dir)
    step = select_snapshot_step(elements, "peak-tension")
    snapshot = axial_snapshot(nodes, elements, step, field)
    vmin, vmax = _field_range(elements, field)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.0, 3.4))
    artist = plot_tension_2d(ax, snapshot, half_width_m, deformation_scale, norm, cmap)
    fig.colorbar(artist, ax=ax, label=snapshot["label"])
    fig.tight_layout()
    fig.savefig(output_dir / "tension_2d_peak.png", dpi=170)
    plt.close(fig)

    fig = plt.figure(figsize=(8.0, 5.2))
    ax3 = fig.add_subplot(111, projection="3d")
    artist3 = plot_tension_3d(
        ax3,
        snapshot,
        half_width_m,
        half_thickness_m,
        deformation_scale,
        norm,
        cmap,
    )
    fig.colorbar(artist3, ax=ax3, shrink=0.75, pad=0.1, label=snapshot["label"])
    fig.tight_layout()
    fig.savefig(output_dir / "tension_3d_peak.png", dpi=170)
    plt.close(fig)

    return {
        "step": int(step),
        "time_s": float(snapshot["time_s"]),
        "field_min": float(np.min(snapshot["scalar"])),
        "field_max": float(np.max(snapshot["scalar"])),
    }


class TensionViewer:
    """Small matplotlib UI with a time slider and 2D/3D display switch."""

    def __init__(
        self,
        input_dir: Path,
        half_width_m: float,
        half_thickness_m: float,
        deformation_scale: float,
    ) -> None:
        self.nodes, self.elements = load_fem_history(input_dir)
        self.steps = np.unique(self.elements["step"]).astype(int)
        self.half_width_m = half_width_m
        self.half_thickness_m = half_thickness_m
        self.deformation_scale = deformation_scale
        self.view = "2D"
        self.field = "stress"
        self.fig = plt.figure(figsize=(10.0, 6.2))
        self.main_ax = None
        self.colorbar = None

        slider_ax = self.fig.add_axes([0.24, 0.05, 0.52, 0.03])
        self.slider = Slider(
            slider_ax,
            "time step",
            0,
            len(self.steps) - 1,
            valinit=0,
            valstep=1,
        )
        view_ax = self.fig.add_axes([0.02, 0.06, 0.12, 0.11])
        field_ax = self.fig.add_axes([0.84, 0.06, 0.13, 0.11])
        self.view_radio = RadioButtons(view_ax, ("2D", "3D"), active=0)
        self.field_radio = RadioButtons(field_ax, ("stress", "strain"), active=0)
        self.slider.on_changed(self._on_slider)
        self.view_radio.on_clicked(self._on_view)
        self.field_radio.on_clicked(self._on_field)
        self.redraw()

    def _on_slider(self, _value) -> None:
        self.redraw()

    def _on_view(self, label: str) -> None:
        self.view = label
        self.redraw()

    def _on_field(self, label: str) -> None:
        self.field = label
        self.redraw()

    def redraw(self) -> None:
        if self.colorbar is not None:
            self.colorbar.remove()
            self.colorbar = None
        if self.main_ax is not None:
            self.fig.delaxes(self.main_ax)

        index = int(round(self.slider.val))
        step = int(self.steps[index])
        snapshot = axial_snapshot(self.nodes, self.elements, step, self.field)
        vmin, vmax = _field_range(self.elements, self.field)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("viridis")

        if self.view == "2D":
            self.main_ax = self.fig.add_axes([0.16, 0.18, 0.68, 0.72])
            artist = plot_tension_2d(
                self.main_ax,
                snapshot,
                self.half_width_m,
                self.deformation_scale,
                norm,
                cmap,
            )
        else:
            self.main_ax = self.fig.add_axes([0.16, 0.18, 0.68, 0.72], projection="3d")
            artist = plot_tension_3d(
                self.main_ax,
                snapshot,
                self.half_width_m,
                self.half_thickness_m,
                self.deformation_scale,
                norm,
                cmap,
            )
        self.colorbar = self.fig.colorbar(artist, ax=self.main_ax, shrink=0.76, pad=0.08)
        self.colorbar.set_label(snapshot["label"])
        self.fig.canvas.draw_idle()

    def show(self) -> None:
        plt.show()


def launch_ui(
    input_dir: Path,
    half_width_m: float = 0.0025,
    half_thickness_m: float = 0.001,
    deformation_scale: float = 1.0,
) -> None:
    """Launch the interactive tensile-only 2D/3D viewer."""
    TensionViewer(input_dir, half_width_m, half_thickness_m, deformation_scale).show()


def main() -> None:
    parser = argparse.ArgumentParser(description="2D/3D UI for strictly 1D tensile FEM results")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures/fem1d_demo"))
    parser.add_argument("--half-width-m", type=float, default=0.0025)
    parser.add_argument("--half-thickness-m", type=float, default=0.001)
    parser.add_argument("--deformation-scale", type=float, default=1.0)
    parser.add_argument("--field", choices=("stress", "strain"), default="stress")
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--save-previews", action="store_true")
    args = parser.parse_args()

    if not args.interactive and not args.save_previews:
        args.save_previews = True

    if args.save_previews:
        summary = save_preview_images(
            args.input_dir,
            args.output_dir,
            args.half_width_m,
            args.half_thickness_m,
            args.deformation_scale,
            args.field,
        )
        print(summary)
    if args.interactive:
        launch_ui(
            args.input_dir,
            args.half_width_m,
            args.half_thickness_m,
            args.deformation_scale,
        )


if __name__ == "__main__":
    main()
