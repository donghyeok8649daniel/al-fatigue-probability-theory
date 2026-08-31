# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: TensionRunConfig, FEMTensionApp
# - 주요 함수/메서드: TensionRunConfig.length_m, TensionRunConfig.width_m, TensionRunConfig.thickness_m
#   TensionRunConfig.area_m2, TensionRunConfig.young_pa, validate_run_config, repository_root
#   solver_executable, _solver_sources, _needs_rebuild, build_fem_solver, solver_command, run_fem_solver
#   FEMTensionApp.__init__, FEMTensionApp._create_parameter_panel, FEMTensionApp._create_result_controls
#   FEMTensionApp._read_config, FEMTensionApp._set_status, FEMTensionApp._on_run, FEMTensionApp._on_save
#   FEMTensionApp._on_slider, FEMTensionApp._on_view, FEMTensionApp._on_field
#   FEMTensionApp._clear_main_axes, FEMTensionApp.redraw, FEMTensionApp.show, run_headless_smoke, main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Integrated GUI for the standalone one-dimensional tensile FEM scaffold.

The application deliberately keeps mechanics one-dimensional.  Width and
thickness define the axial cross-sectional area and the display extrusion,
but they do not introduce transverse degrees of freedom or constitutive laws.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import os
import shutil
import subprocess
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.widgets import Button, RadioButtons, Slider, TextBox
import numpy as np

from simulations.fem_tension_ui import (
    _field_range,
    axial_snapshot,
    load_fem_history,
    plot_tension_2d,
    plot_tension_3d,
    save_preview_images,
)
from simulations.ftgsim_format import create_ftgsim, extract_results, open_ftgsim


@dataclass(frozen=True)
class TensionRunConfig:
    """User-facing inputs for one quasistatic cyclic axial-tension solve."""

    length_mm: float = 50.0
    width_mm: float = 10.0
    thickness_mm: float = 1.0
    young_gpa: float = 69.0
    elements: int = 40
    stress_mean_mpa: float = 50.0
    stress_amplitude_mpa: float = 100.0
    frequency_hz: float = 20.0
    cycles: int = 2
    steps_per_cycle: int = 80
    deformation_scale: float = 1.0

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
    def area_m2(self) -> float:
        return self.width_m * self.thickness_m

    @property
    def young_pa(self) -> float:
        return self.young_gpa * 1.0e9


def validate_run_config(config: TensionRunConfig) -> None:
    """Reject nonphysical or numerically meaningless GUI inputs."""
    positive = {
        "length_mm": config.length_mm,
        "width_mm": config.width_mm,
        "thickness_mm": config.thickness_mm,
        "young_gpa": config.young_gpa,
        "frequency_hz": config.frequency_hz,
        "deformation_scale": config.deformation_scale,
    }
    for name, value in positive.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if config.elements < 1:
        raise ValueError("elements must be at least 1")
    if config.cycles < 1:
        raise ValueError("cycles must be at least 1")
    if config.steps_per_cycle < 2:
        raise ValueError("steps_per_cycle must be at least 2")
    if not np.isfinite(config.stress_mean_mpa):
        raise ValueError("stress_mean_mpa must be finite")
    if not np.isfinite(config.stress_amplitude_mpa) or config.stress_amplitude_mpa < 0.0:
        raise ValueError("stress_amplitude_mpa must be finite and nonnegative")


def config_from_ftgsim(path: Path) -> tuple[TensionRunConfig, dict, dict]:
    """Load and validate the tensile-only setup from an `.ftgsim` bundle."""
    bundle = open_ftgsim(path)
    if bundle.setup.get("physics_model") != "1d_normal_tensile":
        raise ValueError("ftgsim project is not a 1D normal-tensile model")
    values = bundle.setup.get("tension_run")
    if not isinstance(values, dict):
        raise ValueError("ftgsim setup is missing tension_run")
    allowed = set(TensionRunConfig.__dataclass_fields__)
    unknown = set(values) - allowed
    if unknown:
        raise ValueError(f"unknown tension_run fields: {sorted(unknown)}")
    config = TensionRunConfig(**values)
    validate_run_config(config)
    return config, bundle.geometry, bundle.display


def save_tension_ftgsim(path: Path, config: TensionRunConfig, output_dir: Path | None = None,
                        *, view: str = "2D", field: str = "stress") -> Path:
    """Save setup and any existing FEM CSV results as an open `.ftgsim` bundle."""
    validate_run_config(config)
    files: dict[str, Path] = {}
    if output_dir is not None:
        source_dir = Path(output_dir)
        for name in ("nodes.csv", "elements.csv", "metadata.csv", "summary.json"):
            source = source_dir / name
            if source.is_file():
                files[f"results/{name}"] = source
    return create_ftgsim(
        path,
        setup={
            "physics_model": "1d_normal_tensile",
            "material_scope": "pure_single_crystal_aluminum",
            "tension_run": asdict(config),
            "probability_model": {
                "enabled": False,
                "coordinate": "local_homogeneous_spacing",
                "energy": "exact_riemann_zeta_bulk_lattice",
                "initiation_definition": "first_tangent_stiffness_loss",
                "critical_stretch_rule": "((m+1)/(n+1))**(1/(m-n))",
                "note": "Parameters and escape coupling must be supplied before activation.",
            },
            "result_references": sorted(files),
        },
        geometry={
            "mesh_dimension": 1,
            "loading_axis": [1.0, 0.0, 0.0],
            "geometry_kind": "uniform_bar",
            "length_mm": config.length_mm,
            "width_mm": config.width_mm,
            "thickness_mm": config.thickness_mm,
            "elements": config.elements,
        },
        display={"view": view, "field": field, "deformation_scale": config.deformation_scale},
        files=files,
        generator={"application": "fem_tension_app", "format_extension": ".ftgsim"},
    )


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def solver_executable(repo_root: Path | None = None) -> Path:
    """Return the platform-specific expected C solver path."""
    root = repository_root() if repo_root is None else Path(repo_root)
    suffix = ".exe" if os.name == "nt" else ""
    return root / "fem1d" / "bin" / f"fem1d_solver{suffix}"


def _solver_sources(repo_root: Path) -> tuple[Path, ...]:
    return (
        repo_root / "fem1d" / "src" / "main.c",
        repo_root / "fem1d" / "src" / "fem1d.c",
        repo_root / "fem1d" / "include" / "fem1d.h",
    )


def _needs_rebuild(solver: Path, repo_root: Path) -> bool:
    if not solver.exists():
        return True
    solver_time = solver.stat().st_mtime
    return any(source.exists() and source.stat().st_mtime > solver_time for source in _solver_sources(repo_root))


def build_fem_solver(repo_root: Path | None = None, force: bool = False) -> Path:
    """Build the C solver with make when available, otherwise a C compiler directly."""
    root = repository_root() if repo_root is None else Path(repo_root)
    solver = solver_executable(root)
    if not force and not _needs_rebuild(solver, root):
        return solver

    solver.parent.mkdir(parents=True, exist_ok=True)
    build_errors: list[str] = []

    for make_name in ("make", "mingw32-make"):
        make = shutil.which(make_name)
        if make is None:
            continue
        completed = subprocess.run(
            [make, "-C", str(root / "fem1d")],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode == 0 and solver.exists():
            return solver
        build_errors.append(f"{make_name}: {completed.stderr.strip() or completed.stdout.strip()}")

    for compiler_name in ("cc", "gcc", "clang"):
        compiler = shutil.which(compiler_name)
        if compiler is None:
            continue
        command = [
            compiler,
            "-I",
            str(root / "fem1d" / "include"),
            "-O2",
            "-std=c11",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            str(root / "fem1d" / "src" / "main.c"),
            str(root / "fem1d" / "src" / "fem1d.c"),
            "-o",
            str(solver),
            "-lm",
        ]
        completed = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
        if completed.returncode == 0 and solver.exists():
            return solver
        build_errors.append(f"{compiler_name}: {completed.stderr.strip() or completed.stdout.strip()}")

    detail = "\n".join(error for error in build_errors if error)
    raise RuntimeError(
        "Could not build fem1d solver. Install make+gcc/clang or provide --solver."
        + (f"\n{detail}" if detail else "")
    )


def solver_command(config: TensionRunConfig, solver: Path, output_dir: Path) -> list[str]:
    """Translate GUI inputs into the C solver's explicit command-line interface."""
    validate_run_config(config)
    return [
        str(solver),
        "--elements",
        str(config.elements),
        "--length-m",
        f"{config.length_m:.17g}",
        "--area-m2",
        f"{config.area_m2:.17g}",
        "--young-pa",
        f"{config.young_pa:.17g}",
        "--stress-mean-mpa",
        f"{config.stress_mean_mpa:.17g}",
        "--stress-amplitude-mpa",
        f"{config.stress_amplitude_mpa:.17g}",
        "--frequency-hz",
        f"{config.frequency_hz:.17g}",
        "--cycles",
        str(config.cycles),
        "--steps-per-cycle",
        str(config.steps_per_cycle),
        "--outdir",
        str(output_dir),
    ]


def run_fem_solver(
    config: TensionRunConfig,
    output_dir: Path,
    solver: Path | None = None,
    auto_build: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the C solver and return its captured process result."""
    validate_run_config(config)
    root = repository_root()
    executable = Path(solver) if solver is not None else solver_executable(root)
    if not executable.exists():
        if not auto_build:
            raise FileNotFoundError(f"FEM solver not found: {executable}")
        executable = build_fem_solver(root)

    output_dir = Path(output_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    completed = subprocess.run(
        solver_command(config, executable, output_dir),
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "C FEM solve failed.\n"
            + (completed.stdout.strip() + "\n" if completed.stdout.strip() else "")
            + completed.stderr.strip()
        )
    for required in ("nodes.csv", "elements.csv", "metadata.csv"):
        if not (output_dir / required).exists():
            raise RuntimeError(f"C FEM finished without required output: {required}")
    return completed


class FEMTensionApp:
    """Matplotlib desktop app that runs the C solver and visualizes its axial output."""

    _INPUT_SPECS = (
        ("length_mm", "Length [mm]", "50"),
        ("width_mm", "Width [mm]", "10"),
        ("thickness_mm", "Thickness [mm]", "1"),
        ("young_gpa", "E [GPa]", "69"),
        ("elements", "Elements", "40"),
        ("stress_mean_mpa", "Mean stress [MPa]", "50"),
        ("stress_amplitude_mpa", "Amplitude [MPa]", "100"),
        ("frequency_hz", "Frequency [Hz]", "20"),
        ("cycles", "Cycles", "2"),
        ("steps_per_cycle", "Steps/cycle", "80"),
        ("deformation_scale", "Deformation scale", "1"),
    )

    def __init__(
        self,
        output_dir: Path,
        solver: Path | None = None,
        auto_build: bool = True,
        project_path: Path | None = None,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.solver = None if solver is None else Path(solver)
        self.auto_build = auto_build
        self.nodes: np.ndarray | None = None
        self.elements: np.ndarray | None = None
        self.steps = np.array([0], dtype=int)
        self.config = TensionRunConfig()
        self.view = "2D"
        self.field = "stress"
        self.main_ax = None
        self.colorbar = None

        self.fig = plt.figure(figsize=(13.4, 7.5))
        self.fig.canvas.manager.set_window_title("1D Tensile FEM — C Solver + 2D/3D Viewer")
        self.textboxes: dict[str, TextBox] = {}
        self._create_parameter_panel()
        self._create_result_controls()
        self.status = self.fig.text(
            0.245,
            0.955,
            "Ready — mechanics scope: 1D axial tension only",
            ha="left",
            va="top",
        )
        if project_path is not None:
            self._load_project(Path(project_path))
        self.redraw()

    def _apply_config_to_boxes(self) -> None:
        for key, _label, _initial in self._INPUT_SPECS:
            self.textboxes[key].set_val(str(getattr(self.config, key)))

    def _load_project(self, path: Path) -> None:
        self.config, _geometry, display = config_from_ftgsim(path)
        self._apply_config_to_boxes()
        self.view = display.get("view", "2D") if display.get("view") in {"2D", "3D"} else "2D"
        self.field = display.get("field", "stress") if display.get("field") in {"stress", "strain"} else "stress"
        bundle = open_ftgsim(path)
        extracted = extract_results(bundle, self.output_dir)
        if {"nodes.csv", "elements.csv"}.issubset({item.name for item in extracted}):
            self.nodes, self.elements = load_fem_history(self.output_dir)
            self.steps = np.unique(self.elements["step"]).astype(int)
            self.slider.valmax = max(len(self.steps) - 1, 1)
            self.slider.ax.set_xlim(self.slider.valmin, self.slider.valmax)
        self.status.set_text(f"Opened {path.name} (1D normal tension only)")

    def _create_parameter_panel(self) -> None:
        self.fig.text(0.025, 0.955, "Tensile test inputs", ha="left", va="top", weight="bold")
        top = 0.89
        spacing = 0.057
        for i, (key, label, initial) in enumerate(self._INPUT_SPECS):
            ax = self.fig.add_axes([0.025, top - i * spacing, 0.175, 0.035])
            box = TextBox(ax, label, initial=initial, label_pad=0.03)
            self.textboxes[key] = box

        run_ax = self.fig.add_axes([0.025, 0.205, 0.082, 0.048])
        self.run_button = Button(run_ax, "Run FEM")
        self.run_button.on_clicked(self._on_run)

        save_ax = self.fig.add_axes([0.118, 0.205, 0.082, 0.048])
        self.save_button = Button(save_ax, "Save views")
        self.save_button.on_clicked(self._on_save)

        project_ax = self.fig.add_axes([0.025, 0.145, 0.175, 0.040])
        self.project_button = Button(project_ax, "Save .ftgsim")
        self.project_button.on_clicked(self._on_save_project)

        self.fig.text(
            0.025,
            0.13,
            "A = width × thickness\n2D/3D = display only\nNo shear / von-Mises / Poisson model",
            ha="left",
            va="top",
            fontsize=8.5,
        )

    def _create_result_controls(self) -> None:
        slider_ax = self.fig.add_axes([0.33, 0.055, 0.43, 0.03])
        self.slider = Slider(slider_ax, "time step", 0, 1, valinit=0, valstep=1)
        self.slider.on_changed(self._on_slider)

        view_ax = self.fig.add_axes([0.79, 0.035, 0.075, 0.09])
        field_ax = self.fig.add_axes([0.89, 0.035, 0.085, 0.09])
        self.view_radio = RadioButtons(view_ax, ("2D", "3D"), active=0)
        self.field_radio = RadioButtons(field_ax, ("stress", "strain"), active=0)
        self.view_radio.on_clicked(self._on_view)
        self.field_radio.on_clicked(self._on_field)

    def _read_config(self) -> TensionRunConfig:
        def floating(key: str) -> float:
            return float(self.textboxes[key].text.strip())

        def integer(key: str) -> int:
            raw = float(self.textboxes[key].text.strip())
            if not raw.is_integer():
                raise ValueError(f"{key} must be an integer")
            return int(raw)

        config = TensionRunConfig(
            length_mm=floating("length_mm"),
            width_mm=floating("width_mm"),
            thickness_mm=floating("thickness_mm"),
            young_gpa=floating("young_gpa"),
            elements=integer("elements"),
            stress_mean_mpa=floating("stress_mean_mpa"),
            stress_amplitude_mpa=floating("stress_amplitude_mpa"),
            frequency_hz=floating("frequency_hz"),
            cycles=integer("cycles"),
            steps_per_cycle=integer("steps_per_cycle"),
            deformation_scale=floating("deformation_scale"),
        )
        validate_run_config(config)
        return config

    def _set_status(self, message: str) -> None:
        self.status.set_text(message)
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def _on_run(self, _event) -> None:
        try:
            self.config = self._read_config()
            self._set_status("Running C FEM...")
            completed = run_fem_solver(
                self.config,
                self.output_dir,
                solver=self.solver,
                auto_build=self.auto_build,
            )
            self.nodes, self.elements = load_fem_history(self.output_dir)
            self.steps = np.unique(self.elements["step"]).astype(int)
            self.slider.valmax = max(len(self.steps) - 1, 1)
            self.slider.ax.set_xlim(self.slider.valmin, self.slider.valmax)
            self.slider.valstep = 1
            self.slider.set_val(0)
            self._set_status(
                f"Solved: {self.config.elements} elements, A={self.config.area_m2:.4g} m² — "
                + (completed.stdout.strip().splitlines()[0] if completed.stdout.strip() else "C FEM complete")
            )
            self.redraw()
        except Exception as exc:  # GUI boundary: surface exact failure text to the user.
            self._set_status(f"ERROR: {exc}")

    def _on_save(self, _event) -> None:
        if self.nodes is None or self.elements is None:
            self._set_status("Run FEM before saving views.")
            return
        try:
            preview_dir = self.output_dir.parent / f"{self.output_dir.name}_figures"
            result = save_preview_images(
                self.output_dir,
                preview_dir,
                half_width_m=self.config.width_m / 2.0,
                half_thickness_m=self.config.thickness_m / 2.0,
                deformation_scale=self.config.deformation_scale,
                field=self.field,
            )
            self._set_status(f"Saved 2D/3D peak-tension views to {preview_dir} (step {result['step']}).")
        except Exception as exc:
            self._set_status(f"ERROR saving views: {exc}")

    def _on_save_project(self, _event) -> None:
        try:
            self.config = self._read_config()
            target = self.output_dir.parent / f"{self.output_dir.name}.ftgsim"
            saved = save_tension_ftgsim(target, self.config, self.output_dir,
                                        view=self.view, field=self.field)
            self._set_status(f"Saved project: {saved}")
        except Exception as exc:
            self._set_status(f"ERROR saving project: {exc}")

    def _on_slider(self, _value) -> None:
        self.redraw()

    def _on_view(self, label: str) -> None:
        self.view = label
        self.redraw()

    def _on_field(self, label: str) -> None:
        self.field = label
        self.redraw()

    def _clear_main_axes(self) -> None:
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

    def redraw(self) -> None:
        self._clear_main_axes()
        projection = "3d" if self.view == "3D" else None
        self.main_ax = self.fig.add_axes([0.245, 0.15, 0.72, 0.76], projection=projection)

        if self.nodes is None or self.elements is None:
            self.main_ax.text(
                0.5,
                0.5,
                "Enter tensile-test parameters and press Run FEM",
                transform=self.main_ax.transAxes,
                ha="center",
                va="center",
            )
            self.main_ax.set_axis_off()
            self.fig.canvas.draw_idle()
            return

        index = min(int(round(self.slider.val)), len(self.steps) - 1)
        step = int(self.steps[index])
        snapshot = axial_snapshot(self.nodes, self.elements, step, self.field)
        vmin, vmax = _field_range(self.elements, self.field)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("viridis")

        if self.view == "2D":
            artist = plot_tension_2d(
                self.main_ax,
                snapshot,
                self.config.width_m / 2.0,
                self.config.deformation_scale,
                norm,
                cmap,
            )
        else:
            artist = plot_tension_3d(
                self.main_ax,
                snapshot,
                self.config.width_m / 2.0,
                self.config.thickness_m / 2.0,
                self.config.deformation_scale,
                norm,
                cmap,
            )
        self.colorbar = self.fig.colorbar(artist, ax=self.main_ax, shrink=0.76, pad=0.07)
        self.colorbar.set_label(snapshot["label"])
        self.fig.canvas.draw_idle()

    def show(self) -> None:
        plt.show()


def run_headless_smoke(
    output_dir: Path,
    preview_dir: Path,
    solver: Path | None = None,
    auto_build: bool = True,
) -> dict[str, float | int]:
    """CI smoke test: run C FEM, load outputs, and render both display dimensions."""
    config = TensionRunConfig(elements=12, cycles=1, steps_per_cycle=24)
    run_fem_solver(config, output_dir, solver=solver, auto_build=auto_build)
    nodes, elements = load_fem_history(output_dir)
    if len(np.unique(elements["element"])) != config.elements:
        raise RuntimeError("headless smoke: element count mismatch")
    summary = save_preview_images(
        output_dir,
        preview_dir,
        half_width_m=config.width_m / 2.0,
        half_thickness_m=config.thickness_m / 2.0,
        deformation_scale=config.deformation_scale,
        field="stress",
    )
    summary["nodes_rows"] = int(nodes.size)
    summary["elements_rows"] = int(elements.size)
    return summary


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Integrated GUI for strictly 1D tensile FEM")
    parser.add_argument("project", nargs="?", type=Path, help="optional .ftgsim project to open")
    parser.add_argument("--output-dir", type=Path, default=Path("results/data/fem1d_ui_run"))
    parser.add_argument("--preview-dir", type=Path, default=Path("results/figures/fem1d_ui_run"))
    parser.add_argument("--solver", type=Path, default=None)
    parser.add_argument("--no-auto-build", action="store_true")
    parser.add_argument("--headless-smoke", action="store_true")
    parser.add_argument("--save-project", type=Path, default=None,
                        help="write the headless result as a .ftgsim bundle")
    args = parser.parse_args(argv)

    if args.headless_smoke:
        summary = run_headless_smoke(
            args.output_dir,
            args.preview_dir,
            solver=args.solver,
            auto_build=not args.no_auto_build,
        )
        for key, value in summary.items():
            print(f"{key}={value}")
        if args.save_project is not None:
            saved = save_tension_ftgsim(args.save_project,
                TensionRunConfig(elements=12, cycles=1, steps_per_cycle=24), args.output_dir)
            print(f"project={saved}")
        return

    app = FEMTensionApp(
        output_dir=args.output_dir,
        solver=args.solver,
        auto_build=not args.no_auto_build,
        project_path=args.project,
    )
    app.show()


if __name__ == "__main__":
    main()
