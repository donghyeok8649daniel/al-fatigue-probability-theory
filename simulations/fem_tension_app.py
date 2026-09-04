# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다.
# - 주요 클래스: TensionRunConfig, FEMTensionApp
# - 주요 함수/메서드: TensionRunConfig.length_m, TensionRunConfig.width_m, TensionRunConfig.thickness_m
#   TensionRunConfig.area_m2, TensionRunConfig.young_pa, TensionRunConfig.cubic_constants
#   TensionRunConfig.elastic_calibration_mode, validate_run_config, config_from_ftgsim
#   save_tension_ftgsim, initiation_snapshot, repository_root, solver_executable, _solver_sources
#   _needs_rebuild, build_fem_solver, solver_command, run_fem_solver, FEMTensionApp.__init__
#   FEMTensionApp._apply_config_to_boxes, FEMTensionApp._load_project
#   FEMTensionApp._create_parameter_panel, FEMTensionApp._create_result_controls
#   FEMTensionApp._read_config, FEMTensionApp._set_status, FEMTensionApp._on_run, FEMTensionApp._on_save
#   FEMTensionApp._on_save_project, FEMTensionApp._open_geometry, FEMTensionApp._on_open_geometry
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
import csv
from dataclasses import asdict, dataclass
from pathlib import Path
import os
import shutil
import subprocess
import sys
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
from simulations.ftgsim_format import create_ftgsim, extract_geometry, extract_results, open_ftgsim
from simulations.fvm1d_solver import run as run_fvm_backend
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.solver import LoadParams, SolverParams, run_ensemble
from simulations.mesh_viewer import MeshViewport, SUPPORTED_EXTENSIONS, load_mesh
from simulations.visualize_fem1d import load_numeric_csv
from theory.cubic_normal_orientation import (
    CubicElasticConstants,
    directional_young_modulus,
    miller_unit_vector,
)


@dataclass(frozen=True)
class TensionRunConfig:
    """User-facing inputs for one quasistatic cyclic axial-tension solve."""

    length_mm: float = 50.0
    width_mm: float = 10.0
    thickness_mm: float = 1.0
    young_gpa: float = 69.0
    loading_h: int = 1
    loading_k: int = 0
    loading_l: int = 0
    cubic_c11_gpa: float | None = None
    cubic_c12_gpa: float | None = None
    cubic_c44_gpa: float | None = None
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
        constants = self.cubic_constants
        if constants is not None:
            return directional_young_modulus(
                constants, self.loading_h, self.loading_k, self.loading_l)
        return self.young_gpa * 1.0e9

    @property
    def cubic_constants(self) -> CubicElasticConstants | None:
        values = (self.cubic_c11_gpa, self.cubic_c12_gpa, self.cubic_c44_gpa)
        if all(value is None for value in values):
            return None
        if any(value is None for value in values):
            raise ValueError("provide all of C11, C12 and C44, or none")
        return CubicElasticConstants(*(float(value)*1e9 for value in values))

    @property
    def elastic_calibration_mode(self) -> str:
        return "cubic_direction_projection" if self.cubic_constants is not None else "user_supplied_axis_modulus"


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
    if not all(isinstance(value, int) for value in (config.loading_h, config.loading_k, config.loading_l)):
        raise ValueError("Miller direction components must be integers")
    miller_unit_vector(config.loading_h, config.loading_k, config.loading_l)
    if config.cubic_constants is not None:
        config.cubic_constants.validate()
        if not np.isfinite(config.young_pa) or config.young_pa <= 0:
            raise ValueError("directional Young modulus must be positive")
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
                        *, view: str = "2D", field: str = "stress",
                        geometry_source: Path | None = None) -> Path:
    """Save setup and any existing FEM CSV results as an open `.ftgsim` bundle."""
    validate_run_config(config)
    files: dict[str, Path] = {}
    if output_dir is not None:
        source_dir = Path(output_dir)
        for name in ("nodes.csv", "elements.csv", "metadata.csv", "summary.json",
                     "probability_elements.csv", "initiation_elements.csv"):
            source = source_dir / name
            if source.is_file():
                files[f"results/{name}"] = source
    geometry_member = None
    source_dimension = None
    if geometry_source is not None:
        source = Path(geometry_source)
        if source.suffix.lower() not in SUPPORTED_EXTENSIONS or not source.is_file():
            raise ValueError("geometry_source must be an existing OBJ/STL/PLY/VTK file")
        geometry_member = f"geometry/source{source.suffix.lower()}"
        files[geometry_member] = source
        source_dimension = load_mesh(source).dimension
    has_initiation_results = "results/initiation_elements.csv" in files
    return create_ftgsim(
        path,
        setup={
            "physics_model": "1d_normal_tensile",
            "material_scope": "pure_single_crystal_aluminum",
            "tension_run": asdict(config),
            "probability_model": {
                "enabled": has_initiation_results,
                "calibration_status": (
                    "parameters_not_embedded_or_calibrated" if has_initiation_results else "not_solved"
                ),
                "coordinate": "local_homogeneous_spacing",
                "energy": "exact_riemann_zeta_bulk_lattice",
                "initiation_definition": "first_tangent_stiffness_loss",
                "critical_stretch_rule": "((m+1)/(n+1))**(1/(m-n))",
                "note": "Parameters and escape coupling must be supplied before activation.",
            },
            "result_references": sorted(name for name in files if name.startswith("results/")),
        },
        geometry={
            "mesh_dimension": 1,
            "loading_axis": miller_unit_vector(
                config.loading_h, config.loading_k, config.loading_l).tolist(),
            "crystal_loading_direction_hkl": [
                config.loading_h, config.loading_k, config.loading_l],
            "elastic_calibration_mode": config.elastic_calibration_mode,
            "directional_young_modulus_pa": config.young_pa,
            "geometry_kind": "uniform_bar",
            "length_mm": config.length_mm,
            "width_mm": config.width_mm,
            "thickness_mm": config.thickness_mm,
            "elements": config.elements,
            "source_member": geometry_member,
            "source_dimension": source_dimension,
        },
        display={"view": view, "field": field, "deformation_scale": config.deformation_scale},
        files=files,
        generator={"application": "fem_tension_app", "format_extension": ".ftgsim"},
    )


def initiation_snapshot(
    nodes: np.ndarray,
    elements: np.ndarray,
    initiation_elements: np.ndarray,
    step: int,
    field: str,
) -> dict:
    """Map an optional first-passage scalar channel onto the axial FEM elements."""
    columns = {
        "initiation": ("initiation_probability", "cumulative initiation probability"),
        "survival": ("survival", "intact survival probability"),
        "hazard": ("hazard_per_s", "initiation hazard [1/s]"),
    }
    if field not in columns:
        raise ValueError("field must be initiation, survival or hazard")
    column, label = columns[field]
    if column not in (initiation_elements.dtype.names or ()):
        raise ValueError(f"initiation result is missing column: {column}")
    base = axial_snapshot(nodes, elements, step, "stress")
    result_rows = initiation_elements[initiation_elements["step"] == step]
    element_rows = elements[elements["step"] == step]
    result_rows = result_rows[np.argsort(result_rows["element"])]
    element_rows = element_rows[np.argsort(element_rows["element"])]
    if result_rows.size != element_rows.size or not np.array_equal(
        result_rows["element"].astype(int), element_rows["element"].astype(int)
    ):
        raise ValueError("initiation/FEM element identifiers do not align")
    base["scalar"] = np.asarray(result_rows[column], dtype=float)
    base["field"] = field
    base["label"] = label
    return base


def repository_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
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


def run_selected_solver(
    config: TensionRunConfig,
    output_dir: Path,
    backend: str,
    solver: Path | None = None,
    auto_build: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the selected backend while preserving the shared CSV output contract."""
    if backend == "Theory":
        model_p = ModelParams()
        load_p = LoadParams(force_max=2.5 + 0.009 * config.stress_amplitude_mpa, cycles=config.cycles)
        solver_p = SolverParams(
            dt=load_p.period / max(config.steps_per_cycle, 2),
            n_trajectories=32,
            first_passage_stride=5,
            record_stride=1,
        )
        out = run_ensemble(model_p, load_p, solver_p)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        with (output_dir / "nodes.csv").open("w", newline="", encoding="utf-8") as h:
            w = csv.writer(h); w.writerow(["time_s", "step", "node", "x_m", "displacement_m", "applied_stress_pa"])
            for step, (t, force, strain) in enumerate(zip(out["time"], out["force"], out["strain"])):
                w.writerow([t, step, 0, 0.0, 0.0, force])
                w.writerow([t, step, 1, 1.0, strain, force])
        survival = np.asarray(out["survival"], dtype=float)
        time = np.asarray(out["time"], dtype=float)
        probability = 1.0 - survival
        hazard = np.zeros_like(survival)
        if len(survival) > 1:
            safe = np.maximum(survival, 1.0e-12)
            hazard[1:] = np.maximum(0.0, -np.diff(np.log(safe)) / np.maximum(np.diff(time), 1.0e-12))
        with (output_dir / "elements.csv").open("w", newline="", encoding="utf-8") as h:
            w = csv.writer(h); w.writerow(["time_s", "step", "element", "x_mid_m", "strain", "stress_pa", "applied_stress_pa"])
            for step, (t, force, strain) in enumerate(zip(time, out["force"], out["strain"])):
                w.writerow([t, step, 0, 0.5, strain, force, force])
        with (output_dir / "initiation_elements.csv").open("w", newline="", encoding="utf-8") as h:
            w = csv.writer(h); w.writerow(["time_s", "step", "element", "initiation_probability", "survival", "hazard_per_s"])
            for step, (t, p, s, hz) in enumerate(zip(time, probability, survival, hazard)):
                w.writerow([t, step, 0, p, s, hz])
        with (output_dir / "metadata.csv").open("w", newline="", encoding="utf-8") as h:
            csv.writer(h).writerows([["solver", "theory_core_v1_probability"], ["cycles", config.cycles]])
        return subprocess.CompletedProcess(["theory_core_v1"], 0, "Theory Core v1 complete\n", "")
    if backend == "FEM":
        return run_fem_solver(config, output_dir, solver=solver, auto_build=auto_build)
    if backend != "FVM":
        raise ValueError(f"unknown solver backend: {backend}")
    run_fvm_backend(
        elements=config.elements,
        length_m=config.length_m,
        area_m2=config.area_m2,
        young_pa=config.young_pa,
        stress_mean_mpa=config.stress_mean_mpa,
        stress_amplitude_mpa=config.stress_amplitude_mpa,
        frequency_hz=config.frequency_hz,
        cycles=config.cycles,
        steps_per_cycle=config.steps_per_cycle,
        outdir=output_dir,
    )
    return subprocess.CompletedProcess(["fvm1d_solver"], 0, "FVM complete\n", "")


class FEMTensionApp:
    """Matplotlib desktop app that runs the C solver and visualizes its axial output."""

    _INPUT_SPECS = (
        ("length_mm", "L [mm]", "50"),
        ("width_mm", "W [mm]", "10"),
        ("thickness_mm", "T [mm]", "1"),
        ("young_gpa", "E [GPa]", "69"),
        ("loading_direction", "Axis [h k l]", "1 0 0"),
        ("elements", "Cells", "40"),
        ("stress_mean_mpa", "Mean stress MPa", "50"),
        ("stress_amplitude_mpa", "Stress amp MPa", "100"),
        ("frequency_hz", "Freq Hz", "20"),
        ("cycles", "Cycles", "2"),
        ("steps_per_cycle", "Steps/cycle", "80"),
        ("deformation_scale", "Deform scale", "1"),
    )

    def __init__(
        self,
        output_dir: Path,
        solver: Path | None = None,
        auto_build: bool = True,
        project_path: Path | None = None,
        geometry_path: Path | None = None,
        backend: str = "FVM",
    ) -> None:
        self.output_dir = Path(output_dir)
        self.solver = None if solver is None else Path(solver)
        self.auto_build = auto_build
        self.backend = backend if backend in {"FEM", "FVM", "Theory"} else "Theory"
        self.nodes: np.ndarray | None = None
        self.elements: np.ndarray | None = None
        self.initiation_elements: np.ndarray | None = None
        self.steps = np.array([0], dtype=int)
        self.config = TensionRunConfig()
        self.view = "2D"
        self.field = "stress"
        self.main_ax = None
        self.colorbar = None
        self.mesh_viewports: list[MeshViewport] = []
        self.geometry_source_path: Path | None = None

        self.fig = plt.figure(figsize=(13.4, 7.5))
        self._place_window_on_screen()
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
        if geometry_path is not None:
            self._open_geometry(Path(geometry_path))
        self.redraw()

    def _place_window_on_screen(self) -> None:
        """Center the GUI on the current monitor and keep it inside the work area."""
        window = getattr(self.fig.canvas.manager, "window", None)
        if window is None:
            return
        try:
            window.update_idletasks()
            screen_w = int(window.winfo_screenwidth())
            screen_h = int(window.winfo_screenheight())
            width = min(1400, max(960, int(screen_w * 0.88)))
            height = min(900, max(620, int(screen_h * 0.82)))
            x = max(0, (screen_w - width) // 2)
            y = max(0, (screen_h - height) // 2)
            window.geometry(f"{width}x{height}+{x}+{y}")
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return

    def _apply_config_to_boxes(self) -> None:
        for key, _label, _initial in self._INPUT_SPECS:
            if key == "loading_direction":
                value = f"{self.config.loading_h} {self.config.loading_k} {self.config.loading_l}"
            else:
                value = str(getattr(self.config, key))
            self.textboxes[key].set_val(value)

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
        initiation_path = self.output_dir / "initiation_elements.csv"
        if initiation_path.is_file():
            self.initiation_elements = load_numeric_csv(initiation_path)
        geometry_files = extract_geometry(bundle, self.output_dir / "imported_geometry")
        if geometry_files:
            self._open_geometry(geometry_files[0])
        self.status.set_text(f"Opened {path.name} (1D normal tension only)")

    def _create_parameter_panel(self) -> None:
        self.fig.text(0.025, 0.955, "Tensile test inputs", ha="left", va="top", weight="bold")
        top = 0.89
        spacing = 0.057
        for i, (key, label, initial) in enumerate(self._INPUT_SPECS):
            ax = self.fig.add_axes([0.025, top - i * spacing, 0.175, 0.035])
            box = TextBox(ax, label, initial=initial, label_pad=0.03)
            self.textboxes[key] = box

        backend_ax = self.fig.add_axes([0.025, 0.235, 0.175, 0.055])
        options = ("Theory", "FVM", "FEM")
        self.backend_radio = RadioButtons(backend_ax, options, active=options.index(self.backend))
        self.backend_radio.on_clicked(self._on_backend)

        run_ax = self.fig.add_axes([0.025, 0.175, 0.082, 0.048])
        self.run_button = Button(run_ax, "Run solver")
        self.run_button.on_clicked(self._on_run)

        save_ax = self.fig.add_axes([0.118, 0.175, 0.082, 0.048])
        self.save_button = Button(save_ax, "Save views")
        self.save_button.on_clicked(self._on_save)

        project_ax = self.fig.add_axes([0.025, 0.115, 0.175, 0.040])
        self.project_button = Button(project_ax, "Save .ftgsim")
        self.project_button.on_clicked(self._on_save_project)

        geometry_ax = self.fig.add_axes([0.025, 0.065, 0.175, 0.040])
        self.geometry_button = Button(geometry_ax, "Open 1D/2D/3D mesh")
        self.geometry_button.on_clicked(self._on_open_geometry)

        self.fig.text(
            0.025,
            0.045,
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
        field_ax = self.fig.add_axes([0.875, 0.02, 0.115, 0.18])
        self.view_radio = RadioButtons(view_ax, ("2D", "3D"), active=0)
        self.field_radio = RadioButtons(
            field_ax, ("stress", "strain", "initiation", "survival", "hazard"), active=0)
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

        direction = self.textboxes["loading_direction"].text.replace(",", " ").split()
        if len(direction) != 3:
            raise ValueError("crystal axis must contain three Miller integers: h k l")
        try:
            h, k, l = (int(value) for value in direction)
        except ValueError:
            raise ValueError("crystal axis must contain integer Miller indices") from None

        config = TensionRunConfig(
            length_mm=floating("length_mm"),
            width_mm=floating("width_mm"),
            thickness_mm=floating("thickness_mm"),
            young_gpa=floating("young_gpa"),
            loading_h=h,
            loading_k=k,
            loading_l=l,
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

    def _on_backend(self, label: str) -> None:
        self.backend = label
        if self.backend != "Theory" and self.field not in {"stress", "strain"}:
            self.field = "stress"
            self.field_radio.set_active(0)
        self._set_status(f"Selected {label} backend")

    def _on_run(self, _event) -> None:
        try:
            self.config = self._read_config()
            self._set_status(f"Running {self.backend}...")
            self.initiation_elements = None
            completed = run_selected_solver(
                self.config,
                self.output_dir,
                self.backend,
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
                + (completed.stdout.strip().splitlines()[0] if completed.stdout.strip() else f"{self.backend} complete")
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
                                        view=self.view, field=self.field,
                                        geometry_source=self.geometry_source_path)
            self._set_status(f"Saved project: {saved}")
        except Exception as exc:
            self._set_status(f"ERROR saving project: {exc}")

    def _open_geometry(self, path: Path) -> None:
        mesh = load_mesh(path)
        self.geometry_source_path = Path(path)
        loading_axis = miller_unit_vector(
            self.config.loading_h, self.config.loading_k, self.config.loading_l)
        viewport = MeshViewport(mesh, loading_axis=loading_axis)
        self.mesh_viewports.append(viewport)
        viewport.figure.show()
        self._set_status(
            f"Opened {path.name}: {mesh.vertices.shape[0]} nodes, "
            f"{len(mesh.faces)} faces, inferred {mesh.dimension}D, "
            f"crystal axis [{self.config.loading_h} {self.config.loading_k} {self.config.loading_l}]"
        )

    def _on_open_geometry(self, _event) -> None:
        try:
            from tkinter import Tk, filedialog
            root = Tk(); root.withdraw()
            selected = filedialog.askopenfilename(
                title="Open mesh geometry",
                filetypes=[("Supported mesh", "*.obj *.stl *.ply *.vtk"),
                           ("OBJ", "*.obj"), ("STL", "*.stl"),
                           ("PLY", "*.ply"), ("Legacy VTK", "*.vtk")],
            )
            root.destroy()
            if selected:
                self._open_geometry(Path(selected))
        except Exception as exc:
            self._set_status(f"ERROR opening geometry: {exc}")

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
        if self.field in {"stress", "strain"}:
            snapshot = axial_snapshot(self.nodes, self.elements, step, self.field)
            vmin, vmax = _field_range(self.elements, self.field)
        else:
            if self.initiation_elements is None:
                self.main_ax.text(0.5, 0.5,
                    "No initiation result in this project.\n"
                    "Calibrated/declared probability parameters are required before solving.",
                    transform=self.main_ax.transAxes, ha="center", va="center")
                self.main_ax.set_axis_off(); self.fig.canvas.draw_idle(); return
            snapshot = initiation_snapshot(
                self.nodes, self.elements, self.initiation_elements, step, self.field)
            values = np.asarray(self.initiation_elements[
                {"initiation": "initiation_probability", "survival": "survival",
                 "hazard": "hazard_per_s"}[self.field]], dtype=float)
            vmin, vmax = float(np.min(values)), float(np.max(values))
            if np.isclose(vmin, vmax):
                margin = max(abs(vmin), 1.0) * 1e-9; vmin -= margin; vmax += margin
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
    parser = argparse.ArgumentParser(description="Integrated GUI for 1D tensile FEM/FVM")
    parser.add_argument("input", nargs="?", type=Path,
                        help="optional .ftgsim project or OBJ/STL/PLY/VTK geometry")
    parser.add_argument("--geometry", type=Path, default=None,
                        help="mesh geometry to open in the CAD-style viewport")
    parser.add_argument("--output-dir", type=Path, default=Path("results/data/fem1d_ui_run"))
    parser.add_argument("--preview-dir", type=Path, default=Path("results/figures/fem1d_ui_run"))
    parser.add_argument("--solver", type=Path, default=None)
    parser.add_argument("--backend", choices=("Theory", "FVM", "FEM"), default="Theory")
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

    project_path = args.input if args.input and args.input.suffix.lower() == ".ftgsim" else None
    geometry_path = args.geometry
    if args.input and args.input.suffix.lower() in SUPPORTED_EXTENSIONS:
        if geometry_path is not None:
            parser.error("specify geometry either positionally or with --geometry, not both")
        geometry_path = args.input
    elif args.input and project_path is None:
        parser.error("input must be .ftgsim, .obj, .stl, .ply or .vtk")

    app = FEMTensionApp(
        output_dir=args.output_dir,
        solver=args.solver,
        auto_build=not args.no_auto_build,
        project_path=project_path,
        geometry_path=geometry_path,
        backend=args.backend,
    )
    app.show()


if __name__ == "__main__":
    main()
