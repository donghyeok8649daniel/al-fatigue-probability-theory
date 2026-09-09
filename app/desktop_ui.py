"""Lightweight Pre/Solve/Post desktop UI for the probability-PDE solver."""
from __future__ import annotations

import os
import queue
import sys
import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import matplotlib

matplotlib.use("TkAgg")
if os.name == "nt":
    matplotlib.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter
import numpy as np

from .i18n import (
    DEFAULT_LANGUAGE,
    FIELD_TEXT_KEYS,
    LANGUAGE_NAMES,
    Localizer,
    plot_strings,
    tr,
)
from .solver_adapter import (
    UIAnalysisConfig,
    load_interpretation,
    physical_load_conversion,
    run_ui_analysis,
)
from .convergence_check import run_convergence_check
from .specimen_probability import (
    BELOW_RESOLUTION,
    aggregate_specimen_probability,
)
from solver_v1.physical_time import load_time_calibration
from solver_v1.kinetic_calibration_workflow import validate_for_energy_model
from solver_v1.energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    TWO_ROW_LJ_REFERENCE,
    energy_model_metadata,
)


APP_BG = "#eef1f4"
PANEL_BG = "#ffffff"
HEADER_BG = "#253746"
SIDEBAR_BG = "#e4e9ed"
ACCENT = "#1677c8"
TEXT = "#202830"
MUTED = "#65717c"
WINDOW_TITLE = tr("app.title", DEFAULT_LANGUAGE)
_INSTANCE_MUTEX = None


def acquire_single_instance() -> bool:
    """Prevent accidental duplicate desktop windows on Windows."""

    global _INSTANCE_MUTEX
    if os.name != "nt":
        return True
    import ctypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    user32 = ctypes.WinDLL("user32", use_last_error=True)
    kernel32.CreateMutexW.argtypes = (ctypes.c_void_p, ctypes.c_bool, ctypes.c_wchar_p)
    kernel32.CreateMutexW.restype = ctypes.c_void_p
    kernel32.CloseHandle.argtypes = (ctypes.c_void_p,)
    user32.FindWindowW.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p)
    user32.FindWindowW.restype = ctypes.c_void_p
    handle = kernel32.CreateMutexW(None, False, "Local\\AlFatigueProbabilityPDE")
    if not handle:
        return True
    if ctypes.get_last_error() == 183:
        window = user32.FindWindowW(None, WINDOW_TITLE)
        if window:
            user32.ShowWindow(window, 9)
            user32.SetForegroundWindow(window)
        kernel32.CloseHandle(handle)
        return False
    _INSTANCE_MUTEX = handle
    return True


class DesktopApp:
    """Responsive engineering workspace with explicit Pre/Solve/Post stages."""

    PARAMS = (
        ("young_gpa", "field.young_gpa", "69", "unit.gpa"),
        ("stress_mean_mpa", "field.stress_mean", "50", "unit.mpa"),
        ("stress_amplitude_mpa", "field.stress_amplitude", "100", "unit.mpa"),
        ("model_frequency", "field.model_frequency", "25", "unit.cycles_model_time"),
        ("cycles", "field.cycles", "1", "unit.cycles"),
        ("steps_per_cycle", "field.output_resolution", "40", "unit.steps_cycle"),
        ("tensile_direction", "field.axial_direction", "1 0 0", "unit.direction"),
    )
    LOAD_PRESETS = {
        "custom": ("preset.custom", None, None),
        "small_signal": ("preset.small_signal", 0.0, 50.0),
        "moderate": ("preset.moderate", 0.0, 150.0),
        "compression": ("preset.compression", -500.0, 400.0),
        "extreme": ("preset.extreme", 900.0, 2000.0),
    }

    def __init__(self, *, root=None) -> None:
        self.localizer = Localizer(DEFAULT_LANGUAGE)
        self.time_calibration = load_time_calibration(
            Path(__file__).resolve().parents[1]
            / "solver_v1" / "data" / "aluminum_kinetic_calibration.json"
        )
        # One interpreter per desktop process. Tests/embedded workspaces can
        # supply a Toplevel in that interpreter instead of repeatedly loading
        # and destroying independent native Tcl/Tk runtimes.
        self.root = root if root is not None else tk.Tk()
        self.root.title(self._tr("app.title"))
        self.root.configure(bg=APP_BG)
        self.root.minsize(940, 620)
        self._center(1180, 760)

        self.entries: dict[str, ttk.Entry] = {}
        self.language_display = tk.StringVar(value=LANGUAGE_NAMES[DEFAULT_LANGUAGE])
        self.spatial_backend = tk.StringVar(value="FVM")
        self.analysis_quality = tk.StringVar(value=self._tr("quality.preview_option"))
        self.analysis_quality_code = "preview"
        self.probability_scale = tk.StringVar(value=self._tr("option.local"))
        self.probability_scale_code = "local"
        self.load_preset = tk.StringVar(value=self._tr("preset.custom"))
        self.load_preset_code = "custom"
        self.time_basis_code = "model"
        self.time_basis_display = tk.StringVar(value=self._tr("option.model_time"))
        self.energy_model_code = TWO_ROW_LJ_REFERENCE
        self.energy_model_display = tk.StringVar(
            value=self._tr("energy_model.lj_reference")
        )
        self.time_warning_display = tk.StringVar(value="")
        self._frequency_values = {"model": "25", "physical": "1"}
        self.stress_context = tk.StringVar(value="")
        self.specimen_N_eff = tk.StringVar(value="N_eff = —")
        self.local_floor_display = tk.StringVar(value="—")
        self.plastic_floor_display = tk.StringVar(value="—")
        self.configurational_barrier_display = tk.StringVar(value="—")
        self.opening_barrier_display = tk.StringVar(value="—")
        self.plasticity_detail_display = tk.StringVar(value="—")
        self.probability_status_display = tk.StringVar(value="—")
        self.local_probability_display = tk.StringVar(value="—")
        self.specimen_extrapolation_display = tk.StringVar(value="—")
        self.specimen_certified_display = tk.StringVar(value="—")
        self._specimen_status_key = "status.local_result_required"
        self.field = tk.StringVar(value="strain_components")
        self.status = tk.StringVar(value=self._tr("status.ready"))
        self.busy = False
        self.stop_event = threading.Event()
        self.result: dict[str, object] | None = None
        self.live_records: list[dict[str, float]] = []
        self._queue: queue.Queue[tuple[str, object]] = queue.Queue(maxsize=512)
        self._poll_job = None
        self._last_draw = 0.0
        self._last_field: str | None = None
        self._view_limits: dict[str, tuple[tuple[float, float], tuple[float, float]]] = {}
        self._pan_origin = None
        self._text_bindings: list[tuple[object, str, str]] = []
        self._tree_text: list[tuple[str, str]] = []
        self._status_key = "status.ready"
        self._status_values: dict[str, object] = {}
        self._summary_kind = "ready"
        self._summary_payload: dict[str, object] = {}
        self._button_key = "button.run"
        self.last_config: UIAnalysisConfig | None = None

        self._styles()
        self._header()
        self._workspace()
        self._statusbar()
        self._connect_plot_events()
        self.root.protocol("WM_DELETE_WINDOW", self._close)
        self._update_stress_context()
        self._refresh_time_basis_ui()
        self._refresh_specimen_labels()
        self._plot()

    def _tr(self, key: str, **values: object) -> str:
        return self.localizer.text(key, **values)

    def _bind_text(self, widget, key: str, option: str = "text"):
        widget.configure(**{option: self._tr(key)})
        self._text_bindings.append((widget, option, key))
        return widget

    def _center(self, width: int, height: int) -> None:
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        width = min(width, max(900, sw - 60))
        height = min(height, max(600, sh - 90))
        x, y = max(0, (sw - width) // 2), max(0, (sh - height) // 2)
        self.root.geometry(f"{width}x{height}+{x}+{y}")

    def _styles(self) -> None:
        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("App.TFrame", background=APP_BG)
        style.configure("Panel.TFrame", background=PANEL_BG)
        style.configure("Side.TFrame", background=SIDEBAR_BG)
        style.configure(
            "Header.TLabel",
            background=HEADER_BG,
            foreground="white",
            font=("Segoe UI", 13, "bold"),
        )
        style.configure(
            "SubHeader.TLabel",
            background=HEADER_BG,
            foreground="#cbd5dd",
            font=("Segoe UI", 9),
        )
        style.configure(
            "Section.TLabel", background=PANEL_BG, foreground=TEXT,
            font=("Segoe UI", 10, "bold")
        )
        style.configure("Property.TLabel", background=PANEL_BG, foreground=TEXT)
        style.configure("Unit.TLabel", background=PANEL_BG, foreground=MUTED)
        style.configure(
            "Accent.TButton", background=ACCENT, foreground="white",
            font=("Segoe UI", 9, "bold"), padding=(14, 7)
        )
        style.configure(
            "Field.Toolbutton", background="#e7edf2", foreground=TEXT,
            padding=(8, 5), relief="flat"
        )
        style.map(
            "Field.Toolbutton",
            background=[("selected", ACCENT), ("active", "#cfe4f5")],
            foreground=[("selected", "white")],
        )
        style.configure("TNotebook.Tab", font=("Segoe UI", 10, "bold"), padding=(18, 7))

    def _header(self) -> None:
        header = tk.Frame(self.root, bg=HEADER_BG, height=56)
        header.pack(fill="x")
        header.pack_propagate(False)
        brand = self._bind_text(
            ttk.Label(header, style="Header.TLabel"), "app.brand"
        )
        brand.pack(
            side="left", padx=(18, 10), pady=15
        )
        subtitle = self._bind_text(
            ttk.Label(header, style="SubHeader.TLabel"), "app.subtitle"
        )
        subtitle.pack(side="left", pady=18)
        language_frame = tk.Frame(header, bg=HEADER_BG)
        language_frame.pack(side="right", padx=16, pady=12)
        language_label = self._bind_text(
            ttk.Label(language_frame, style="SubHeader.TLabel"), "language.label"
        )
        language_label.pack(side="left", padx=(0, 7))
        self.language_selector = ttk.Combobox(
            language_frame,
            textvariable=self.language_display,
            values=tuple(LANGUAGE_NAMES.values()),
            state="readonly",
            width=10,
        )
        self.language_selector.pack(side="left")
        self.language_selector.bind("<<ComboboxSelected>>", self._on_language_selected)

    def _workspace(self) -> None:
        body = ttk.Frame(self.root, style="App.TFrame")
        body.pack(fill="both", expand=True)
        sidebar = ttk.Frame(body, style="Side.TFrame", width=205)
        sidebar.pack(side="left", fill="y")
        sidebar.pack_propagate(False)
        tree_header = tk.Label(
            sidebar, bg=SIDEBAR_BG, fg=MUTED,
            font=("Segoe UI", 8, "bold")
        )
        self._bind_text(tree_header, "nav.study_tree").pack(
            anchor="w", padx=14, pady=(14, 5)
        )
        self.tree = ttk.Treeview(sidebar, show="tree", height=20)

        def tree_item(parent: str, key: str, *, open: bool = False) -> str:
            item = self.tree.insert(
                parent, "end", text=f"  {self._tr(key)}", open=open
            )
            self._tree_text.append((item, key))
            return item

        project = tree_item("", "tree.study", open=True)
        pre = tree_item(project, "tree.pre", open=True)
        tree_item(pre, "tree.material_load")
        tree_item(pre, "tree.axial_direction")
        solve = tree_item(project, "tree.solve", open=True)
        tree_item(solve, "tree.probability_pde")
        tree_item(solve, "tree.spatial")
        post = tree_item(project, "tree.post", open=True)
        for key in FIELD_TEXT_KEYS.values():
            tree_item(post, key)
        self.tree.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        main = ttk.Frame(body, style="App.TFrame")
        main.pack(side="left", fill="both", expand=True)
        self.notebook = ttk.Notebook(main)
        self.notebook.pack(fill="both", expand=True, padx=8, pady=8)
        self.pre_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.solve_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.post_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.notebook.add(self.pre_tab, text=self._tr("tab.pre"))
        self.notebook.add(self.solve_tab, text=self._tr("tab.solve"))
        self.notebook.add(self.post_tab, text=self._tr("tab.post"))
        self._pre_tab()
        self._solve_tab()
        self._post_tab()

    def _pre_tab(self) -> None:
        section = self._bind_text(
            ttk.Label(self.pre_tab, style="Section.TLabel"), "section.axial_inputs"
        )
        section.grid(
            row=0, column=0, columnspan=3, sticky="w", padx=20, pady=(18, 10)
        )
        energy_label = self._bind_text(
            ttk.Label(self.pre_tab, style="Property.TLabel"), "field.energy_model"
        )
        energy_label.grid(row=1, column=0, sticky="w", padx=(20, 8), pady=7)
        self.energy_model_selector = ttk.Combobox(
            self.pre_tab,
            textvariable=self.energy_model_display,
            values=self._energy_model_values(),
            state="readonly",
            width=36,
        )
        self.energy_model_selector.grid(
            row=1, column=1, columnspan=2, sticky="ew", padx=(4, 20), pady=7
        )
        self.energy_model_selector.bind(
            "<<ComboboxSelected>>", self._on_energy_model_selected
        )
        time_label = self._bind_text(
            ttk.Label(self.pre_tab, style="Property.TLabel"), "field.time_basis"
        )
        time_label.grid(row=2, column=0, sticky="w", padx=(20, 8), pady=7)
        self.time_basis_selector = ttk.Combobox(
            self.pre_tab,
            textvariable=self.time_basis_display,
            values=self._time_basis_values(),
            state="readonly",
            width=20,
        )
        self.time_basis_selector.grid(row=2, column=1, sticky="ew", padx=4, pady=7)
        self.time_basis_selector.bind(
            "<<ComboboxSelected>>", self._on_time_basis_selected
        )
        self.load_kinetics_button = self._bind_text(
            ttk.Button(self.pre_tab, command=self._choose_kinetic_calibration),
            "button.load_kinetics",
        )
        self.load_kinetics_button.grid(row=3, column=0, sticky="w", padx=(20, 8), pady=7)
        self.time_warning_label = ttk.Label(
            self.pre_tab,
            textvariable=self.time_warning_display,
            style="Unit.TLabel",
            wraplength=360,
            justify="left",
        )
        self.time_warning_label.grid(row=2, column=2, rowspan=2, sticky="w", padx=(6, 20), pady=7)
        for row, (key, label_key, default, unit_key) in enumerate(self.PARAMS, start=4):
            label = self._bind_text(
                ttk.Label(self.pre_tab, style="Property.TLabel"), label_key
            )
            label.grid(
                row=row, column=0, sticky="w", padx=(20, 8), pady=7
            )
            entry = ttk.Entry(self.pre_tab, width=20)
            entry.insert(0, default)
            entry.grid(row=row, column=1, sticky="ew", padx=4, pady=7)
            unit = self._bind_text(
                ttk.Label(self.pre_tab, style="Unit.TLabel"), unit_key
            )
            unit.grid(
                row=row, column=2, sticky="w", padx=(6, 20), pady=7
            )
            self.entries[key] = entry
            if key == "tensile_direction":
                entry.configure(state="disabled")
            if key == "model_frequency":
                self.frequency_label = label
                self.frequency_unit = unit
        self.pre_tab.columnconfigure(1, weight=1)
        preset_frame = self._bind_text(
            ttk.LabelFrame(self.pre_tab, padding=10), "section.load_presets"
        )
        preset_frame.grid(
            row=len(self.PARAMS) + 4,
            column=0,
            columnspan=3,
            sticky="ew",
            padx=20,
            pady=(10, 4),
        )
        preset_label = self._bind_text(
            ttk.Label(preset_frame, style="Property.TLabel"), "field.load_preset"
        )
        preset_label.grid(row=0, column=0, sticky="w", padx=(0, 8))
        self.preset_selector = ttk.Combobox(
            preset_frame,
            textvariable=self.load_preset,
            values=self._preset_values(),
            state="readonly",
            width=31,
        )
        self.preset_selector.grid(row=0, column=1, sticky="ew")
        self.preset_selector.bind("<<ComboboxSelected>>", self._on_preset_selected)
        self.apply_preset_button = self._bind_text(
            ttk.Button(preset_frame, command=self._apply_load_preset),
            "button.apply_preset",
        )
        self.apply_preset_button.grid(row=0, column=2, padx=(8, 0))
        preset_frame.columnconfigure(1, weight=1)
        self.stress_context_label = ttk.Label(
            preset_frame,
            textvariable=self.stress_context,
            style="Unit.TLabel",
            justify="left",
        )
        self.stress_context_label.grid(
            row=1, column=0, columnspan=3, sticky="ew", pady=(8, 0)
        )
        for key in ("young_gpa", "stress_mean_mpa", "stress_amplitude_mpa"):
            self.entries[key].bind("<FocusOut>", self._update_stress_context)
            self.entries[key].bind("<Return>", self._update_stress_context)
        note = self._bind_text(
            ttk.LabelFrame(self.pre_tab, padding=12), "section.scope"
        )
        note.grid(row=len(self.PARAMS) + 5, column=0, columnspan=3, sticky="ew", padx=20, pady=8)
        scope = self._bind_text(ttk.Label(note, justify="left"), "scope.text")
        scope.pack(anchor="w")

    def _solve_tab(self) -> None:
        left = ttk.Frame(self.solve_tab, style="Panel.TFrame")
        left.pack(side="left", fill="y", padx=20, pady=18)
        self._bind_text(
            ttk.Label(left, style="Section.TLabel"), "section.probability_pde"
        ).pack(anchor="w")
        self._bind_text(ttk.Label(
            left, foreground=ACCENT,
            background=PANEL_BG, font=("Segoe UI", 9, "bold")
        ), "status.always_active").pack(anchor="w", pady=(3, 15))
        self._bind_text(
            ttk.Label(left, style="Section.TLabel"), "section.spatial_backend"
        ).pack(anchor="w")
        ttk.Combobox(
            left, textvariable=self.spatial_backend, values=("FVM", "FEM"),
            state="readonly", width=22
        ).pack(anchor="w", pady=(6, 15))
        self._bind_text(
            ttk.Label(left, style="Section.TLabel"), "section.grid_quality"
        ).pack(
            anchor="w"
        )
        self.quality_selector = ttk.Combobox(
            left,
            textvariable=self.analysis_quality,
            values=self._quality_values(),
            state="readonly",
            width=29,
        )
        self.quality_selector.pack(anchor="w", pady=(6, 15))
        self.quality_selector.bind("<<ComboboxSelected>>", self._on_quality_selected)
        explanation = self._bind_text(
            ttk.Label(left, style="Property.TLabel", justify="left"),
            "solve.explanation",
        )
        explanation.pack(anchor="w", pady=(0, 10))
        specimen = self._bind_text(
            ttk.LabelFrame(left, padding=8), "section.specimen_probability"
        )
        specimen.pack(fill="x", pady=(0, 10))
        scale_label = self._bind_text(
            ttk.Label(specimen), "field.probability_scale"
        )
        scale_label.grid(row=0, column=0, columnspan=2, sticky="w")
        self.probability_scale_selector = ttk.Combobox(
            specimen,
            textvariable=self.probability_scale,
            values=self._probability_scale_values(),
            state="readonly",
            width=20,
        )
        self.probability_scale_selector.grid(
            row=1, column=0, columnspan=3, sticky="ew", pady=(3, 6)
        )
        self.probability_scale_selector.bind(
            "<<ComboboxSelected>>", self._on_probability_scale_selected
        )
        for row, (key, text_key) in enumerate((
            ("correlation_area_mm2", "field.correlation_area"),
            ("stressed_area_mm2", "field.stressed_area"),
        ), start=2):
            label = self._bind_text(ttk.Label(specimen), text_key)
            label.grid(row=row, column=0, sticky="w", pady=2)
            entry = ttk.Entry(specimen, width=10)
            entry.grid(row=row, column=1, sticky="ew", padx=4, pady=2)
            entry.bind("<FocusOut>", self._update_specimen_probability)
            entry.bind("<Return>", self._update_specimen_probability)
            self.entries[key] = entry
            unit = self._bind_text(ttk.Label(specimen), "unit.mm2")
            unit.grid(row=row, column=2, sticky="w", pady=2)
        specimen.columnconfigure(1, weight=1)
        ttk.Label(specimen, textvariable=self.specimen_N_eff).grid(
            row=4, column=0, columnspan=3, sticky="w", pady=(5, 0)
        )
        ttk.Label(specimen, textvariable=self.local_floor_display).grid(
            row=5, column=0, columnspan=3, sticky="w"
        )
        ttk.Label(specimen, textvariable=self.plastic_floor_display).grid(
            row=6, column=0, columnspan=3, sticky="w"
        )
        ttk.Label(specimen, textvariable=self.local_probability_display).grid(
            row=7, column=0, columnspan=3, sticky="w"
        )
        ttk.Label(specimen, textvariable=self.specimen_extrapolation_display).grid(
            row=8, column=0, columnspan=3, sticky="w"
        )
        ttk.Label(specimen, textvariable=self.specimen_certified_display).grid(
            row=9, column=0, columnspan=3, sticky="w"
        )
        ttk.Label(specimen, textvariable=self.probability_status_display).grid(
            row=10, column=0, columnspan=3, sticky="w"
        )
        mechanism = self._bind_text(
            ttk.LabelFrame(left, padding=8), "section.mechanism_diagnostics"
        )
        mechanism.pack(fill="x", pady=(0, 10))
        ttk.Label(
            mechanism, textvariable=self.configurational_barrier_display
        ).pack(anchor="w")
        ttk.Label(
            mechanism, textvariable=self.opening_barrier_display
        ).pack(anchor="w")
        ttk.Label(
            mechanism,
            textvariable=self.plasticity_detail_display,
            justify="left",
            wraplength=235,
        ).pack(anchor="w", pady=(3, 0))
        self._bind_text(
            ttk.Label(mechanism, wraplength=235, justify="left"),
            "diagnostic.barrier_scope",
        ).pack(anchor="w", pady=(3, 0))
        self.run_button = ttk.Button(
            left, text=self._tr("button.run"), style="Accent.TButton",
            command=self._start_solve
        )
        self.run_button.pack(fill="x")
        self.convergence_button = ttk.Button(
            left,
            text=self._tr("button.run_convergence"),
            command=self._start_convergence_check,
            state="disabled",
        )
        self.convergence_button.pack(fill="x", pady=(6, 0))
        self._text_bindings.append(
            (self.convergence_button, "text", "button.run_convergence")
        )
        self.progress = ttk.Progressbar(left, mode="indeterminate", length=235)
        self.progress.pack(fill="x", pady=12)
        self.summary = tk.Text(
            self.solve_tab, wrap="word", relief="flat", bg="#f7f9fa",
            fg=TEXT, font=("Consolas", 9), padx=14, pady=12
        )
        self.summary.pack(side="left", fill="both", expand=True, padx=(0, 20), pady=18)
        self._render_summary()

    def _post_tab(self) -> None:
        toolbar = ttk.Frame(self.post_tab, style="Panel.TFrame")
        toolbar.pack(fill="x", padx=12, pady=(10, 3))
        self.field_buttons: dict[str, ttk.Radiobutton] = {}
        for index, (key, text_key) in enumerate(FIELD_TEXT_KEYS.items()):
            row, column = divmod(index, 5)
            button = ttk.Radiobutton(
                toolbar, text=self._tr(text_key), value=key, variable=self.field,
                style="Field.Toolbutton", command=self._plot
            )
            button.grid(row=row, column=column, sticky="ew", padx=2, pady=2)
            self._text_bindings.append((button, "text", text_key))
            self.field_buttons[key] = button
        for column in range(5):
            toolbar.columnconfigure(column, weight=1, uniform="field")
        self.figure = Figure(figsize=(8, 5), dpi=96, facecolor=PANEL_BG)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.post_tab)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=(4, 0))
        nav = ttk.Frame(self.post_tab, style="Panel.TFrame")
        nav.pack(fill="x", padx=10, pady=(0, 7))
        NavigationToolbar2Tk(self.canvas, nav, pack_toolbar=False).pack(side="left")
        nav_help = self._bind_text(
            ttk.Label(nav, foreground=MUTED), "nav.instructions"
        )
        nav_help.pack(side="right", padx=8)

    def _statusbar(self) -> None:
        tk.Label(
            self.root, textvariable=self.status, anchor="w", bg="#dbe2e7",
            fg=TEXT, font=("Segoe UI", 9), padx=10, pady=4
        ).pack(fill="x", side="bottom")

    def _quality_values(self) -> tuple[str, str]:
        return (
            self._tr("quality.preview_option"),
            self._tr("quality.resolved_option"),
        )

    def _time_basis_values(self) -> tuple[str, ...]:
        values = [self._tr("option.model_time")]
        if self._kinetics_available():
            values.append(self._tr("option.physical_time"))
        return tuple(values)

    def _kinetics_available(self) -> bool:
        try:
            validate_for_energy_model(self.time_calibration, self.energy_model_code)
            return True
        except ValueError:
            return False

    def load_kinetic_calibration_path(self, path) -> None:
        """Validate before modifying setup; keep all existing result arrays/views."""
        if self.busy:
            raise ValueError(self._tr("error.kinetics_busy"))
        calibration = load_time_calibration(path)
        validate_for_energy_model(calibration, self.energy_model_code)
        if self.time_basis_code == "physical":
            # Return to the same model frequency before replacing the clock.
            self.time_basis_display.set(self._tr("option.model_time"))
            self._on_time_basis_selected()
            if self.time_basis_code != "model":
                raise ValueError(self._tr("error.time_frequency"))
        self.time_calibration = calibration
        self._refresh_time_basis_ui()
        self._set_status("status.kinetics_loaded")

    def _choose_kinetic_calibration(self) -> None:
        if self.busy:
            self._set_status("error.kinetics_busy")
            return
        path = filedialog.askopenfilename(title=self._tr("dialog.kinetics_title"),
            filetypes=[(self._tr("filetype.kinetics_json"), "*.json")])
        if not path:
            return
        try:
            self.load_kinetic_calibration_path(path)
        except (ValueError, TypeError, OSError) as exc:
            messagebox.showerror(self._tr("dialog.kinetics_title"),
                self._tr("error.kinetics_invalid", detail=str(exc)))

    def _energy_model_values(self) -> tuple[str, ...]:
        return tuple(
            self._tr(energy_model_metadata(model_id).display_key)
            for model_id in (
                TWO_ROW_LJ_REFERENCE,
                ANALYTIC_LJ_EAM_HYPOTHETICAL,
                AL_TARGET_BEST_FEASIBLE,
            )
        )

    def _on_energy_model_selected(self, _event=None) -> None:
        if self.time_basis_code == "physical":
            self.time_basis_display.set(self._tr("option.model_time"))
            self._on_time_basis_selected()
            if self.time_basis_code != "model":
                self.energy_model_display.set(self._tr(
                    energy_model_metadata(self.energy_model_code).display_key))
                return
        selected = self.energy_model_display.get()
        for model_id in (
            TWO_ROW_LJ_REFERENCE,
            ANALYTIC_LJ_EAM_HYPOTHETICAL,
            AL_TARGET_BEST_FEASIBLE,
        ):
            if selected == self._tr(energy_model_metadata(model_id).display_key):
                self.energy_model_code = model_id
                break
        self.energy_model_display.set(
            self._tr(energy_model_metadata(self.energy_model_code).display_key)
        )
        self._refresh_time_basis_ui()

    def _on_time_basis_selected(self, _event=None) -> None:
        selected = self.time_basis_display.get()
        requested = (
            "physical" if selected == self._tr("option.physical_time") else "model"
        )
        if requested == "physical" and not self._kinetics_available():
            self.time_basis_code = "model"
            self.time_basis_display.set(self._tr("option.model_time"))
            self._set_status("status.physical_time_unavailable")
            self._refresh_time_basis_ui()
            return
        if hasattr(self, "entries") and "model_frequency" in self.entries:
            if requested != self.time_basis_code:
                try:
                    old = float(self.entries["model_frequency"].get())
                    if not np.isfinite(old) or old <= 0:
                        raise ValueError("frequency")
                    t0 = float(self.time_calibration.t0_seconds)
                    converted = old / t0 if requested == "physical" else old * t0
                    self._frequency_values[requested] = f"{converted:.16g}"
                except (ValueError, TypeError):
                    self.time_basis_display.set(self._tr("option.physical_time"
                        if self.time_basis_code == "physical" else "option.model_time"))
                    self._set_status("error.time_frequency")
                    return
            self._frequency_values[self.time_basis_code] = self.entries[
                "model_frequency"
            ].get()
            self.time_basis_code = requested
            entry = self.entries["model_frequency"]
            entry.delete(0, "end")
            entry.insert(0, self._frequency_values[requested])
        else:
            self.time_basis_code = requested
        self._refresh_time_basis_ui()

    def _refresh_time_basis_ui(self) -> None:
        if not hasattr(self, "time_basis_selector"):
            return
        self.time_basis_selector.configure(values=self._time_basis_values())
        self.time_basis_display.set(
            self._tr(
                "option.physical_time"
                if self.time_basis_code == "physical"
                else "option.model_time"
            )
        )
        if hasattr(self, "frequency_label"):
            self.frequency_label.configure(
                text=self._tr(
                    "field.frequency"
                    if self.time_basis_code == "physical"
                    else "field.model_frequency"
                )
            )
            self.frequency_unit.configure(
                text=self._tr(
                    "unit.hz"
                    if self.time_basis_code == "physical"
                    else "unit.cycles_model_time"
                )
            )
        warning = (
            "status.kinetic_calibrated"
            if self._kinetics_available()
            else "status.kinetic_uncalibrated"
        )
        self.time_warning_display.set(self._tr(warning))
        if self._kinetics_available():
            c = self.time_calibration
            self.time_warning_display.set(self._tr("status.kinetics_details",
                t0=c.t0_seconds, temperature=c.temperature_K, source=c.source))

    def _preset_values(self) -> tuple[str, ...]:
        return tuple(self._tr(values[0]) for values in self.LOAD_PRESETS.values())

    def _on_preset_selected(self, _event=None) -> None:
        selected = self.load_preset.get()
        for code, (key, _mean, _amplitude) in self.LOAD_PRESETS.items():
            if selected == self._tr(key):
                self.load_preset_code = code
                return

    def _apply_load_preset(self) -> None:
        self._on_preset_selected()
        _key, mean, amplitude = self.LOAD_PRESETS[self.load_preset_code]
        if mean is None or amplitude is None:
            return
        for field, value in (
            ("stress_mean_mpa", mean),
            ("stress_amplitude_mpa", amplitude),
        ):
            self.entries[field].delete(0, "end")
            self.entries[field].insert(0, f"{value:g}")
        self._update_stress_context()

    def _update_stress_context(self, _event=None) -> None:
        """Update a non-blocking load interpretation without changing inputs."""

        try:
            values = load_interpretation(
                young_gpa=float(self.entries["young_gpa"].get()),
                stress_mean_mpa=float(self.entries["stress_mean_mpa"].get()),
                stress_amplitude_mpa=float(
                    self.entries["stress_amplitude_mpa"].get()
                ),
            )
        except ValueError:
            self.stress_context.set("—")
            return
        regime_key = f"stress.regime.{values['regime']}"
        self.stress_context.set(
            self._tr(
                "stress.context",
                sigma_min=values["stress_min_mpa"],
                sigma_max=values["stress_max_mpa"],
                reduced_min=values["reduced_stress_min"],
                reduced_max=values["reduced_stress_max"],
                regime=self._tr(regime_key),
            )
        )

    def _probability_scale_values(self) -> tuple[str, str]:
        return self._tr("option.local"), self._tr("option.specimen")

    def _on_probability_scale_selected(self, _event=None) -> None:
        self.probability_scale_code = (
            "specimen"
            if self.probability_scale.get() == self._tr("option.specimen")
            else "local"
        )
        counterpart = {
            ("local", "specimen_survival_probability"): "local_survival_probability",
            ("local", "specimen_initiation_probability"): "local_initiation_probability",
            ("specimen", "local_survival_probability"): "specimen_survival_probability",
            ("specimen", "local_initiation_probability"): "specimen_initiation_probability",
        }.get((self.probability_scale_code, self.field.get()))
        if counterpart is not None:
            self.field.set(counterpart)
        self._update_specimen_probability()

    def _update_specimen_probability(self, _event=None) -> None:
        """Recompute cheap specimen post-processing; never invoke the PDE."""

        if self.result is None:
            self._specimen_status_key = "status.local_result_required"
            self._refresh_specimen_labels()
            return
        try:
            correlation = float(self.entries["correlation_area_mm2"].get())
            stressed = float(self.entries["stressed_area_mm2"].get())
        except ValueError:
            self._specimen_status_key = "status.area_required"
            self.result.pop("specimen_initiation_probability", None)
            self.result.pop("specimen_survival_probability", None)
            self.result.pop("specimen_probability_extrapolation", None)
            self._refresh_specimen_labels()
            self._plot()
            return
        aggregated = aggregate_specimen_probability(
            self.result["local_initiation_probability"],
            correlation_area_mm2=correlation,
            stressed_area_mm2=stressed,
            local_numerical_floor=float(self.result["local_rare_event_floor"]),
            resolution_certified=bool(
                self.result.get("probability_resolution_certified", False)
            ),
        )
        self.result.update(
            {
                "specimen_initiation_probability": (
                    aggregated.specimen_initiation_probability
                ),
                "specimen_survival_probability": (
                    aggregated.specimen_survival_probability
                ),
                "specimen_probability_extrapolation": (
                    aggregated.mathematical_extrapolation
                ),
                "specimen_probability_resolved_mask": aggregated.resolved_mask,
                "N_eff": aggregated.N_eff,
                "specimen_probability_status": aggregated.status,
            }
        )
        if self.result.get("probability_resolution_status") == BELOW_RESOLUTION:
            self._specimen_status_key = "status.below_numerical_resolution"
        elif not bool(self.result.get("probability_resolution_certified", False)):
            self._specimen_status_key = "status.requires_convergence"
        elif aggregated.status == BELOW_RESOLUTION:
            self._specimen_status_key = "status.below_numerical_resolution"
        else:
            self._specimen_status_key = "status.resolved"
        self._refresh_specimen_labels()
        self._plot()

    def _refresh_specimen_labels(self) -> None:
        N_eff = self.result.get("N_eff") if self.result is not None else None
        floor = (
            self.result.get("local_rare_event_floor")
            if self.result is not None
            else None
        )
        plastic_floor = (
            self.result.get("plastic_signal_floor")
            if self.result is not None
            else None
        )
        configurational_barrier = (
            self.result.get("configurational_barrier")
            if self.result is not None
            else None
        )
        opening_barrier = (
            self.result.get("opening_barrier_at_configurational_saddle")
            if self.result is not None
            else None
        )
        local_probability = None
        extrapolation = None
        certified = None
        if self.result is not None:
            local_values = np.asarray(
                self.result.get("local_initiation_probability", []), dtype=float
            )
            if local_values.size:
                local_probability = float(local_values[-1])
            extrapolated_values = np.asarray(
                self.result.get("specimen_probability_extrapolation", []), dtype=float
            )
            if extrapolated_values.size:
                extrapolation = float(extrapolated_values[-1])
            certified_values = np.asarray(
                self.result.get("specimen_initiation_probability", []), dtype=float
            )
            if certified_values.size and np.isfinite(certified_values[-1]):
                certified = float(certified_values[-1])
        self.specimen_N_eff.set(
            f"{self._tr('diagnostic.N_eff')}: "
            + (f"{float(N_eff):.8g}" if N_eff is not None else "—")
        )
        self.local_floor_display.set(
            f"{self._tr('diagnostic.rare_event_floor')}: "
            + (f"{float(floor):.3e}" if floor is not None else "—")
        )
        self.plastic_floor_display.set(
            f"{self._tr('diagnostic.plastic_floor')}: "
            + (f"{float(plastic_floor):.3e}" if plastic_floor is not None else "—")
            + (
                " ("
                + self._tr(
                    "status.plasticity_resolved"
                    if self.result is not None
                    and self.result.get("plastic_resolution_status")
                    == "resolved_interwell_transfer"
                    else "status.plasticity_unresolved"
                )
                + ")"
                if plastic_floor is not None
                else ""
            )
        )
        self.configurational_barrier_display.set(
            f"{self._tr('diagnostic.configurational_barrier')}: "
            + (
                f"{float(configurational_barrier):.3e}"
                if configurational_barrier is not None
                and np.isfinite(float(configurational_barrier))
                else "—"
            )
        )
        self.opening_barrier_display.set(
            f"{self._tr('diagnostic.opening_barrier')}: "
            + (
                f"{float(opening_barrier):.3e}"
                if opening_barrier is not None
                and np.isfinite(float(opening_barrier))
                else "—"
            )
        )
        if self.result is None or not np.asarray(
            self.result.get("plastic_strain", []), dtype=float
        ).size:
            self.plasticity_detail_display.set("—")
        else:
            plastic = np.asarray(self.result.get("plastic_strain", []), dtype=float)
            net = np.asarray(
                self.result.get("accumulated_net_registry_transfer", []), dtype=float
            )
            gross = np.asarray(
                self.result.get("cumulative_gross_registry_activity", []), dtype=float
            )
            status = str(
                self.result.get("plastic_resolution_status", "requires_convergence")
            )
            status_key = (
                "status.plasticity_resolved"
                if status == "resolved_interwell_transfer"
                else "status.plasticity_unresolved"
            )
            self.plasticity_detail_display.set(
                f"{self._tr('diagnostic.max_plastic_strain')}: "
                f"{float(np.max(np.abs(plastic))):.3e}\n"
                f"{self._tr('diagnostic.final_plastic_strain')}: "
                f"{float(plastic[-1]):.3e}\n"
                f"{self._tr('diagnostic.cumulative_net_registry')}: "
                f"{float(net[-1]):.3e}\n"
                f"{self._tr('diagnostic.cumulative_gross_registry')}: "
                f"{float(gross[-1]):.3e}\n"
                f"{self._tr('diagnostic.plasticity_status')}: "
                f"{self._tr(status_key)}"
            )
        self.local_probability_display.set(
            f"{self._tr('diagnostic.local_probability')}: "
            + (f"{local_probability:.8g}" if local_probability is not None else "—")
        )
        self.specimen_extrapolation_display.set(
            f"{self._tr('diagnostic.specimen_extrapolation')}: "
            + (f"{extrapolation:.8g}" if extrapolation is not None else "—")
        )
        self.specimen_certified_display.set(
            f"{self._tr('diagnostic.specimen_certified')}: "
            + (f"{certified:.8g}" if certified is not None else "—")
        )
        self.probability_status_display.set(
            f"{self._tr(self._specimen_status_key)} · "
            f"{self._tr('status.independence_assumption')}"
        )

    def _quality_text(self, code: str) -> str:
        return self._tr(f"quality.{code}")

    def _on_quality_selected(self, _event=None) -> None:
        selected = self.analysis_quality.get()
        self.analysis_quality_code = (
            "resolved"
            if selected == self._tr("quality.resolved_option")
            else "preview"
        )

    def _on_language_selected(self, _event=None) -> None:
        selected = self.language_display.get()
        language = next(
            code for code, display in LANGUAGE_NAMES.items() if display == selected
        )
        if self.localizer.set_language(language):
            self._refresh_language()

    def _refresh_language(self) -> None:
        """Redraw text only; solver inputs and numerical results are untouched."""

        if hasattr(self, "ax") and self.ax.has_data():
            self._remember_view()
        self.root.title(self._tr("app.title"))
        for widget, option, key in self._text_bindings:
            widget.configure(**{option: self._tr(key)})
        for item, key in self._tree_text:
            self.tree.item(item, text=f"  {self._tr(key)}")
        for tab, key in (
            (self.pre_tab, "tab.pre"),
            (self.solve_tab, "tab.solve"),
            (self.post_tab, "tab.post"),
        ):
            self.notebook.tab(tab, text=self._tr(key))
        self.quality_selector.configure(values=self._quality_values())
        self.analysis_quality.set(
            self._tr(f"quality.{self.analysis_quality_code}_option")
        )
        self.preset_selector.configure(values=self._preset_values())
        self.load_preset.set(
            self._tr(self.LOAD_PRESETS[self.load_preset_code][0])
        )
        self.probability_scale_selector.configure(
            values=self._probability_scale_values()
        )
        self.probability_scale.set(self._tr(f"option.{self.probability_scale_code}"))
        self.energy_model_selector.configure(values=self._energy_model_values())
        self.energy_model_display.set(
            self._tr(energy_model_metadata(self.energy_model_code).display_key)
        )
        self._refresh_time_basis_ui()
        self._set_status(self._status_key, **self._status_values)
        self.run_button.configure(text=self._tr(self._button_key))
        self._render_summary()
        self._update_stress_context()
        self._refresh_specimen_labels()
        self._plot()

    def _set_status(self, key: str, **values: object) -> None:
        self._status_key = key
        self._status_values = dict(values)
        self.status.set(self._tr(key, **values))

    def _set_button(self, key: str, *, state: str = "normal") -> None:
        self._button_key = key
        self.run_button.configure(text=self._tr(key), state=state)

    def _show_summary(self, kind: str, payload: dict[str, object] | None = None) -> None:
        self._summary_kind = kind
        self._summary_payload = payload or {}
        self._render_summary()

    def _render_summary(self) -> None:
        if not hasattr(self, "summary"):
            return
        if self._summary_kind == "ready":
            text = self._tr("summary.ready")
        elif self._summary_kind == "solving":
            values = self._summary_payload
            config = values["config"]
            conversion = values["conversion"]
            text = self._tr(
                "summary.solving",
                backend=values["backend"],
                quality=self._quality_text(config.analysis_quality),
                na=config.grid_n_a,
                ns=config.grid_n_s,
                integrator=self._tr(f"integrator.{config.integration_method}"),
                stress_min=config.stress_min_mpa,
                stress_max=config.stress_max_mpa,
                reduced_min=conversion["reduced_stress_min"],
                reduced_max=conversion["reduced_stress_max"],
                period=conversion["model_period"],
                de_fast=conversion["de_fast"],
                de_slow=conversion["de_slow"],
                kappa=conversion["relaxed_axial_kappa"],
                preload=conversion["preload_force"],
            )
        elif self._summary_kind == "complete":
            result = self._summary_payload
            text = self._tr(
                "summary.complete",
                source=self._tr("model.probability_source"),
                quality=self._quality_text(str(result["analysis_quality"])),
                grid=result["grid_shape"],
                integrator=self._tr(f"integrator.{result['integration_method']}"),
                initial_condition=self._tr("model.initial_condition"),
                initial_stress=result["initial_stress_mpa"],
                frequency=result["model_frequency"],
                period=result["model_period"],
                tau_fast=result["tau_fast"],
                tau_slow=result["tau_slow"],
                de_fast=result["de_fast"],
                de_slow=result["de_slow"],
                transfer=result["linear_transfer_magnitude"],
                phase=result["linear_transfer_phase_degrees"],
                kappa=result["relaxed_axial_kappa"],
                survival=float(np.asarray(result["survival"])[-1]),
                initiation=float(np.asarray(result["initiation_probability"])[-1]),
                residual=float(np.max(np.abs(result["mass_balance_residual"]))),
                time_status=self._tr(
                    "status.kinetic_calibrated"
                    if result.get("time_basis") == "physical"
                    else "status.kinetic_uncalibrated"
                ),
            )
        else:
            text = str(self._summary_payload.get("text", ""))
        if self._summary_kind == "solving":
            config = self._summary_payload["config"]
            conversion = self._summary_payload["conversion"]
            text += "\n" + self._energy_model_summary(conversion)
            text += "\n" + self._time_summary_text(config.time_basis, conversion)
        elif self._summary_kind == "complete":
            text += "\n" + self._energy_model_summary(self._summary_payload)
            text += "\n" + self._time_summary_text(
                str(self._summary_payload.get("time_basis", "model")),
                self._summary_payload,
            )
            if self._summary_payload.get("analysis_quality") == "preview":
                text += "\n" + self._tr("status.preview_uncertified") + "\n"
            rows = self._summary_payload.get("per_cycle_diagnostics", [])
            if rows:
                text += "\n" + self._tr("section.cycle_diagnostics") + "\n"
                text += self._tr("cycle.table_header") + "\n"
                for row in rows[:12]:
                    text += (
                        f"{int(row['cycle']):>5d} | {row['mean_total_strain']:.3e} | "
                        f"{row['total_strain_amplitude']:.3e} | "
                        f"{row['mean_normal_strain']:.3e} | "
                        f"{row['mean_intrawell_strain']:.3e} | "
                        f"{row['plastic_strain_increment']:.3e} | "
                        f"{row['cumulative_plastic_strain']:.3e} | "
                        f"{row['net_registry_transfer']:.3e} | "
                        f"{row['gross_registry_activity']:.3e} | "
                        f"{row['absorbed_mass']:.3e} | "
                        f"{row['minimum_configurational_barrier']:.3e} | "
                        f"{row['minimum_opening_barrier']:.3e} | "
                        f"{row['peak_first_passage_flux']:.3e} | "
                        f"{row['survival_at_cycle_end']:.10g}\n"
                    )
        self._set_summary(text)

    def _energy_model_summary(self, values: dict[str, object]) -> str:
        return self._tr(
            "summary.energy_model",
            name=self._tr(str(values["energy_model_display_key"])),
            python_class=values["energy_model_python_class"],
            status=values["energy_model_calibration_status"],
            source=values["energy_model_parameter_source"],
            a0=values["energy_model_a0"],
            b=values["energy_model_b"],
            chi=values["energy_model_chi"],
            kappa=values["energy_model_kappa_axial"],
        )

    def _time_summary_text(
        self, time_basis: str, values: dict[str, object]
    ) -> str:
        if time_basis != "physical":
            return self._tr("summary.model_time_basis")
        return self._tr(
            "summary.physical_time_basis",
            frequency_hz=values["frequency_hz"],
            period_seconds=values["physical_period_seconds"],
            duration_seconds=values["physical_duration_seconds"],
            t0=values["t0_seconds"],
            tau_fast=(
                values.get("tau_fast_seconds")
                or float(values["tau_fast"]) * float(values["t0_seconds"])
            ),
            tau_slow=(
                values.get("tau_slow_seconds")
                or float(values["tau_slow"]) * float(values["t0_seconds"])
            ),
            mobility_a=(
                values.get("physical_M_a")
                or self.time_calibration.M_a_phys_m2_per_J_s
            ),
            mobility_s=(
                values.get("physical_M_s")
                or self.time_calibration.M_s_phys_m2_per_J_s
            ),
            source=(
                values.get("kinetic_calibration_source")
                or self.time_calibration.source
            ),
        )

    def _config(self) -> UIAnalysisConfig:
        self._on_quality_selected()
        resolved = self.analysis_quality_code == "resolved"
        entered_frequency = float(self.entries["model_frequency"].get())
        return UIAnalysisConfig(
            young_gpa=float(self.entries["young_gpa"].get()),
            stress_mean_mpa=float(self.entries["stress_mean_mpa"].get()),
            stress_amplitude_mpa=float(self.entries["stress_amplitude_mpa"].get()),
            model_frequency=(
                entered_frequency if self.time_basis_code == "model" else 1.0
            ),
            physical_frequency_hz=(
                entered_frequency if self.time_basis_code == "physical" else None
            ),
            time_basis=self.time_basis_code,
            time_calibration=(
                self.time_calibration if self.time_basis_code == "physical" else None
            ),
            energy_model=self.energy_model_code,
            cycles=float(self.entries["cycles"].get()),
            steps_per_cycle=int(self.entries["steps_per_cycle"].get()),
            grid_n_a=81 if resolved else 21,
            grid_n_s=91 if resolved else 31,
            integration_method="implicit" if resolved else "explicit",
            analysis_quality="resolved" if resolved else "preview",
        )

    def _start_solve(self) -> None:
        if self.busy:
            self.stop_event.set()
            self._set_status("status.stopping")
            self._set_button("button.stopping", state="disabled")
            return
        try:
            config = self._config()
            config.validate()
            conversion = physical_load_conversion(config)
        except Exception as exc:
            messagebox.showerror(
                self._tr("dialog.invalid_title"),
                self._tr("dialog.invalid_detail", detail=str(exc)),
                parent=self.root,
            )
            return
        self.busy = True
        self.last_config = config
        self.stop_event.clear()
        self.result = None
        self._specimen_status_key = "status.local_result_required"
        self._refresh_specimen_labels()
        self.live_records.clear()
        self._view_limits.clear()
        self._set_button("button.stop")
        self.convergence_button.configure(state="disabled")
        self.progress.start(12)
        self._set_status("status.solving")
        self._show_summary(
            "solving",
            {
                "backend": self.spatial_backend.get(),
                "config": config,
                "conversion": conversion,
            },
        )
        threading.Thread(target=self._solve_worker, args=(config,), daemon=True).start()
        if self._poll_job is None:
            self._poll_job = self.root.after(30, self._drain_queue)

    def _solve_worker(self, config: UIAnalysisConfig) -> None:
        try:
            def emit(record: dict[str, float]) -> None:
                while not self.stop_event.is_set():
                    try:
                        self._queue.put(("record", record), timeout=0.05)
                        return
                    except queue.Full:
                        continue

            result = run_ui_analysis(
                config, record_callback=emit, stop_requested=self.stop_event.is_set
            )
            self._queue.put(("done", result))
        except Exception as exc:
            self._queue.put(("error", str(exc)))

    def _start_convergence_check(self) -> None:
        if self.busy:
            return
        if (
            self.result is None
            or self.last_config is None
            or self.result.get("analysis_quality") != "resolved"
        ):
            self._set_status("status.preview_uncertified")
            return
        self.busy = True
        self.stop_event.clear()
        self.run_button.configure(state="disabled")
        self.convergence_button.configure(
            text=self._tr("button.convergence_running"), state="disabled"
        )
        self.progress.start(12)
        self._set_status("button.convergence_running")
        threading.Thread(target=self._convergence_worker, daemon=True).start()
        if self._poll_job is None:
            self._poll_job = self.root.after(30, self._drain_queue)

    def _convergence_worker(self) -> None:
        try:
            report = run_convergence_check(
                self.last_config,
                reference_result=self.result,
                stop_requested=self.stop_event.is_set,
            )
            self._queue.put(("convergence_done", report.certified_result))
        except Exception as exc:
            self._queue.put(("error", str(exc)))

    def _drain_queue(self) -> None:
        self._poll_job = None
        for _ in range(128):
            try:
                kind, payload = self._queue.get_nowait()
            except queue.Empty:
                break
            if kind == "record":
                self.live_records.append(payload)
                now = time.monotonic()
                if now - self._last_draw > 0.20:
                    self._last_draw = now
                    self._plot()
                self._set_status(
                    (
                        "status.live_physical"
                        if payload.get("physical_time_seconds") is not None
                        else "status.live"
                    ),
                    time=(
                        payload["physical_time_seconds"]
                        if payload.get("physical_time_seconds") is not None
                        else payload["model_time"]
                    ),
                    survival=payload["survival"],
                    strain=payload["strain"],
                )
            elif kind == "done":
                self._solve_done(payload)
                return
            elif kind == "convergence_done":
                self._convergence_done(payload)
                return
            elif kind == "error":
                self._solve_failed(str(payload))
                return
        if self.busy or not self._queue.empty():
            self._poll_job = self.root.after(30, self._drain_queue)

    def _solve_done(self, result: dict[str, object]) -> None:
        self.result = result
        self.busy = False
        self.progress.stop()
        self._set_button("button.run")
        self.convergence_button.configure(
            state=("normal" if result["analysis_quality"] == "resolved" else "disabled")
        )
        self.notebook.select(self.post_tab)
        self._set_status(
            "status.complete",
            records=len(result["model_time"]),
            survival=float(np.asarray(result["survival"])[-1]),
        )
        self._show_summary("complete", result)
        self._update_specimen_probability()
        self._plot(force_auto=True)

    def _convergence_done(self, result: dict[str, object]) -> None:
        self.result = result
        self.busy = False
        self.progress.stop()
        self._set_button("button.run")
        self.convergence_button.configure(
            text=self._tr("button.run_convergence"), state="normal"
        )
        self._set_status("status.convergence_complete")
        self._show_summary("complete", result)
        self._update_specimen_probability()
        self._plot(force_auto=True)

    def _solve_failed(self, detail: str) -> None:
        self.busy = False
        self.progress.stop()
        self._set_button("button.run")
        self.convergence_button.configure(
            state=(
                "normal"
                if self.result is not None
                and self.result.get("analysis_quality") == "resolved"
                else "disabled"
            )
        )
        self.convergence_button.configure(text=self._tr("button.run_convergence"))
        self._set_status("status.failed")
        messagebox.showerror(
            self._tr("dialog.solver_title"),
            self._tr("dialog.solver_detail", detail=detail),
            parent=self.root,
        )

    def _set_summary(self, text: str) -> None:
        self.summary.configure(state="normal")
        self.summary.delete("1.0", "end")
        self.summary.insert("end", text)
        self.summary.configure(state="disabled")

    def _plot_data(self) -> dict[str, np.ndarray] | None:
        if self.result is not None:
            return {key: np.asarray(value) for key, value in self.result.items() if isinstance(value, np.ndarray)}
        if not self.live_records:
            return None
        keys = self.live_records[0].keys()
        return {key: np.asarray([row[key] for row in self.live_records], dtype=float) for key in keys}

    def _plot(self, force_auto: bool = False) -> None:
        field = self.field.get()
        time_basis = (
            str(self.result.get("time_basis", "model"))
            if self.result is not None
            else self.time_basis_code
        )
        text = plot_strings(
            field, self.localizer.language, time_basis=time_basis
        )
        if self._last_field is not None and self._last_field != field and self.ax.has_data():
            self._view_limits[self._last_field] = (self.ax.get_xlim(), self.ax.get_ylim())
        self._last_field = field
        data = self._plot_data()
        self.ax.clear()
        if data is None or "model_time" not in data:
            self.ax.text(
                0.5, 0.5, self._tr("plot.empty"), ha="center", va="center",
                transform=self.ax.transAxes, color=MUTED
            )
            self.ax.set_axis_off()
            self.canvas.draw_idle()
            return
        self.ax.set_axis_on()
        x = np.asarray(data.get("plot_time", data["model_time"]), dtype=float)
        rate_scale = 1.0
        if time_basis == "physical" and self.result is not None:
            rate_scale = 1.0 / float(self.result["t0_seconds"])

        def y_values(key: str) -> np.ndarray:
            values = np.asarray(data[key], dtype=float)
            if key in {
                "first_passage_flux",
                "interwell_net_flux",
                "interwell_gross_flux",
                "net_plastic_flow_rate",
                "gross_configurational_slip_activity",
                "selective_opening_plastic_rate",
            }:
                return values * rate_scale
            return values
        if field.startswith("specimen_") and (
            field not in data or not np.any(np.isfinite(data[field]))
        ):
            self.ax.text(
                0.5, 0.5, self._tr("plot.specimen_unavailable"),
                ha="center", va="center", transform=self.ax.transAxes, color=MUTED,
            )
            self.ax.set_axis_off()
            self.canvas.draw_idle()
            return
        if field == "strain_components":
            for (key, color), label in zip((
                ("normal_strain", "#4e79a7"),
                ("intrawell_strain", "#59a14f"),
                ("plastic_strain", "#e15759"),
                ("strain", "#111827"),
            ), text["legend"]):
                self.ax.plot(x, data[key], label=label, color=color, linewidth=1.8)
            self.ax.legend(loc="best", fontsize=8, frameon=False)
        elif field == "mass_diagnostics":
            for key, label in zip((
                "raw_intact_mass",
                "raw_one_minus_survival",
                "cumulative_absorbed_mass",
                "mass_balance_residual",
                "negative_mass_correction",
                "cumulative_negative_mass_correction",
                "flux_consistency_residual",
            ), text["legend"]):
                self.ax.plot(x, data[key], label=label, linewidth=1.5)
            self.ax.legend(loc="best", fontsize=8, frameon=False)
        elif field == "well_populations":
            if "well_populations" not in data:
                self.ax.text(
                    0.5, 0.5, self._tr("status.local_result_required"),
                    ha="center", va="center", transform=self.ax.transAxes,
                    color=MUTED,
                )
            else:
                for index, well in enumerate(data["well_indices"]):
                    self.ax.plot(
                        x, data["well_populations"][:, index],
                        label=f"P_{int(well):+d}", linewidth=1.5,
                    )
                self.ax.legend(loc="best", fontsize=8, frameon=False)
        elif field == "interwell_flux":
            if "interwell_net_flux" not in data:
                self.ax.text(
                    0.5, 0.5, self._tr("status.local_result_required"),
                    ha="center", va="center", transform=self.ax.transAxes,
                    color=MUTED,
                )
            else:
                boundaries = data["interwell_boundary_lower_index"]
                for index, lower in enumerate(boundaries):
                    tag = f"{int(lower):+d}→{int(lower) + 1:+d}"
                    self.ax.plot(
                        x, y_values("interwell_net_flux")[:, index],
                        label=f"{tag} {self._tr('legend.interwell_net')}",
                        linewidth=1.5,
                    )
                    self.ax.plot(
                        x, y_values("interwell_gross_flux")[:, index],
                        label=f"{tag} {self._tr('legend.interwell_gross')}",
                        linewidth=1.0, linestyle="--",
                    )
                self.ax.legend(loc="best", fontsize=7, frameon=False)
        elif field == "plastic_flow":
            for key, label, style in zip(
                (
                    "net_plastic_flow_rate",
                    "gross_configurational_slip_activity",
                    "selective_opening_plastic_rate",
                ),
                text["legend"],
                ("-", "--", ":"),
            ):
                self.ax.plot(
                    x, y_values(key), label=label, linewidth=1.5, linestyle=style
                )
            self.ax.legend(loc="best", fontsize=7, frameon=False)
        elif field == "registry_transfer":
            for key, label, style in zip(
                ("accumulated_net_registry_transfer", "absorbed_registry_moment"),
                text["legend"],
                ("-", "--"),
            ):
                self.ax.plot(
                    x, data[key], label=label, linewidth=1.5, linestyle=style
                )
            self.ax.legend(loc="best", fontsize=7, frameon=False)
        elif field == "gross_registry_activity":
            for key, label, style in zip(
                (
                    "cumulative_forward_registry_activity",
                    "cumulative_backward_registry_activity",
                    "cumulative_gross_registry_activity",
                ),
                text["legend"],
                ("-", "--", ":"),
            ):
                self.ax.plot(
                    x, data[key], label=label, linewidth=1.5, linestyle=style
                )
            self.ax.legend(loc="best", fontsize=7, frameon=False)
        else:
            key = "applied_stress_mpa" if field == "stress" else field
            self.ax.plot(x, y_values(key), color=ACCENT, linewidth=1.8)
        if field == "specimen_probability_extrapolation":
            self.ax.text(
                0.01,
                0.02,
                self._tr("plot.specimen_extrapolation_note"),
                transform=self.ax.transAxes,
                fontsize=8,
                color=MUTED,
            )
        self.ax.set_xlabel(text["xlabel"])
        self.ax.set_ylabel(text["ylabel"])
        self.ax.set_title(text["title"], loc="left", fontsize=12, fontweight="bold")
        self.ax.grid(True, color="#d5dde3", linewidth=0.7, alpha=0.8)
        formatter = ScalarFormatter(useMathText=True, useOffset=True)
        formatter.set_powerlimits((-3, 3))
        self.ax.yaxis.set_major_formatter(formatter)
        if force_auto:
            self._view_limits.pop(field, None)
        if field in self._view_limits:
            self.ax.set_xlim(*self._view_limits[field][0])
            self.ax.set_ylim(*self._view_limits[field][1])
        else:
            self.ax.relim()
            self.ax.autoscale_view()
        self.figure.tight_layout(pad=1.1)
        self.canvas.draw_idle()

    def _connect_plot_events(self) -> None:
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("button_press_event", self._on_press)
        self.canvas.mpl_connect("button_release_event", self._on_release)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("key_press_event", self._on_key)

    def _remember_view(self) -> None:
        self._view_limits[self.field.get()] = (self.ax.get_xlim(), self.ax.get_ylim())

    def _on_scroll(self, event) -> None:
        if event.inaxes is not self.ax or event.xdata is None or event.ydata is None:
            return
        factor = 0.80 if event.button == "up" else 1.25
        key = (event.key or "").lower()
        zoom_x = "control" not in key and "ctrl" not in key
        zoom_y = "shift" not in key
        if zoom_x:
            left, right = self.ax.get_xlim()
            self.ax.set_xlim(
                event.xdata - (event.xdata - left) * factor,
                event.xdata + (right - event.xdata) * factor,
            )
        if zoom_y:
            low, high = self.ax.get_ylim()
            self.ax.set_ylim(
                event.ydata - (event.ydata - low) * factor,
                event.ydata + (high - event.ydata) * factor,
            )
        self._remember_view()
        self.canvas.draw_idle()

    def _on_press(self, event) -> None:
        if event.dblclick and event.inaxes is self.ax:
            self._reset_view()
            return
        if event.button == 1 and event.inaxes is self.ax and event.xdata is not None:
            self._pan_origin = (
                event.xdata, event.ydata, self.ax.get_xlim(), self.ax.get_ylim()
            )

    def _on_motion(self, event) -> None:
        if self._pan_origin is None or event.inaxes is not self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        x0, y0, xlim, ylim = self._pan_origin
        dx, dy = event.xdata - x0, event.ydata - y0
        self.ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
        self.ax.set_ylim(ylim[0] - dy, ylim[1] - dy)
        self._remember_view()
        self.canvas.draw_idle()

    def _on_release(self, _event) -> None:
        self._pan_origin = None

    def _on_key(self, event) -> None:
        if (event.key or "").lower() == "home":
            self._reset_view()

    def _reset_view(self) -> None:
        self._view_limits.pop(self.field.get(), None)
        self._plot(force_auto=True)

    def _close(self) -> None:
        self.stop_event.set()
        self.root.destroy()

    def run(self) -> None:
        self.root.mainloop()


def startup_smoke() -> dict[str, float]:
    """Exercise imports and canonical conversion without opening a window."""

    conversion = physical_load_conversion(UIAnalysisConfig())
    return {
        "a0": conversion["a0"],
        "relaxed_axial_kappa": conversion["relaxed_axial_kappa"],
        "model_frequency": conversion["model_frequency"],
    }


def main() -> None:
    if "--smoke" in sys.argv:
        values = startup_smoke()
        print(
            f"desktop-ui smoke OK: a0={values['a0']:.16g}, "
            f"kappa={values['relaxed_axial_kappa']:.16g}"
        )
        return
    if not acquire_single_instance():
        return
    DesktopApp().run()


if __name__ == "__main__":
    main()
