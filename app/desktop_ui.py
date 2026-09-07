"""Lightweight Pre/Solve/Post desktop UI for the probability-PDE solver."""
from __future__ import annotations

import os
import queue
import sys
import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

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
from .solver_adapter import UIAnalysisConfig, physical_load_conversion, run_ui_analysis


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

    def __init__(self) -> None:
        self.localizer = Localizer(DEFAULT_LANGUAGE)
        self.root = tk.Tk()
        self.root.title(self._tr("app.title"))
        self.root.configure(bg=APP_BG)
        self.root.minsize(940, 620)
        self._center(1180, 760)

        self.entries: dict[str, ttk.Entry] = {}
        self.language_display = tk.StringVar(value=LANGUAGE_NAMES[DEFAULT_LANGUAGE])
        self.spatial_backend = tk.StringVar(value="FVM")
        self.analysis_quality = tk.StringVar(value=self._tr("quality.preview_option"))
        self.analysis_quality_code = "preview"
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

        self._styles()
        self._header()
        self._workspace()
        self._statusbar()
        self._connect_plot_events()
        self.root.protocol("WM_DELETE_WINDOW", self._close)
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
        for row, (key, label_key, default, unit_key) in enumerate(self.PARAMS, start=1):
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
        self.pre_tab.columnconfigure(1, weight=1)
        note = self._bind_text(
            ttk.LabelFrame(self.pre_tab, padding=12), "section.scope"
        )
        note.grid(row=len(self.PARAMS) + 2, column=0, columnspan=3, sticky="ew", padx=20, pady=14)
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
        explanation.pack(anchor="w", pady=(0, 16))
        self.run_button = ttk.Button(
            left, text=self._tr("button.run"), style="Accent.TButton",
            command=self._start_solve
        )
        self.run_button.pack(fill="x")
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
        self._set_status(self._status_key, **self._status_values)
        self.run_button.configure(text=self._tr(self._button_key))
        self._render_summary()
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
                time_status=self._tr("model.time_warning"),
            )
        else:
            text = str(self._summary_payload.get("text", ""))
        self._set_summary(text)

    def _config(self) -> UIAnalysisConfig:
        direction = self.entries["tensile_direction"].get().replace(",", " ").split()
        if len(direction) != 3 or np.linalg.norm([float(v) for v in direction]) == 0.0:
            raise ValueError(self._tr("error.direction"))
        self._on_quality_selected()
        resolved = self.analysis_quality_code == "resolved"
        return UIAnalysisConfig(
            young_gpa=float(self.entries["young_gpa"].get()),
            stress_mean_mpa=float(self.entries["stress_mean_mpa"].get()),
            stress_amplitude_mpa=float(self.entries["stress_amplitude_mpa"].get()),
            model_frequency=float(self.entries["model_frequency"].get()),
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
        self.stop_event.clear()
        self.result = None
        self.live_records.clear()
        self._view_limits.clear()
        self._set_button("button.stop")
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
                    "status.live",
                    time=payload["model_time"],
                    survival=payload["survival"],
                    strain=payload["strain"],
                )
            elif kind == "done":
                self._solve_done(payload)
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
        self.notebook.select(self.post_tab)
        self._set_status(
            "status.complete",
            records=len(result["model_time"]),
            survival=float(np.asarray(result["survival"])[-1]),
        )
        self._show_summary("complete", result)
        self._plot(force_auto=True)

    def _solve_failed(self, detail: str) -> None:
        self.busy = False
        self.progress.stop()
        self._set_button("button.run")
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
        text = plot_strings(field, self.localizer.language)
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
        x = np.asarray(data["model_time"], dtype=float)
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
            ), text["legend"]):
                self.ax.plot(x, data[key], label=label, linewidth=1.5)
            self.ax.legend(loc="best", fontsize=8, frameon=False)
        else:
            key = "applied_stress_mpa" if field == "stress" else field
            self.ax.plot(x, data[key], color=ACCENT, linewidth=1.8)
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
