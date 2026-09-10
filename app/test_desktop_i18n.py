import tkinter as tk
import sys

import numpy as np
import pytest

import app.desktop_ui as desktop_ui
from solver_v1.energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE,
    TWO_ROW_LJ_REFERENCE,
    energy_model_metadata,
)


@pytest.fixture(scope="module")
def tk_root():
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        # Headless Unix is optional; broken Tcl installation/access is NOT a
        # successful skip. On Windows run approved GUI checks outside a sandbox
        # that prevents Tcl scripts from loading. Never swallow application bugs.
        message = str(exc).lower()
        if sys.platform != "win32" and ("no display name" in message or "couldn't connect to display" in message):
            pytest.skip(f"headless Tk display: {exc}")
        raise
    root.withdraw()
    yield root
    root.destroy()


def desktop_for_test(tk_root):
    # Isolate widgets/state per test while sharing the production-like SINGLE
    # Tcl interpreter. No retry or skip conceals native/runtime failures.
    root = tk.Toplevel(tk_root)
    root.withdraw()
    return desktop_ui.DesktopApp(root=root)


def test_desktop_language_switch_preserves_state_and_does_not_run_solver(
    monkeypatch: pytest.MonkeyPatch, tk_root,
) -> None:
    calls = 0

    def forbidden_analysis(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise AssertionError("language switching must not run the PDE")

    monkeypatch.setattr(desktop_ui, "run_ui_analysis", forbidden_analysis)
    app = desktop_for_test(tk_root)
    try:
        app.entries["young_gpa"].delete(0, "end")
        app.entries["young_gpa"].insert(0, "70.5")
        app.analysis_quality_code = "resolved"
        app.analysis_quality.set(app._tr("quality.resolved_option"))
        app.field.set("local_initiation_probability")
        result = {
            "model_time": np.array([0.0, 0.01, 0.02]),
            "local_initiation_probability": np.array([0.0, 1.0e-12, 2.0e-12]),
        }
        app.result = result
        assert app.energy_model_code == TWO_ROW_LJ_REFERENCE
        app.energy_model_display.set(
            app._tr(energy_model_metadata(AL_TARGET_BEST_FEASIBLE).display_key)
        )
        app._on_energy_model_selected()
        assert app.energy_model_code == AL_TARGET_BEST_FEASIBLE
        identities = {key: id(value) for key, value in result.items()}
        app._plot()
        app.ax.set_xlim(0.004, 0.016)
        app.ax.set_ylim(0.4e-12, 1.6e-12)
        app._remember_view()

        app.language_display.set("English")
        app._on_language_selected()

        assert calls == 0
        assert app.entries["young_gpa"].get() == "70.5"
        assert app.analysis_quality_code == "resolved"
        assert app.field.get() == "local_initiation_probability"
        assert app.energy_model_code == AL_TARGET_BEST_FEASIBLE
        assert "best-feasible" in app.energy_model_display.get()
        assert {key: id(value) for key, value in app.result.items()} == identities
        np.testing.assert_allclose(app.ax.get_xlim(), (0.004, 0.016))
        np.testing.assert_allclose(app.ax.get_ylim(), (0.4e-12, 1.6e-12))
        assert app.ax.get_title(loc="left") == "Local initiation probability"
        assert app.ax.get_xlabel() == "Dimensionless solver time"
        assert app.ax.get_ylabel() == "Initiation probability [-]"

        app.language_display.set("한국어")
        app._on_language_selected()
        assert calls == 0
        assert app.ax.get_title(loc="left") == "국소 균열 개시 확률"
        assert app.ax.get_xlabel() == "무차원 솔버 시간"
    finally:
        app.root.destroy()


def test_specimen_area_postprocessing_and_stress_warning_do_not_run_solver(
    monkeypatch: pytest.MonkeyPatch, tk_root,
) -> None:
    calls = 0

    def forbidden_analysis(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise AssertionError("post-processing controls must not run the PDE")

    monkeypatch.setattr(desktop_ui, "run_ui_analysis", forbidden_analysis)
    app = desktop_for_test(tk_root)
    try:
        app.result = {
            "model_time": np.array([0.0, 0.1]),
            "local_initiation_probability": np.array([0.0, 1.0e-5]),
            "local_rare_event_floor": 1.0e-8,
            "probability_resolution_certified": True,
            "plastic_signal_floor": 1.0e-9,
        }
        app.field.set("local_initiation_probability")
        for key, value in (
            ("correlation_area_mm2", "0.5"),
            ("stressed_area_mm2", "5.0"),
        ):
            app.entries[key].delete(0, "end")
            app.entries[key].insert(0, value)
        app._update_specimen_probability()
        assert calls == 0
        assert app.result["N_eff"] == 10.0
        assert np.isfinite(app.result["specimen_initiation_probability"][-1])

        original = {
            key: app.entries[key].get()
            for key in ("young_gpa", "stress_mean_mpa", "stress_amplitude_mpa")
        }
        app._update_stress_context()
        assert calls == 0
        assert {
            key: app.entries[key].get()
            for key in original
        } == original
        assert app.stress_context.get()
    finally:
        app.root.destroy()


def test_load_bound_kinetics_switch_units_and_preserve_results(tmp_path, monkeypatch, tk_root):
    from app.test_physical_time_ui import _calibration
    from solver_v1.physical_time import save_time_calibration
    from dataclasses import replace
    def forbidden(*_a, **_kw):
        raise AssertionError("loading calibration/units/language must not rerun PDE")
    monkeypatch.setattr(desktop_ui, "run_ui_analysis", forbidden)
    path = tmp_path / "hypothetical.json"
    calibration = _calibration(.2)
    save_time_calibration(calibration, path)
    app = desktop_for_test(tk_root)
    try:
        assert app._time_basis_values() == (app._tr("option.model_time"),)
        assert app.entries["tensile_direction"].instate(["disabled"])
        result = {"model_time": np.array([0., 1.]), "strain": np.array([0., .001])}
        app.result = result
        app.field.set("strain")
        app._plot()
        app.ax.set_xlim(.2, .7)
        app._remember_view()
        app.load_kinetic_calibration_path(path)
        assert app.result is result
        assert len(app._time_basis_values()) == 2
        app.time_basis_display.set(app._tr("option.physical_time"))
        app._on_time_basis_selected()
        assert float(app.entries["model_frequency"].get()) == pytest.approx(125.)
        assert app._config().effective_model_frequency == pytest.approx(25.)
        assert app.frequency_unit.cget("text") == "Hz"
        app.language_display.set("English")
        app._on_language_selected()
        assert app.time_basis_code == "physical"
        # The old result is STILL model time, not silently relabeled by setup.
        assert app.ax.get_xlabel() == "Dimensionless solver time"
        np.testing.assert_allclose(app.ax.get_xlim(), [.2, .7])
        assert app.result is result
        bad = tmp_path / "wrong.json"
        save_time_calibration(replace(calibration, energy_model_id="wrong"), bad)
        with pytest.raises(ValueError, match="energy_model_id"):
            app.load_kinetic_calibration_path(bad)
        assert app.time_calibration == calibration
        app.energy_model_display.set(app._tr(energy_model_metadata(AL_TARGET_BEST_FEASIBLE).display_key))
        app._on_energy_model_selected()
        assert app.time_basis_code == "model"
        assert float(app.entries["model_frequency"].get()) == pytest.approx(25.)
        assert len(app._time_basis_values()) == 1
    finally:
        app.root.destroy()


@pytest.mark.parametrize("height", [620, 480])
def test_small_window_scroll_and_fixed_solve_actions(monkeypatch, tk_root, height):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("scrolling/resize/language must not run the PDE")
    monkeypatch.setattr(desktop_ui, "run_ui_analysis", forbidden)
    app = desktop_for_test(tk_root)
    try:
        app.root.geometry(f"940x{height}+0+0")
        app.root.deiconify()
        app.notebook.select(app.solve_tab)
        app.root.update()
        values = {key: entry.get() for key, entry in app.entries.items()}
        result = {"model_time": np.array([0., 1.]), "strain": np.array([0., .001])}
        app.result = result
        app.field.set("strain")
        app._plot()
        app.ax.set_xlim(.2, .7)
        app._remember_view()

        def visible(widget):
            assert widget.winfo_ismapped()
            top = widget.winfo_rooty()
            assert app.solve_tab.winfo_rooty() <= top
            assert top + widget.winfo_height() <= (
                app.solve_tab.winfo_rooty() + app.solve_tab.winfo_height()
            )

        for language in ("English", "한국어"):
            app.language_display.set(language)
            app._on_language_selected()
            app.root.update()
            visible(app.run_button)
            visible(app.convergence_button)
            assert app.solve_scroll.canvas.yview()[1] < 1
            quality = app.analysis_quality.get()
            # Scoped handler runs BEFORE Combobox's native wheel-selection.
            app.quality_selector.event_generate("<MouseWheel>", delta=-120)
            app.root.update()
            assert app.solve_scroll.canvas.yview()[0] > 0
            assert app.analysis_quality.get() == quality
            app.solve_scroll.canvas.yview_moveto(1)
            app.root.update()
            visible(app.run_button)
            visible(app.convergence_button)
            app.solve_scroll.canvas.yview_moveto(0)

        assert {key: entry.get() for key, entry in app.entries.items()} == values
        assert app.result is result
        np.testing.assert_array_equal(app.result["strain"], [0., .001])
        np.testing.assert_allclose(app.ax.get_xlim(), [.2, .7])
        # The form bindtag must not intercept Matplotlib zoom events.
        assert app.solve_scroll._tag not in app.canvas.get_tk_widget().bindtags()
        assert app.summary.cget("yscrollcommand")

        app.notebook.select(app.pre_tab)
        app.root.update()
        assert app.pre_scroll.canvas.yview()[1] < 1
        app.pre_scroll.canvas.focus_force()
        app.root.update()
        app.pre_scroll.canvas.event_generate("<Next>")
        app.root.update()
        assert app.pre_scroll.canvas.yview()[0] > 0
    finally:
        app.root.destroy()


def test_scroll_tab_focus_reveals_fields_and_cleans_private_bindings(tk_root):
    from app.scrollable_panel import ScrollablePanel
    from tkinter import ttk
    window = tk.Toplevel(tk_root)
    window.geometry("340x240+0+0")
    try:
        panel = ScrollablePanel(window)
        panel.pack(fill="both", expand=True)
        fields = []
        for index in range(30):
            field = ttk.Entry(panel.content)
            field.insert(0, str(index))
            field.pack(pady=4)
            fields.append(field)
        window.update()
        fields[-1].event_generate("<FocusIn>")
        window.update()
        assert panel.canvas.yview()[0] > .5
        assert fields[-1].winfo_rooty() >= panel.canvas.winfo_rooty()
        assert fields[-1].winfo_rooty() + fields[-1].winfo_height() <= (
            panel.canvas.winfo_rooty() + panel.canvas.winfo_height() + 18
        )
        assert [field.get() for field in fields] == [str(i) for i in range(30)]
        tag = panel._tag
        assert window.bind_class(tag, "<MouseWheel>")
        panel.destroy()
        assert not window.bind_class(tag, "<MouseWheel>")
    finally:
        window.destroy()
