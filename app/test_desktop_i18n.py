import tkinter as tk

import numpy as np
import pytest

import app.desktop_ui as desktop_ui
from solver_v1.energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE,
    TWO_ROW_LJ_REFERENCE,
    energy_model_metadata,
)


def test_desktop_language_switch_preserves_state_and_does_not_run_solver(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = 0

    def forbidden_analysis(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise AssertionError("language switching must not run the PDE")

    monkeypatch.setattr(desktop_ui, "run_ui_analysis", forbidden_analysis)
    try:
        app = desktop_ui.DesktopApp()
    except tk.TclError as exc:
        pytest.skip(f"Tk display unavailable: {exc}")
    app.root.withdraw()
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
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = 0

    def forbidden_analysis(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise AssertionError("post-processing controls must not run the PDE")

    monkeypatch.setattr(desktop_ui, "run_ui_analysis", forbidden_analysis)
    try:
        app = desktop_ui.DesktopApp()
    except tk.TclError as exc:
        pytest.skip(f"Tk display unavailable: {exc}")
    app.root.withdraw()
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
