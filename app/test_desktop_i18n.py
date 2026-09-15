import tkinter as tk
import sys

import numpy as np
import pytest

import app.desktop_ui as desktop_ui

def test_real_result_and_setup_restore_without_pde(tk_root,tmp_path,monkeypatch):
    from app.project_file import save_bundle, load_bundle, capture, restore
    from app.test_project_file import exact
    from app.solver_adapter import UIAnalysisConfig, run_ui_analysis
    app=desktop_for_test(tk_root)
    try:
        config=UIAnalysisConfig(model_frequency=1000.,cycles=.2,steps_per_cycle=10)
        result=run_ui_analysis(config)
        app.last_config=config
        app.geometry_workflow.generate()
        app.load_workflow.normal_mean.set('125')
        app.load_workflow.apply()
        app.local_only.set(True)
        app.entries['correlation_area_mm2'].insert(0,'.1')
        app.result=result; app.field.set('strain_components')
        app._view_limits={'strain_components':((0.,.001),(-.01,.02))}
        before=capture(app)
        expected=before['result']
        path=tmp_path/'run.ftgsim'; save_bundle(path,before)
        exact(before,load_bundle(path))
        monkeypatch.setattr('app.desktop_ui.run_ui_analysis',lambda *_a,**_k:pytest.fail('restore reran solver'))
        app.result=None; app.load_workflow.loads.clear(); app.load_workflow.normal_mean.set('999')
        restore(app,load_bundle(path))
        exact(expected,app.result)
        assert app.entries['stress_mean_mpa'].get()=='125'
        assert app.last_config==config
        assert len(app.load_workflow.loads)==1 and app.local_only.get()
        np.testing.assert_array_equal(app.geometry_workflow.mesh.vertices,before['mesh']['vertices'])
        assert app._view_limits==before['views']
        app.load_workflow.show_map()
        assert app.load_workflow.maps[-1].colorbar is not None
        exact(expected,app.result)
        app.load_workflow.maps[-1].window.destroy()
        # Save/open dialog path also operates without a new solve.
        monkeypatch.setattr('tkinter.filedialog.asksaveasfilename',lambda **k:str(path))
        monkeypatch.setattr('tkinter.messagebox.askyesno',lambda *a,**k:True)
        monkeypatch.setattr('tkinter.messagebox.showerror',lambda *a,**k:pytest.fail(str(a)))
        app.project_files.save(as_new=True)
        app.project_files.open(path)
        exact(expected,app.result)
    finally: app.root.destroy()

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


def test_setup_editor_and_ai_threads_are_inert_until_review(tk_root, monkeypatch):
    import time
    from app.setup_language_view import SetupLanguageView
    from app.ai_chat_view import AIChatView
    from app.test_aft_ai_controller import MockTransport
    app = desktop_for_test(tk_root)
    editor = chat = None
    try:
        monkeypatch.setattr(desktop_ui, 'run_ui_analysis', lambda *_: pytest.fail('setup/chat ran PDE'))
        app.geometry_workflow.generate()
        app.load_workflow.apply()
        original = list(app.load_workflow.loads)
        editor = SetupLanguageView(app)
        assert editor.validate() is not None
        assert app.load_workflow.loads == original
        transport = MockTransport('AFT 1 unvalidated text only')
        chat = AIChatView(app, lambda: editor.text.get('1.0','end'), transport=transport)
        chat.model.set('mock'); chat.new_session()
        chat.message.insert('1.0','Explain single crystal limitations')
        chat.reviewed.set(True); chat.attach.set(True); chat.prepare()
        assert not transport.calls
        chat.attach.set(False); chat.send()
        assert not transport.calls
        app.field.set('strain')
        app.result = {'model_time':np.array([0.,1.]), 'strain':np.array([0.,.001]), 'cumulative_absorbed_mass':np.array([0.,1e-16]),
                      'local_rare_event_floor':1e-14, 'unselected_private_note':'do not transmit'}
        chat.attach_result.set(True)
        assert 'do not transmit' not in chat.result_summary()
        chat.prepare(); chat.send()
        deadline = time.monotonic()+2
        controller, sid = chat.controllers[chat.current]
        while len(controller.history(sid)) < 2 and time.monotonic() < deadline:
            app.root.update(); time.sleep(.01)
        assert len(controller.history(sid)) == 2
        chat.new_session()
        other, sid2 = chat.controllers[chat.current]
        assert other.history(sid2) == ()
        assert app.load_workflow.loads == original
        assert len(transport.calls) == 1
        import json
        payload = json.loads(transport.calls[0][1]['input'][-1]['content'])
        assert 'result_summary' in payload['selected_non_secret_summaries']
        assert 'setup_summary' not in payload['selected_non_secret_summaries']
        app.language_display.set('English'); app._on_language_selected()
        assert len(transport.calls) == 1
        assert chat.window.title() == app._tr('ai.title')
        assert editor.window.title() == app._tr('setup.title')
        chat.close(); chat = None
        editor.window.destroy(); editor = None
        app.language_display.set('한국어'); app._on_language_selected()
    finally:
        if chat is not None: chat.close()
        if editor is not None: editor.window.destroy()
        app.root.destroy()


def test_matrix_multiple_loads_and_one_balance_confirmation(tk_root, monkeypatch):
    from app.load_balance import resultant
    app = desktop_for_test(tk_root)
    try:
        monkeypatch.setattr(desktop_ui, 'run_ui_analysis', lambda *_: pytest.fail('setup ran PDE'))
        app.geometry_workflow.generate()
        load = app.load_workflow
        assert app.entries['model_frequency'].master is app.load_inputs
        load.tensor_text.set('0,0,0;0,0,0;0,0,100+20*sin(2*pi*f*t)')
        assert load.matrix_vars[2][2].get() == '100+20*sin(2*pi*f*t)'
        load.matrix_vars[2][2].set('100+10*sin(2*pi*f*t)')
        assert '100+10*sin' in load.tensor_text.get()
        load.region_code = 'top'; load.apply()
        load.apply()
        assert len(load.loads) == 2
        load.region_code = 'bottom'
        calls = []
        monkeypatch.setattr('app.load_workflow.messagebox.askyesno', lambda *a, **k: False)
        load.balance()
        assert load.correction is None
        monkeypatch.setattr('app.load_workflow.messagebox.askyesno', lambda *a, **k: calls.append(1) or True)
        load.balance(); load.balance()
        assert len(calls) == 1
        for t in (0., .0123, .08):
            np.testing.assert_allclose(resultant(app.geometry_workflow.mesh, load.traction_at(t)), 0, atol=1e-8)
        correction = load.correction
        app.language_display.set('English'); app._on_language_selected()
        assert load.correction is correction and len(load.loads) == 2
        with pytest.raises(ValueError): app._config()
        app.local_only.set(True)
        config = app._config()
        assert config.stress_mean_mpa == float(load.normal_mean.get())
        assert config.stress_amplitude_mpa == float(load.normal_amplitude.get())
        assert load.correction is correction and len(load.loads) == 2
        app.local_only.set(False)
        with pytest.raises(ValueError): app._config()
        load.load_list.selection_set(0); load.remove_load()
        assert load.correction is None and len(load.loads) == 1
    finally:
        app.root.destroy()


@pytest.mark.parametrize('format_name', ['ascii_stl', 'binary_stl', 'obj'])
def test_imported_surface_is_visible_and_pickable_without_generate(
    tk_root, tmp_path, monkeypatch, format_name,
):
    import struct
    from types import SimpleNamespace
    from mpl_toolkits.mplot3d import proj3d

    vertices = np.array([[10, 20, 30], [14, 20, 30], [14, 24, 30], [10, 24, 30],
                         [10, 20, 32], [14, 20, 32], [14, 24, 32], [10, 24, 32]], dtype=float)
    faces = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
                      [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                      [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    path = tmp_path / ('translated.obj' if format_name == 'obj' else 'translated.stl')
    if format_name == 'binary_stl':
        path.write_bytes(b' '*80 + struct.pack('<I', len(faces)) + b''.join(
            struct.pack('<12fH', 0, 0, 0, *triangle.ravel(), 0) for triangle in vertices[faces]))
    elif format_name == 'ascii_stl':
        path.write_text('solid test\n' + ''.join(
            'facet normal 0 0 0\nouter loop\n' + ''.join(
                'vertex ' + ' '.join(map(str, vertex)) + '\n' for vertex in triangle)
            + 'endloop\nendfacet\n' for triangle in vertices[faces]) + 'endsolid test\n')
    else:
        path.write_text(''.join('v ' + ' '.join(map(str, v)) + '\n' for v in vertices)
                        + ''.join('f ' + ' '.join(map(str, f+1)) + '\n' for f in faces))
    monkeypatch.setattr('app.geometry_workflow.filedialog.askopenfilename', lambda **_: str(path))
    monkeypatch.setattr(desktop_ui, 'run_ui_analysis', lambda *_a, **_k: pytest.fail('import ran PDE'))
    app = desktop_for_test(tk_root)
    try:
        flow, load = app.geometry_workflow, app.load_workflow
        flow.generate()
        load.apply()
        assert load.loads
        old_mesh = flow.mesh
        values = {k: entry.get() for k, entry in app.entries.items()}
        result = {'model_time': np.array([0., 1.]), 'strain': np.array([0., .001])}
        app.result = result; app.field.set('strain')
        flow.unit_scale.set('0.5')
        flow.import_model()
        mesh = flow.mesh
        assert mesh is not None  # Imported triangles are usable before Generate.
        assert mesh is flow.geometry and mesh is not old_mesh
        assert mesh.source == path.name and flow.cylinder_dimensions is None
        np.testing.assert_array_equal(mesh.vertices[mesh.faces], .5*vertices[faces])
        assert load.picker.mesh is mesh and len(load.picker.ax.collections) == 1
        assert not load.loads and load.correction is None
        assert len(load._indices()) == 2

        app.root.geometry('940x480+0+0')
        app.root.deiconify()
        app.notebook.select(app.load_tab)
        app.root.update()
        picker = load.picker
        canvas = picker.canvas.get_tk_widget()
        assert canvas.winfo_ismapped() and canvas.winfo_width() > 100
        picker.ax.view_init(elev=75, azim=15)
        picker.canvas.draw()
        center = mesh.vertices[mesh.faces[load._indices()[0]]].mean(axis=0)
        x, y, _ = proj3d.proj_transform(*center, picker.ax.get_proj())
        picker._click(SimpleNamespace(button=1, inaxes=picker.ax, xdata=x, ydata=y, key=None))
        assert load.region_code == 'custom' and len(load._indices()) == 2
        assert np.all(mesh.vertices[mesh.faces[load._indices()], 2] == 16.)
        load.apply()
        assert load.applied.area_mm2 == pytest.approx(4.)
        stored_load = load.applied
        view = (picker.ax.elev, picker.ax.azim, picker.ax.get_xlim())
        app.language_display.set('English'); app._on_language_selected()
        assert load.applied is stored_load and flow.mesh is mesh
        assert view == (picker.ax.elev, picker.ax.azim, picker.ax.get_xlim())
        assert app.result is result and values == {k: entry.get() for k, entry in app.entries.items()}
        # Explicit refinement still works and invalidates old triangle assignments.
        flow.target.set('1.5'); flow.generate()
        assert len(flow.mesh.faces) > len(mesh.faces)
        assert load.picker.mesh is flow.mesh and not load.loads
    finally:
        app.root.destroy()


def test_geometry_mesh_workflow_preserves_pde_and_language(monkeypatch, tk_root):
    def forbidden(*args, **kwargs):
        raise AssertionError('Geometry must not run probability analysis')
    monkeypatch.setattr(desktop_ui, 'run_ui_analysis', forbidden)
    app = desktop_for_test(tk_root)
    try:
        assert len(app.notebook.tabs()) == 6
        assert app.notebook.select() == str(app.model_tab)
        inputs = {k: v.get() for k, v in app.entries.items()}
        result = {'model_time': np.array([0., 1.]), 'strain': np.array([0., .001])}
        app.result = result; app.field.set('strain')
        flow = app.geometry_workflow
        flow.generate()
        mesh = flow.mesh
        assert mesh.closed
        assert len(flow.ax.collections) == 1
        load = app.load_workflow
        load.region_code = "top"
        load.refresh()
        load.normal_mean.set("120")
        load.shear_amplitude.set("30")
        load.apply()
        assert load.applied.face_indices
        assert load.applied.normal_mean_mpa == 120
        assert load.applied.shear_amplitude_mpa == 30
        load.preview_tensor()
        assert 'MPa' in load.tensor_preview.get()
        assert 'sin' in load.applied.tensor_expression
        assert app.result is result
        app.language_display.set('English'); app._on_language_selected()
        assert app.notebook.tab(app.mesh_tab, 'text') == '2  MESH'
        assert 'Triangles' in flow.mesh_description.get()
        assert flow.mesh is mesh and app.result is result
        flow.radius.set('6'); flow.use_cylinder()
        assert flow.mesh is None
        assert app.result is result
        inputs['stress_mean_mpa'] = '120'
        assert inputs == {k: v.get() for k, v in app.entries.items()}
    finally:
        app.root.destroy()


def test_mesh_map_is_uncertified_and_preserves_result(tk_root, monkeypatch):
    app = desktop_for_test(tk_root)
    try:
        monkeypatch.setattr(desktop_ui, 'run_ui_analysis', lambda *_: pytest.fail('map reran PDE'))
        app.geometry_workflow.generate()
        app.entries['correlation_area_mm2'].insert(0,'.1')
        result = {'model_time':np.array([0.,1.]), 'local_initiation_probability':np.array([0.,1e-16])}
        app.result = result
        app.field.set('local_initiation_probability')
        app.load_workflow.show_map()
        view = app.load_workflow.maps[-1]
        view.index.set(1); view.update()
        assert '1.00000000e-16' in view.info.get()
        assert len(view.picker.ax.collections) == 1
        assert len(view.picker.figure.axes) == 2
        view.picker.ax.view_init(30,50)
        app.language_display.set('English'); app._on_language_selected()
        assert 'Uncertified' in view.info.get()
        assert view.picker.ax.azim == 50
        bounds = view.picker.ax.get_position(original=True).bounds
        from types import SimpleNamespace
        limits = np.array([view.picker.ax.get_xlim3d(), view.picker.ax.get_ylim3d(), view.picker.ax.get_zlim3d()])
        for _ in range(20):
            view.index.set(0); view.update()
            view.index.set(1); view.update()
            view.picker._zoom(SimpleNamespace(inaxes=view.picker.ax, button='up'))
            assert np.ptp(view.picker.ax.get_xlim3d()) < np.ptp(limits[0])
            view.picker._zoom(SimpleNamespace(inaxes=view.picker.ax, button='down'))
            np.testing.assert_allclose(view.picker.ax.get_position(original=True).bounds, bounds)
            np.testing.assert_allclose([view.picker.ax.get_xlim3d(), view.picker.ax.get_ylim3d(), view.picker.ax.get_zlim3d()], limits)
        assert len(view.picker.figure.axes) == 2
        assert app.result is result
        assert 'probability_resolution_certified' not in result
        app.entries['correlation_area_mm2'].delete(0,'end')
        view.update()
        assert view.colorbar is None
        assert len(view.picker.figure.axes) == 1
        view.window.destroy()
    finally: app.root.destroy()


def test_single_load_input_and_projected_click(tk_root):
    from types import SimpleNamespace
    from mpl_toolkits.mplot3d import proj3d
    app = desktop_for_test(tk_root)
    try:
        load = app.load_workflow
        load.normal_mean.set('125')
        load.normal_amplitude.set('25')
        assert app._config().stress_mean_mpa == 125
        assert app._config().stress_amplitude_mpa == 25
        assert app.entries['stress_mean_mpa'].master is app.load_inputs
        load.shear_mean.set('1')
        with pytest.raises(ValueError): app._config()
        load.shear_mean.set('0')
        app.geometry_workflow.generate()
        picker = load.picker
        picker.ax.view_init(elev=75, azim=15)
        picker.canvas.draw()
        x, y, _ = proj3d.proj_transform(1., 0., 30., picker.ax.get_proj())
        picker._click(SimpleNamespace(button=1, inaxes=picker.ax, xdata=x, ydata=y, key=None))
        assert load.region_code == 'custom'
        mesh = app.geometry_workflow.mesh
        assert len(load._indices())
        assert np.all(mesh.vertices[mesh.faces[load._indices()], 2] == 30)
        view = (picker.ax.elev, picker.ax.azim, picker.ax.get_xlim())
        app.language_display.set('English'); app._on_language_selected()
        assert view == (picker.ax.elev, picker.ax.azim, picker.ax.get_xlim())
        load.apply()
        assert load.applied is not None
        app.geometry_workflow.generate()
        assert load.applied is None
    finally:
        app.root.destroy()


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

        app.notebook.select(app.load_tab)
        app.root.update()
        assert app.load_scroll.canvas.yview()[1] < 1
        app.load_scroll.canvas.focus_force()
        app.root.update()
        app.load_scroll.canvas.event_generate("<Next>")
        app.root.update()
        assert app.load_scroll.canvas.yview()[0] > 0
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
