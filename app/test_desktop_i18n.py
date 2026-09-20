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
        # A saved nominal local history is not a spatial field.
        assert app.load_workflow.maps[-1].colorbar is None
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


def test_material_controls_language_and_small_window_preserve_al_result(tk_root, monkeypatch):
    from app.materials import ALUMINUM, SILICON_WAFER, DOPANT_SPECIES
    from app.solver_adapter import UIAnalysisConfig

    def forbidden(*_a, **_k):
        pytest.fail('Material or language selection reran the solver')

    monkeypatch.setattr(desktop_ui, 'run_ui_analysis', forbidden)
    app = desktop_for_test(tk_root)
    try:
        app.root.geometry('940x480+0+0')
        app.root.deiconify()
        app.notebook.select(app.pre_tab)
        app.root.update()
        before = {k: v.get() for k, v in app.entries.items()}
        assert app.material_id == ALUMINUM
        assert app.aluminum_form.winfo_ismapped() and not app.wafer_panel.winfo_ismapped()
        app.last_config = UIAnalysisConfig()
        result = dict(model_time=np.array([0., 1.]), strain=np.array([0., .001]), material_id=ALUMINUM)
        app.result = result
        app.field.set('strain')
        app._plot()
        app.ax.set_xlim(.1, .7)
        app._remember_view()
        app.material_selector.current(1)
        app.material_selector.event_generate('<<ComboboxSelected>>')
        app.root.update()
        assert app.material_id == SILICON_WAFER
        assert app.wafer_panel.winfo_ismapped() and not app.aluminum_form.winfo_ismapped()
        assert app.dopant_concentration_entry.instate(['disabled'])
        app.doping_toggle.invoke()
        app.dopant_concentration.set('2.50E16')
        for index, species in enumerate(DOPANT_SPECIES):
            app.dopant_selector.current(index)
            app.dopant_selector.event_generate('<<ComboboxSelected>>')
            app.root.update()
            config = app._config(allow_face_loads=True)
            assert config.dopant_species == species
            assert config.dopant_concentration_cm3 == 2.5e16
        for language in ('English', '한국어'):
            app.language_display.set(language)
            app._on_language_selected()
            app.root.update()
            assert app.material_selector.get() == app._tr('material.silicon_wafer')
            assert app.dopant_selector.get() == app._tr('dopant.Sb')
            assert app.dopant_concentration.get() == '2.50E16'
            assert app.dopant_concentration_entry.winfo_rootx() + app.dopant_concentration_entry.winfo_width() <= (
                app.pre_scroll.canvas.winfo_rootx() + app.pre_scroll.canvas.winfo_width())
            for button in (app.run_button, app.axial_run_button, app.solid_run_button, app.convergence_button):
                assert button.instate(['disabled'])
            assert app.result is result
            assert app.figure._suptitle.get_text() == app._result_material_text(ALUMINUM)
            np.testing.assert_allclose(app.ax.get_xlim(), [.1, .7])
        app.dopant_concentration.set('NaN')
        assert app.dopant_feedback.get() == app._tr('error.dopant_concentration')
        app.doping_toggle.invoke()
        assert app._config(allow_face_loads=True).dopant_concentration_cm3 is None
        app.material_selector.current(0)
        app.material_selector.event_generate('<<ComboboxSelected>>')
        app.root.update()
        assert app.run_button.instate(['!disabled'])
        assert app.aluminum_form.winfo_ismapped()
        assert {k: v.get() for k, v in app.entries.items()} == before
        assert app._config(allow_face_loads=True).material_id == ALUMINUM
    finally:
        app.root.destroy()


def test_wafer_project_roundtrip_retains_historical_al_and_legacy_defaults(tk_root, tmp_path, monkeypatch):
    from copy import deepcopy
    from app.materials import ALUMINUM, SILICON_WAFER, MaterialInputError
    from app.project_file import capture, restore, save_bundle, load_bundle
    from app.test_project_file import exact
    from app.solver_adapter import UIAnalysisConfig, run_ui_analysis

    app = desktop_for_test(tk_root)
    try:
        app.geometry_workflow.generate()
        config = UIAnalysisConfig(model_frequency=1000., cycles=.2, steps_per_cycle=10)
        app.last_config = config
        app.result = run_ui_analysis(config)
        assert app.result['material_id'] == ALUMINUM and app.result['dopant_species'] is None
        app.material_selector.current(1)
        app._on_material_selected()
        app.doping_toggle.invoke()
        app.dopant_species_code = 'As'
        app.dopant_concentration.set(' 3.125E18 ')
        app._refresh_material_ui()
        app._view_limits = {'strain_components': ((0., .001), (-.01, .02))}
        state = capture(app)
        path = tmp_path / 'wafer.ftgsim'
        save_bundle(path, state)
        exact(state, load_bundle(path))
        monkeypatch.setattr(desktop_ui, 'run_ui_analysis', lambda *_a, **_k: pytest.fail('restore ran PDE'))
        app.material_id = ALUMINUM
        app.dopant_concentration.set('')
        restore(app, load_bundle(path))
        assert app.material_id == SILICON_WAFER and app.doping_enabled.get()
        assert app.dopant_species_code == 'As' and app.dopant_concentration.get() == ' 3.125E18 '
        exact(state['result'], app.result)
        assert app.last_config == config and app.last_config.material_id == ALUMINUM
        assert app._view_limits == state['views']
        assert app._result_material_text(ALUMINUM) in app.summary.get('1.0', 'end')
        assert app.figure._suptitle.get_text() == app._result_material_text(ALUMINUM)
        app.load_workflow.show_map()
        view = app.load_workflow.maps[-1]
        assert app._result_material_text(ALUMINUM) in view.info.get()
        assert view.colorbar is None  # A local history still cannot color a mesh.
        view.window.destroy()
        # Unfinished input is a saveable draft, not a valid numerical configuration.
        state['material']['dopant_concentration_cm3'] = '1e'
        restore(app, state)
        assert capture(app)['material']['dopant_concentration_cm3'] == '1e'
        state['material']['doping_enabled'] = False
        restore(app, state)
        assert not app.doping_enabled.get() and app.dopant_concentration_entry.instate(['disabled'])
        assert app.dopant_species_code == 'As' and app.dopant_concentration.get() == '1e'
        # Corrupt material data must fail before changing any input or result.
        original_result = app.result
        invalid = deepcopy(state)
        invalid['material']['dopant_species'] = 'unknown'
        with pytest.raises(MaterialInputError):
            restore(app, invalid)
        assert app.result is original_result and app.material_id == SILICON_WAFER
        invalid = deepcopy(state)
        invalid['result']['material_id'] = SILICON_WAFER
        with pytest.raises(ValueError, match='Unsupported result material'):
            restore(app, invalid)
        assert app.result is original_result
        # Existing version-2 projects have no material fields anywhere.
        legacy = deepcopy(state)
        del legacy['material']
        for key in ('material_id', 'doping_enabled', 'dopant_species', 'dopant_concentration_cm3'):
            legacy['last_config'].pop(key)
            legacy['result'].pop(key)
        legacy['result'].pop('material_display_key')
        legacy['result'].pop('dopant_concentration_basis')
        restore(app, legacy)
        assert app.material_id == ALUMINUM and app.last_config.material_id == ALUMINUM
        assert not app.doping_enabled.get()
        assert app.run_button.instate(['!disabled'])
    finally:
        app.root.destroy()


def test_silicon_actions_stop_before_load_preflight_and_keep_result(tk_root, monkeypatch):
    app = desktop_for_test(tk_root)
    try:
        app.material_selector.current(1)
        app._on_material_selected()
        original = {'saved': True}
        app.result = original
        messages = []
        monkeypatch.setattr(desktop_ui.messagebox, 'showerror', lambda *a, **_k: messages.append(a))
        monkeypatch.setattr(app, '_config', lambda **_k: pytest.fail('Si reached Al preflight'))
        for kwargs in ({}, {'spatial': True}, {'solid': True}):
            app._start_solve(**kwargs)
        assert len(messages) == 3
        assert all(app._tr('error.silicon_backend_unavailable') in text[1] for text in messages)
        app._start_convergence_check()
        assert not app.busy and app.result is original
        assert app._poll_job is None
    finally:
        app.root.destroy()


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


def test_local_only_result_never_colors_mesh_by_triangle_area(tk_root, monkeypatch):
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
        assert app._tr('map.need_result') in view.info.get()
        assert len(view.picker.ax.collections) == 1
        assert len(view.picker.figure.axes) == 1
        view.picker.ax.view_init(30,50)
        app.language_display.set('English'); app._on_language_selected()
        assert 'No saved spatial specimen result' in view.info.get()
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
        assert len(view.picker.figure.axes) == 1
        assert app.result is result
        assert 'probability_resolution_certified' not in result
        app.entries['correlation_area_mm2'].delete(0,'end')
        view.update()
        assert view.colorbar is None
        assert len(view.picker.figure.axes) == 1
        view.window.destroy()
    finally: app.root.destroy()


def test_actual_spatial_map_round_trip_language_view_and_area_independence(tk_root, tmp_path, monkeypatch):
    from app.axial_specimen import prepare_axial_sections, run_axial_probability
    from app.solver_adapter import UIAnalysisConfig, run_ui_analysis
    from app.test_axial_specimen import necked_specimen, bottom
    from app.project_file import capture, restore, save_bundle, load_bundle
    from app.test_project_file import exact
    from app.specimen_mesh import refine_selected
    from types import SimpleNamespace
    mesh = necked_specimen()
    config = UIAnalysisConfig(model_frequency=1000, cycles=.2, steps_per_cycle=10,
                               stress_mean_mpa=500, stress_amplitude_mpa=100)
    nominal = run_ui_analysis(config)
    nominal['axial_specimen'] = run_axial_probability(config,
        prepare_axial_sections(mesh, bottom(mesh), 6), nominal_result=nominal)
    app = desktop_for_test(tk_root)
    try:
        app.geometry_workflow.geometry = app.geometry_workflow.mesh = mesh
        app.load_workflow.refresh()
        app.result = nominal; app.last_config = config; app.axial_count.set('6')
        app.load_workflow.show_map(); view = app.load_workflow.maps[-1]
        view.window.withdraw()
        monkeypatch.setattr(desktop_ui, 'run_ui_analysis', lambda *_a, **_k: pytest.fail('view ran PDE'))
        monkeypatch.setattr('app.axial_specimen.run_ui_analysis', lambda *_a, **_k: pytest.fail('view ran PDE'))
        assert view.colorbar is not None
        assert '6' in view.info.get()
        assert len(view.picker.figure.axes) == 2
        before = capture(app)
        view.picker.ax.view_init(30, 50)
        bounds = view.picker.ax.get_position(original=True).bounds
        original_colors = view.picker.ax.collections[0].get_facecolor().copy()
        # A_c and field/language/time/zoom are post-processing only.
        app.entries['correlation_area_mm2'].insert(0, '1e-20'); view.update()
        np.testing.assert_array_equal(view.picker.ax.collections[0].get_facecolor(), original_colors)
        for _ in range(3):
            for field in ('stress', 'survival', 'initiation'):
                view.field.set(field); view.update()
                view.index.set(0); view.update()
                view.index.set(len(nominal['axial_specimen']['model_time'])-1); view.update()
            view.picker._zoom(SimpleNamespace(inaxes=view.picker.ax, button='up'))
            view.picker._zoom(SimpleNamespace(inaxes=view.picker.ax, button='down'))
            np.testing.assert_allclose(view.picker.ax.get_position(original=True).bounds, bounds)
        app.language_display.set('English'); app._on_language_selected()
        assert 'uncertified' in view.info.get()
        assert view.picker.ax.azim == 50
        exact(before['result']['axial_specimen'], app.result['axial_specimen'])
        save_bundle(tmp_path/'spatial.ftgsim', capture(app))
        restore(app, load_bundle(tmp_path/'spatial.ftgsim'))
        assert app.axial_count.get() == '6'
        exact(before['result']['axial_specimen'], app.result['axial_specimen'])
        assert view.colorbar is not None
        # A changed mesh never borrows the old face indices or colors.
        app.geometry_workflow.mesh, _ = refine_selected(mesh, [0, 1])
        app.load_workflow.refresh()
        assert view.colorbar is None
        assert app._tr('axial.stale') in view.info.get()
        view.window.destroy()
    finally:
        app.root.destroy()


def test_balanced_multiple_face_loads_run_from_main_button_and_restore(tk_root, tmp_path, monkeypatch):
    import time
    from app.specimen_mesh import cylinder
    from app.solid_mechanics import solid_field
    from app.project_file import capture, restore, save_bundle, load_bundle
    from app.test_project_file import exact
    app = desktop_for_test(tk_root)
    errors = []
    monkeypatch.setattr('tkinter.messagebox.showerror', lambda *a, **k: errors.append(a))
    try:
        mesh = cylinder(radius_mm=2, length_mm=8, target_mm=3)
        app.geometry_workflow.geometry = app.geometry_workflow.mesh = mesh
        app.load_workflow.refresh()
        for key, value in dict(model_frequency='1000', cycles='.2', steps_per_cycle='10',
                               stress_mean_mpa='100', stress_amplitude_mpa='30',
                               tensile_direction='0 0 1').items():
            app.entries[key].delete(0, 'end'); app.entries[key].insert(0, value)
        app.axial_count.set('2')
        app.solid_target.set('2')
        app.entries['correlation_area_mm2'].insert(0, '1')
        app.entries['stressed_area_mm2'].insert(0, '100')
        app.field.set('specimen_probability_extrapolation')
        live_extrapolation = []
        app.load_workflow.tensor_text.set('0,0,0;0,0,0;0,0,normal_mean+normal_amp*sin(2*pi*f*t)')
        for region in ('top', 'bottom'):
            app.load_workflow.region_code = region
            app.load_workflow.apply()
        assert len(app.load_workflow.loads) == 2
        app._start_solve()  # stored balanced loads select 3D, no local-only bypass
        deadline = time.monotonic()+150
        while app.busy and time.monotonic() < deadline:
            tk_root.update(); time.sleep(.01)
            if app.result is None and app.live_records:
                data = app._plot_data()
                live_extrapolation.append(float(data['specimen_probability_extrapolation'][-1]))
        assert not app.busy and not errors
        assert live_extrapolation, '3D solve must publish the reference PDE before completion'
        assert app.result is not None and 'solid_specimen' in app.result
        expected = -np.expm1(100*np.log1p(-app.result['local_initiation_probability']))
        np.testing.assert_array_equal(app.result['specimen_probability_extrapolation'], expected)
        np.testing.assert_array_equal(app.ax.lines[0].get_ydata(), expected)
        assert '—' not in app.specimen_extrapolation_display.get()
        assert np.isnan(app.result['specimen_initiation_probability']).all()
        solid = app.result['solid_specimen']
        assert solid['target_mm'] == 2.
        assert app.geometry_workflow.target.get() == '3'
        assert len(app.load_workflow.loads) == 2
        np.testing.assert_allclose(solid_field(solid, 'zz', 0), 100., atol=1e-8)
        assert len(solid['sample_cells']) == 1
        view = app.load_workflow.maps[-1]; view.window.withdraw()
        assert view.colorbar is not None
        before = capture(app)
        monkeypatch.setattr('app.solid_mechanics.run_solid_probability', lambda *_a, **_k: pytest.fail('display reran solver'))
        for field in ('mises', 'xx', 'yy', 'zz', 'xy', 'xz', 'yz', 'stress', 'initiation'):
            view.field.set(field); view.update()
            assert view.colorbar is not None
        view.field.set('hazard_density'); view.update()
        assert view.colorbar is None
        assert app._tr('risk.need_volume') in view.info.get()
        app.correlation_volume.set('1')
        assert view.colorbar is not None
        for field in ('hazard_density', 'cell_risk', 'unit_risk'):
            view.field.set(field); view.update()
            assert view.colorbar is not None
        bounds = view.picker.ax.get_position(original=True).bounds
        limits = view.picker.ax.get_xlim3d(), view.picker.ax.get_ylim3d(), view.picker.ax.get_zlim3d()
        view.slice_fraction.set(.5); view.update()
        assert len(np.unique(view._surface_cells)) < len(solid['tetrahedra'])
        np.testing.assert_allclose(view.picker.ax.get_position(original=True).bounds, bounds)
        np.testing.assert_allclose([view.picker.ax.get_xlim3d(), view.picker.ax.get_ylim3d(),
                                    view.picker.ax.get_zlim3d()], limits)
        view.slice_fraction.set(1); view.update()
        app.language_display.set('English'); app._on_language_selected()
        assert '3D' in view.window.title()
        state = capture(app)
        state['result'] = dict(state['result'])
        # Extrapolation is cheap post-processing of saved raw data and areas.
        # A missing or stale display cache must not hide area extrapolation.
        for key in ('specimen_initiation_probability', 'specimen_survival_probability',
                    'specimen_probability_extrapolation', 'specimen_probability_resolved_mask',
                    'N_eff', 'specimen_probability_status'):
            state['result'].pop(key, None)
        for cached in (False, True):
            if cached:
                state['result']['specimen_probability_extrapolation'] = np.full_like(expected, .5)
                state['result']['N_eff'] = -1.
            save_bundle(tmp_path/'actual-solid.ftgsim', state)
            restore(app, load_bundle(tmp_path/'actual-solid.ftgsim'))
            np.testing.assert_array_equal(app.result['specimen_probability_extrapolation'], expected)
            np.testing.assert_array_equal(app.ax.lines[0].get_ydata(), expected)
            assert app.result['N_eff'] == 100.
            assert '—' not in app.specimen_extrapolation_display.get()
            assert app._specimen_status_key == 'status.below_numerical_resolution'
        exact(before['result']['solid_specimen'], app.result['solid_specimen'])
        assert app.correlation_volume.get() == '1'
        assert app.solid_target.get() == '2'
        assert view.colorbar is not None
        assert str(app.convergence_button['state']) == 'disabled'
        view.window.destroy()
    finally:
        app.stop_event.set()
        app.root.destroy()


def test_volume_size_error_keeps_counts_through_worker_and_localizes(tk_root, monkeypatch):
    from app.solid_mechanics import SolidMeshLimitError, MAX_SOLID_NODES
    from app.solver_adapter import UIAnalysisConfig
    from app.project_file import capture, restore
    app = desktop_for_test(tk_root)
    messages = []
    try:
        def too_large(*_a, **_k):
            raise SolidMeshLimitError(MAX_SOLID_NODES, 1245123, .3)
        monkeypatch.setattr('app.solid_mechanics.run_solid_probability', too_large)
        monkeypatch.setattr('tkinter.messagebox.showerror', lambda *a, **_k: messages.append(a))
        app._solve_worker(UIAnalysisConfig(), solid_request={})
        kind, payload = app._queue.get_nowait()
        assert kind == 'error' and payload['key'] == 'solid.size_detail'
        app._queue.put((kind, payload)); app._drain_queue()
        assert messages and f'{MAX_SOLID_NODES:,}' in messages[-1][1] and '1,245,123' in messages[-1][1]
        assert '체적 메시 목표 길이' in messages[-1][1]
        state = capture(app)
        state.pop('solid_target_mm')  # old projects used the surface target for both
        state['geometry_inputs']['target'] = '7'
        restore(app, state)
        assert app.solid_target.get() == '7'
    finally:
        app.root.destroy()


def test_unassigned_stress_stays_draft_and_basic_run_requests_unloaded_3d(tk_root, monkeypatch):
    from types import SimpleNamespace
    from .specimen_mesh import cylinder
    app = desktop_for_test(tk_root)
    calls = []
    try:
        app.geometry_workflow.mesh = cylinder(radius_mm=2.,length_mm=8.,target_mm=3.)
        app.load_workflow.refresh()
        app.load_workflow.normal_mean.set('1200')
        app.load_workflow.normal_amplitude.set('1500')
        app.load_workflow.refresh_balance()
        assert '무하중' in app.load_workflow.balance_summary.get()
        monkeypatch.setattr('tkinter.messagebox.showerror', lambda *a,**k: pytest.fail(str(a)))
        monkeypatch.setattr(desktop_ui.threading,'Thread',lambda **k:SimpleNamespace(start=lambda:calls.append(k['args'])))
        app._start_solve()
        config, sections, request = calls[0]
        assert config.stress_mean_mpa == config.stress_amplitude_mpa == 0.
        assert sections is None and request['loads'] == ()
        assert app.entries['stress_mean_mpa'].get() == '1200'
        assert app.entries['stress_amplitude_mpa'].get() == '1500'
    finally:
        app.root.destroy()


def test_initial_ensemble_choice_survives_language_and_project_restore(tk_root, monkeypatch):
    from app.project_file import capture, restore
    app = desktop_for_test(tk_root)
    try:
        assert app._config().initialization == 'loaded_gibbs'
        app.initialization.set('zero_load_gibbs')
        assert app._config().initialization == 'zero_load_gibbs'
        state = capture(app)
        monkeypatch.setattr('app.desktop_ui.run_ui_analysis', lambda *_a, **_k: pytest.fail('display ran PDE'))
        app.language_display.set('English'); app._on_language_selected()
        assert app.initialization.get() == 'zero_load_gibbs'
        app.initialization.set('loaded_gibbs')
        restore(app, state)
        assert app.initialization.get() == 'zero_load_gibbs'
        del state['initialization']
        restore(app, state)
        assert app.initialization.get() == 'loaded_gibbs'
    finally:
        app.root.destroy()


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
        # The 3D axial projection is now editable; it does not unlock the
        # independent physical-time/kinetic calibration gate.
        assert not app.entries["tensile_direction"].instate(["disabled"])
        app.entries['tensile_direction'].delete(0, 'end')
        app.entries['tensile_direction'].insert(0, '0 0 1')
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
        assert app.entries['tensile_direction'].get() == '0 0 1'
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
            visible(app.axial_run_button)
            visible(app.solid_run_button)
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
