"""Version 2 .ftgsim: portable setup and exact numerical results, no pickle."""
from dataclasses import asdict
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
import zipfile
import numpy as np

LIMIT = 256 * 1024 * 1024
SCHEMA = 'aft.probability-project/2'


def save_bundle(path, state):
    path = Path(path)
    payloads = {}
    def pack(x):
        if isinstance(x, np.ndarray):
            if x.dtype.hasobject: raise ValueError('Object arrays cannot be saved')
            name = f'arrays/{len(payloads)}.npy'
            buf = io.BytesIO(); np.save(buf, x, allow_pickle=False)
            payloads[name] = buf.getvalue()
            return {'array': name}
        if isinstance(x, np.generic): return pack(x.item())
        if isinstance(x, dict): return {'dict': {str(k): pack(v) for k,v in x.items()}}
        if isinstance(x, (tuple, list)): return {'tuple' if isinstance(x, tuple) else 'list': [pack(v) for v in x]}
        if isinstance(x, float) and not math.isfinite(x): return {'float': repr(x)}
        if x is None or type(x) in (bool, int, float, str): return x
        raise ValueError(f'Unsupported saved value: {type(x).__name__}')
    tree = pack(state)
    payloads['state.json'] = json.dumps(tree, ensure_ascii=False, allow_nan=False).encode('utf-8')
    manifest = dict(format='ftgsim', schema_version='2.0.0', application_schema=SCHEMA,
                    generator='Al Fatigue probability PDE desktop',
                    checksums_sha256={k:hashlib.sha256(v).hexdigest() for k,v in payloads.items()})
    payloads['ftgsim-manifest.json'] = json.dumps(manifest).encode()
    if len(payloads)>2048 or sum(map(len,payloads.values()))>LIMIT: raise ValueError('Project exceeds size limit')
    handle, tmp = tempfile.mkstemp(prefix=path.stem+'-', suffix='.tmp', dir=path.parent)
    os.close(handle)
    try:
        with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for name,data in payloads.items(): z.writestr(name,data)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
    return path


def load_bundle(path):
    with zipfile.ZipFile(path) as z:
        infos=z.infolist(); names=[i.filename for i in infos]
        if len(names)>2048 or len(set(names))!=len(names) or sum(i.file_size for i in infos)>LIMIT:
            raise ValueError('Invalid project size or duplicate members')
        manifest=json.loads(z.read('ftgsim-manifest.json'))
        if manifest.get('application_schema')!=SCHEMA:
            raise ValueError('Legacy or unsupported .ftgsim model. Open it with its original AlFatigue executable; it is not a probability-PDE project.')
        hashes=manifest['checksums_sha256']
        if set(names)!=set(hashes)|{'ftgsim-manifest.json'}: raise ValueError('Project member mismatch')
        data={}
        for name,digest in hashes.items():
            raw=z.read(name)
            if hashlib.sha256(raw).hexdigest()!=digest: raise ValueError('Project checksum mismatch')
            data[name]=raw
    def unpack(x):
        if not isinstance(x,dict): return x
        if len(x)!=1: raise ValueError('Invalid project value')
        kind,value=next(iter(x.items()))
        if kind=='dict': return {k:unpack(v) for k,v in value.items()}
        if kind in ('list','tuple'):
            v=[unpack(y) for y in value]; return tuple(v) if kind=='tuple' else v
        if kind=='float' and value in ('nan','inf','-inf'): return float(value)
        if kind=='array':
            raw=data[value]; stream=io.BytesIO(raw)
            version=np.lib.format.read_magic(stream)
            if version==(1,0): shape,order,dtype=np.lib.format.read_array_header_1_0(stream)
            elif version==(2,0): shape,order,dtype=np.lib.format.read_array_header_2_0(stream)
            else: raise ValueError('Unsupported array version')
            if dtype.hasobject or math.prod(shape)*dtype.itemsize!=len(raw)-stream.tell():
                raise ValueError('Invalid numerical array')
            return np.load(io.BytesIO(raw),allow_pickle=False)
        raise ValueError('Invalid project value tag')
    return unpack(json.loads(data['state.json']))


def capture(app):
    from .surface_setup import encode
    g,l=app.geometry_workflow,app.load_workflow
    def mesh(m):
        return None if m is None else dict(vertices=m.vertices,faces=m.faces,source=m.source)
    result = None if app.result is None else dict(app.result)
    if result is not None:
        from solver_v1.model import TwoRowLJ
        from solver_v1.probability_pde_2d import Grid2D
        if isinstance(result.get('grid'), Grid2D):
            result['grid'] = asdict(result['grid'])
        if isinstance(result.get('model'), TwoRowLJ):
            model = result['model']
            result['model'] = dict(python_class=type(model).__module__+'.'+type(model).__name__,
                parameters=asdict(model.p), representation='metadata only; no executable object')
            if hasattr(model,'embedding'):
                result['model'].update(embedding=asdict(model.embedding), density=asdict(model.density_params))
    return dict(schema=SCHEMA, entries={k:v.get() for k,v in app.entries.items()},
        geometry=mesh(g.geometry), mesh=mesh(g.mesh), cylinder_dimensions=g.cylinder_dimensions,
        geometry_inputs={k:getattr(g,k).get() for k in ('radius','length','unit_scale','target')},
        loads=None if g.mesh is None else encode(g.mesh,l.loads,l.correction),
        load_editor=dict(tensor=l.tensor_text.get(),shear_mean=l.shear_mean.get(),shear_amplitude=l.shear_amplitude.get(),
                         region=l.region_code,selected=l.custom_indices),
        local_only=app.local_only.get(), energy_model=app.energy_model_code,
        quality=app.analysis_quality_code, probability_scale=app.probability_scale_code,
        time_basis=app.time_basis_code, calibration=asdict(app.time_calibration),
        frequency_values=dict(app._frequency_values), field=app.field.get(), views=dict(app._view_limits),
        result=result, last_config=None if app.last_config is None else asdict(app.last_config))


def restore(app, state):
    from .specimen_mesh import SurfaceMesh
    from .surface_setup import decode
    from .load_balance import correction_operator
    from .i18n import FIELD_TEXT_KEYS
    from .solver_adapter import UIAnalysisConfig
    from solver_v1.energy_model_registry import energy_model_metadata
    from solver_v1.physical_time import PhysicalTimeCalibration
    from solver_v1.kinetic_calibration_workflow import validate_for_energy_model
    if state.get('schema')!=SCHEMA: raise ValueError('Invalid project schema')
    def mesh(row):
        if row is None: return None
        faces=np.asarray(row['faces'])
        if faces.dtype.kind not in 'iu': raise ValueError('Invalid mesh indices')
        return SurfaceMesh(row['vertices'],faces,row['source'])
    geometry=mesh(state['geometry']); grid=mesh(state['mesh'])
    if geometry is None: raise ValueError('Missing geometry')
    loads,faces=([],None) if state['loads'] is None else decode(state['loads'],grid)
    correction=None if faces is None else correction_operator(grid,faces)
    model=state['energy_model']; energy_model_metadata(model)
    calibration=PhysicalTimeCalibration(**state['calibration']); calibration.validate()
    if state['time_basis']=='physical': validate_for_energy_model(calibration,model)
    elif state['time_basis']!='model': raise ValueError('Invalid time basis')
    if state['quality'] not in ('preview','resolved') or state['probability_scale'] not in ('local','specimen'):
        raise ValueError('Invalid analysis options')
    if state['field'] not in FIELD_TEXT_KEYS: raise ValueError('Invalid plot field')
    if set(state['entries'])!=set(app.entries) or not all(isinstance(v,str) for v in state['entries'].values()):
        raise ValueError('Unsupported setup fields')
    config=state['last_config']
    if config is not None:
        config=dict(config)
        if config['time_calibration'] is not None: config['time_calibration']=PhysicalTimeCalibration(**config['time_calibration'])
        config=UIAnalysisConfig(**config); config.validate()
    result=state['result']
    if result is not None:
        if not isinstance(result,dict) or not isinstance(result.get('model_time'),np.ndarray) or not len(result['model_time']):
            raise ValueError('Invalid result history')
        for key in ('strain','survival','local_initiation_probability'):
            if np.shape(result.get(key))!=result['model_time'].shape: raise ValueError('Result history length mismatch')
    editor=state['load_editor']
    from .tensor_load import compile_tensor_matrix
    compile_tensor_matrix(editor['tensor'])
    if editor['region'] not in app.load_workflow.REGIONS: raise ValueError('Invalid face selection')
    selected=np.asarray(editor['selected'])
    if selected.dtype.kind not in 'iu' or (len(selected) and (grid is None or selected.min()<0 or selected.max()>=len(grid.faces))):
        raise ValueError('Invalid selected faces')
    # Validated data is applied without running a solver or opening external CAD.
    g,l=app.geometry_workflow,app.load_workflow
    app.result=None; app.live_records.clear()
    g.geometry,g.mesh=geometry,grid; g.cylinder_dimensions=state['cylinder_dimensions']
    for k,v in state['geometry_inputs'].items(): getattr(g,k).set(v)
    for k,v in state['entries'].items(): app.entries[k].delete(0,'end'); app.entries[k].insert(0,v)
    l.refresh()
    l.loads=loads; l.correction=correction; l.applied=loads[-1] if loads else None
    l.tensor_text.set(editor['tensor']); l.shear_mean.set(editor['shear_mean']); l.shear_amplitude.set(editor['shear_amplitude'])
    l.region_code=editor['region']; l.custom_indices=selected
    app.local_only.set(state['local_only']); app.energy_model_code=model
    app.analysis_quality_code=state['quality']; app.probability_scale_code=state['probability_scale']
    app.time_calibration=calibration; app.time_basis_code=state['time_basis']; app._frequency_values=state['frequency_values']
    app.last_config=config; app.field.set(state['field']); app._view_limits=state['views']
    app.load_preset_code='custom'
    app._summary_kind='ready'; app._summary_payload={}
    app.ax.clear(); app._refresh_language(); g.draw()
    app.result=result
    if result is not None:
        app._show_summary('complete',result)
        app.notebook.select(app.post_tab)
    app.convergence_button.configure(state='normal' if result is not None and config is not None and result.get('analysis_quality')=='resolved' else 'disabled')
    app._refresh_specimen_labels(); app._plot()
    for view in l.maps: view.update()


class ProjectFiles:
    def __init__(self,app): self.app=app; self.path=None
    def save(self,as_new=False):
        from tkinter import filedialog,messagebox
        app=self.app
        if app.busy:
            messagebox.showinfo(app._tr('project.save'),app._tr('project.busy'),parent=app.root); return
        path=None if as_new else self.path
        if path is None: path=filedialog.asksaveasfilename(parent=app.root,defaultextension='.ftgsim',filetypes=[('Al Fatigue project','*.ftgsim')])
        if not path: return
        try: self.path=save_bundle(path,capture(app)); app._set_status('project.saved',path=str(self.path))
        except Exception as exc: messagebox.showerror(app._tr('project.save'),str(exc),parent=app.root)
    def open(self,path=None):
        from tkinter import filedialog,messagebox
        app=self.app
        if app.busy:
            messagebox.showinfo(app._tr('project.open'),app._tr('project.busy'),parent=app.root); return
        if path is None: path=filedialog.askopenfilename(parent=app.root,filetypes=[('Al Fatigue project','*.ftgsim')])
        if not path: return
        if not messagebox.askyesno(app._tr('project.open'),app._tr('project.replace'),parent=app.root): return
        try: restore(app,load_bundle(path)); self.path=Path(path); app._set_status('project.loaded',path=str(self.path))
        except Exception as exc: messagebox.showerror(app._tr('project.open'),str(exc),parent=app.root)
