"""Version 2 .ftgsim: portable setup and exact numerical results, no pickle."""
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import zipfile
import numpy as np

LIMIT = 16 * 1024 * 1024 * 1024
STREAM_CHUNK = 1024 * 1024
METADATA_LIMIT = 256 * 1024 * 1024
MANIFEST_LIMIT = 1024 * 1024
SCHEMA = 'aft.probability-project/2'


def save_bundle(path, state):
    path = Path(path)
    payloads = {}
    def pack(x):
        if isinstance(x, np.ndarray):
            if x.dtype.hasobject: raise ValueError('Object arrays cannot be saved')
            name = f'arrays/{len(payloads)}.npy'
            payloads[name] = x
            return {'array': name}
        if isinstance(x, np.generic): return pack(x.item())
        if isinstance(x, dict): return {'dict': {str(k): pack(v) for k,v in x.items()}}
        if isinstance(x, (tuple, list)): return {'tuple' if isinstance(x, tuple) else 'list': [pack(v) for v in x]}
        if isinstance(x, float) and not math.isfinite(x): return {'float': repr(x)}
        if x is None or type(x) in (bool, int, float, str): return x
        raise ValueError(f'Unsupported saved value: {type(x).__name__}')
    tree = pack(state)
    payloads['state.json'] = json.dumps(tree, ensure_ascii=False, allow_nan=False).encode('utf-8')
    if len(payloads['state.json']) > METADATA_LIMIT:
        raise ValueError('Project metadata exceeds size limit')
    if len(payloads)+1>2048 or sum(v.nbytes if isinstance(v,np.ndarray) else len(v)
                                  for v in payloads.values())>LIMIT:
        raise ValueError('Project exceeds size limit')
    total = 0
    class CheckedWriter:
        def __init__(self, stream): self.stream=stream; self.digest=hashlib.sha256()
        def write(self, data):
            nonlocal total
            total += len(data)
            if total > LIMIT: raise ValueError('Project exceeds size limit')
            self.digest.update(data)
            return self.stream.write(data)
        def flush(self): self.stream.flush()
    handle, tmp = tempfile.mkstemp(prefix=path.stem+'-', suffix='.tmp', dir=path.parent)
    os.close(handle)
    try:
        with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED) as z:
            hashes = {}
            for name,data in payloads.items():
                with z.open(name,'w',force_zip64=True) as member:
                    writer = CheckedWriter(member)
                    if isinstance(data,np.ndarray): np.save(writer,data,allow_pickle=False)
                    else: writer.write(data)
                    hashes[name] = writer.digest.hexdigest()
            manifest = dict(format='ftgsim', schema_version='2.0.0', application_schema=SCHEMA,
                            generator='Al Fatigue probability PDE desktop', checksums_sha256=hashes)
            manifest_bytes = json.dumps(manifest).encode()
            if len(manifest_bytes) > MANIFEST_LIMIT:
                raise ValueError('Project manifest exceeds size limit')
            with z.open('ftgsim-manifest.json','w') as member:
                CheckedWriter(member).write(manifest_bytes)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
    return path


def load_bundle(path):
    with zipfile.ZipFile(path) as z:
        infos=z.infolist(); names=[i.filename for i in infos]
        if len(names)>2048 or len(set(names))!=len(names) or sum(i.file_size for i in infos)>LIMIT:
            raise ValueError('Invalid project size or duplicate members')
        if z.getinfo('ftgsim-manifest.json').file_size > MANIFEST_LIMIT:
            raise ValueError('Invalid project manifest size')
        manifest=json.loads(z.read('ftgsim-manifest.json'))
        if manifest.get('application_schema')!=SCHEMA:
            raise ValueError('Legacy or unsupported .ftgsim model. Open it with its original AlFatigue executable; it is not a probability-PDE project.')
        hashes=manifest['checksums_sha256']
        if set(names)!=set(hashes)|{'ftgsim-manifest.json'}: raise ValueError('Project member mismatch')
        data={}
        for name,digest in hashes.items():
            if name == 'state.json':
                if z.getinfo(name).file_size > METADATA_LIMIT:
                    raise ValueError('Invalid project metadata size')
                raw=z.read(name)
                if hashlib.sha256(raw).hexdigest()!=digest: raise ValueError('Project checksum mismatch')
                data[name]=raw
            elif name.startswith('arrays/') and name.endswith('.npy'):
                data[name] = _read_array(z, name, digest)
            else: raise ValueError('Invalid project member')
    def unpack(x):
        if not isinstance(x,dict): return x
        if len(x)!=1: raise ValueError('Invalid project value')
        kind,value=next(iter(x.items()))
        if kind=='dict': return {k:unpack(v) for k,v in value.items()}
        if kind in ('list','tuple'):
            v=[unpack(y) for y in value]; return tuple(v) if kind=='tuple' else v
        if kind=='float' and value in ('nan','inf','-inf'): return float(value)
        if kind=='array':
            array=data[value]
            if not isinstance(array,np.ndarray): raise ValueError('Invalid numerical array')
            return array
        raise ValueError('Invalid project value tag')
    return unpack(json.loads(data['state.json']))


def _read_array(archive, name, expected_digest):
    """Validate the NPY header before allocation, then read into the final array.

    The compressed bytes, a complete uncompressed NPY buffer, and a second copy
    of every array no longer coexist. Object/pickle payloads remain forbidden.
    """
    digest = hashlib.sha256()
    with archive.open(name) as member:
        class HeaderReader:
            def read(self, size=-1):
                raw=member.read(size); digest.update(raw); return raw
        reader = HeaderReader()
        version=np.lib.format.read_magic(reader)
        if version==(1,0): shape,order,dtype=np.lib.format.read_array_header_1_0(reader)
        elif version==(2,0): shape,order,dtype=np.lib.format.read_array_header_2_0(reader)
        else: raise ValueError('Unsupported array version')
        if dtype.hasobject or math.prod(shape)*dtype.itemsize != archive.getinfo(name).file_size-member.tell():
            raise ValueError('Invalid numerical array')
        array=np.empty(shape,dtype=dtype,order='F' if order else 'C')
        buffer=memoryview(array.ravel(order='K').view(np.uint8))
        for start in range(0,len(buffer),STREAM_CHUNK):
            chunk=buffer[start:start+STREAM_CHUNK]
            if member.readinto(chunk) != len(chunk): raise ValueError('Invalid numerical array')
            digest.update(chunk)
    if digest.hexdigest()!=expected_digest: raise ValueError('Project checksum mismatch')
    return array


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
        local_only=app.local_only.get(), axial_count=app.axial_count.get(),
        solid_poisson=app.solid_poisson.get(), solid_target_mm=app.solid_target.get(),
        initialization=app.initialization.get(),
        material=app._material_draft(),
        energy_model=app.energy_model_code,
        correlation_volume_mm3=app.correlation_volume.get(),
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
    from .materials import ALUMINUM, default_material_draft, validate_material_draft
    from solver_v1.energy_model_registry import energy_model_metadata
    from solver_v1.physical_time import PhysicalTimeCalibration
    from solver_v1.kinetic_calibration_workflow import validate_for_energy_model
    if state.get('schema')!=SCHEMA: raise ValueError('Invalid project schema')
    material = validate_material_draft(state.get('material', default_material_draft()))
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
    if state.get('initialization', 'loaded_gibbs') not in ('loaded_gibbs', 'zero_load_gibbs'):
        raise ValueError('Invalid initial ensemble')
    if set(state['entries'])!=set(app.entries) or not all(isinstance(v,str) for v in state['entries'].values()):
        raise ValueError('Unsupported setup fields')
    config=state['last_config']
    if config is not None:
        config=dict(config)
        if config['time_calibration'] is not None: config['time_calibration']=PhysicalTimeCalibration(**config['time_calibration'])
        config=UIAnalysisConfig(**config); config.validate()
    result=state['result']
    if result is not None:
        if not isinstance(result,dict):
            raise ValueError('Invalid result history')
        # Persisted historical results cannot be relabelled by a new draft.
        # Only the Al backend has ever produced supported app results.
        if result.get('material_id', ALUMINUM) != ALUMINUM or (config is not None and config.material_id != ALUMINUM):
            raise ValueError('Unsupported result material')
        if not isinstance(result.get('model_time'),np.ndarray) or not len(result['model_time']):
            raise ValueError('Invalid result history')
        for key in ('strain','survival','local_initiation_probability'):
            if np.shape(result.get(key))!=result['model_time'].shape: raise ValueError('Result history length mismatch')
        if 'axial_specimen' in result:
            from .axial_specimen import validate_spatial_result
            validate_spatial_result(result['axial_specimen'])
        if 'solid_specimen' in result:
            from .solid_mechanics import validate_solid_result
            validate_solid_result(result['solid_specimen'])
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
    app.material_id=material['material_id']; app.doping_enabled.set(material['doping_enabled'])
    app.dopant_species_code=material['dopant_species']
    app.dopant_concentration.set(material['dopant_concentration_cm3'])
    app.axial_count.set(state.get('axial_count', '16'))
    app.solid_poisson.set(state.get('solid_poisson', '0.33'))
    app.solid_target.set(state.get('solid_target_mm', state['geometry_inputs'].get('target', '3')))
    app.initialization.set(state.get('initialization', 'loaded_gibbs'))
    app.correlation_volume.set(state.get('correlation_volume_mm3', ''))
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
    app._refresh_material_ui()
    # Derived specimen values may be absent/stale in a saved result. Rebuild
    # them from its original local history and restored areas without a solve.
    app._update_specimen_probability()
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
