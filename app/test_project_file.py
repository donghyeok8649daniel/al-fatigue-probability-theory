"""Project round-trip with a real local PDE history, no solver on restore."""
import io
import json
import zipfile
import numpy as np
import pytest
from app.project_file import save_bundle, load_bundle, capture, restore
from app.solver_adapter import UIAnalysisConfig, run_ui_analysis


def exact(a,b):
    if isinstance(a,np.ndarray):
        assert isinstance(b,np.ndarray) and a.dtype==b.dtype and a.shape==b.shape
        assert a.tobytes()==b.tobytes()
    elif isinstance(a,dict):
        assert set(a)==set(b)
        for k in a: exact(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert type(a)==type(b) and len(a)==len(b)
        for x,y in zip(a,b): exact(x,y)
    elif isinstance(a,float) and np.isnan(a): assert np.isnan(b)
    else: assert a==b


def test_bundle_exact_and_corruption(tmp_path):
    p=tmp_path/'한글 project.ftgsim'
    state={'v':np.array([0.,1e-20,np.nan,np.inf,-0.]),'t':(1,None,False),'grid':np.ones((3,2),dtype=np.int64)}
    save_bundle(p,state); exact(state,load_bundle(p))
    with zipfile.ZipFile(p) as z: files={n:z.read(n) for n in z.namelist()}
    files['state.json']+=b' '
    with zipfile.ZipFile(p,'w') as z:
        for n,v in files.items(): z.writestr(n,v)
    with pytest.raises(ValueError,match='checksum'): load_bundle(p)
    with pytest.raises(ValueError): save_bundle(p,{'bad':np.array([object()])})


def test_legacy_project_is_not_reinterpreted(tmp_path):
    p=tmp_path/'legacy.ftgsim'
    with zipfile.ZipFile(p,'w') as z: z.writestr('ftgsim-manifest.json',json.dumps({'schema_version':'1.0.0'}))
    with pytest.raises(ValueError,match='Legacy'): load_bundle(p)


def test_streamed_arrays_preserve_layout_endianness_and_atomic_size_failure(tmp_path, monkeypatch):
    from app import project_file
    p = tmp_path/'streamed.ftgsim'
    base = np.arange(2000, dtype='>f8').reshape(40,50)
    state = {'strided': base[::3, ::2], 'fortran': np.asfortranarray(base),
             'scalar': np.array(-0.), 'empty': np.empty((0,3)),
             'dates': np.array(['2026-09-16', 'NaT'], dtype='datetime64[D]'),
             'text': np.array(['형상', 'mesh']),
             'structured': np.array([(1, 2.5)], dtype=[('index', '>i4'), ('stress', '<f8')])}
    save_bundle(p, state)
    original_read = zipfile.ZipFile.read
    def metadata_only(self, name, *args, **kwargs):
        assert not str(name).endswith('.npy'), 'Array payload must be streamed'
        return original_read(self, name, *args, **kwargs)
    monkeypatch.setattr(zipfile.ZipFile, 'read', metadata_only)
    exact(state, load_bundle(p))
    original = p.read_bytes()
    # Includes streaming headers/manifest, beyond the array-only preflight.
    monkeypatch.setattr(project_file, 'LIMIT', 1000)
    with pytest.raises(ValueError, match='size limit'):
        save_bundle(p, {'a': np.arange(100, dtype=np.float64)})
    assert p.read_bytes() == original
    monkeypatch.setattr(project_file, 'METADATA_LIMIT', 100)
    with pytest.raises(ValueError, match='metadata exceeds size limit'):
        save_bundle(p, {'label': 'x'*100})
    assert p.read_bytes() == original
