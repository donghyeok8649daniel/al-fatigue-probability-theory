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
