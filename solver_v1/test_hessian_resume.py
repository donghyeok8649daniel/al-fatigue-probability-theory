import json
import numpy as np
import pytest
from solver_v1.hessian_resume import load_prefix,check_replayed_row,sha,IDENTITY

@pytest.fixture
def prefix(tmp_path):
    dim=4; basis=np.eye(dim);gradient=np.arange(dim)*1e-5
    h=np.array([[4.,1.,2.,.5],[1.,3.,.2,.3],[2.,.2,5.,.1],[.5,.3,.1,2.]])
    p={k:'same' for k in IDENTITY};p['dimension']=dim
    np.savez(tmp_path/'partial.npz',hessian_rows=h[:2],basis=basis,gradient=gradient)
    (tmp_path/'protocol.json').write_text(json.dumps(p))
    (tmp_path/'summary.json').write_text(json.dumps(dict(complete=False,completed_rows=2,dimension=dim)))
    return tmp_path,p,basis,gradient,h

def test_valid_prefix_preserves_off_block_columns(prefix):
    folder,p,b,g,h=prefix
    rows,record=load_prefix(folder,p,b,g,sha(folder/'partial.npz'))
    np.testing.assert_array_equal(rows,h[:2])
    assert record['reused_rows']==2
    check_replayed_row(rows[0],h[0])

@pytest.mark.parametrize('field',IDENTITY)
def test_context_mismatch_rejected(prefix,field):
    folder,p,b,g,h=prefix;p={**p,field:'changed'}
    with pytest.raises(ValueError,match='context changed'):
        load_prefix(folder,p,b,g,sha(folder/'partial.npz'))

@pytest.mark.parametrize('kind',['hash','gradient','basis','count','nan','asymmetry'])
def test_corruption_rejected(prefix,kind):
    folder,p,b,g,h=prefix;expected=sha(folder/'partial.npz')
    if kind=='hash':expected='0'*64
    if kind=='gradient':g=g+.01
    if kind=='basis':b=b[:,::-1]
    if kind=='count':
        (folder/'summary.json').write_text(json.dumps(dict(complete=False,completed_rows=3,dimension=4)))
    if kind in ('nan','asymmetry'):
        rows=h[:2].copy();rows[0,1]=np.nan if kind=='nan' else 2.
        np.savez(folder/'partial.npz',hessian_rows=rows,basis=b,gradient=g);expected=sha(folder/'partial.npz')
    with pytest.raises(ValueError):load_prefix(folder,p,b,g,expected)

def test_wrong_recomputed_row_rejected(prefix):
    *_,h=prefix
    with pytest.raises(ValueError,match='does not reproduce'):check_replayed_row(h[0],h[0]+.01)
    with pytest.raises(ValueError,match='invalid'):check_replayed_row(h[0],h[0]*np.nan)
