"""Meaningful restart tests using a coupled nonlinear target and rejections."""
from copy import deepcopy
import json
import numpy as np
import pytest
from solver_v1.silicon_chain_continuation import advance,save_checkpoint,load_checkpoint,parent_state

def initial():
    u=np.array([.5,-.8,.2]);phi,g,obs=evaluate(u)
    return dict(u=u,phi=phi,gradient=g,observation=obs,rng_state=np.random.default_rng(559).bit_generator.state,
                completed=0,retained_evaluation=0)
def evaluate(u):
    v=np.array([.8,.6,0.]);q=float(v@u)
    if np.linalg.norm(u)>4.:return np.inf,np.full(3,np.nan),np.full(2,np.nan)
    return 1.3*q**4,5.2*q**3*v,np.array([q,np.linalg.norm(u)])
def run(state,count):
    records=[]
    for _ in range(count):
        state,tx=advance(evaluate,state,step=.7,steps=(2,4))
        state['retained_evaluation']=state['completed'] if tx['accepted'] else state.get('retained_evaluation',0)
        records.append(tx)
    return state,records
def same(a,b):
    assert a.keys()==b.keys()
    for key in a:
        if isinstance(a[key],np.ndarray):np.testing.assert_array_equal(a[key],b[key])
        else:assert a[key]==b[key]
def test_uninterrupted_equals_saved_restored_continuation(tmp_path):
    start=initial();whole,records=run(deepcopy(start),48)
    partial,first=run(deepcopy(start),17)
    save_checkpoint(tmp_path/'state.npz',partial,first[-1]);restored,tx=load_checkpoint(tmp_path/'state.npz')
    same(partial,restored);assert tx==first[-1]
    final,last=run(restored,31);same(whole,final);assert records==first+last
    assert any(t['accepted'] for t in records) and not all(t['accepted'] for t in records)
def test_interruption_does_not_consume_rng_or_change_caller():
    start=initial();before=deepcopy(start);calls=0
    def broken(u):
        nonlocal calls
        calls+=1
        if calls==2:raise TimeoutError('budget')
        return evaluate(u)
    with pytest.raises(TimeoutError):advance(broken,start,step=.3,steps=(4,4))
    same(start,before)
    result,_=advance(evaluate,start,step=.3,steps=(4,4))
    expected,_=advance(evaluate,before,step=.3,steps=(4,4));same(result,expected)
def test_rejected_state_repeat_is_preserved():
    start=initial()
    def outside(u):return np.inf,np.full(3,np.nan),np.full(2,np.nan)
    result,tx=advance(outside,start,step=.3,steps=(2,4))
    assert tx['domain_rejected'] and not tx['accepted'] and tx['uniform'] is None
    for key in ('u','phi','gradient','observation'):np.testing.assert_array_equal(result[key],start[key])
def test_parent_uses_retained_endpoint_not_rejected_last_evaluation(tmp_path):
    state=initial();record=dict(u=state['u'],observation=state['observation'],correction=state['phi'],
        gradient=state['gradient'],inside=True)
    np.savez(tmp_path/'evaluation_0007.npz',**record)
    np.savez(tmp_path/'evaluation_0009.npz',**dict(record,u=state['u']+1))
    np.savez(tmp_path/'chain.npz',samples=[state['u']],observations=[state['observation']],
        potentials=[state['phi']],retained_evaluation=[7])
    (tmp_path/'summary.json').write_text(json.dumps(dict(complete=True,completed_proposals=1,final_rng_state=state['rng_state'])))
    (tmp_path/'transactions.json').write_text(json.dumps([dict(retained_evaluation=7,last_evaluation=9,accepted=False)]))
    restored=parent_state(tmp_path);np.testing.assert_array_equal(restored['u'],state['u'])
    assert restored['retained_evaluation']==7
    record['u']=state['u']+1;np.savez(tmp_path/'evaluation_0007.npz',**record)
    with pytest.raises(ValueError,match='bytes mismatch'):parent_state(tmp_path)
def test_invalid_checkpoint_rejected(tmp_path):
    state=initial();state['gradient'][0]=np.nan
    save_checkpoint(tmp_path/'bad.npz',state)
    with pytest.raises(ValueError,match='invalid committed'):load_checkpoint(tmp_path/'bad.npz')
