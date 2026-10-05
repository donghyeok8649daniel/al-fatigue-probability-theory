"""Independent continuation replay detects corrupted RNG, paths and repeats."""
import json
import numpy as np
import pytest
from solver_v1.silicon_chain_continuation import advance,parent_state,save_checkpoint
from solver_v1.silicon_sampling_score import nonlinear_correction_gradient
from results.silicon_wafer_feasibility.test_force_thermal_replay import fixture as parent_fixture
from results.silicon_wafer_feasibility.audit_thermal_continuation import replay_extension

def fixture(tmp_path):
    parent_root=tmp_path/'parent';parent_root.mkdir()
    source,target,_,common,parent=parent_fixture(parent_root)
    metadata=json.loads((parent/'summary.json').read_text(encoding='utf-8'))
    metadata['completed_proposals']=len(json.loads((parent/'transactions.json').read_text(encoding='utf-8')))
    (parent/'summary.json').write_text(json.dumps(metadata),encoding='utf-8')
    folder=tmp_path/'extension';folder.mkdir()
    domain=dict(center=common,direction=np.array([1.,0.,0.,0.,0.,0.,0.,0.,0.]))
    protocol=dict(additional_draws_per_chain=32,temperature_K=300.,step=.3,steps=[2,4])
    state=parent_state(parent);count=0;iteration=-1
    def evaluate(u):
        nonlocal count
        r=target.coordinates(u);inside=target.inside(r)
        if inside:
            x=r[0];v=np.array([.8,.6,0.]);q=v@x
            energy=.5*x@source['hessian']@x+100*q**4
            forces=np.zeros_like(r);forces[0]=-source['hessian']@x-400*q**3*v
            phi=target.correction(u,energy,0.);score=nonlinear_correction_gradient(target,u,forces)
            obs=np.array([energy,r[0,0],np.linalg.norm(r[0])])
        else:
            energy=np.nan;forces=np.full_like(r,np.nan);phi=np.inf;score=np.full_like(u,np.nan);obs=np.full(3,np.nan)
        np.savez_compressed(folder/f'evaluation_{count:04d}.npz',u=u,positions=r,energy=energy,forces=forces,
            correction=phi,gradient=score,inside=inside,observation=obs,
            context='restart_verification' if iteration<0 else 'sampling',iteration=iteration)
        count+=1
        return phi,score,obs
    evaluate(state['u']);state['retained_evaluation']=0;save_checkpoint(folder/'checkpoint_0000.npz',state)
    for iteration in range(32):
        old=state;first=count
        state,tx=advance(evaluate,old,step=protocol['step'],steps=protocol['steps'])
        state['retained_evaluation']=count-1 if tx['accepted'] else old['retained_evaluation']
        tx.update(iteration=iteration,global_proposal=56+iteration,start_evaluation=old['retained_evaluation'],
            first_evaluation=first,last_evaluation=count-1,retained_evaluation=state['retained_evaluation'])
        save_checkpoint(folder/f'checkpoint_{iteration+1:04d}.npz',state,tx)
    return folder,parent,source,target,domain,protocol

def test_independent_extension_replay(tmp_path):
    arguments=fixture(tmp_path)
    result,records=replay_extension(*arguments)
    assert len(records)==32 and result['all_rng_states_replayed']
    assert max(result['maximum_errors'].values())<1e-11

@pytest.mark.parametrize('kind',['uniform','retention','rng','path'])
def test_independent_extension_rejects_corruption(tmp_path,kind):
    arguments=fixture(tmp_path);folder=arguments[0]
    path=folder/('evaluation_0001.npz' if kind=='path' else 'checkpoint_0001.npz')
    with np.load(path,allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
    if kind=='path':data['u']+=.1
    elif kind=='rng':
        rng=json.loads(str(data['rng_json']));rng['state']['state']+=1;data['rng_json']=np.array(json.dumps(rng))
    else:
        tx=json.loads(str(data['transaction_json']))
        if kind=='uniform':tx['uniform']=.123456
        else:tx['retained_evaluation']+=1
        data['transaction_json']=np.array(json.dumps(tx))
    np.savez_compressed(path,**data)
    with pytest.raises((ValueError,AssertionError)):replay_extension(*arguments)
