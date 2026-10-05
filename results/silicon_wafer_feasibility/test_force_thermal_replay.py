"""Synthetic record integrity tests, not Si material/equilibrium validation."""
import json
from pathlib import Path
import numpy as np
import pytest
from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_sampling_score import nonlinear_correction_gradient
from solver_v1.silicon_thermal_research import harmonic_split_proposal
from results.silicon_wafer_feasibility.audit_force_thermal import replay

def dump(path,value):path.write_text(json.dumps(value,allow_nan=False),encoding='utf-8')

def fixture(tmp_path):
    folder=tmp_path/'chain_0';folder.mkdir()
    positions=np.array([[0.,0.,0.],[-5.,0.,0.],[5.,0.,0.]])
    basis=np.zeros((9,3));basis[:3]=np.eye(3)
    hessian=np.diag([2.5,4.,3.])
    source=dict(positions=positions,basis=basis,hessian=hessian,gradient=np.zeros(3),energy=np.array(0.),forces=np.zeros((3,3)))
    target=FixedGripTarget(positions,basis,hessian,np.zeros(3),300.,positions)
    protocol=dict(warmup=4,draws=32,seeds=[991],temperature_K=300.,step=1.1,steps=[2,4])
    np.savez_compressed(tmp_path/'common_domain.npz',direction=np.array([1.,0.,0.,0.,0.,0.,0.,0.,0.]))
    count=0;iteration=-1;last=None
    def evaluate(u):
        nonlocal count,last
        r=target.coordinates(u);inside=target.inside(r)
        if inside:
            x=r[0];v=np.array([.8,.6,0.]);q=v@x
            energy=.5*x@hessian@x+100*q**4
            forces=np.zeros_like(r);forces[0]=-hessian@x-400*q**3*v
            phi=target.correction(u,energy,0.)
            gradient=nonlinear_correction_gradient(target,u,forces)
            observation=np.array([energy,r[0,0],np.linalg.norm(r[0])])
        else:
            energy=np.nan;forces=np.full_like(r,np.nan);phi=np.inf
            gradient=np.full_like(u,np.nan);observation=np.full(3,np.nan)
        last=count
        np.savez_compressed(folder/f'evaluation_{count:04d}.npz',u=u,positions=r,energy=energy,forces=forces,
            correction=phi,gradient=gradient,inside=inside,observation=observation,
            context='sampling' if iteration>=0 else 'preflight',iteration=iteration)
        count+=1
        return phi,gradient,observation
    u=np.zeros(3);current=evaluate(u);current_index=0;rng=np.random.default_rng(991)
    records=[];samples=[];obs=[];phis=[];accepts=[];retained=[];errors=[]
    for iteration in range(36):
        p=rng.standard_normal(3);length=int(rng.integers(2,5));first=count;start=current_index
        h0=.5*(u@u+p@p)+current[0]
        proposal=harmonic_split_proposal(evaluate,u,p,1.1,length,current)
        uniform=None;error=None;accept=False
        if proposal is not None:
            error=float(.5*(proposal[0]@proposal[0]+proposal[1]@proposal[1])+proposal[2]-h0)
            uniform=float(rng.random());accept=bool(np.log(uniform)<min(0.,-error))
            if accept:u=proposal[0].copy();current=proposal[2:];current_index=last
        records.append(dict(iteration=iteration,steps=length,start_evaluation=start,first_evaluation=first,
            last_evaluation=last,retained_evaluation=current_index,delta_h=error,uniform=uniform,
            accepted=accept,domain_rejected=proposal is None))
        if iteration>=4:
            samples.append(u.copy());obs.append(current[2]);phis.append(current[0]);accepts.append(accept)
            retained.append(current_index);errors.append(np.inf if error is None else error)
    np.savez_compressed(folder/'chain.npz',samples=samples,observations=obs,potentials=phis,
        accepted_at_saved_step=accepts,retained_evaluation=retained,hamiltonian_errors=errors)
    dump(folder/'transactions.json',records)
    dump(folder/'summary.json',dict(complete=True,final_rng_state=rng.bit_generator.state,acceptance_fraction=float(np.mean(accepts))))
    assert any(accepts) and not all(accepts)
    return source,target,protocol,positions,folder

def test_independent_force_path_acceptance_and_rejected_repeats(tmp_path):
    source,target,protocol,common,folder=fixture(tmp_path)
    result,thermal,features=replay(tmp_path,tmp_path,0,source,target,protocol,common)
    assert result['maximum_score_error']<1e-11
    assert result['maximum_path_error']<1e-11
    assert result['distinct_retained_evaluations']<result['correlated_draws']
    assert thermal.shape==(32,3) and features.shape==(32,3)

def test_replay_rejects_corrupted_retention_history(tmp_path):
    source,target,protocol,common,folder=fixture(tmp_path)
    path=folder/'transactions.json';records=json.loads(path.read_text())
    records[7]['retained_evaluation']+=1;dump(path,records)
    with pytest.raises(ValueError,match='repeat lost'):
        replay(tmp_path,tmp_path,0,source,target,protocol,common)

def test_replay_rejects_corrupted_saved_score(tmp_path):
    source,target,protocol,common,folder=fixture(tmp_path)
    path=folder/'evaluation_0001.npz'
    with np.load(path) as z:values={k:z[k].copy() for k in z.files}
    values['gradient']+=1.
    np.savez_compressed(path,**values)
    with pytest.raises((ValueError,AssertionError)):
        replay(tmp_path,tmp_path,0,source,target,protocol,common)
