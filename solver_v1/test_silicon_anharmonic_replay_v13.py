"""Manufactured cubic/quartic energy and deliberate corruption checks.

These are implementation tests, not Si calculations or finite-T validation.
"""
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.constants import Boltzmann, electron_volt
from results.silicon_wafer_feasibility.replay_anharmonic_v13 import audit, sha


def dump(path,obj):
    path.write_text(json.dumps(obj),encoding='utf-8')


@pytest.fixture
def data(tmp_path):
    root=tmp_path/'synthetic';root.mkdir();probes=root/'probes';probes.mkdir()
    b=np.eye(9)[:,:6];g=np.full(6,.003);h=np.diag(np.arange(2.,8.))
    r=np.arange(9.,dtype=float).reshape(3,3);z=np.full(3,14);cell=np.eye(3)*20
    kb=Boltzmann/electron_volt; kt=300*kb; c=.04; d=.01
    points=[];sources=[]
    for name,e0 in [('loading8',-12.),('return8',-13.)]:
        folder=root/('dense_'+name);folder.mkdir()
        base_f=(-b@g).reshape(3,3)
        np.savez(folder/'raw_hessian.npz',hessian=h,basis=b,gradient=g,positions=r,
            numbers=z,cell=cell,energy=e0,forces=base_f)
        np.savez(folder/'spectrum.npz',eigenvalues=np.diag(h),eigenvectors=np.eye(6))
        sources.append(dict(state=name,raw_sha256=sha(folder/'raw_hessian.npz'),spectrum_sha256=sha(folder/'spectrum.npz')))
        energies=[];forces=[];evaluated=[];directions=[]
        for mode in range(3):
            lam=h[mode,mode];v=np.eye(6)[:,mode];directions.append((b@v).reshape(3,3))
            amplitudes=[('max_atom_'+str(a),a) for a in (.01,.05,.1,.2)]
            amplitudes += [('quadratic_sigma_'+str(a),a*np.sqrt(kt/lam)) for a in (1.,2.)]
            for label,q in amplitudes:
                for sign in (-1,1):
                    x=sign*q;displacement=x*v
                    # Independent analytic polynomial evaluated in all six coordinates.
                    energy=e0+g@displacement+.5*displacement@h@displacement+c*np.sum(displacement**3)+d*np.sum(displacement**4)
                    derivative=g+h@displacement+3*c*displacement**2+4*d*displacement**3
                    forces.append((-b@derivative).reshape(3,3));energies.append(energy);evaluated.append((mode,x))
                    points.append(dict(state=name,mode=mode,label=label,sign=sign,coordinate_A=x,
                        maximum_atom_displacement_A=abs(x),curvature_eV_A2=lam,
                        baseline_directional_gradient_eV_A=.003,actual_delta_energy_eV=energy-e0,
                        linear_quadratic_delta_energy_eV=.003*x+.5*lam*x*x,
                        nonlinear_energy_residual_eV=c*x**3+d*x**4,
                        nonlinear_energy_residual_over_kBT300=(c*x**3+d*x**4)/kt,
                        actual_directional_gradient_eV_A=derivative[mode],
                        nonlinear_full_gradient_norm_eV_A=abs(3*c*x*x+4*d*x**3),
                        nonlinear_directional_gradient_eV_A=3*c*x*x+4*d*x**3,
                        symmetric_nonlinear_energy_eV=d*q**4,
                        antisymmetric_nonlinear_energy_eV=c*q**3))
        np.savez(probes/(name+'_raw.npz'),positions=r,numbers=z,cell=cell,basis=b,
            directions=np.asarray(directions),evaluated=np.asarray(evaluated),
            energies=np.asarray(energies),forces=np.asarray(forces),baseline_energy=e0,baseline_forces=base_f)
    dump(probes/'points.json',points)
    dump(probes/'protocol.json',dict(modes=[0,1,2],fixed_max_atom_displacements_A=[.01,.05,.1,.2],
        quadratic_coordinate_amplitudes_at_300K=[1.,2.],maximum_atom_bound_A=.3,sources=sources))
    dump(probes/'summary.json',dict(complete=True,first_crack_certified=False,physical_clock=None,
        actual_points=72,new_model_calls=74,completed_source_replays=[dict(state='loading8'),dict(state='return8')],skipped=[]))
    return root,probes


def test_cubic_quartic_raw_replay(data):
    result=audit(*data)
    assert result['complete'] and result['raw_replay_passed']
    assert result['verified_points']==72 and result['complete_signed_pairs']==36
    assert result['maximum_scalar_replay_error']<1e-10
    assert result['new_model_calls']==0


@pytest.mark.parametrize('field',[
    'actual_delta_energy_eV','coordinate_A','nonlinear_full_gradient_norm_eV_A',
    'symmetric_nonlinear_energy_eV','antisymmetric_nonlinear_energy_eV'])
def test_rejects_corrupted_observations(data,field):
    _,probes=data;path=probes/'points.json';rows=json.loads(path.read_text());rows[0][field]+=.1;dump(path,rows)
    with pytest.raises(ValueError,match='mismatch'):audit(*data)


def test_rejects_nonfinite_report(data):
    _,probes=data;path=probes/'points.json';rows=json.loads(path.read_text());rows[0]['actual_delta_energy_eV']=float('nan');dump(path,rows)
    with pytest.raises(ValueError,match='nonfinite'):audit(*data)


def test_rejects_partial_data_marked_complete(data):
    _,probes=data;path=probes/'points.json';rows=json.loads(path.read_text());dump(path,rows[:-2])
    rawpath=probes/'return8_raw.npz'
    with np.load(rawpath) as d:raw={k:d[k].copy() for k in d.files}
    for key in ('energies','forces','evaluated'):raw[key]=raw[key][:-2]
    np.savez(rawpath,**raw)
    path=probes/'summary.json';summary=json.loads(path.read_text());summary.update(actual_points=70,new_model_calls=72);dump(path,summary)
    with pytest.raises(ValueError,match='incomplete schedule'):audit(*data)
    summary['complete']=False;dump(path,summary)
    result=audit(*data)
    assert result['raw_replay_passed'] and not result['complete'] and result['verified_points']==70


def test_rejects_changed_source(data):
    _,probes=data;path=probes/'protocol.json';protocol=json.loads(path.read_text());protocol['sources'][0]['raw_sha256']='bad';dump(path,protocol)
    with pytest.raises(ValueError,match='source hash mismatch'):audit(*data)


def test_rejects_call_double_count(data):
    _,probes=data;path=probes/'summary.json';summary=json.loads(path.read_text());summary['new_model_calls']=75;dump(path,summary)
    with pytest.raises(ValueError,match='call accounting'):audit(*data)
