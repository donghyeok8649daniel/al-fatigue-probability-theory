import numpy as np
import pytest
from .collective_generator_bridge import (plane_zero_sum_basis,covariance_mobility,
    harmonic_generator,energy_rescale,projected_susceptibility)
from .run_weak_replica_v31 import definition
from .run_finite_q_compatibility_v31 import points


def test_weak_campaign_null_and_driven_work_share_csv_schema(tmp_path):
    import csv
    from .run_weak_replica_v31 import internal_work_diagnostics
    from .run_low_stress_cyclic_diagnostic import write_csv
    thermo=np.zeros((3,4));thermo[:,3]=[10.,11.,12.]
    null=internal_work_diagnostics(thermo)
    driven=internal_work_diagnostics(thermo,np.array([0.,.9,1.8]))
    assert list(null)==list(driven) and all(v is None for v in null.values())
    np.testing.assert_allclose(driven['max_work_residual_eV'],.2)
    np.testing.assert_allclose(driven['rms_work_residual_eV'],np.sqrt(.05/3))
    path=tmp_path/'checks.csv';write_csv(path,[null,driven])
    with path.open(newline='') as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==2 and rows[0]['work_eV']=='' and float(rows[1]['work_eV'])==1.8
    with pytest.raises(ValueError):internal_work_diagnostics(thermo,[0.,1.])


def test_completed_weak_campaign_export_end_to_end(tmp_path,monkeypatch,capsys):
    import csv
    import json
    from . import run_weak_replica_v31 as runner
    study=tmp_path/'study';study.mkdir();out=tmp_path/'out'
    protocol=runner.definition(dict(force_scale_eV_A=[3.,1.],potential_sha256='test-only'))
    (study/'protocol.json').write_text(json.dumps(protocol))
    cases=[]
    for seed in protocol['seeds']:
        for axis,sign in [(None,0),(0,-1),(0,1),(1,-1),(1,1)]:
            case=dict(name=f'{seed}_{axis}_{sign}',seed=seed,axis=axis,sign=sign,restart_sha256='bound')
            drive=None
            if axis is not None:
                case['force_eV_A']=sign*protocol['force_scale_eV_A'][axis]
                drive=dict(axis=axis,force_eV_A=case['force_eV_A'],frequency_per_ps=.2)
            cases.append(case);record=study/case['name'];record.mkdir()
            meta=dict(completed=True,restart_sha256='bound',ensemble='nve',conjugate_drive=drive,
                      **{k:protocol[k] for k in ('dt_ps','frame_ps','duration_ps')})
            (record/'summary.json').write_text(json.dumps(meta))
            t=np.array([0.,100.,500.,1100.,1500.,1950.])*1e-12
            thermo=np.zeros((len(t),4));thermo[:,0]=300.
            np.savez(record/'plane_coordinates.npz',time_seconds=t,
                     coordinates_m=np.zeros((len(t),3,3)),thermo=thermo,internal_step_work_eV=np.zeros(len(t)))
    (study/'cases.json').write_text(json.dumps(cases))
    # Numerical estimators are tested elsewhere; this fixture tests report plumbing.
    monkeypatch.setattr(runner,'harmonic_response',lambda *a,**k:dict(real_A2_eV=.01,imag_A2_eV=-.001))
    monkeypatch.setattr(runner,'null_lockin_windows',lambda *a,**k:[dict(real_A=.0001,imag_A=.0001)])
    monkeypatch.setattr(runner,'multitaper_spectrum',lambda *a,**k:None)
    monkeypatch.setattr(runner,'band_integral_proxy',lambda *a,**k:dict(integral_proxy_m2_seconds=np.eye(3)*1e-40))
    runner.analyze(study,out)
    result=json.loads((out/'summary.json').read_bytes())
    assert result['completed'] and not result['production_clock_calibrated']
    assert len(json.loads(capsys.readouterr().out))==2
    with (out/'checks.csv').open(newline='') as stream:checks=list(csv.DictReader(stream))
    assert len(checks)==10 and sum(r['work_eV']=='' for r in checks)==2


def test_existing_shape_signed_screening_is_not_logarithmized():
    from .run_joint_shape_v31 import shape_coordinates,decode_shape
    shape=np.array([4.,6.,7.,100.,8.,-.5])
    for vary in (False,True):
        x=shape_coordinates(shape,vary)
        np.testing.assert_allclose(decode_shape(x,shape[-1],vary),shape)
    assert shape_coordinates(shape,True)[-1]==-.5
    with pytest.raises(ValueError):shape_coordinates([4.,6.,7.,100.,8.,1.01],True)
    # Reject a nondecaying isolated-neighbor limit, not a numerical penalty repair.
    with pytest.raises(ValueError):decode_shape(np.r_[np.log([8.,6.,7.,100.,2.]),-1.],0,True)


@pytest.mark.parametrize('planes',[3,6,12])
def test_zero_sum_projection_keeps_all_nontranslation_modes(planes):
    B=plane_zero_sum_basis(planes)
    np.testing.assert_allclose(B.T@B,np.eye(3*(planes-1)),atol=1e-15)
    P=np.kron(np.eye(planes)-np.ones((planes,planes))/planes,np.eye(3))
    np.testing.assert_allclose(B@B.T,P,atol=1e-15)


def test_coupled_covariance_recovers_full_mobility_not_scalar_proxy():
    H=np.array([[3.,.8],[.8,2.]]);M=np.array([[2.,.7],[.7,1.]])
    kBT=.03;C=kBT*np.linalg.inv(H);K=np.linalg.solve(M@H,C)
    gotH,gotM=covariance_mobility(C,K,kBT)
    np.testing.assert_allclose(gotH,H,rtol=1e-13)
    np.testing.assert_allclose(gotM,M,rtol=1e-13)
    scalar=C[0,0]**2/(kBT*K[0,0])
    assert abs(scalar-M[0,0])>.1


def test_per_site_energy_requires_per_site_thermal_term():
    H=np.array([[3.,.2],[.2,1.]]);M=np.diag([1.,.05]);T=.02
    new=energy_rescale(H,M,T,144)
    before=harmonic_generator(H,M,T);after=harmonic_generator(*new)
    for a,b in zip(before,after):np.testing.assert_allclose(a,b)
    wrong=harmonic_generator(new[0],new[1],T)
    np.testing.assert_allclose(wrong[1],144*before[1])


def test_si_nonsymmetric_covariance_not_hidden_by_absolute_tolerance():
    C=1e-24*np.array([[2.,.3],[0.,1.]])
    with pytest.raises(ValueError):covariance_mobility(C,np.eye(2)*1e-36,4e-21)
    with pytest.raises(ValueError):covariance_mobility(np.eye(2)*1e-24,np.diag([1e-36,1e-54]),4e-21)


def test_replica_protocol_is_weak_bound_and_not_clock_certification():
    d=definition(dict(force_scale_eV_A=[3.,1.],potential_sha256='bound'))
    assert len(set(d['seeds']))==2
    assert d['duration_ps']==2000 and d['force_fraction']==1
    assert d['dt_ps']==.00125 and not d['production_clock_calibrated']
    pp=points();fit={q for role,q in pp if role=='fit'};excluded={q for role,q in pp if role=='excluded'}
    assert len(fit)==4 and len(excluded)==4 and not fit&excluded


def test_eliminated_memory_has_distinct_low_and_high_frequency_mobility():
    H=np.array([[3.,.8],[.8,2.]])
    M=np.array([[2.,.7],[.7,1.]])
    v=np.array([1.,0.]);thermal=.03
    C=thermal*np.linalg.inv(H);K=np.linalg.solve(M@H,C)
    static=float(v@C@v/thermal)
    low=projected_susceptibility(H,M,v,1e-6)
    np.testing.assert_allclose(low.real,static,rtol=1e-11)
    np.testing.assert_allclose(-low.imag/1e-6,v@K@v/thermal,rtol=1e-11)
    scalar=static**2/(v@K@v/thermal)
    assert scalar < v@M@v
    high=1j*1e8*projected_susceptibility(H,M,v,1e8)
    np.testing.assert_allclose(high.real,v@M@v,rtol=1e-11)
    # Matching one scalar constant cannot match both asymptotic coefficients.
    scalar_response=1/(1/static+1j*1./scalar)
    assert abs(scalar_response-projected_susceptibility(H,M,v,1.))>1e-3


def test_reversible_overdamped_storage_does_not_exceed_static_compliance():
    H=np.array([[3.,.8],[.8,2.]])
    M=np.array([[2.,.7],[.7,1.]])
    v=np.array([1.,-.3])
    static=float(v@np.linalg.solve(H,v))
    chi=np.array([projected_susceptibility(H,M,v,w) for w in np.logspace(-7,7,81)])
    assert np.all(chi.real<=static*(1+1e-14))
    assert np.all(chi.real>0) and np.all(chi.imag<0)
    # Consistent SI and Angstrom/eV units give the same measured compliance.
    L=1e-10;E=1.602176634e-19;t=1e-12
    physical=projected_susceptibility(H*E/L**2,M*L**2/(E*t),v,2/t)
    np.testing.assert_allclose(physical*16.02176634,projected_susceptibility(H,M,v,2),rtol=1e-13)


def test_completion_workflow_never_analyzes_partial_campaign(tmp_path,monkeypatch):
    import json
    from . import finish_weak_replica_v31 as workflow
    study=tmp_path/'study';study.mkdir()
    (study/'cases.json').write_text(json.dumps([dict(name='unfinished',seed=1)]))
    ticks=iter([0.,4000.,4001.,4002.])
    monkeypatch.setattr(workflow.time,'monotonic',lambda:next(ticks))
    monkeypatch.setattr(workflow,'analyze',lambda *args:pytest.fail('partial analysis forbidden'))
    with pytest.raises(TimeoutError):
        workflow.finish(study,tmp_path/'out',tmp_path/'bridges',tmp_path/'watch',1.)
    status=json.loads((tmp_path/'watch/status.json').read_text())
    assert not status['completed'] and not status['production_clock_calibrated']


def test_completion_analysis_does_not_approve_calibration(tmp_path,monkeypatch):
    import json
    from . import finish_weak_replica_v31 as workflow
    study=tmp_path/'study';study.mkdir();record=study/'seed1_null';record.mkdir()
    (study/'cases.json').write_text(json.dumps([dict(name='seed1_null',seed=1)]))
    (record/'summary.json').write_text('{}')
    calls=[]
    monkeypatch.setattr(workflow,'analyze',lambda *args:calls.append('analyze'))
    monkeypatch.setattr(workflow,'bridge',lambda *args:calls.append('bridge'))
    workflow.finish(study,tmp_path/'out',tmp_path/'bridges',tmp_path/'watch')
    assert calls==['analyze','bridge']
    status=json.loads((tmp_path/'watch/status.json').read_text())
    assert status['completed'] and not status['production_clock_calibrated']
    assert not status['material_accepted']
