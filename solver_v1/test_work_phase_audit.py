import numpy as np
import pytest
from solver_v1.work_phase_audit import work_phase, duration_sensitivity


@pytest.mark.parametrize('force', [2., -2.])
@pytest.mark.parametrize('loss', [.03, -.02, 0.])
def test_signed_phase_and_endpoint(force, loss):
    t = np.linspace(1.25, 21.25, 40001)
    omega = 2*np.pi*.2
    storage, drift = .4, .001
    q = .8+drift*t+force*(storage*np.cos(omega*t)+loss*np.sin(omega*t))
    # Exact primitive of F dq, including nonperiodic reversible/thermal drift.
    w = (force*drift*np.sin(omega*t)/omega
         +force**2*storage*np.cos(2*omega*t)/4
         +force**2*loss*(omega*t/2+np.sin(2*omega*t)/4))
    out = work_phase(t, q, w, force, .2)
    expected = loss+2*drift*(-t[-1]*np.cos(omega*t[-1])
                           +t[0]*np.cos(omega*t[0]))/(omega*force*(t[-1]-t[0]))
    assert out['loss_native_A2_eV'] == pytest.approx(expected, abs=1e-12)
    assert out['loss_coordinate_A2_eV'] == pytest.approx(expected, abs=1e-9)
    assert not out['production_clock_calibrated']


def test_endpoint_is_not_dissipation():
    t = np.linspace(0, 20, 20001)
    q = .003*t
    w = 2*.003*np.sin(2*np.pi*.2*t)/(2*np.pi*.2)
    out = work_phase(t, q, w, 2, .2)
    assert abs(out['loss_uncorrected_work_A2_eV']) < 1e-15
    assert out['loss_native_A2_eV'] < 0
    # A finite thermal trend may exchange energy with the forcing; no abs repair.
    assert out['endpoint_work_eV'] == pytest.approx(.12)


def test_validation_and_planning():
    t = np.linspace(0, 20, 101)
    with pytest.raises(ValueError):
        work_phase(t[:-1], t[:-1], t[:-1], 1, .2)
    with pytest.raises(ValueError):
        work_phase(t, t, t, 0, .2)
    with pytest.raises(ValueError):
        duration_sensitivity(100, .02, -.01)
    assert duration_sensitivity(100, .02, .01) == 3600


def test_report_source_binding_and_no_calibration(tmp_path):
    import json
    from solver_v1.run_work_phase_audit import run
    from solver_v1.run_low_frequency_forcing_v29 import sha
    from solver_v1.run_low_stress_cyclic_diagnostic import write_csv
    study, previous, out = [tmp_path/name for name in ('raw','previous','out')]
    study.mkdir(); previous.mkdir()
    protocol=dict(frequency_per_ps=.2,duration_ps=20.,exclude_ps=0.,
                  frame_ps=.025,force_scale_eV_A=[1.,1.],seeds=[123])
    cases=[dict(name='null',axis=None,seed=123,sign=0,restart_sha256='test')]
    cases += [dict(name=f'{axis}_{sign}',axis=axis,seed=123,sign=sign,
                   force_eV_A=float(sign),restart_sha256='test')
              for axis in (0,1) for sign in (-1,1)]
    checks=[]
    t=np.arange(801)*.025; omega=2*np.pi*.2
    for case in cases:
        folder=study/case['name']; folder.mkdir()
        q=np.zeros((len(t),6,3))
        q[:]=(.0001*np.sin(omega*t))[:,None,None]
        metadata=dict(completed=True,restart_sha256='test',elapsed_seconds=1.,duration_ps=20.)
        extra={}
        if case['axis'] is not None:
            force=case['force_eV_A']; axis=case['axis']
            q[:,0,axis]=force*(.4*np.cos(omega*t)+.03*np.sin(omega*t))
            extra['internal_step_work_eV']=force**2*(.4*np.cos(2*omega*t)/4
                +.03*(omega*t/2+np.sin(2*omega*t)/4))
            metadata['conjugate_drive']=dict(axis=axis,force_eV_A=force,frequency_per_ps=.2)
        np.savez(folder/'plane_coordinates.npz',time_seconds=t*1e-12,
                 coordinates_m=q*1e-10,thermo=np.full((len(t),5),300.),**extra)
        (folder/'summary.json').write_text(json.dumps(metadata))
        checks.append(dict(case=case['name'],trajectory_sha256=sha(folder/'plane_coordinates.npz')))
    (study/'protocol.json').write_text(json.dumps(protocol))
    (study/'cases.json').write_text(json.dumps(cases))
    (previous/'summary.json').write_text(json.dumps(dict(completed=True,protocol=protocol,checks=checks)))
    write_csv(previous/'fdt.csv',[dict(axis=axis,parts=1,imag_A2_eV=-.03) for axis in (0,1)])
    run(study,previous,out)
    result=json.loads((out/'summary.json').read_text())
    assert result['completed'] and not result['new_MD_run']
    assert not result['production_clock_calibrated']
    assert result['max_abs_native_coordinate_loss_difference_A2_eV'] < 1e-14
    with pytest.raises(FileExistsError):run(study,previous,out)
    checks[0]['trajectory_sha256']='wrong'
    (previous/'summary.json').write_text(json.dumps(dict(completed=True,protocol=protocol,checks=checks)))
    with pytest.raises(ValueError,match='binding'):run(study,previous,tmp_path/'bad')


def test_thermal_alias_can_fake_low_frequency_loss():
    # Not an Al simulation: exact counterexample to checking only drive Nyquist.
    t=np.arange(8001)*.0025
    omega, thermal=2*np.pi*.2, 2*np.pi*9.8
    q=.01*np.sin(thermal*t)
    w=.01*thermal/2*(np.sin((thermal-omega)*t)/(thermal-omega)
                      +np.sin((thermal+omega)*t)/(thermal+omega))
    fine=work_phase(t[::10],q[::10],w[::10],1.,.2)
    coarse=work_phase(t[::40],q[::40],w[::40],1.,.2)
    assert abs(fine['loss_coordinate_A2_eV']) < 1e-14
    assert abs(coarse['loss_native_A2_eV']) < 1e-14
    assert coarse['loss_coordinate_A2_eV'] == pytest.approx(-.01,abs=1e-13)
