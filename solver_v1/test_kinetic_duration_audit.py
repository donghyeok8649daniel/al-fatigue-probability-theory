import json

import numpy as np
import pytest
from scipy.linalg import toeplitz

from solver_v1.kinetic_duration_audit import (KB_EV_K, cycle_windows, sine_weights,
    stationary_quadrature_variance, conditional_precision_budget, thermal_work_summary)
from solver_v1.run_kinetic_duration_audit import run, THERMO_COLUMNS
from solver_v1.run_low_frequency_forcing_v29 import sha
from solver_v1.run_low_stress_cyclic_diagnostic import write_csv
from solver_v1.work_phase_audit import duration_sensitivity


def test_signed_phase_partition_preserves_all_samples():
    t = np.linspace(1.25, 21.25, 4001)
    q = 7.+2*(.4*np.cos(2*np.pi*.2*t)-.03*np.sin(2*np.pi*.2*t))+.001*t
    windows = cycle_windows(t, .2, 5., 1.25)
    estimates = [sine_weights(t[i:j+1], .2, 2) @ q[i:j+1] for i, j in windows]
    full = sine_weights(t, .2, 2) @ q
    assert np.mean(estimates) == pytest.approx(full, abs=1e-14)
    assert sine_weights(t, .2, -2) @ q == pytest.approx(-full, abs=1e-14)
    assert full < 0  # Signed loss and thermal trend are never clipped.
    assert windows == [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4000)]


@pytest.mark.parametrize('kind', ['partial', 'nonuniform', 'offgrid', 'nyquist', 'empty'])
def test_window_rejects_invalid_sampling(kind):
    t = np.linspace(0, 20, 4001)
    with pytest.raises(ValueError):
        if kind == 'partial':
            cycle_windows(t[:-1], .2, 5.)
        elif kind == 'nonuniform':
            t[5] += .001
            cycle_windows(t, .2, 5.)
        elif kind == 'offgrid':
            cycle_windows(t, .2, 5., .001)
        elif kind == 'nyquist':
            cycle_windows(t, 110., 5.)
        else:
            sine_weights([], .2, 2.)


def test_stationary_variance_matches_independent_dense_covariance():
    t = np.linspace(.4, 20.4, 401)
    weights = sine_weights(t, .2, 2.)
    covariance = .003*np.exp(-2*(t-t[0]))
    exact = weights @ toeplitz(covariance) @ weights
    assert stationary_quadrature_variance(weights, covariance) == pytest.approx(exact, rel=2e-14)
    with pytest.raises(ValueError, match='negative'):
        stationary_quadrature_variance([1., -1.], [1., 2.])


def test_ou_finite_window_and_two_sided_fdt_factor():
    # Analytic OU covariance integral, independent of the FFT implementation.
    variance, decay, force, temperature, frequency = .003, 2., 2., 300., .2
    omega = 2*np.pi*frequency
    denominator = decay**2+omega**2
    loss = variance*decay*omega/(KB_EV_K*temperature*denominator)
    duration = 20.
    asymptotic = 4*variance*decay/(force**2*duration*denominator)
    boundary = 8*variance*omega**2*(-np.expm1(-decay*duration))/(force**2*duration**2*denominator**2)
    errors = []
    for dt in (.01, .005):
        t = np.arange(round(duration/dt)+1)*dt
        actual = stationary_quadrature_variance(sine_weights(t, frequency, force),
                                               variance*np.exp(-decay*t))
        errors.append(abs(actual-asymptotic-boundary))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert errors[1] < 2e-5*asymptotic
    budget = conditional_precision_budget(force, frequency, loss, temperature, duration)
    assert budget['predicted_loss_std_A2_eV']**2 == pytest.approx(asymptotic)
    assert budget['conditional_work_at_snr_eV'] == pytest.approx(2*KB_EV_K*temperature*9)


def test_weaker_force_changes_duration_but_not_work_at_fixed_snr():
    first = conditional_precision_budget(2., .2, .03, 300., 20.)
    weaker = conditional_precision_budget(1., .2, .03, 300., 20.)
    assert weaker['conditional_duration_ps'] == pytest.approx(4*first['conditional_duration_ps'])
    assert weaker['conditional_work_at_snr_eV'] == pytest.approx(first['conditional_work_at_snr_eV'])
    assert not first['production_clock_calibrated'] and not first['confidence_interval']
    with pytest.raises(ValueError):
        conditional_precision_budget(2., .2, -.03, 300., 20.)


def test_temperature_trend_and_native_energy_residual_use_sliced_origin():
    t = np.linspace(100., 200., 101)
    w = 7.+.002*t
    h = np.zeros((len(t), 5))
    h[:, 0] = 300.+.001*t
    h[:, 3] = -10.+w+1e-5*t
    out = thermal_work_summary(t, h, w)
    assert out['temperature_linear_trend_K_ns'] == pytest.approx(1.)
    assert out['work_eV'] == pytest.approx(.2)
    assert out['max_work_residual_eV'] == pytest.approx(.001)
    assert out['endpoint_work_residual_eV'] == pytest.approx(.001)
    null = thermal_work_summary(t, h)
    assert null['work_eV'] is None and null['max_work_residual_eV'] is None


def _synthetic_campaign(tmp_path):
    study, previous, work_audit = [tmp_path/name for name in ('raw', 'previous', 'work')]
    for path in (study, previous, work_audit):
        path.mkdir()
    restart = study/'seed123_init'/'equilibrated.restart'
    restart.parent.mkdir()
    restart.write_bytes(b'synthetic restart, not MD')
    restart_hash = sha(restart)
    protocol = dict(seeds=[123], duration_ps=20., exclude_ps=0., dt_ps=.00125,
                    frame_ps=.025, frequency_per_ps=.2, force_scale_eV_A=[2., 1.])
    cases = [dict(name='null', seed=123, axis=None, sign=0, restart_sha256=restart_hash)]
    cases += [dict(name=f'axis{axis}_sign{sign}', seed=123, axis=axis, sign=sign,
                   force_eV_A=sign*protocol['force_scale_eV_A'][axis], restart_sha256=restart_hash)
              for axis in (0, 1) for sign in (-1, 1)]
    checks = []
    t = np.arange(801)*.025
    for case in cases:
        folder = study/case['name']
        folder.mkdir()
        q = .0001*np.sin(2*np.pi*.2*t)[:, None, None]*np.array([1., -1.])[None, :, None]*np.ones((1, 1, 3))
        h = np.zeros((len(t), 5))
        h[:, 0] = 300.+.001*t
        h[:, 3] = -10.+.002*t
        extra = {} if case['axis'] is None else dict(internal_step_work_eV=.002*t)
        drive = None if case['axis'] is None else dict(axis=case['axis'], force_eV_A=case['force_eV_A'], frequency_per_ps=.2)
        np.savez(folder/'plane_coordinates.npz', time_seconds=t*1e-12, coordinates_m=q*1e-10, thermo=h, **extra)
        metadata = dict(completed=True, ensemble='nve', thermostat_damping_ps=None,
            dt_ps=.00125, frame_ps=.025, duration_ps=20., conjugate_drive=drive,
            restart_sha256=restart_hash, internal_step_work=drive is not None,
            thermo_columns=THERMO_COLUMNS, atom_count=8, repeats=2, atoms_per_plane=4,
            elapsed_seconds=1.)
        (folder/'summary.json').write_text(json.dumps(metadata))
        checks.append(dict(case=case['name'], trajectory_sha256=sha(folder/'plane_coordinates.npz')))
    (study/'protocol.json').write_text(json.dumps(protocol))
    (study/'cases.json').write_text(json.dumps(cases))
    (previous/'summary.json').write_text(json.dumps(dict(completed=True, protocol=protocol, checks=checks)))
    write_csv(previous/'fdt.csv', [dict(seed=123, axis=axis, parts=1, nw=2., half_band_cycles_ps=.01,
                                     imag_A2_eV=-.03) for axis in (0, 1)])
    (work_audit/'summary.json').write_text(json.dumps(dict(completed=True, new_MD_run=False,
        raw_hashes={r['case']: r['trajectory_sha256'] for r in checks}, planning=[dict(axis=axis,
            conditional_record_duration_ps=duration_sensitivity(20., .0001/protocol['force_scale_eV_A'][axis], .03),
            independent_fdt_loss_A2_eV=.03,
            observed_imaginary_null_A2_eV=.0001/protocol['force_scale_eV_A'][axis],
            desired_ratio=3., estimated_single_record_wall_hours=1.)
            for axis in (0, 1)])))
    return study, previous, work_audit


def test_report_binding_negative_results_and_replay(tmp_path):
    paths = _synthetic_campaign(tmp_path)
    out = tmp_path/'report'
    summary = run(*paths, out, parts=(1, 2, 4), thermal_window_ps=5.)
    assert summary['completed'] and not summary['new_MD_run']
    assert summary['max_null_plane_sum_A'] == 0.
    assert not summary['stationarity_certified'] and not summary['production_clock_calibrated']
    assert summary['full_record_null'][0]['max_abs_mean_over_planes_A2_eV'] == 0.
    assert summary['full_record_null'][0]['rms_loss_A2_eV'] > 0.
    replay = tmp_path/'replay'
    run(*paths, replay, parts=(1, 2, 4), thermal_window_ps=5.)
    for path in out.glob('*.csv'):
        assert path.read_bytes() == (replay/path.name).read_bytes()
    with pytest.raises(FileExistsError):
        run(*paths, out, parts=(1, 2, 4), thermal_window_ps=5.)
    # Temperature is taken from raw thermo, so altered bytes must be rejected.
    raw = paths[0]/'null'/'plane_coordinates.npz'
    raw.write_bytes(raw.read_bytes()+b'changed')
    with pytest.raises(ValueError, match='binding'):
        run(*paths, tmp_path/'bad', parts=(1, 2, 4), thermal_window_ps=5.)
    assert not (tmp_path/'bad').exists()


@pytest.mark.parametrize('kind', ['duplicate', 'force', 'restart', 'thermostat', 'planning'])
def test_report_rejects_mismatched_campaign_before_writing(tmp_path, kind):
    paths = _synthetic_campaign(tmp_path)
    study = paths[0]
    if kind in ('duplicate', 'force'):
        path = study/'cases.json'
        cases = json.loads(path.read_bytes())
        if kind == 'duplicate':
            cases[-1] = cases[-2]
        else:
            cases[-1]['force_eV_A'] *= 2
        path.write_text(json.dumps(cases))
    elif kind == 'restart':
        (study/'seed123_init'/'equilibrated.restart').write_bytes(b'changed')
    elif kind == 'planning':
        path = paths[2]/'summary.json'
        metadata = json.loads(path.read_bytes())
        metadata['planning'][0]['conditional_record_duration_ps'] *= 2
        path.write_text(json.dumps(metadata))
    else:
        path = study/'null'/'summary.json'
        metadata = json.loads(path.read_bytes())
        metadata['thermostat_damping_ps'] = .1
        path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError):
        run(*paths, tmp_path/'bad', parts=(1, 2, 4), thermal_window_ps=5.)
    assert not (tmp_path/'bad').exists()
