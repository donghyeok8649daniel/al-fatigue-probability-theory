import numpy as np
import pytest
from .impedance_mobility import (
    MOBILITY_SI_PER_A2_EV_PS, scalar_impedance, inverse_disk_bounds,
    scalar_low_band_drag, fit_constant_drag,
)


def test_exact_inertial_response_does_not_turn_mass_into_overdamped_clock():
    H, mass, drag, f = 4., .2, .3, .4
    w = 2*np.pi*f
    chi = 1/(H-mass*w*w+1j*w*drag)
    r = scalar_impedance(chi, f, static_chi=1/H)
    assert r['drag_eV_ps_A2'] == pytest.approx(drag)
    assert r['dispersion_coefficient_eV_ps2_A2'] == pytest.approx(mass)
    assert r['mobility_m2_per_J_s'] == pytest.approx(MOBILITY_SI_PER_A2_EV_PS/drag)
    assert not r['production_clock_calibrated']


def test_exact_inverse_disk_encloses_boundary_and_attains_extrema():
    chi, radius, f = .7-.2j, .03, .6
    theta = np.linspace(0, 2*np.pi, 100001)
    gamma = (1/(chi+radius*np.exp(1j*theta))).imag/(2*np.pi*f)
    r = inverse_disk_bounds(chi, radius, f)
    assert gamma.min() >= r['drag_lower_eV_ps_A2']-1e-14
    assert gamma.max() <= r['drag_upper_eV_ps_A2']+1e-14
    assert gamma.min() == pytest.approx(r['drag_lower_eV_ps_A2'], abs=1e-10)
    assert gamma.max() == pytest.approx(r['drag_upper_eV_ps_A2'], abs=1e-10)


def test_unresolved_loss_has_no_finite_mobility_upper_bound():
    r = inverse_disk_bounds(1-.001j, .002, 1.)
    assert not r['positive_drag_resolved']
    assert r['mobility_upper_m2_per_J_s'] is None
    assert r['mobility_lower_m2_per_J_s'] > 0
    assert not inverse_disk_bounds(1j, .1, 1.)['passive_drag_feasible']
    assert not inverse_disk_bounds(.001j, .002, 1.)['inverse_set_bounded']
    assert not inverse_disk_bounds(1j, 1., 1.)['passive_drag_feasible']
    assert scalar_impedance(1+.01j, 1.)['mobility_m2_per_J_s'] is None


def test_zero_radius_exact_and_sign_not_repaired():
    chi, f = 1-.4j, .3
    r = inverse_disk_bounds(chi, 0., f)
    point = scalar_impedance(chi, f)
    assert r['mobility_lower_m2_per_J_s'] == pytest.approx(point['mobility_m2_per_J_s'])
    assert r['mobility_upper_m2_per_J_s'] == pytest.approx(point['mobility_m2_per_J_s'])


def test_low_band_scalar_is_not_matrix_inverse_diagonal():
    C = np.eye(2)*2e-24
    K = np.array([[2., 1.], [1., 3.]])*1e-39
    kT = 1.380649e-23*300
    scalar_M = 1/scalar_low_band_drag(C[0,0], K[0,0], 300)
    matrix_M = C@np.linalg.inv(K)@C/kT
    assert scalar_M != pytest.approx(matrix_M[0,0])


def test_deterministic_drag_fit_scaling_and_rejection():
    y = np.array([2., 3., 4.]); s = np.array([1., 2., 1.])
    r = fit_constant_drag(y, s)
    assert r['drag'] == pytest.approx(3.)
    assert fit_constant_drag(y*1e-11, s*1e-11)['drag'] == pytest.approx(3e-11)
    assert not fit_constant_drag(-y, s)['admissible']
    assert not r['production_clock_calibrated']
    assert MOBILITY_SI_PER_A2_EV_PS == pytest.approx(6.241509074460763e10)


@pytest.mark.parametrize('value', [0, -1, np.nan])
def test_invalid_inputs_rejected(value):
    with pytest.raises(ValueError): scalar_impedance(1j, value)
    with pytest.raises(ValueError): scalar_low_band_drag(1., value, 300)


def test_completed_record_calibration_reproducible_and_clock_closed(tmp_path):
    import csv
    import json
    from .run_impedance_calibration import run
    first, second = tmp_path/'first', tmp_path/'second'
    run(first); run(second)
    assert {p.name:p.read_bytes() for p in first.iterdir()} == {
        p.name:p.read_bytes() for p in second.iterdir()}
    with pytest.raises(FileExistsError): run(first)
    candidate=json.loads((first/'calibration_candidate.json').read_bytes())
    assert not candidate['production_clock_calibrated']
    assert candidate['t0_seconds'] is None
    assert not candidate['zero_frequency_limit_certified']
    assert len(candidate['fits']) == 2
    with (first/'low_band_residuals.csv').open() as stream:
        rows=list(csv.DictReader(stream))
    assert all((r['fitted_band']=='True') == (float(r['lower_per_ps']) >= .02) for r in rows)
    assert all(r['normalized_residual']=='' for r in rows if r['fitted_band']=='False')
    with (first/'forced_impedance.csv').open() as stream:
        rows=list(csv.DictReader(stream))
    low=[r for r in rows if float(r['frequency_per_ps'])==.05 and int(r['parts'])==1]
    assert all(r['mobility_upper_m2_per_J_s']=='' for r in low)
