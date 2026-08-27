import numpy as np

from theory.hamiltonian_slip_bath import (
    SlipBathParameters,
    cycle_work,
    simulate_slip_bath,
)


def test_energy_balance_running_case():
    result = simulate_slip_bath(
        SlipBathParameters(force_amplitude=0.50),
        bath_sites=500,
        dt=0.015,
        cycles=8,
        record_stride=10,
    )
    assert result.final_energy_balance_relative_error < 2.0e-4


def test_bounded_intrawell_case_has_no_secular_drift():
    result = simulate_slip_bath(
        SlipBathParameters(force_amplitude=0.34),
        bath_sites=500,
        dt=0.015,
        cycles=8,
        record_stride=10,
    )
    increments = np.diff(result.cycle_slip[-4:])
    assert np.max(np.abs(increments)) < 5.0e-3


def test_running_case_changes_cycle_state():
    result = simulate_slip_bath(
        SlipBathParameters(force_amplitude=0.50),
        bath_sites=500,
        dt=0.015,
        cycles=8,
        record_stride=10,
    )
    increments = np.diff(result.cycle_slip[-4:])
    assert np.mean(increments) < -0.8
    assert np.mean(increments) > -1.2


def test_hysteresis_work_is_positive_in_running_case():
    result = simulate_slip_bath(
        SlipBathParameters(force_amplitude=0.50),
        bath_sites=500,
        dt=0.015,
        cycles=8,
        record_stride=5,
    )
    area = cycle_work(result, 6)
    assert area > 0.0
