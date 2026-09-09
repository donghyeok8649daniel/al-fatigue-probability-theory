"""Numerical/semantic tests, not evidence of an Al fatigue calibration."""
import numpy as np
import pytest

from .interface_static_scenarios import InterfaceUnits
from .low_stress_cyclic_diagnostic import (
    CyclicTractions, aligned_grid, registry_partition, gibbs_cell_mass,
    interface_traffic, periodic_energy_samples, run_reflecting_cycles,
)
from .probability_pde_2d import _sg_generator_2d, _sg_rates_and_rhs
from types import SimpleNamespace


def toy():
    grid = aligned_grid(1., 1., n_a=9, cells_per_well=9, lower_over_h=.7, upper_over_h=1.3)
    energy = 2*(grid.a[:, None]-1)**2+.12*(1-np.cos(2*np.pi*grid.s[None, :]))
    units = InterfaceUnits(3e-10, 8e-20)
    return grid, energy, units


def test_low_stress_projection_and_units_do_not_create_hz():
    load = CyclicTractions.axial_45(20., period_model=4.)
    np.testing.assert_allclose(load.phase_value_mpa(np.pi/2), [10, 10])
    u = InterfaceUnits(3e-10, 8e-20)
    np.testing.assert_allclose(u.force_to_traction(u.traction_to_force(.02)), .02)
    assert not hasattr(load, "frequency_hz")
    with pytest.raises(ValueError):
        CyclicTractions(period_model=0.)


def test_aligned_negative_registry_and_periodic_samples():
    grid, _, _ = toy()
    n, wells, faces, err = registry_partition(grid, 1.)
    assert tuple(wells) == (-1, 0, 1) and err < 1e-14
    np.testing.assert_array_equal(n[faces+1]-n[faces], 1)
    samples = periodic_energy_samples(lambda a, s: [a*a+np.cos(2*np.pi*s)], grid, 1.)
    np.testing.assert_allclose(samples, grid.a[:, None]**2+np.cos(2*np.pi*grid.s[None, :]), atol=2e-15)
    reflected = periodic_energy_samples(lambda a, s: [a*a+np.cos(2*np.pi*s)], grid, 1., mirror_registry=True)
    np.testing.assert_allclose(reflected, samples, atol=2e-15)


def test_gibbs_stationarity_and_exact_sg_generator_identity():
    grid, energy, _ = toy()
    model = SimpleNamespace(p=SimpleNamespace(kT=.025, mobility_a=1., mobility_s=.05))
    n, _, _, _ = registry_partition(grid, 1.)
    mass = gibbs_cell_mass(energy, .025, n, central_only=False)
    generator = _sg_generator_2d(energy, model, grid)
    np.testing.assert_allclose(generator @ mass.ravel(), 0., atol=2e-15)
    tilted = mass*(1+.1*grid.s[None, :])
    rhs, _ = _sg_rates_and_rhs(tilted, energy, model, grid)
    np.testing.assert_allclose(generator @ tilted.ravel(), rhs.ravel(), atol=2e-15)


def test_reflecting_cycles_conserve_population_moment_and_have_no_crack_field():
    grid, energy, units = toy()
    result = run_reflecting_cycles(energy, grid, h=1., period=1., units=units,
        load=CyclicTractions.axial_45(30., period_model=2.), cycles=3,
        hold_cycles=2, steps_per_cycle=32)
    audit = result["audit"]
    assert audit["max_mass_residual"] < 1e-12
    assert audit["max_well_balance_residual"] < 2e-14
    assert audit["max_registry_balance_residual"] < 2e-14
    assert audit["max_strain_decomposition_residual"] < 1e-14
    assert audit["minimum_cell_mass"] >= 0 and audit["numerical_repair"] == 0
    assert audit["opening_probability"] is None
    assert not audit["physical_hz_available"] and not audit["production_model_changed"]
    end = result["history"][-1]
    assert abs(end["mean_well_index"]-end["cumulative_net_transfer"]) < 1e-12
    assert end["cumulative_gross_grid_traffic"] >= abs(end["cumulative_net_transfer"])
    assert result["cycles"][-1]["segment"] == "zero_stress_hold"


def test_unloaded_stationary_ensemble_is_not_fatigue_despite_gross_traffic():
    grid, energy, units = toy()
    result = run_reflecting_cycles(energy, grid, h=1., period=1., units=units,
        load=CyclicTractions(period_model=1.), cycles=2, hold_cycles=0,
        steps_per_cycle=16, central_only=False)
    np.testing.assert_allclose(result["initial_mass"], result["final_mass"], atol=2e-15)
    assert abs(result["history"][-1]["cumulative_net_transfer"]) < 1e-14
    assert result["history"][-1]["cumulative_gross_grid_traffic"] > 0


def test_gross_sg_traffic_has_grid_recrossing_divergence_not_physical_hop_count():
    gross, rescaled = [], []
    for cells in (9, 19, 39):
        grid = aligned_grid(1., 1., n_a=5, cells_per_well=cells)
        n, _, faces, _ = registry_partition(grid, 1.)
        energy = np.zeros((len(grid.a), len(grid.s)))
        mass = gibbs_cell_mass(energy, .025, n, central_only=False)
        plus, minus = interface_traffic(mass, energy, grid, .025, .05, faces)
        np.testing.assert_allclose(plus-minus, 0., atol=1e-16)
        gross.append(float(np.sum(plus+minus)))
        rescaled.append(gross[-1]*grid.ds)
    assert gross[2] > gross[1] > gross[0]
    # Integer face spacing is represented by floating-point linspace; allow
    # its accumulated arithmetic error, not a physical traffic tolerance.
    np.testing.assert_allclose(rescaled, rescaled[0], rtol=64*np.finfo(float).eps)


def test_backward_euler_timestep_convergence_of_driven_density():
    grid, energy, units = toy()
    values = []
    for steps in (16, 32, 64, 128):
        r = run_reflecting_cycles(energy, grid, h=1., period=1., units=units,
            load=CyclicTractions.axial_45(40., period_model=2.), cycles=2,
            hold_cycles=0, steps_per_cycle=steps)
        values.append(r["final_mass"])
    errors = [np.sum(abs(v-values[-1])) for v in values[:-1]]
    assert errors[2] < errors[1] < errors[0]


def test_different_mpa_does_not_change_static_energy_or_mobility():
    grid, energy, units = toy()
    original = energy.copy()
    values = []
    for amplitude in (0., 50.):
        r = run_reflecting_cycles(energy, grid, h=1., period=1., units=units,
            load=CyclicTractions.axial_45(amplitude, period_model=2.), cycles=1,
            hold_cycles=0, steps_per_cycle=32)
        values.append(r["cycles"][0]["registry_shear_amplitude"])
        assert r["audit"]["mobility_s_model"] == .05
    np.testing.assert_array_equal(energy, original)
    assert values[1] > values[0]


def test_explicit_implicit_same_generator_converges_and_keeps_flux_balance():
    grid, energy, units = toy()
    implicit, explicit = [], []
    for steps in (16, 64):
        values = []
        for integrator in ("implicit", "explicit"):
            r = run_reflecting_cycles(energy, grid, h=1., period=1., units=units,
                load=CyclicTractions.axial_45(30., period_model=.05), cycles=1,
                hold_cycles=0, steps_per_cycle=steps, integrator=integrator)
            assert r["audit"]["max_well_balance_residual"] < 1e-14
            values.append(r["final_mass"])
        implicit.append(values[0]); explicit.append(values[1])
    assert np.sum(abs(implicit[1]-explicit[1])) < np.sum(abs(implicit[0]-explicit[0]))


def test_signed_shear_reflection_reverses_registry_not_gross_activity():
    grid, energy, units = toy()
    outputs = []
    for phase in (0., np.pi):
        outputs.append(run_reflecting_cycles(energy, grid, h=1., period=1., units=units,
            load=CyclicTractions(normal_amplitude_mpa=15., shear_amplitude_mpa=15.,
                                 phase_s=phase, period_model=2.),
            cycles=3, hold_cycles=2, steps_per_cycle=32))
    np.testing.assert_allclose(outputs[0]["final_mass"], outputs[1]["final_mass"][:, ::-1], atol=5e-15)
    x, y = [r["history"][-1] for r in outputs]
    assert abs(x["mean_well_index"]+y["mean_well_index"]) < 1e-14
    np.testing.assert_allclose(x["cumulative_gross_grid_traffic"], y["cumulative_gross_grid_traffic"], atol=1e-14)


def test_lossless_output_and_explicit_snapshot_sampling(tmp_path):
    from .run_low_stress_cyclic_diagnostic import write_csv
    from .summarize_low_stress_cycles import rows, snapshot
    file = tmp_path/"history.csv.gz"
    source = [{"cycle": 1, "x": float(np.nextafter(.123, 1.))}, {"cycle": 2, "x": 1e-18}]
    write_csv(file, source)
    assert rows(tmp_path/"history.csv") == source
    full = np.arange(10*3*5).reshape(10, 3, 5)
    compact = dict(cycle_indices=np.array([0, 1, 7, 8, 9]), same_phase=full[[0, 1, 7, 8, 9]])
    np.testing.assert_array_equal(snapshot(compact, 8), full[8])
    np.testing.assert_array_equal(snapshot(compact, -1), full[-1])
    with pytest.raises(ValueError):
        snapshot(compact, 6)
