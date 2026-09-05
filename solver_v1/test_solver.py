from __future__ import annotations
import numpy as np
from solver_v1.deterministic_normal import (
    DeterministicRunParams,
    InitialMeasureAtom,
    NormalChainParams,
    empirical_phase_space_support,
    phi_second,
    run_deterministic_pushforward,
    spacing_acceleration,
)
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.solver import LoadParams, SolverParams, run_ensemble

def test_batch_gradient_matches_scalar():
    p = ModelParams(n_cells=3)
    m = TwoRowLJ(p)
    a = np.array([m.a0*1.01, m.a0*0.99, m.a0*1.02])
    s = np.array([0.02, -0.01, 0.03])
    _, ga, gs = m.energy_gradient(a, s, 1.2)
    _, gab, gsb = m.energy_gradient_batch(a[None,:], s[None,:], 1.2)
    assert np.allclose(ga, gab[0], rtol=1e-11, atol=1e-11)
    assert np.allclose(gs, gsb[0], rtol=1e-11, atol=1e-11)

def test_full_hessian_matches_finite_difference_of_full_gradient():
    p = ModelParams(n_cells=3)
    m = TwoRowLJ(p)
    a = np.array([m.a0*1.01, m.a0*0.99, m.a0*1.02])
    s = np.array([0.02, -0.01, 0.03])
    haa, has, hss = m.hessian_blocks(a, s)
    analytical = np.block([[haa, has], [has.T, hss]])
    state = np.concatenate((a, s))
    numerical = np.empty_like(analytical)
    delta = 1.0e-6
    for column in range(len(state)):
        plus = state.copy(); plus[column] += delta
        minus = state.copy(); minus[column] -= delta
        _, ga_plus, gs_plus = m.energy_gradient(plus[:3], plus[3:], 1.2)
        _, ga_minus, gs_minus = m.energy_gradient(minus[:3], minus[3:], 1.2)
        numerical[:, column] = (
            np.concatenate((ga_plus, gs_plus)) - np.concatenate((ga_minus, gs_minus))
        ) / (2.0*delta)
    assert np.allclose(analytical, numerical, rtol=2e-5, atol=2e-5)

def test_strain_components_sum_to_declared_total_strain():
    p = ModelParams(n_cells=3)
    m = TwoRowLJ(p)
    a = np.array([[m.a0*1.01, m.a0*0.99, m.a0*1.02]])
    s = np.array([[1.12, -0.06, 2.03]])
    total, normal, intrawell, plastic = m.strain_components_batch(a, s)
    assert np.isclose(total[0], m.strain(a[0], s[0]))
    assert np.isclose(total[0], normal[0] + intrawell[0] + plastic[0])

def test_reference_normal_tangent_matches_uniform_finite_difference():
    m = TwoRowLJ(ModelParams(n_cells=3))
    tangent = m.reference_normal_tangent_force_per_strain()
    delta = 1.0e-6
    a_plus = np.full(3, m.a0*(1.0 + delta))
    a_minus = np.full(3, m.a0*(1.0 - delta))
    s = np.zeros(3)
    _, ga_plus, _ = m.energy_gradient(a_plus, s, 0.0)
    _, ga_minus, _ = m.energy_gradient(a_minus, s, 0.0)
    numerical = float(np.mean(ga_plus - ga_minus) / (2.0*delta))

    assert tangent > 0.0
    assert np.isclose(tangent, numerical, rtol=2.0e-7, atol=2.0e-7)

def test_coupled_opening_escape_uses_outward_full_normal_mode():
    m = TwoRowLJ(ModelParams(n_cells=1))
    force = 3.0
    minimum, saddle = m.opening_stationary_points(0.0, force)
    assert minimum is not None and saddle is not None
    inside = np.array([[0.5*(minimum + saddle)]])
    outside = np.array([[saddle*1.01]])
    s = np.zeros((1, 1))
    inside_escape, *_ = m.coupled_opening_escape_batch(inside, s, force)
    outside_escape, *_ = m.coupled_opening_escape_batch(outside, s, force)
    assert not inside_escape[0]
    assert outside_escape[0]

def test_opening_barrier_decreases_with_tension():
    m = TwoRowLJ(ModelParams())
    m._build_opening_table()
    b1 = m.opening_barrier(0.0, 2.0)
    b2 = m.opening_barrier(0.0, 3.0)
    assert b2 < b1
    assert b2 > 0.0

def test_periodic_registry_energy():
    p = ModelParams()
    m = TwoRowLJ(p)
    a = m.a0
    e0 = m.local_energy(a, 0.13)
    e1 = m.local_energy(a, 0.13 + p.b)
    assert np.isclose(e0, e1, rtol=1e-11, atol=1e-11)

def test_solver_returns_per_trajectory_first_passage_or_censoring():
    result = run_ensemble(
        ModelParams(),
        LoadParams(force_min=0.0, force_max=3.6, period=1.0, cycles=1),
        SolverParams(dt=0.02, n_trajectories=4, seed=7, record_stride=10),
    )
    passage = result["first_passage_time"]
    assert passage.shape == (4,)
    assert result["invalid_trajectory"].shape == (4,)
    assert not np.any(result["invalid_trajectory"] & np.isfinite(passage))
    assert np.all(np.isnan(passage) | ((passage > 0.0) & (passage <= 1.0)))
    assert result["observation_end_time"] == 1.0

def test_solver_stream_can_be_stopped_without_retaining_history():
    records = []

    result = run_ensemble(
        ModelParams(),
        LoadParams(force_min=0.0, force_max=1.0, period=1.0, cycles=100),
        SolverParams(dt=0.02, n_trajectories=4, seed=7, record_stride=1),
        record_callback=records.append,
        stop_requested=lambda: len(records) >= 3,
        retain_history=False,
    )

    assert len(records) == 3
    assert result["time"].size == 0
    assert result["observation_end_time"] < 100.0

def test_active_normal_potential_has_unit_reference_tangent_and_declared_boundary():
    params = NormalChainParams(n_cells=3)
    assert np.isclose(phi_second(np.array([1.0]), params)[0], 1.0)
    assert np.isclose(phi_second(np.array([params.lambda_c]), params)[0], 0.0, atol=1e-12)

def test_deterministic_spacing_acceleration_matches_declared_boundaries():
    params = NormalChainParams(n_cells=3)
    spacing = np.array([[0.98, 1.01, 1.03]])
    traction = 0.002
    acceleration = spacing_acceleration(spacing, traction, params)[0]
    from solver_v1.deterministic_normal import phi_prime
    gradient = phi_prime(spacing[0], params)
    expected = np.array([
        gradient[1] - gradient[0],
        gradient[2] - 2.0*gradient[1] + gradient[0],
        traction + gradient[1] - 2.0*gradient[2],
    ])
    np.testing.assert_allclose(acceleration, expected)

def test_discrete_initial_measure_is_pushed_forward_without_sampling():
    params = NormalChainParams(n_cells=2)
    measure = (
        InitialMeasureAtom(0.75, np.ones(2), np.zeros(2)),
        InitialMeasureAtom(0.25, np.full(2, params.lambda_c*1.001), np.zeros(2)),
    )
    run = DeterministicRunParams(dt=0.01, duration=0.0, record_stride=1)
    first = run_deterministic_pushforward(params, run, lambda _t: 0.0, initial_measure=measure)
    second = run_deterministic_pushforward(params, run, lambda _t: 0.0, initial_measure=measure)

    assert first["survival"][0] == 0.75
    assert first["specimen_survival"][0] == 0.75
    assert np.array_equal(first["survival"], second["survival"])
    assert np.array_equal(first["strain"], second["strain"])

    lam, velocity, mass = empirical_phase_space_support(
        first["spacing_support"][0],
        first["spacing_rate_support"][0],
        first["measure_weights"],
    )
    assert lam.shape == velocity.shape == mass.shape == (4,)
    assert np.isclose(np.sum(mass), 1.0)
    np.testing.assert_allclose(mass, [0.375, 0.375, 0.125, 0.125])

def test_unforced_deterministic_chain_conserves_mechanical_energy_numerically():
    params = NormalChainParams(n_cells=3)
    measure = (
        InitialMeasureAtom(
            1.0,
            np.array([1.01, 0.995, 1.005]),
            np.array([0.002, -0.001, 0.0005]),
        ),
    )
    result = run_deterministic_pushforward(
        params,
        DeterministicRunParams(dt=0.001, duration=1.0, record_stride=10),
        lambda _t: 0.0,
        initial_measure=measure,
    )
    energy = result["mechanical_energy"]
    assert np.max(np.abs(energy - energy[0])) < 1.0e-9
