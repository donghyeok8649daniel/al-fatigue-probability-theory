from __future__ import annotations
import numpy as np
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
