from __future__ import annotations
import numpy as np
from solver_v1.model import ModelParams, TwoRowLJ

def test_batch_gradient_matches_scalar():
    p = ModelParams(n_cells=3)
    m = TwoRowLJ(p)
    a = np.array([m.a0*1.01, m.a0*0.99, m.a0*1.02])
    s = np.array([0.02, -0.01, 0.03])
    _, ga, gs = m.energy_gradient(a, s, 1.2)
    _, gab, gsb = m.energy_gradient_batch(a[None,:], s[None,:], 1.2)
    assert np.allclose(ga, gab[0], rtol=1e-11, atol=1e-11)
    assert np.allclose(gs, gsb[0], rtol=1e-11, atol=1e-11)

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
