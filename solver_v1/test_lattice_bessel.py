import numpy as np

from solver_v1.lattice_bessel import (
    FourierLatticeConfig,
    two_row_lj_direct_reference,
    two_row_lj_infinite_energy_gradient,
)
from solver_v1.model import ModelParams, TwoRowLJ


CFG = FourierLatticeConfig(tol=1.0e-13, max_modes=64)


def _analytic(a, s):
    return two_row_lj_infinite_energy_gradient(
        a,
        s,
        epsilon=1.0,
        sigma_lj=0.82,
        b=1.0,
        config=CFG,
    )


def test_bessel_kernel_matches_large_direct_lattice_sum():
    samples = [
        (0.70, 0.00),
        (0.70, 0.27),
        (0.85, -0.41),
        (1.00, 0.10),
        (1.40, 0.49),
    ]
    for a, s in samples:
        energy, deda, deds, modes = _analytic(a, s)
        ref_energy, ref_deda, ref_deds = two_row_lj_direct_reference(
            a,
            s,
            epsilon=1.0,
            sigma_lj=0.82,
            b=1.0,
            images=4096,
        )
        assert abs(float(energy - ref_energy)) < 2.0e-11
        assert abs(float(deda - ref_deda)) < 2.0e-10
        assert abs(float(deds - ref_deds)) < 2.0e-10
        assert modes < CFG.max_modes


def test_bessel_kernel_is_periodic_in_registry_coordinate():
    a = np.array([0.72, 0.86, 1.10])
    s = np.array([-0.37, 0.11, 0.43])
    e0, ga0, gs0, _ = _analytic(a, s)
    e1, ga1, gs1, _ = _analytic(a, s + 1.0)
    np.testing.assert_allclose(e0, e1, rtol=0.0, atol=2.0e-12)
    np.testing.assert_allclose(ga0, ga1, rtol=0.0, atol=2.0e-11)
    np.testing.assert_allclose(gs0, gs1, rtol=0.0, atol=2.0e-11)


def test_bessel_analytic_gradients_match_finite_difference():
    a = 0.91
    s = 0.23
    h = 2.0e-6
    energy, deda, deds, _ = _analytic(a, s)
    ep, _, _, _ = _analytic(a + h, s)
    em, _, _, _ = _analytic(a - h, s)
    sp, _, _, _ = _analytic(a, s + h)
    sm, _, _, _ = _analytic(a, s - h)
    fd_a = float((ep - em) / (2.0 * h))
    fd_s = float((sp - sm) / (2.0 * h))
    assert abs(float(deda) - fd_a) < 2.0e-7
    assert abs(float(deds) - fd_s) < 2.0e-7
    assert np.isfinite(float(energy))


def test_model_scalar_and_batch_energy_gradients_agree_with_bessel_kernel():
    params = ModelParams(n_cells=1, lattice_fourier_tol=1.0e-13)
    model = TwoRowLJ(params)
    a = np.array([0.93])
    s = np.array([0.18])
    force = 0.37

    energy, ga, gs = model.energy_gradient(a, s, force)
    benergy, bga, bgs = model.energy_gradient_batch(
        a.reshape(1, 1), s.reshape(1, 1), force
    )

    np.testing.assert_allclose([energy], benergy, rtol=0.0, atol=2.0e-12)
    np.testing.assert_allclose(ga.reshape(1, 1), bga, rtol=0.0, atol=2.0e-12)
    np.testing.assert_allclose(gs.reshape(1, 1), bgs, rtol=0.0, atol=2.0e-12)
