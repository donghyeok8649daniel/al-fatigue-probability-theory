import math

import numpy as np
import pytest

from solver_v1.analytic_lj_eam import SquareRootEmbedding
from solver_v1.energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
)
from solver_v1.fcc111_full_energy import FullFCC111StackEnergy
from solver_v1.fcc111_geometry import DIRECT_110, SHOCKLEY_112, registry_path


@pytest.fixture(scope="module")
def full_lj():
    return FullFCC111StackEnergy.from_reduced_model(
        build_energy_model(TWO_ROW_LJ_REFERENCE)
    )


@pytest.fixture(scope="module")
def full_hypothetical():
    return FullFCC111StackEnergy.from_reduced_model(
        build_energy_model(ANALYTIC_LJ_EAM_HYPOTHETICAL)
    )


def test_existing_reduced_models_remain_numerically_unchanged():
    lj = build_energy_model(TWO_ROW_LJ_REFERENCE)
    hybrid = build_energy_model(ANALYTIC_LJ_EAM_HYPOTHETICAL)
    calibrated = build_energy_model(AL_TARGET_BEST_FEASIBLE)
    assert lj.a0 == pytest.approx(0.7713438268704838, abs=2e-15)
    assert lj.sigma_over_E_force_scale() == pytest.approx(86.29296488740997)
    assert hybrid.a0 == pytest.approx(0.9677866043190319)
    assert calibrated.a0 == pytest.approx(0.8164957221838465)


def test_full_lj_equilibrium_hessian_and_kappa_are_stable(full_lj):
    assert full_lj.a0 == pytest.approx(0.6845871380591213, rel=2e-11)
    assert abs(full_lj.local_deda(full_lj.a0, 0.0)) < 1e-10
    assert abs(full_lj.local_deds(full_lj.a0, 0.0)) < 1e-11
    hessian = full_lj.local_hessian(full_lj.a0, 0.0)
    np.testing.assert_allclose(hessian, hessian.T, atol=1e-14)
    assert np.min(np.linalg.eigvalsh(hessian)) > 0.0
    assert full_lj.sigma_over_E_force_scale() == pytest.approx(
        135.17737716036478, rel=2e-10
    )


def test_full_hybrid_applies_embedding_once_after_total_density(full_hypothetical):
    result = full_hypothetical.evaluate(full_hypothetical.a0, 0.0)
    assert result.density.value > 0.0
    expected = full_hypothetical.embedding.value(result.density.value)
    assert result.embedding_energy == pytest.approx(float(expected), rel=2e-14)
    assert result.embedding_energy != pytest.approx(2.0 * float(expected))
    assert result.energy == pytest.approx(
        result.pair.value + result.embedding_energy, rel=2e-14
    )


def test_zero_embedding_recovers_full_lj_pair_exactly(full_lj):
    reduced = build_energy_model(ANALYTIC_LJ_EAM_HYPOTHETICAL)
    zero = FullFCC111StackEnergy(
        reduced.p,
        density_params=reduced.density_params,
        embedding=SquareRootEmbedding(amplitude=0.0, rho_ref=1.0),
    )
    pair_only = FullFCC111StackEnergy(reduced.p)
    point = (0.91, 0.137)
    assert zero.energy(*point) == pytest.approx(pair_only.energy(*point), abs=2e-14)
    np.testing.assert_allclose(zero.grad(*point), pair_only.grad(*point), atol=2e-13)
    np.testing.assert_allclose(
        zero.hessian(*point), pair_only.hessian(*point), atol=3e-12
    )


@pytest.mark.parametrize("path_id", (DIRECT_110, SHOCKLEY_112))
def test_full_energy_is_periodic_along_declared_scalar_path(path_id):
    model = FullFCC111StackEnergy.from_reduced_model(
        build_energy_model(TWO_ROW_LJ_REFERENCE), path_id=path_id
    )
    period = registry_path(path_id).period_over_b * model.p.b
    a = 1.03 * model.a0
    left = model.evaluate(a, 0.123)
    right = model.evaluate(a, 0.123 + period)
    assert left.energy == pytest.approx(right.energy, abs=2e-12)
    assert left.d_da == pytest.approx(right.d_da, abs=2e-11)
    assert left.d_ds == pytest.approx(right.d_ds, abs=2e-11)


@pytest.mark.parametrize("model_name", (TWO_ROW_LJ_REFERENCE, ANALYTIC_LJ_EAM_HYPOTHETICAL))
def test_full_analytic_gradient_and_hessian_match_finite_difference(model_name):
    model = FullFCC111StackEnergy.from_reduced_model(build_energy_model(model_name))
    a = 1.03 * model.a0
    s = 0.137
    step = 1.0e-4
    value = model.energy
    result = model.evaluate(a, s)
    grad_a = (value(a + step, s) - value(a - step, s)) / (2.0 * step)
    grad_s = (value(a, s + step) - value(a, s - step)) / (2.0 * step)
    h_aa = (value(a + step, s) - 2.0 * value(a, s) + value(a - step, s)) / step**2
    h_ss = (value(a, s + step) - 2.0 * value(a, s) + value(a, s - step)) / step**2
    h_as = (
        value(a + step, s + step)
        - value(a + step, s - step)
        - value(a - step, s + step)
        + value(a - step, s - step)
    ) / (4.0 * step**2)
    assert result.d_da == pytest.approx(grad_a, rel=3e-6, abs=2e-7)
    assert result.d_ds == pytest.approx(grad_s, rel=3e-6, abs=2e-7)
    assert result.d2_daa == pytest.approx(h_aa, rel=3e-6, abs=1e-4)
    assert result.d2_das == pytest.approx(h_as, rel=3e-6, abs=1e-4)
    assert result.d2_dss == pytest.approx(h_ss, rel=3e-6, abs=1e-4)


def test_full_reciprocal_energy_agrees_with_independent_direct_sum(full_hypothetical):
    a = 1.03 * full_hypothetical.a0
    s = 0.137
    analytic = full_hypothetical.evaluate(a, s)
    direct = full_hypothetical.direct_reference(
        a, s, radial_index=60, layers=60
    )
    assert abs(analytic.energy - direct.energy) <= 1.2 * direct.estimated_pair_tail
    assert analytic.d_da == pytest.approx(direct.d_da, rel=5e-5, abs=3e-6)
    assert analytic.d_ds == pytest.approx(direct.d_ds, rel=5e-7, abs=3e-9)
    assert analytic.d2_daa == pytest.approx(direct.d2_daa, rel=5e-5, abs=1e-5)
    assert analytic.d2_das == pytest.approx(direct.d2_das, rel=2e-6, abs=3e-9)
    assert analytic.d2_dss == pytest.approx(direct.d2_dss, rel=5e-6, abs=2e-6)


def test_best_feasible_reduced_fit_is_unbound_in_full_stack_geometry():
    model = FullFCC111StackEnergy.from_reduced_model(
        build_energy_model(AL_TARGET_BEST_FEASIBLE),
        require_stable_equilibrium=False,
    )
    assert math.isnan(model.a0)
    assert model.equilibrium_error == (
        "full FCC surface has no stable zero-load normal root"
    )
    probes = np.array([0.5, 0.8, 1.2, 2.0, 4.0])
    assert np.all([model.local_deda(float(a), 0.0) < 0.0 for a in probes])
