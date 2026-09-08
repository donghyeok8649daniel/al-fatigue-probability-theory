import math

import numpy as np
import pytest

from solver_v1.aluminum_full_fcc_calibration import (
    FullFCCParameters,
    build_full_fcc_calibration_model,
    calibrate_full_fcc,
    evaluate_full_fcc_bulk_observables,
    full_fcc_bulk_targets,
    full_fcc_equilibrium_near_target,
    sensitivity_jacobian,
)
from solver_v1.energy_model_registry import (
    AL_BEST_FEASIBLE_PARAMETERS,
    AL_TARGET_BEST_FEASIBLE,
    build_energy_model,
)


CANDIDATE = FullFCCParameters(
    0.051463163185591194, 0.9278241602685452, 4.6853176139362835,
    5.679412614403284, 2.760290609784449,
)


def test_full_fcc_elastic_target_mapping_from_cubic_constants():
    target = full_fcc_bulk_targets()
    assert target.a0_over_b == pytest.approx(math.sqrt(2/3))
    assert target.cohesive_ev_atom == pytest.approx(3.36)
    assert target.strain_h_aa_ev == pytest.approx(12.646038954154417)
    assert target.strain_h_bb_ev == pytest.approx(37.31618052045565)
    assert target.strain_h_ab_ev == pytest.approx(12.024102612146821)


def test_full_fcc_candidate_matches_bulk_targets_and_is_stable():
    target = full_fcc_bulk_targets()
    value = evaluate_full_fcc_bulk_observables(CANDIDATE)
    residual = (value.values()-target.values())/target.residual_scales
    assert np.max(np.abs(residual)) < 4e-5
    hessian = np.array([[value.strain_h_aa_ev, value.strain_h_ab_ev],
                        [value.strain_h_ab_ev, value.strain_h_bb_ev]])
    assert np.min(np.linalg.eigvalsh(hessian)) > 0
    model, a0 = full_fcc_equilibrium_near_target(CANDIDATE)
    assert a0 == pytest.approx(target.a0_over_b, abs=1e-6)
    assert abs(model.local_deda(a0, 0.0)) < 1e-9
    atomized_energy = -CANDIDATE.embedding_linear_ev
    assert atomized_energy-model.energy(a0, 0.0) == pytest.approx(3.36, abs=2e-5)


def test_full_fcc_density_gauge_sets_target_environment_to_one():
    model = build_full_fcc_calibration_model(CANDIDATE, find_equilibrium=False)
    rho = model.full_environment_density(math.sqrt(2/3), 0).value
    assert model.embedding is not None
    assert rho/model.embedding.rho_ref == pytest.approx(1.0, rel=2e-10)


def test_legacy_elastic_observations_are_structurally_dependent():
    target = full_fcc_bulk_targets()
    assert target.strain_h_bb_ev == pytest.approx(
        2*target.strain_h_aa_ev + target.strain_h_ab_ev)
    value = evaluate_full_fcc_bulk_observables(CANDIDATE, strain_step=2e-4)
    assert value.lateral_force == pytest.approx(2*value.normal_force, abs=2e-5)
    assert value.strain_h_bb_ev == pytest.approx(
        2*value.strain_h_aa_ev + value.strain_h_ab_ev, abs=5e-5)
    jacobian = sensitivity_jacobian(CANDIDATE, extended=True)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    # A fifth floating-point singular value does NOT create a fifth observable.
    assert singular[0]/singular[-1] > 1e4


def test_old_reduced_calibration_is_unchanged():
    assert AL_BEST_FEASIBLE_PARAMETERS.epsilon_lj_ev == pytest.approx(0.1551326727)
    model = build_energy_model(AL_TARGET_BEST_FEASIBLE)
    assert model.a0 == pytest.approx(0.8164957221838465)
    assert model.sigma_over_E_force_scale() == pytest.approx(13.855051447120799)


def test_deterministic_calibration_restarts_reproducibly_at_solution():
    fit = calibrate_full_fcc(
        extended=True, starts=(CANDIDATE,), max_nfev=2
    )[0]
    assert fit.cost < 2e-8
    assert fit.rank <= 4
    assert math.isinf(fit.condition_number)
    np.testing.assert_allclose(
        fit.parameters.as_array(extended=True),
        CANDIDATE.as_array(extended=True), rtol=2e-3,
    )
