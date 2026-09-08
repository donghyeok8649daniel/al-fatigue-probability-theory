import numpy as np
import pytest

from solver_v1.analytic_lj_eam import AnalyticLJEAM
from solver_v1.energy_model_registry import (
    AL_BEST_FEASIBLE_PARAMETERS,
    AL_TARGET_BEST_FEASIBLE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    ENERGY_MODEL_IDS,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
    energy_model_metadata,
    energy_model_result_metadata,
)
from solver_v1.model import TwoRowLJ


def test_energy_models_are_explicitly_distinguishable() -> None:
    models = {model_id: build_energy_model(model_id) for model_id in ENERGY_MODEL_IDS}
    assert type(models[TWO_ROW_LJ_REFERENCE]) is TwoRowLJ
    assert isinstance(models[ANALYTIC_LJ_EAM_HYPOTHETICAL], AnalyticLJEAM)
    assert isinstance(models[AL_TARGET_BEST_FEASIBLE], AnalyticLJEAM)
    assert energy_model_metadata(TWO_ROW_LJ_REFERENCE).calibration_status.startswith(
        "dimensionless reference"
    )
    assert "hypothetical" in energy_model_metadata(
        ANALYTIC_LJ_EAM_HYPOTHETICAL
    ).calibration_status
    assert "best-feasible" in energy_model_metadata(
        AL_TARGET_BEST_FEASIBLE
    ).calibration_status


def test_lj_reference_static_calibration_is_unchanged() -> None:
    model = build_energy_model(TWO_ROW_LJ_REFERENCE)
    assert model.a0 == pytest.approx(0.7713438268704838, abs=2.0e-15)
    assert model.sigma_over_E_force_scale() == pytest.approx(
        86.29296488740997, rel=2.0e-14
    )
    assert model.p.chi_axial_projection == pytest.approx(0.20)


def test_best_feasible_factory_uses_stored_parameter_set() -> None:
    model = build_energy_model(AL_TARGET_BEST_FEASIBLE)
    assert model.a0 == pytest.approx(0.8164957221, rel=3.0e-8)
    assert model.p.epsilon == pytest.approx(AL_BEST_FEASIBLE_PARAMETERS.epsilon_lj_ev)
    assert model.p.sigma_lj == pytest.approx(
        AL_BEST_FEASIBLE_PARAMETERS.sigma_lj_over_b
    )
    assert model.embedding.linear == pytest.approx(AL_BEST_FEASIBLE_PARAMETERS.linear_ev)
    assert abs(model.local_deda(model.a0, 0.0)) < 1.0e-9


def test_result_metadata_cannot_label_lj_as_calibrated_hybrid() -> None:
    model = build_energy_model(TWO_ROW_LJ_REFERENCE)
    metadata = energy_model_result_metadata(TWO_ROW_LJ_REFERENCE, model)
    assert metadata["energy_model_python_class"].endswith("TwoRowLJ")
    assert "not calibrated Al" in metadata["energy_model_calibration_status"]
    assert metadata["energy_model_id"] != AL_TARGET_BEST_FEASIBLE
    assert np.isfinite(metadata["energy_model_kappa_axial"])


def test_unknown_energy_model_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown energy model"):
        build_energy_model("calibrated_al_by_filename_only")
