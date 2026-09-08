from types import SimpleNamespace

import numpy as np
import pytest

from solver_v1.plasticity_audit import (
    final_well_populations,
    outside_central_well_mass,
    summarize_plasticity_run,
)


def _audit_result() -> dict[str, object]:
    time_rows = 3
    well_indices = np.array([-1, 0, 1])
    populations = np.array(
        [
            [0.0, 1.0, 0.0],
            [0.01, 0.97, 0.01],
            [0.02, 0.94, 0.03],
        ]
    )
    cumulative_net = np.array(
        [[0.0, 0.0], [-0.01, 0.01], [-0.02, 0.03]]
    )
    cumulative_gross = np.array(
        [[0.0, 0.0], [0.01, 0.01], [0.02, 0.03]]
    )
    return {
        "model": SimpleNamespace(p=SimpleNamespace(chi_axial_projection=0.2)),
        "grid": SimpleNamespace(a=np.arange(5), s=np.arange(7)),
        "well_indices": well_indices,
        "well_populations": populations,
        "plastic_strain": np.array([0.0, 1.0e-4, 2.0e-4]),
        "intrawell_strain": np.array([0.0, 2.0e-4, 3.0e-4]),
        "strain_decomposition_residual": np.array([0.0, 2.0e-17, -1.0e-17]),
        "unnormalized_registry_moment": populations @ well_indices,
        "cumulative_interwell_net_transfer": cumulative_net,
        "cumulative_interwell_gross_transfer": cumulative_gross,
        "cumulative_interwell_forward_transfer": 0.5
        * (cumulative_gross + cumulative_net),
        "cumulative_interwell_backward_transfer": 0.5
        * (cumulative_gross - cumulative_net),
        "accumulated_net_registry_transfer": np.sum(cumulative_net, axis=1),
        "cumulative_absorbed_mass": np.array([0.0, 0.01, 0.01]),
        "mass_balance_residual": np.array([0.0, 3.0e-15, -2.0e-15]),
        "well_population_balance_residual": np.zeros((time_rows, 3)),
        "cumulative_negative_mass_correction": np.zeros(time_rows),
        "flux_consistency_residual": np.zeros(time_rows),
        "interwell_boundary_alignment_error": np.array([0.0, 1.0e-16]),
    }


def test_audit_summary_preserves_net_gross_and_numerical_residuals():
    result = _audit_result()
    summary = summarize_plasticity_run(
        result,
        label="reference",
        dt=1.0e-3,
        integrator="implicit",
        energy_model="two_row_lj_reference",
    )

    assert summary.grid_a == 5
    assert summary.grid_s == 7
    assert summary.max_epsilon_p == pytest.approx(2.0e-4)
    assert summary.max_decomposition_residual == pytest.approx(2.0e-17)
    assert summary.max_outside_central_well_mass == pytest.approx(0.05)
    assert summary.cumulative_forward_transfer == pytest.approx(0.03)
    assert summary.cumulative_backward_transfer == pytest.approx(0.02)
    assert summary.cumulative_net_registry_transfer == pytest.approx(0.01)
    assert summary.cumulative_gross_registry_activity == pytest.approx(0.05)
    assert summary.mass_residual == pytest.approx(3.0e-15)
    assert summary.well_balance_residual == 0.0


def test_audit_well_helpers_do_not_renormalize_survivor_mass():
    result = _audit_result()

    np.testing.assert_allclose(
        outside_central_well_mass(result), np.array([0.0, 0.02, 0.05])
    )
    assert final_well_populations(result) == pytest.approx(
        {-1: 0.02, 0: 0.94, 1: 0.03}
    )
    assert sum(final_well_populations(result).values()) == pytest.approx(0.99)


def test_saved_audit_metadata_identifies_prior_ui_model_and_time_status():
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    metadata = json.loads(
        (root / "results" / "plasticity_audit" / "active_model_metadata.json")
        .read_text(encoding="utf-8")
    )

    assert metadata["previous_ui_default"] == "two_row_lj_reference"
    assert metadata["two_row_lj_reference"]["python_class"].endswith("TwoRowLJ")
    assert "not calibrated Al" in metadata["two_row_lj_reference"][
        "calibration_status"
    ]
    assert metadata["time_basis"] == "model"
    assert metadata["physical_seconds_available"] is False
    assert metadata["physical_hz_available"] is False


def test_saved_aligned_refinement_reports_real_values_without_rescaling():
    import csv
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    with (
        root / "results" / "plasticity_audit" / "refinement_summary.csv"
    ).open(encoding="utf-8", newline="") as stream:
        rows = {row["label"]: row for row in csv.DictReader(stream)}

    finest = rows["G_aligned60_dt"]
    assert int(finest["grid_a"]) == 41
    assert int(finest["grid_s"]) == 60
    assert float(finest["well_boundary_alignment_error"]) < 1.0e-12
    assert float(finest["max_epsilon_p"]) == pytest.approx(
        2.61234844327e-12, rel=1.0e-12
    )
    assert float(finest["cumulative_gross_registry_activity"]) > float(
        finest["cumulative_net_registry_transfer"]
    )
    assert float(finest["mass_residual"]) < 1.0e-13
    assert finest["plasticity_resolution_status"] == "resolved-tiny-model-flow"


def test_saved_well_populations_and_resolution_classifications_are_consistent():
    import csv
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    result_dir = root / "results" / "plasticity_audit"
    with (result_dir / "well_population_summary.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        populations = {
            int(row["well_index"]): float(row["final_absolute_intact_mass"])
            for row in csv.DictReader(stream)
        }
    classification = json.loads(
        (result_dir / "classification_summary.json").read_text(encoding="utf-8")
    )

    assert set(populations) == {-1, 0, 1}
    assert sum(populations.values()) == pytest.approx(
        1.0
        - classification["finest_run"]["opening_absorbed_mass"]
        + 3.21964677141e-15,
        abs=2.0e-15,
    )
    plasticity = classification["plasticity_resolution"]
    opening = classification["opening_resolution"]
    assert plasticity["signal_to_floor_ratio"] > 1.0
    assert "resolved_extremely_small" in plasticity["classification"]
    assert opening["signal_to_floor_ratio"] < 1.0
    assert opening["classification"].startswith("unresolved")
    assert classification["unload_hold"]["classification"].endswith(
        "not_converged"
    )
