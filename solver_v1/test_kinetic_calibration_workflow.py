"""All kinetic numbers in this file are SYNTHETIC verification, never Al data."""
from dataclasses import replace
import json
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from .energy_model_registry import TWO_ROW_LJ_REFERENCE, AL_TARGET_BEST_FEASIBLE, build_energy_model
from .kinetic_calibration_workflow import (
    bind_to_energy_model, build_time_basis_model, calibration_from_collective_correlations,
    infer_from_correlation_matrices, validate_for_energy_model,
)
from .physical_time import (
    BOLTZMANN_J_PER_K as KB, calibration_from_mobilities, load_time_calibration,
    model_time_to_seconds,
)


def synthetic_correlation(H, ma=2., ms=.3, temperature=300.):
    t = np.array([0., .015, .03, .06])
    M = np.diag([ma, ms])
    C0 = KB*temperature*np.linalg.inv(H)
    C = np.array([expm(-M@H*lag)@C0 for lag in t])
    return t,C


def test_coupled_covariance_identifies_mobility_without_swapped_rate_branch():
    H = np.array([[12.,2.],[2.,5.]])
    t,C = synthetic_correlation(H)
    fitted = infer_from_correlation_matrices(t,C,H,temperature_K=300.,relative_tolerance=1e-8)
    assert fitted["M_a_phys_m2_per_J_s"] == pytest.approx(2.,rel=2e-13)
    assert fitted["M_s_phys_m2_per_J_s"] == pytest.approx(.3,rel=2e-13)
    assert fitted["lag_rate_relative_error"] < 1e-12


@pytest.mark.parametrize("failure", ["temperature", "lag", "offdiagonal", "oscillatory"])
def test_incompatible_correlation_cannot_unlock_clock(failure):
    H = np.array([[12.,2.],[2.,5.]])
    t,C = synthetic_correlation(H)
    if failure == "temperature":
        C *= 1.5
    if failure == "lag":
        C[-1] *= .9
    if failure == "offdiagonal":
        M = np.array([[2.,.4],[.4,.3]])
        C = np.array([expm(-M@H*lag)@C[0] for lag in t])
    if failure == "oscillatory":
        C[1] *= -1.
    with pytest.raises(ValueError):
        infer_from_correlation_matrices(t,C,H,temperature_K=300.,relative_tolerance=.01)


def test_clock_pipeline_preserves_independently_measured_mobility_ratio():
    model = build_energy_model(TWO_ROW_LJ_REFERENCE)
    L,E = 2.8e-10,1.6e-19
    H = model.hessian(model.a0,0.)*E/L**2
    t,C = synthetic_correlation(H)
    c,fit = calibration_from_collective_correlations(t,C,model_id=TWO_ROW_LJ_REFERENCE,
        length_scale_m=L,energy_scale_J=E,temperature_K=300.,
        source="SYNTHETIC exact matrix regression only",relative_tolerance=1e-8)
    assert c.model_mobility_s == pytest.approx(.15)
    assert c.t0_seconds == pytest.approx(L**2/(2.*E))
    physical = build_time_basis_model(TWO_ROW_LJ_REFERENCE,time_basis="physical",calibration=c)
    assert physical.p.mobility_s == pytest.approx(.15)
    assert physical.p.kT == pytest.approx(KB*300/E)
    np.testing.assert_array_equal(physical.hessian(model.a0,0),model.hessian(model.a0,0))
    reference = build_time_basis_model(TWO_ROW_LJ_REFERENCE,calibration=c)
    assert reference.p.mobility_s == .05
    assert reference.p.kT == .02


def test_uncertain_or_stale_model_and_invalid_units_refused():
    c = calibration_from_mobilities(length_scale_m=2.8e-10,energy_scale_J=1.6e-19,
        M_a_phys_m2_per_J_s=2.,M_s_phys_m2_per_J_s=.1,
        source="hypothetical unit test",temperature_K=300.)
    with pytest.raises(ValueError,match="energy_model_id"):
        validate_for_energy_model(c,TWO_ROW_LJ_REFERENCE)
    c = bind_to_energy_model(c,TWO_ROW_LJ_REFERENCE)
    for broken in (replace(c,energy_model_fingerprint="stale"),
                   replace(c,coordinate_definition="dislocation centre X"),
                   replace(c,temperature_K=None)):
        with pytest.raises(ValueError):
            validate_for_energy_model(broken,TWO_ROW_LJ_REFERENCE)
    with pytest.raises(ValueError):
        bind_to_energy_model(c,AL_TARGET_BEST_FEASIBLE)
    with pytest.raises(ValueError):
        model_time_to_seconds(1.,replace(c,model_mobility_s=float("nan")))


def test_actual_cli_import_reproducible_and_no_overwrite(tmp_path):
    model = build_energy_model(TWO_ROW_LJ_REFERENCE)
    L,E = 2.8e-10,1.6e-19
    t,C = synthetic_correlation(model.hessian(model.a0,0)*E/L**2)
    meta = dict(model_id=TWO_ROW_LJ_REFERENCE,length_scale_m=L,energy_scale_J=E,
        temperature_K=300.,source="SYNTHETIC CLI test",relative_tolerance=1e-8)
    (tmp_path/"meta.json").write_text(json.dumps(meta),encoding="utf-8")
    np.savetxt(tmp_path/"C.csv",np.column_stack([t,C.reshape(-1,4)]),delimiter=",",
        header="time_seconds,C_aa_m2,C_as_m2,C_sa_m2,C_ss_m2",comments="")
    command = [sys.executable,"-m","solver_v1.import_kinetic_calibration",
        str(tmp_path/"meta.json"),str(tmp_path/"C.csv"),"--output",str(tmp_path/"out.json"),
        "--diagnostics",str(tmp_path/"fit.json")]
    subprocess.run(command,check=True,capture_output=True)
    calibration = load_time_calibration(tmp_path/"out.json")
    assert calibration.calibrated
    assert calibration.model_mobility_s == pytest.approx(.15)
    before = (tmp_path/"out.json").read_bytes()
    assert subprocess.run(command,capture_output=True).returncode != 0
    assert (tmp_path/"out.json").read_bytes() == before
