"""ARCHIVAL replay of the initial, underidentified six-observation experiment.

This is NOT the current calibration command.  Use
python -m solver_v1.run_full_fcc_calibration_audit to actually rerun fits.
The saved vectors omitted an independent shear measurement.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

from .aluminum_calibration import cell_ev_to_energy_density, energy_density_to_cell_ev
from .aluminum_full_fcc_calibration import (
    FullFCCParameters,
    evaluate_full_fcc_bulk_observables,
    full_fcc_bulk_targets,
    full_fcc_equilibrium_near_target,
    sensitivity_jacobian,
)
from .energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
)
from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_full_energy import FullFCC111StackEnergy
from .fcc111_geometry import DIRECT_110, SHOCKLEY_112, registry_path


SQUARE_STAGE1 = FullFCCParameters(
    0.10725542771610395, 0.9545342926175423, 0.7913535325526708,
    2.503712142830092, 0.0,
)
SQUARE_JOINT = FullFCCParameters(
    0.0038583102832933933, 1.1681262990758785, 0.5122315833354125,
    3.68017690186849, 0.0,
)
LINEAR_STAGE1 = FullFCCParameters(
    0.09852552212667592, 0.9288532190445209, 1.4007942968283327,
    4.4887389681264835, 1.971890325238724,
)
LINEAR_CANDIDATES = (
    FullFCCParameters(
        0.051463163185591194, 0.9278241602685452, 4.6853176139362835,
        5.679412614403284, 2.760290609784449,
    ),
    FullFCCParameters(
        0.013645377854871132, 1.02668545888573, 4.711004774121304,
        5.655895464976975, 2.30373517720161,
    ),
)


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _bulk_fit_rows(name: str, parameters: FullFCCParameters) -> tuple[dict, list[dict], list[dict], list[dict]]:
    targets = full_fcc_bulk_targets()
    obs = evaluate_full_fcc_bulk_observables(parameters)
    residual = (obs.values()-targets.values())/targets.residual_scales
    extended = parameters.embedding_linear_ev != 0.0
    # Two cubic identities remove duplicate measurements. Keep the old
    # six-observation residual table as historical data, not evidence of rank 5.
    jac = sensitivity_jacobian(parameters, extended=extended)[[0,2,3,5]]
    _, singular, right = np.linalg.svd(jac, full_matrices=False)
    rank = int(np.linalg.matrix_rank(jac))
    condition = (float(singular[0]/singular[-1])
                 if rank==jac.shape[1] else float("inf"))
    model, a0 = full_fcc_equilibrium_near_target(parameters)
    parameter_row = {
        "set": name, "epsilon_LJ_eV": parameters.epsilon_lj_ev,
        "sigma_LJ_over_b": parameters.sigma_lj_over_b,
        "beta_b_over_re": parameters.density_decay_b,
        "A_eV": parameters.embedding_sqrt_ev,
        "B_eV": parameters.embedding_linear_ev,
        "rho0_gauge": 1.0, "rho_ref_gauge": "rho_full_at_ideal_FCC",
        "a0_over_b": a0, "objective": 0.5*float(residual@residual),
        "rank": rank, "condition_number": condition,
    }
    residual_rows = [
        {"set": name, "observable": target_name, "target": target,
         "prediction": prediction, "scale": scale, "normalized_residual": value}
        for target_name, target, prediction, scale, value in zip(
            targets.names, targets.values(), obs.values(), targets.residual_scales, residual
        )
    ]
    parameter_names = ["log_epsilon", "log_sigma", "log_decay", "log_A"]
    if extended:
        parameter_names.append("log_B")
    svd_rows = []
    for index, value in enumerate(singular):
        row = {"set": name, "mode": index+1, "singular_value": value,
               "rank": rank, "condition_number": condition}
        row.update({key: right[index,j] for j,key in enumerate(parameter_names)})
        svd_rows.append(row)
    covariance = np.linalg.pinv(jac.T@jac)
    scale = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    correlation = covariance/np.outer(scale,scale)
    correlation_rows = [
        {"set": name, "parameter_i": parameter_names[i],
         "parameter_j": parameter_names[j], "correlation": correlation[i,j]}
        for i in range(len(parameter_names)) for j in range(i+1,len(parameter_names))
    ]
    return parameter_row, residual_rows, svd_rows, correlation_rows


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    calibration_dir = root/"results"/"aluminum_full_fcc_calibration"/"legacy_replay"
    interface_dir = root/"results"/"fcc111_active_interface"/"legacy_replay"
    calibration_dir.mkdir(parents=True, exist_ok=True)
    interface_dir.mkdir(parents=True, exist_ok=True)

    source = root/"solver_v1"/"data"/"aluminum_reference_targets.csv"
    with source.open(encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle))
    roles = {
        "mishin_lattice": ("fit", "bulk equilibrium; b and h111"),
        "mishin_cohesion": ("fit", "full bulk atomization energy per atom"),
        "mishin_c11": ("fit", "combined into full-FCC strain Hessian"),
        "mishin_c12": ("fit", "combined into full-FCC strain Hessian"),
        "mishin_c44": ("fit", "combined into full-FCC strain Hessian"),
        "lu_usf110": ("held_out", "active-interface direct_110 GSF"),
        "lu_usf112": ("held_out", "active-interface Shockley path"),
        "lu_isf": ("held_out", "active-interface Shockley endpoint"),
        "wei_surface111": ("held_out", "active-interface work of separation"),
        "mishin_surface111": ("held_out", "alternate consistent-EAM surface target"),
        "mishin_vacancy": ("unrepresentable", "no vacancy coordinate/state"),
    }
    target_rows = []
    for row in source_rows:
        if row["target_id"] not in roles:
            continue
        role, mapping = roles[row["target_id"]]
        target_rows.append({**row, "full_fcc_role": role, "full_fcc_mapping": mapping})
    _write(calibration_dir/"calibration_targets.csv", target_rows)

    parameter_rows = [
        {"set": "square_stage1_equilibrium_cohesion", **{
            "epsilon_LJ_eV": SQUARE_STAGE1.epsilon_lj_ev,
            "sigma_LJ_over_b": SQUARE_STAGE1.sigma_lj_over_b,
            "beta_b_over_re": SQUARE_STAGE1.density_decay_b,
            "A_eV": SQUARE_STAGE1.embedding_sqrt_ev,
            "B_eV": 0.0, "rho0_gauge": 1.0,
            "rho_ref_gauge": "rho_full_at_ideal_FCC", "a0_over_b": "",
            "objective": 2.5383352883284335e-9, "rank": 2,
            "condition_number": float("inf")}},
        {"set": "linear_stage1_equilibrium_cohesion", **{
            "epsilon_LJ_eV": LINEAR_STAGE1.epsilon_lj_ev,
            "sigma_LJ_over_b": LINEAR_STAGE1.sigma_lj_over_b,
            "beta_b_over_re": LINEAR_STAGE1.density_decay_b,
            "A_eV": LINEAR_STAGE1.embedding_sqrt_ev,
            "B_eV": LINEAR_STAGE1.embedding_linear_ev, "rho0_gauge": 1.0,
            "rho_ref_gauge": "rho_full_at_ideal_FCC", "a0_over_b": "",
            "objective": 1.0836113400686468e-9, "rank": 2,
            "condition_number": float("inf")}},
    ]
    residual_rows: list[dict] = []
    svd_rows: list[dict] = []
    correlation_rows: list[dict] = []
    targets = full_fcc_bulk_targets()
    for name, parameters, singular in (
        ("square_stage1_equilibrium_cohesion", SQUARE_STAGE1,
         (971.320448, 38.3349062, 1.06281225e-4)),
        ("linear_stage1_equilibrium_cohesion", LINEAR_STAGE1,
         (572.203138, 72.0645187, 9.57308944e-4)),
    ):
        prediction = evaluate_full_fcc_bulk_observables(parameters)
        residual = (prediction.values()-targets.values())/targets.residual_scales
        residual_rows.extend(
            {"set": name, "observable": observable, "target": target,
             "prediction": value, "scale": scale, "normalized_residual": normalized}
            for observable,target,value,scale,normalized in zip(
                targets.names,targets.values(),prediction.values(),targets.residual_scales,residual
            )
        )
        svd_rows.extend(
            {"set": name, "mode": index+1, "singular_value": value,
             "rank": 2, "condition_number": float("inf"),
             "note":"historical numerical singular values include derivative noise"}
            for index,value in enumerate(singular)
        )
    for name, parameters in (
        ("square_joint", SQUARE_JOINT),
        ("linear_joint_start0", LINEAR_CANDIDATES[0]),
        ("linear_joint_start1", LINEAR_CANDIDATES[1]),
    ):
        p_row, r_rows, s_rows, c_rows = _bulk_fit_rows(name, parameters)
        parameter_rows.append(p_row)
        residual_rows.extend(r_rows)
        svd_rows.extend(s_rows)
        correlation_rows.extend(c_rows)
    _write(calibration_dir/"parameter_sets.csv", parameter_rows)
    _write(calibration_dir/"normalized_residuals.csv", residual_rows)
    _write(calibration_dir/"identifiability_svd.csv", svd_rows)
    _write(calibration_dir/"parameter_correlations.csv", correlation_rows)

    area = 7.102490842787127
    heldout_rows = []
    model_comparison = []
    energy_grid = []
    gsf_rows = []
    opening_rows = []
    hessian_rows = []
    barrier_rows = []
    convergence_rows = []
    for candidate_index, parameters in enumerate(LINEAR_CANDIDATES):
        name = f"linear_joint_start{candidate_index}"
        bulk, a0 = full_fcc_equilibrium_near_target(parameters)
        interfaces = {
            path: FCC111ActiveInterface(bulk, path_id=path, tolerance=2e-8)
            for path in (DIRECT_110, SHOCKLEY_112)
        }
        direct = interfaces[DIRECT_110]
        shockley = interfaces[SHOCKLEY_112]
        direct_x = np.linspace(0.0, direct.period, 61)
        direct_e = np.array([direct.energy(direct.h, float(s)) for s in direct_x])
        partial = registry_path(SHOCKLEY_112).partial_increment_over_b
        assert partial is not None
        shockley_x = np.linspace(0.0, partial*bulk.p.b, 61)
        shockley_e = np.array([shockley.energy(shockley.h, float(s)) for s in shockley_x])
        opening_x = np.linspace(direct.h, 30.0*direct.h, 61)
        opening_e = np.array([direct.energy(float(a), 0.0) for a in opening_x])
        opening_force = np.array([direct.evaluate(float(a), 0.0).d_da for a in opening_x])
        direct_usf = float(np.max(direct_e))
        shockley_usf = float(np.max(shockley_e))
        isf = float(shockley_e[-1])
        work_sep = direct.separated_limit()
        traction_peak = minimize_scalar(
            lambda a: -direct.evaluate(float(a), 0.0).d_da,
            bounds=(direct.h, 3.0*direct.h), method="bounded",
            options={"xatol": 2e-8},
        )
        predictions = (
            ("direct_110_usf", direct_usf, energy_density_to_cell_ev(0.250, area)),
            ("shockley_112_usf", shockley_usf, energy_density_to_cell_ev(0.224, area)),
            ("shockley_isf", isf, energy_density_to_cell_ev(0.164, area)),
            ("work_separation_DFT", work_sep, energy_density_to_cell_ev(2.12, area)),
            ("work_separation_Mishin", work_sep, energy_density_to_cell_ev(1.74, area)),
        )
        for observable, prediction, target in predictions:
            heldout_rows.append(
                {"set": name, "observable": observable, "prediction_eV_cell": prediction,
                 "target_eV_cell": target, "prediction_J_m2": cell_ev_to_energy_density(prediction,area),
                 "target_J_m2": cell_ev_to_energy_density(target,area),
                 "relative_error": (prediction-target)/target, "used_for_fit": "no"}
            )
        hessian = direct.hessian(direct.h, 0.0)
        hessian_rows.append(
            {"set": name, "W_aa": hessian[0,0], "W_as": hessian[0,1],
             "W_ss": hessian[1,1], "eig_min": np.min(np.linalg.eigvalsh(hessian)),
             "eig_max": np.max(np.linalg.eigvalsh(hessian))}
        )
        barrier_rows.append(
            {"set": name, "direct_110_barrier_eV_cell": direct_usf,
             "shockley_112_barrier_eV_cell": shockley_usf,
             "shockley_endpoint_eV_cell": isf,
             "opening_work_eV_cell": work_sep,
             "ideal_opening_traction_eV_per_b": float(-traction_peak.fun),
             "ideal_opening_traction_a_over_b": float(traction_peak.x),
             "status": "load-free static; dynamic ordering not inferred"}
        )
        for path, coordinates, energies in (
            (DIRECT_110, direct_x, direct_e), (SHOCKLEY_112, shockley_x, shockley_e)
        ):
            for s, energy in zip(coordinates, energies):
                gsf_rows.append(
                    {"set": name, "path": path, "s_over_b": s/bulk.p.b,
                     "energy_eV_cell": energy,
                     "energy_J_m2": cell_ev_to_energy_density(energy, area)}
                )
        for a, energy, traction in zip(opening_x, opening_e, opening_force):
            opening_rows.append(
                {"set": name, "a_over_b": a, "opening_eV_cell": energy,
                 "opening_J_m2": cell_ev_to_energy_density(energy, area),
                 "traction_eV_per_b": traction}
            )
        test_a, test_s = 1.05*direct.h, 0.17
        reciprocal = direct.evaluate(test_a, test_s)
        real = direct.direct_reference(test_a, test_s, radial_index=50, layers=40)
        for field in ("energy","d_da","d_ds","d2_daa","d2_das","d2_dss"):
            rv, dv = float(getattr(reciprocal,field)), float(getattr(real,field))
            convergence_rows.append(
                {"set": name, "quantity": field, "reciprocal": rv, "direct": dv,
                 "absolute_error": abs(rv-dv),
                 "relative_error": abs(rv-dv)/max(abs(dv),1e-300),
                 "reciprocal_layers": reciprocal.layers_used,
                 "direct_layers": real.layers_used, "direct_radius": 50}
            )
        model_comparison.append(
            {"model": name, "geometry": "newly_calibrated_full_FCC_bulk",
             "a0_over_b": a0, "bulk_cohesion_eV_atom":
                 evaluate_full_fcc_bulk_observables(parameters).cohesive_ev_atom,
             "direct_GSF_eV_cell": "not a homogeneous-bulk observable",
             "shockley_GSF_eV_cell": "not a homogeneous-bulk observable",
             "work_separation_eV_cell": "not a homogeneous-bulk observable",
             "status": "archival underidentified fit; independent cubic shear omitted"}
        )
        model_comparison.append(
            {"model": name, "geometry": "full_FCC_active_interface",
             "a0_over_b": a0, "bulk_cohesion_eV_atom": "inherited bulk candidate",
             "direct_GSF_eV_cell": direct_usf, "shockley_GSF_eV_cell": shockley_usf,
             "work_separation_eV_cell": work_sep,
             "status": "static interface held-out targets fail; not production-valid"}
        )
        if candidate_index == 0:
            a_values = np.linspace(0.92*direct.h, 2.0*direct.h, 23)
            s_values = np.linspace(0.0, direct.period, 31)
            for a in a_values:
                for s in s_values:
                    value = direct.evaluate(float(a), float(s))
                    energy_grid.append(
                        {"set": name, "a_over_b": a, "s_over_b": s,
                         "energy_eV_cell": value.energy, "dW_da": value.d_da,
                         "dW_ds": value.d_ds}
                    )

    for model_id in (TWO_ROW_LJ_REFERENCE, ANALYTIC_LJ_EAM_HYPOTHETICAL, AL_TARGET_BEST_FEASIBLE):
        reduced = build_energy_model(model_id)
        transferred = FullFCC111StackEnergy.from_reduced_model(
            reduced, require_stable_equilibrium=False
        )
        model_comparison.append(
            {"model": model_id, "geometry": "existing_reduced_model",
             "a0_over_b": reduced.a0,
             "bulk_cohesion_eV_atom": "not comparable from reduced normalization",
             "direct_GSF_eV_cell": "homogeneous/reduced barrier only",
             "shockley_GSF_eV_cell": "not represented",
             "work_separation_eV_cell": "homogeneous/reduced opening only",
             "status": "existing model preserved"}
        )
        model_comparison.append(
            {"model": model_id, "geometry": "same_parameters_full_FCC_transfer",
             "a0_over_b": transferred.a0 if np.isfinite(transferred.a0) else "",
             "bulk_cohesion_eV_atom": "not refitted",
             "direct_GSF_eV_cell": "homogeneous shear is not GSF",
             "shockley_GSF_eV_cell": "not evaluated",
             "work_separation_eV_cell": "uniform dilation is not cleavage",
             "status": ("stable but unrefitted transfer"
                        if np.isfinite(transferred.a0) else transferred.equilibrium_error)}
        )

    _write(calibration_dir/"heldout_validation.csv", heldout_rows)
    _write(calibration_dir/"model_comparison.csv", model_comparison)
    _write(interface_dir/"active_interface_energy_grid.csv", energy_grid)
    _write(interface_dir/"gsf_curve.csv", gsf_rows)
    _write(interface_dir/"opening_curve.csv", opening_rows)
    _write(interface_dir/"hessian.csv", hessian_rows)
    _write(interface_dir/"barrier_summary.csv", barrier_rows)
    _write(interface_dir/"convergence_summary.csv", convergence_rows)

    figure, axes = plt.subplots(1,2,figsize=(10,4))
    for name in ("linear_joint_start0","linear_joint_start1"):
        rows=[r for r in gsf_rows if r["set"]==name and r["path"]==DIRECT_110]
        axes[0].plot([r["s_over_b"] for r in rows],[r["energy_J_m2"] for r in rows],label=name)
        rows=[r for r in gsf_rows if r["set"]==name and r["path"]==SHOCKLEY_112]
        axes[1].plot([r["s_over_b"] for r in rows],[r["energy_J_m2"] for r in rows],label=name)
    axes[0].axhline(.250,color="k",ls="--",label="DFT direct USF")
    axes[1].axhline(.224,color="k",ls="--",label="DFT Shockley USF")
    axes[0].set(title="Direct <110> GSF",xlabel="s/b",ylabel="J/m^2")
    axes[1].set(title="Shockley segment",xlabel="s/b",ylabel="J/m^2")
    for ax in axes: ax.grid(alpha=.25); ax.legend(fontsize=7)
    figure.tight_layout();figure.savefig(interface_dir/"gsf_validation.png",dpi=180);plt.close(figure)

    rows=[r for r in opening_rows if r["set"]=="linear_joint_start0"]
    figure,ax=plt.subplots(figsize=(6,4));ax.plot([r["a_over_b"] for r in rows],[r["opening_J_m2"] for r in rows])
    ax.axhline(2.12,color="k",ls="--",label="DFT 2 gamma(111)")
    ax.set(xlabel="active spacing a/b",ylabel="W_int [J/m^2]",title="Active-interface opening")
    ax.grid(alpha=.25);ax.legend();figure.tight_layout();figure.savefig(interface_dir/"opening_validation.png",dpi=180);plt.close(figure)

    grid=np.asarray([[r["a_over_b"],r["s_over_b"],r["energy_eV_cell"]] for r in energy_grid])
    aa=np.unique(grid[:,0]);ss=np.unique(grid[:,1]);zz=grid[:,2].reshape(aa.size,ss.size)
    figure,ax=plt.subplots(figsize=(6,4));cont=ax.contourf(ss,aa,zz,levels=24);figure.colorbar(cont,ax=ax,label="eV/interface cell")
    ax.set(xlabel="s/b",ylabel="a/b",title="Full-FCC active-interface W_int")
    figure.tight_layout();figure.savefig(interface_dir/"active_interface_contour.png",dpi=180);plt.close(figure)

    summary={
        "bulk_calibration":"archival underidentified fit; shear omitted; superseded by audited_v2",
        "active_interface_validation":"fails GSF and decohesion held-out targets",
        "production_PDE_connected":False,
        "physical_time":"uncalibrated",
        "A_c_used":False,
    }
    (calibration_dir/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (interface_dir/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
