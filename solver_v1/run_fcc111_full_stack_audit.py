"""Generate the deterministic FCC(111) full-stack reference audit tables."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar

from .energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
    energy_model_metadata,
)
from .fcc111_full_energy import FullFCC111StackEnergy
from .fcc111_geometry import DIRECT_110, SHOCKLEY_112, fcc111_geometry_from_b, registry_path
from .fcc111_landscape import (
    configurational_spinodal,
    opening_barrier,
    opening_spinodal,
    relaxed_registry_barrier,
)
from .fcc111_lattice_sum import (
    ReciprocalSumConfig,
    plane_exponential_sum_direct,
    plane_exponential_sum_reciprocal,
    plane_power_sum_direct,
    plane_power_sum_reciprocal,
)


MODEL_IDS = (
    TWO_ROW_LJ_REFERENCE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    AL_TARGET_BEST_FEASIBLE,
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _reduced_relaxed_barrier(model) -> float:
    values = []
    for s in (0.0, 0.5 * model.p.b):
        result = minimize_scalar(
            lambda a: model.local_energy(float(a), s),
            bounds=(max(model.p.a_min, 0.45 * model.a0), min(model.p.a_max, 2.5)),
            method="bounded",
            options={"xatol": 2e-11},
        )
        values.append(float(result.fun))
    return values[1] - values[0]


def run(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    geometry = fcc111_geometry_from_b(1.0)
    geometry_rows = [
        {
            "quantity": key,
            "value": value,
            "units": units,
            "note": note,
        }
        for key, value, units, note in (
            ("a_lat", geometry.lattice_constant, "b units", "sqrt(2)*b"),
            ("b", geometry.b, "b units", "in-plane nearest-neighbor spacing"),
            ("h111", geometry.h111, "b units", "a_lat/sqrt(3)"),
            ("A_atomic_cell", geometry.atomic_cell_area, "b^2", "not statistical A_c"),
            ("tau_x", geometry.tau[0], "b units", "ABC shift component"),
            ("tau_y", geometry.tau[1], "b units", "ABC shift component"),
            ("direct_110_period", registry_path(DIRECT_110).period_over_b, "b", "current calibrated scalar path"),
            ("shockley_112_period", registry_path(SHOCKLEY_112).period_over_b, "b", "held-out configurable path"),
        )
    ]
    _write_csv(output_dir / "geometry_summary.csv", geometry_rows)

    validation_rows: list[dict[str, object]] = []
    delta = np.array([0.137, -0.081])
    config = ReciprocalSumConfig(tol=2e-13, max_shell_index=32)
    for p in (3.0, 6.0):
        reciprocal = plane_power_sum_reciprocal(
            0.91, delta, p=p, geometry=geometry, config=config
        )
        direct = plane_power_sum_direct(
            0.91, delta, p=p, geometry=geometry, radius_index=300
        )
        for field in ("value", "d_d", "d2_dd"):
            left, right = float(getattr(reciprocal, field)), float(getattr(direct, field))
            validation_rows.append(
                {"kernel": f"power_p{int(p)}", "quantity": field,
                 "reciprocal": left, "direct": right, "absolute_error": abs(left-right),
                 "relative_error": abs(left-right)/max(abs(right), 1e-300),
                 "direct_tail_estimate": direct.continuum_tail_estimate,
                 "reciprocal_shells": reciprocal.shells_used, "direct_radius": 300}
            )
    reciprocal = plane_exponential_sum_reciprocal(
        0.91, delta, kappa=3.0, amplitude=np.exp(3.0), geometry=geometry, config=config
    )
    direct = plane_exponential_sum_direct(
        0.91, delta, kappa=3.0, amplitude=np.exp(3.0), geometry=geometry, radius_index=40
    )
    for field in ("value", "d_d", "d2_dd"):
        left, right = float(getattr(reciprocal, field)), float(getattr(direct, field))
        validation_rows.append(
            {"kernel": "exponential", "quantity": field,
             "reciprocal": left, "direct": right, "absolute_error": abs(left-right),
             "relative_error": abs(left-right)/max(abs(right), 1e-300),
             "direct_tail_estimate": direct.continuum_tail_estimate,
             "reciprocal_shells": reciprocal.shells_used, "direct_radius": 40}
        )

    equilibrium_rows: list[dict[str, object]] = []
    hessian_rows: list[dict[str, object]] = []
    full_models: dict[str, FullFCC111StackEnergy] = {}
    for model_id in MODEL_IDS:
        reduced = build_energy_model(model_id)
        meta = energy_model_metadata(model_id)
        reduced_h = reduced.local_hessian(reduced.a0, 0.0)
        equilibrium_rows.append(
            {"model": model_id, "geometry": "reduced_two_row", "status": meta.calibration_status,
             "a0": reduced.a0, "W_a": reduced.local_deda(reduced.a0, 0.0),
             "W_s": reduced.local_deds(reduced.a0, 0.0),
             "minimum_hessian_eigenvalue": np.min(np.linalg.eigvalsh(reduced_h)),
             "kappa_axial_chi_0p2": reduced.sigma_over_E_force_scale(),
             "normalization": "energy per reduced two-row cell"}
        )
        hessian_rows.append(
            {"model": model_id, "geometry": "reduced_two_row", "W_aa": reduced_h[0, 0],
             "W_as": reduced_h[0, 1], "W_ss": reduced_h[1, 1],
             "minimum_eigenvalue": np.min(np.linalg.eigvalsh(reduced_h))}
        )
        full = FullFCC111StackEnergy.from_reduced_model(
            reduced, require_stable_equilibrium=False
        )
        full_models[model_id] = full
        if np.isfinite(full.a0):
            h = full.hessian(full.a0, 0.0)
            equilibrium_rows.append(
                {"model": model_id, "geometry": "full_fcc111_homogeneous_stack",
                 "status": "stable_same_parameters_unrefitted", "a0": full.a0,
                 "W_a": full.local_deda(full.a0, 0.0), "W_s": full.local_deds(full.a0, 0.0),
                 "minimum_hessian_eigenvalue": np.min(np.linalg.eigvalsh(h)),
                 "kappa_axial_chi_0p2": full.sigma_over_E_force_scale(),
                 "normalization": "energy per atom"}
            )
            hessian_rows.append(
                {"model": model_id, "geometry": "full_fcc111_homogeneous_stack",
                 "W_aa": h[0, 0], "W_as": h[0, 1], "W_ss": h[1, 1],
                 "minimum_eigenvalue": np.min(np.linalg.eigvalsh(h))}
            )
        else:
            equilibrium_rows.append(
                {"model": model_id, "geometry": "full_fcc111_homogeneous_stack",
                 "status": full.equilibrium_error, "a0": "", "W_a": "", "W_s": "",
                 "minimum_hessian_eigenvalue": "", "kappa_axial_chi_0p2": "",
                 "normalization": "energy per atom"}
            )
    _write_csv(output_dir / "equilibrium_comparison.csv", equilibrium_rows)
    _write_csv(output_dir / "hessian_comparison.csv", hessian_rows)

    barrier_rows: list[dict[str, object]] = []
    for model_id in MODEL_IDS:
        reduced = build_energy_model(model_id)
        barrier_rows.append(
            {"model": model_id, "geometry": "reduced_two_row", "force": 0.0,
             "force_fraction_opening_spinodal": 0.0,
             "configurational_barrier": _reduced_relaxed_barrier(reduced),
             "opening_barrier": "", "barrier_ratio": "",
             "configurational_spinodal": "", "opening_spinodal": "",
             "status": "reduced reference; opening metric uses different geometry"}
        )
        full = full_models[model_id]
        if not np.isfinite(full.a0):
            barrier_rows.append(
                {"model": model_id, "geometry": "full_fcc111_homogeneous_stack", "force": "",
                 "force_fraction_opening_spinodal": "", "configurational_barrier": "",
                 "opening_barrier": "", "barrier_ratio": "",
                 "configurational_spinodal": "", "opening_spinodal": "",
                 "status": "unavailable: no stable full-stack equilibrium"}
            )
            continue
        opening_force, _ = opening_spinodal(full, 0.0)
        config_spin = configurational_spinodal(full)
        for fraction in (0.0, 0.25, 0.50, 0.75):
            force = fraction * opening_force
            registry = relaxed_registry_barrier(full, force, samples=25)
            opening = opening_barrier(full, 0.0, force)
            reg_value = None if registry is None else registry.barrier
            open_value = None if opening is None else opening.barrier
            ratio = (
                None if reg_value is None or open_value is None or open_value <= 0
                else reg_value / open_value
            )
            barrier_rows.append(
                {"model": model_id, "geometry": "full_fcc111_homogeneous_stack",
                 "force": force, "force_fraction_opening_spinodal": fraction,
                 "configurational_barrier": "" if reg_value is None else reg_value,
                 "opening_barrier": "" if open_value is None else open_value,
                 "barrier_ratio": "" if ratio is None else ratio,
                 "configurational_spinodal": "" if config_spin is None else config_spin.force,
                 "opening_spinodal": opening_force,
                 "status": "static barriers; no dynamic ordering claim"}
            )
    _write_csv(output_dir / "barrier_comparison.csv", barrier_rows)

    convergence_rows: list[dict[str, object]] = []
    for model_id in (TWO_ROW_LJ_REFERENCE, ANALYTIC_LJ_EAM_HYPOTHETICAL):
        full = full_models[model_id]
        a, s = 1.03 * full.a0, 0.137
        started = time.perf_counter()
        analytic = full.evaluate(a, s)
        elapsed_ms = 1e3 * (time.perf_counter() - started)
        direct = full.direct_reference(a, s, radial_index=60, layers=60)
        for field, direct_field in (
            ("energy", "energy"), ("d_da", "d_da"), ("d_ds", "d_ds"),
            ("d2_daa", "d2_daa"), ("d2_das", "d2_das"), ("d2_dss", "d2_dss"),
        ):
            left, right = float(getattr(analytic, field)), float(getattr(direct, direct_field))
            convergence_rows.append(
                {"model": model_id, "quantity": field, "reciprocal": left, "direct": right,
                 "absolute_error": abs(left-right),
                 "relative_error": abs(left-right)/max(abs(right), 1e-300),
                 "pair_tail_estimate": direct.estimated_pair_tail,
                 "density_tail_estimate": direct.estimated_density_tail,
                 "reciprocal_pair_layers": analytic.pair.layers_used,
                 "reciprocal_density_layers": analytic.density.layers_used,
                 "maximum_plane_shells": max(analytic.pair.maximum_plane_shells,
                                               analytic.density.maximum_plane_shells),
                 "cached_evaluation_ms": elapsed_ms}
            )
        validation_rows.append(
            {"kernel": f"full_{model_id}", "quantity": "energy", "reciprocal": analytic.energy,
             "direct": direct.energy, "absolute_error": abs(analytic.energy-direct.energy),
             "relative_error": abs(analytic.energy-direct.energy)/max(abs(direct.energy),1e-300),
             "direct_tail_estimate": direct.estimated_pair_tail + direct.estimated_density_tail,
             "reciprocal_shells": max(analytic.pair.maximum_plane_shells,
                                      analytic.density.maximum_plane_shells), "direct_radius": 60}
        )
    _write_csv(output_dir / "lattice_sum_validation.csv", validation_rows)
    _write_csv(output_dir / "convergence_summary.csv", convergence_rows)

    summary = {
        "status": "research_reference_not_production_default",
        "scalar_path": DIRECT_110,
        "alternative_path": SHOCKLEY_112,
        "physical_time": "uncalibrated; model time only",
        "A_c_used": False,
        "best_feasible_full_stack": full_models[AL_TARGET_BEST_FEASIBLE].equilibrium_error,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = run(root / "results" / "fcc111_full_stack")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
