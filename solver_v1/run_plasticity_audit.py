"""Reproducible deterministic audit used by PLASTICITY_AND_STRAIN_AUDIT.md."""
from __future__ import annotations

import argparse
import json

import numpy as np

from .configurational_landscape import bound_configurational_barrier
from .energy_model_registry import TWO_ROW_LJ_REFERENCE, build_energy_model
from .plasticity_audit import (
    certify_plasticity_refinement,
    final_well_populations,
    summarize_plasticity_run,
)
from .probability_pde_2d import (
    CyclicLoad2D,
    Grid2DParams,
    PDETimeParams,
    cyclic_load_from_sigma_over_E,
    run_probability_pde_2d,
)


def _run(
    model,
    *,
    grid_a: int,
    grid_s: int,
    dt: float,
    integrator: str,
    cycles: float,
    period: float = 0.04,
):
    load = cyclic_load_from_sigma_over_E(
        model,
        sigma_over_E_min=-1100.0 / 69000.0,
        sigma_over_E_max=2900.0 / 69000.0,
        period=period,
        cycles=cycles,
    )
    preload = float(model.force_from_sigma_over_E(900.0 / 69000.0))
    return run_probability_pde_2d(
        prepared_model=model,
        grid_params=Grid2DParams(
            n_a=grid_a, n_s=grid_s, s_wells=3, a_upper=1.60
        ),
        time_params=PDETimeParams(
            max_dt=dt,
            cfl=0.40,
            record_interval=period / 20.0,
            integrator=integrator,
        ),
        load=load,
        preload_force=preload,
    )


def _cycle_rows(result, period: float) -> list[dict[str, float]]:
    time = np.asarray(result["time"], dtype=float)
    cycle = time / period
    rows: list[dict[str, float]] = []
    for number in range(1, int(np.ceil(cycle[-1] - 1.0e-12)) + 1):
        start, end = number - 1.0, min(float(number), float(cycle[-1]))
        mask = (cycle >= start - 1.0e-12) & (cycle <= end + 1.0e-12)
        if not np.any(mask):
            continue
        def interp(key: str, where: float) -> float:
            return float(np.interp(where, cycle, np.asarray(result[key], dtype=float)))
        total = np.asarray(result["strain"], dtype=float)[mask]
        rows.append(
            {
                "cycle": number,
                "mean_epsilon_total": float(np.mean(total)),
                "epsilon_total_amplitude": float(0.5 * np.ptp(total)),
                "mean_epsilon_a": float(
                    np.mean(np.asarray(result["normal_strain"])[mask])
                ),
                "mean_epsilon_xi": float(
                    np.mean(np.asarray(result["intrawell_strain"])[mask])
                ),
                "epsilon_p_increment": interp("plastic_strain", end)
                - interp("plastic_strain", start),
                "epsilon_p_end": interp("plastic_strain", end),
                "net_registry_transfer": interp(
                    "accumulated_net_registry_transfer", end
                )
                - interp("accumulated_net_registry_transfer", start),
                "gross_registry_activity": interp(
                    "cumulative_gross_registry_activity", end
                )
                - interp("cumulative_gross_registry_activity", start),
                "absorbed_mass": interp("cumulative_absorbed_mass", end)
                - interp("cumulative_absorbed_mass", start),
            }
        )
    return rows


def _unload_hold(
    model, *, grid_a: int, grid_s: int, dt: float, hold: float = 0.40
):
    period = 0.04
    ramp = 0.02
    duration = period + ramp + hold
    mean = 900.0 / 69000.0
    amplitude = 2000.0 / 69000.0

    def reduced_stress(time: float) -> float:
        if time < period:
            return float(mean + amplitude * np.sin(2.0 * np.pi * time / period))
        if time < period + ramp:
            return float(mean * (1.0 - (time - period) / ramp))
        return 0.0

    load = cyclic_load_from_sigma_over_E(
        model,
        sigma_over_E_min=mean - amplitude,
        sigma_over_E_max=mean + amplitude,
        period=1.0,
        cycles=duration,
        value_function=reduced_stress,
    )
    result = run_probability_pde_2d(
        prepared_model=model,
        grid_params=Grid2DParams(
            n_a=grid_a, n_s=grid_s, s_wells=3, a_upper=1.60
        ),
        time_params=PDETimeParams(
            max_dt=dt,
            cfl=0.40,
            record_interval=0.01,
            integrator="implicit",
        ),
        load=load,
        preload_force=float(model.force_from_sigma_over_E(mean)),
    )
    unload_index = int(np.argmin(np.abs(np.asarray(result["time"]) - period)))
    return result, {
        "grid_a": grid_a,
        "grid_s": grid_s,
        "dt": dt,
        "epsilon_p_end_loading": float(result["plastic_strain"][unload_index]),
        "epsilon_p_end_hold": float(result["plastic_strain"][-1]),
        "mean_well_index_end_hold": float(result["mean_well_index"][-1]),
        "epsilon_a_end_hold": float(result["normal_strain"][-1]),
        "epsilon_xi_end_hold": float(result["intrawell_strain"][-1]),
        "opening_absorbed_mass": float(result["cumulative_absorbed_mass"][-1]),
        "well_populations_end": final_well_populations(result),
    }


def run_audit(*, cycles: float = 3.0) -> dict[str, object]:
    model_id = TWO_ROW_LJ_REFERENCE
    model = build_energy_model(model_id)
    model._build_opening_table()
    specifications = (
        ("A_current", 21, 31, 1.0e-3, "explicit"),
        ("A_implicit", 21, 31, 1.0e-3, "implicit"),
        ("B_finer_aligned_s", 21, 42, 1.0e-3, "implicit"),
        ("C_finer_a", 31, 31, 1.0e-3, "implicit"),
        ("D_dt_half", 21, 31, 5.0e-4, "implicit"),
        ("E_combined", 31, 42, 5.0e-4, "implicit"),
        ("F_aligned60", 31, 60, 5.0e-4, "implicit"),
        ("G_aligned60_dt", 41, 60, 2.5e-4, "implicit"),
    )
    results = []
    summaries = []
    for label, grid_a, grid_s, dt, integrator in specifications:
        result = _run(
            model,
            grid_a=grid_a,
            grid_s=grid_s,
            dt=dt,
            integrator=integrator,
            cycles=cycles,
        )
        result["cumulative_gross_registry_activity"] = np.sum(
            result["cumulative_interwell_gross_transfer"], axis=1
        )
        results.append(result)
        summaries.append(
            summarize_plasticity_run(
                result,
                label=label,
                dt=dt,
                integrator=integrator,
                energy_model=model_id,
            ).as_dict()
        )
    # Only aligned grids enter the final plastic-flow resolution statement.
    # Unaligned runs remain in the table as evidence of partition sensitivity.
    assessment = certify_plasticity_refinement(
        [results[index] for index in (2, 5, 6, 7)]
    )
    peak_force = float(model.force_from_sigma_over_E(2900.0 / 69000.0))
    barrier = bound_configurational_barrier(
        model, peak_force, n_s=81, quadrature_order=32
    )
    hold_coarse_result, hold_coarse = _unload_hold(
        model, grid_a=21, grid_s=42, dt=2.0e-3, hold=2.0
    )
    hold_fine_result, hold_fine = _unload_hold(
        model, grid_a=31, grid_s=60, dt=1.0e-3, hold=2.0
    )
    hold_assessment = certify_plasticity_refinement(
        [hold_coarse_result, hold_fine_result]
    )
    return {
        "refinement": summaries,
        "plasticity_assessment": {
            **assessment.__dict__,
        },
        "per_cycle": _cycle_rows(results[-1], 0.04),
        "final_well_populations": final_well_populations(results[-1]),
        "barrier": {
            "stress_mpa": 2900.0,
            "force": peak_force,
            "configurational_barrier": None if barrier is None else barrier.barrier,
            "opening_barrier_at_minimum": (
                None if barrier is None else barrier.opening_barrier_at_minimum
            ),
            "opening_barrier_at_saddle": (
                None if barrier is None else barrier.opening_barrier_at_saddle
            ),
        },
        "unload_hold": [hold_coarse, hold_fine],
        "unload_hold_assessment": hold_assessment.__dict__,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cycles", type=float, default=3.0)
    args = parser.parse_args()
    print(json.dumps(run_audit(cycles=args.cycles), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
