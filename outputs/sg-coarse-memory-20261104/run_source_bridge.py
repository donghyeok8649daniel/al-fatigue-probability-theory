"""Audit new coarse-memory claims against the existing N=1 SG entry point."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    CyclicLoad2D, Grid2DParams, PDETimeParams, _sg_generator_2d,
    _sg_rates_and_rhs, build_grid, energy_grid, initial_gibbs_density,
    opening_intact_mask, run_probability_pde_2d,
)
from coarse_memory import aggregation, cdf_bracket, closure_witness, exact_memory, markov_audit


def main():
    started = time.perf_counter()
    checks = []

    def check(name, condition, value=None):
        checks.append({"name": name, "passed": bool(condition), "value": value})
        if not condition:
            raise AssertionError(f"{name}: {value}")

    params = ModelParams(n_cells=1, kT=0.02, mobility_a=1.0,
                         mobility_s=0.05, chi_axial_projection=0.2)
    model = TwoRowLJ(params)
    print("Building the original opening table once", flush=True)
    model._build_opening_table()
    grid_params = Grid2DParams(n_a=9, n_s=15, s_wells=3, a_upper=1.8)
    grid = build_grid(model, grid_params)
    forces = [1.0, 1.4, 0.8]
    energies = [energy_grid(model, grid, force) for force in forces]
    masks = [opening_intact_mask(model, grid, force).ravel() for force in forces]
    generators = [_sg_generator_2d(energy, model, grid).toarray() for energy in energies]
    max_rate = max(float(np.max(-np.diag(a))) for a in generators)
    dt = 0.4 / max_rate
    fine_size = grid.a.size * grid.s.size

    def transition(source, target):
        u = np.eye(fine_size) + dt * generators[source]
        k = masks[target][:, None] * u * masks[source][None, :]
        t = np.zeros((fine_size + 1, fine_size + 1))
        t[:-1, :-1] = k
        t[-1, :-1] = 1.0 - np.sum(k, axis=0)
        t[-1, -1] = 1.0
        return t

    initial_density = initial_gibbs_density(model, grid, preload_force=0.0)
    before_mask = initial_density.ravel() * grid.cell_volume
    initial = np.r_[before_mask * masks[0], float(before_mask[~masks[0]].sum())]
    schedule = [0] * 8 + [1] * 8 + [2] * 8 + [0]
    steps = [transition(schedule[n], schedule[n + 1]) for n in range(24)]
    valid = [np.r_[masks[index], True] for index in schedule]
    check("transition_nonnegative", min(float(t.min()) for t in steps) >= -2e-15)
    check("transition_column_mass", max(float(np.max(abs(t.sum(axis=0) - 1))) for t in steps) < 2e-14)
    check("cemetery_absorbing", all(np.array_equal(t[:, -1], np.eye(fine_size + 1)[:, -1]) for t in steps))

    density_probe = np.linspace(0.1, 0.9, fine_size).reshape(energies[0].shape)
    rhs, _ = _sg_rates_and_rhs(density_probe, energies[0], model, grid)
    check("generator_matches_original_rhs", np.max(abs(generators[0] @ density_probe.ravel() - rhs.ravel())) < 1e-10)

    print("Checking the ordinary solver entry point", flush=True)
    constant_steps = [transition(0, 0)] * 24
    constant = initial.copy()
    for t in constant_steps:
        constant = t @ constant
    result = run_probability_pde_2d(
        prepared_model=model, grid_params=grid_params,
        time_params=PDETimeParams(max_dt=dt, cfl=0.45, record_interval=dt),
        load=CyclicLoad2D(force_min=forces[0], force_max=forces[0], period=24 * dt, cycles=1),
    )
    entry_gap = float(np.max(abs(result["density"].ravel() * grid.cell_volume - constant[:-1])))
    cdf_gap = abs(float(result["initiation_probability"][-1]) - constant[-1])
    check("normal_entry_density_match", entry_gap < 2e-12, entry_gap)
    check("normal_entry_cdf_match", cdf_gap < 2e-12, float(cdf_gap))
    check("normal_entry_negative_repair_zero", result["cumulative_negative_mass_correction"][-1] == 0)

    fine_path = [initial]
    for t in steps:
        fine_path.append(t @ fine_path[-1])
    fine_path = np.asarray(fine_path)
    well_labels = np.tile(model.well_index(grid.s), grid.a.size)
    common = np.logical_and.reduce(masks)
    # A fixed, declared numerical lifting supported in every scheduled basin.
    weights = np.r_[common.astype(float), 1.0]
    levels = []
    memory_result = None
    for bands in [1, 2, 3]:
        a_bands = np.minimum(np.arange(grid.a.size) * bands // grid.a.size, bands - 1)
        raw = np.c_[np.repeat(a_bands, grid.s.size), well_labels]
        # Empty common-support geometric groups cannot define a lifting.
        # Merge their labels with the nearest supported normal band in that well.
        for well in np.unique(well_labels):
            active_bands = np.unique(raw[common & (well_labels == well), 0])
            for band in np.unique(raw[well_labels == well, 0]):
                if band not in active_bands:
                    nearest = active_bands[np.argmin(abs(active_bands - band))]
                    raw[(well_labels == well) & (raw[:, 0] == band), 0] = nearest
        _, labels = np.unique(raw, axis=0, return_inverse=True)
        labels = np.r_[labels, labels.max() + 1]
        p, r = aggregation(labels, weights)
        check(f"right_inverse_bands{bands}", np.max(abs(p @ r - np.eye(p.shape[0]))) < 1e-14)
        audit = markov_audit(steps, fine_path, p, r, valid)
        lower, upper = cdf_bracket(steps, p, labels, p @ initial, valid)
        actual = float(fine_path[-1, -1])
        check(f"cdf_envelope_bands{bands}", lower - 2e-13 <= actual <= upper + 2e-13, [lower, actual, upper])
        check(f"posterior_bound_bands{bands}", all(row["aggregate_tv_error"] <= row["posterior_tv_bound"] + 2e-13 for row in audit))
        check(f"uniform_bound_bands{bands}", all(row["aggregate_tv_error"] <= row["uniform_tv_bound"] + 2e-13 for row in audit))
        levels.append({"requested_normal_bands": bands, "actual_macrostates_with_cemetery": p.shape[0],
                       "cdf_lower": lower, "cdf_upper": upper, "audit": audit})
        if bands == 1:
            memory, initial_terms, memory_terms = exact_memory(steps, p, r, initial)
            macro_truth = fine_path @ p.T
            error = float(np.max(abs(memory - macro_truth)))
            check("ordered_memory_matches_fine_aggregation", error < 2e-12, error)
            omitted_initial, _, _ = exact_memory(steps, p, r, initial, retain_initial=False)
            omission_gap = float(np.max(abs(omitted_initial - macro_truth)))
            witness = closure_witness(steps[7], p, labels, valid[7])
            check("admissible_equal_macro_witness", witness["diameter_tv"] > 1e-8, witness["diameter_tv"])
            check("initial_hidden_law_matters", omission_gap > 1e-10, omission_gap)
            memory_result = {"maximum_aggregate_error": error, "omitting_initial_hidden_term_error": omission_gap,
                             "maximum_initial_term_l1": float(np.max(abs(initial_terms).sum(axis=1))),
                             "maximum_memory_term_l1": float(np.max(abs(memory_terms).sum(axis=1))),
                             "closure_witness": witness}

    # Exact identity partition is the limiting envelope and closure reference.
    identity_labels = np.arange(fine_size + 1)
    p_identity, r_identity = aggregation(identity_labels, np.ones(fine_size + 1))
    identity_valid = [np.ones(fine_size + 1, dtype=bool)] * len(valid)
    lower, upper = cdf_bracket(steps, p_identity, identity_labels, initial, identity_valid)
    check("identity_envelope_collapses", max(abs(lower - fine_path[-1, -1]), abs(upper - lower)) < 2e-13)

    receipt = {
        "status": "PASS_SOURCE_BOUND_FINITE_GRID_AGGREGATION_AUDIT",
        "dateKST": "2026-11-04", "pid": os.getpid(),
        "wall_seconds": time.perf_counter() - started,
        "source": "solver_v1/probability_pde_2d.py",
        "model_params": vars(params), "grid_params": vars(grid_params),
        "fine_states_with_cemetery": fine_size + 1, "dt_model": dt,
        "forces_model": forces, "schedule_indices": schedule,
        "initial_absorbed_mass": float(initial[-1]),
        "final_absorbed_mass": float(fine_path[-1, -1]),
        "normal_entry_density_max_error": entry_gap, "normal_entry_cdf_error": float(cdf_gap),
        "memory": memory_result, "partition_levels": levels,
        "checks": checks, "passed_checks": len(checks),
        "interval_roundoff_certificate": False,
        "grid_or_continuum_convergence_claim": False,
        "physical_time_or_material_calibration": None,
        "lifting": "Uniform over the common intact support, chosen numerical lifting",
        "production_source_edits": False, "independent_review": "NOT_PERFORMED",
    }
    output = Path(__file__).with_name("SOURCE_BRIDGE_RESULT.json")
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "checks": len(checks),
                      "seconds": receipt["wall_seconds"], "output": str(output)}), flush=True)


if __name__ == "__main__":
    main()
