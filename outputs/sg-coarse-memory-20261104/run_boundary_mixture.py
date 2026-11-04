"""A declared equal-aggregate preparation pair on original SG operators."""
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
    Grid2DParams, _sg_generator_2d, build_grid, energy_grid,
    initial_gibbs_density, opening_intact_mask,
)
from coarse_memory import aggregation, cdf_bracket, exact_memory, markov_audit


def main():
    started = time.perf_counter()
    model = TwoRowLJ(ModelParams(n_cells=1, kT=0.02, mobility_a=1.0,
                                mobility_s=0.05, chi_axial_projection=0.2))
    model._build_opening_table()
    grid = build_grid(model, Grid2DParams(n_a=9, n_s=15, s_wells=3, a_upper=1.8))
    forces = [1.0, 1.4, 0.8]
    masks = [opening_intact_mask(model, grid, f).ravel() for f in forces]
    generators = [_sg_generator_2d(energy_grid(model, grid, f), model, grid).toarray() for f in forces]
    dt = 0.4 / max(float(np.max(-np.diag(a))) for a in generators)
    size = generators[0].shape[0]
    schedule = [0] * 8 + [1] * 8 + [2] * 8 + [0]
    steps = []
    for source, target in zip(schedule[:-1], schedule[1:]):
        u = np.eye(size) + dt * generators[source]
        k = masks[target][:, None] * u * masks[source][None, :]
        t = np.zeros((size + 1, size + 1))
        t[:-1, :-1] = k
        t[-1, :-1] = 1.0 - k.sum(axis=0)
        t[-1, -1] = 1.0
        steps.append(t)
    valid = [np.r_[masks[j], True] for j in schedule]
    well = np.tile(model.well_index(grid.s), grid.a.size)
    candidate = np.flatnonzero(masks[0] & ~masks[1] & (well == 0))
    if candidate.size == 0:
        raise ValueError("This declared pair requires a principal-well cell crossed by the load mask")
    boundary = int(candidate[np.argmax(steps[7][-1, candidate])])
    raw = initial_gibbs_density(model, grid).ravel() * grid.cell_volume
    base = np.r_[raw * masks[0], raw[~masks[0]].sum()]
    mixed = 0.6 * base
    mixed[boundary] += 0.4 * float(base[:-1].sum())
    mixed[-1] = base[-1]
    paths = []
    for initial in (base, mixed):
        path = [initial]
        for t in steps:
            path.append(t @ path[-1])
        paths.append(np.asarray(path))

    common = np.logical_and.reduce(masks)
    weights = np.r_[common.astype(float), 1.0]
    levels = []
    checked = []

    def check(name, condition, value=None):
        checked.append({"name": name, "passed": bool(condition), "value": value})
        if not condition:
            raise AssertionError(name)

    for normal_resolution in [1, 3, 9]:
        raw_labels = np.c_[np.repeat(np.arange(grid.a.size) * normal_resolution // grid.a.size, grid.s.size), well]
        for group_well in np.unique(well):
            active = np.unique(raw_labels[common & (well == group_well), 0])
            for band in np.unique(raw_labels[well == group_well, 0]):
                if band not in active:
                    target = active[np.argmin(abs(active - band))]
                    raw_labels[(well == group_well) & (raw_labels[:, 0] == band), 0] = target
        _, labels = np.unique(raw_labels, axis=0, return_inverse=True)
        labels = np.r_[labels, labels.max() + 1]
        p, r = aggregation(labels, weights)
        records = []
        for initial, path in zip((base, mixed), paths):
            audit = markov_audit(steps, path, p, r, valid)
            lower, upper = cdf_bracket(steps, p, labels, p @ initial, valid)
            truth = float(path[-1, -1])
            check(f"bracket_{normal_resolution}_{len(records)}", lower - 2e-13 <= truth <= upper + 2e-13)
            check(f"tv_bound_{normal_resolution}_{len(records)}", all(row['cdf_error'] <= row['posterior_tv_bound'] + 2e-13 for row in audit))
            records.append({"truth": truth, "lower": lower, "upper": upper, "final_audit": audit[-1]})
        if normal_resolution == 1:
            check("same_raw_macro_initial_law", np.max(abs(p @ base - p @ mixed)) < 1e-14)
            check("same_initial_absorbed_mass", base[-1] == mixed[-1])
            check("distinct_finite_horizon_cdf", abs(paths[0][-1, -1] - paths[1][-1, -1]) > 0.1)
            memory, _, _ = exact_memory(steps, p, r, mixed)
            memory_gap = float(np.max(abs(memory - paths[1] @ p.T)))
            check("boundary_mixture_exact_memory", memory_gap < 2e-12, memory_gap)
            check("same_macro_predictor_for_both", abs(records[0]['final_audit']['cdf_markov'] - records[1]['final_audit']['cdf_markov']) < 1e-14)
        levels.append({"normal_resolution": normal_resolution, "macrostates": p.shape[0], "cases": records})

    receipt = {
        "status": "PASS_EQUAL_AGGREGATE_SOURCE_OPERATOR_COUNTEREXAMPLE",
        "dateKST": "2026-11-04", "pid": os.getpid(), "wall_seconds": time.perf_counter() - started,
        "declared_preparation": "60% original post-mask Gibbs survivor law plus40% surviving mass in one principal-well cell crossed by the next load increase; same A0",
        "actual_physical_initial_law": None,
        "boundary_cell_index": boundary, "boundary_a": float(grid.a[boundary // grid.s.size]),
        "boundary_s": float(grid.s[boundary % grid.s.size]),
        "initial_macro_equality": True, "levels": levels,
        "checks": checked, "passed_checks": len(checked),
        "new_normal_entrypoint_runs": 0,
        "normal_entrypoint_reference": "SOURCE_BRIDGE_RESULT.json",
        "independent_review": "NOT_PERFORMED", "outward_roundoff_certificate": False,
    }
    Path(__file__).with_name('BOUNDARY_MIXTURE_RESULT.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding='utf-8')
    np.savez_compressed(Path(__file__).with_name('SOURCE_OPERATOR_PAIR.npz'),
                        steps=np.asarray(steps), valid=np.asarray(valid), base=base, mixed=mixed,
                        fine_paths=np.asarray(paths), a=grid.a, s=grid.s)
    print(json.dumps({"status": receipt['status'], "checks": len(checked),
                      "cdf_base": float(paths[0][-1, -1]), "cdf_mixed": float(paths[1][-1, -1])}), flush=True)


if __name__ == '__main__':
    main()
