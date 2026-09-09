"""Read ACTUALLY executed runs; compare driven density with matched zero load.

No optimizer, PDE re-execution, invented crack probability, or fatigue fitting.
Floors distinguish raw grid differences from paired driven-minus-control
differences. An uncertified signal is not made physical by a large ratio.
"""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import numpy as np

from .run_low_stress_cyclic_diagnostic import OUT, write_csv, save_json


def rows(path):
    if not path.exists() and path.suffix == ".csv":
        path = path.with_suffix(".csv.gz")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        result = list(csv.DictReader(stream))
    for r in result:
        for k, v in r.items():
            try: r[k] = float(v)
            except (ValueError, TypeError): pass
    return result


def snapshot(data, cycle):
    """Read either complete legacy output or explicitly sampled cycle indices."""
    if "cycle_indices" not in data:
        return data["same_phase"][cycle]
    indices = data["cycle_indices"]
    if cycle == -1:
        return data["same_phase"][-1]
    selected = np.flatnonzero(indices == cycle)
    if len(selected) != 1:
        raise ValueError(f"cycle {cycle} was not serialized")
    return data["same_phase"][selected[0]]


def harmonic(history, cycle, period, field, phase=0.):
    selected = [r for r in history if r["cycle"] == cycle]
    angle = 2*np.pi*np.array([r["time_model"] for r in selected])/period
    y = np.array([r[field] for r in selected])
    fit = np.linalg.lstsq(np.column_stack((np.ones(len(angle)), np.sin(angle), np.cos(angle))), y, rcond=None)[0]
    theta = np.arctan2(fit[2], fit[1])-phase
    return float(np.hypot(fit[1], fit[2])), float(np.degrees(np.arctan2(np.sin(theta), np.cos(theta))))


def load_runs():
    result = []
    for path in sorted((OUT/"cycles").glob("*/audit.json")):
        audit = json.loads(path.read_text())
        cycle_rows = rows(path.with_name("per_cycle.csv"))
        history = rows(path.with_name("history.csv"))
        active = [r for r in cycle_rows if r["segment"] == "cyclic"]
        key = (audit["period_model"], audit["grid_a"], audit["grid_s"],
               audit["steps_per_cycle"], tuple(audit["domain_a_edges"]),
               audit["initial_central_only"], len(active), len(cycle_rows)-len(active))
        key = (*key, audit.get("integrator", "implicit"))
        result.append(dict(path=path.parent, audit=audit, cycles=cycle_rows,
                           history=history, active=active, key=key))
    return result


def main():
    runs = load_runs(); zero = {r["key"]: r for r in runs if r["audit"]["case"] == "zero"}
    summary = []; paired_cycles = []
    for run in runs:
        audit = run["audit"]; active = run["active"]; end = active[-1]; held = run["cycles"][-1]
        control = zero.get(run["key"])
        pair = control["active"][-1] if control else None
        control_end = control["cycles"][-1] if control else None
        diff = lambda key: end[key]-pair[key] if pair else None
        amp, phase = harmonic(run["history"], end["cycle"], audit["period_model"],
            "registry_shear", audit["load"]["phase_s"])
        lamp, lphase = harmonic(run["history"], end["cycle"], audit["period_model"], "normal_strain")
        row = dict(run=run["path"].name, case=audit["case"], period_model=audit["period_model"],
            grid_a=audit["grid_a"], grid_s=audit["grid_s"], grid_ds=audit["grid_ds"],
            steps_per_cycle=audit["steps_per_cycle"], wells=audit["wells"],
            integrator=audit.get("integrator", "implicit"),
            upper_a_over_h=(audit["domain_a_edges"][1]/np.sqrt(2/3)),
            driven_cycles=len(active), hold_cycles=len(run["cycles"])-len(active),
            central_initial=audit["initial_central_only"], matched_zero_available=control is not None,
            end_outside_mass=end["outside_central_mass"], zero_end_outside_mass=pair["outside_central_mass"] if pair else None,
            excess_outside_mass=diff("outside_central_mass"), end_mean_n=end["mean_well_index"],
            excess_mean_n=diff("mean_well_index"), final_hold_mean_n=held["mean_well_index"],
            hold_excess_mean_n=held["mean_well_index"]-control_end["mean_well_index"] if control_end else None,
            hold_intrawell_shear=held["intrawell_shear"],
            cumulative_net_transfer=end["cumulative_net_transfer"],
            cumulative_gross_grid_traffic=end["cumulative_gross_grid_traffic"],
            last_cycle_net_transfer=end["cycle_net_transfer"],
            last_cycle_gross_grid_traffic=end["cycle_gross_grid_traffic"],
            last_cycle_external_work_ev=end["external_work_ev_cell"],
            cycle1_external_work_ev=active[0]["external_work_ev_cell"],
            last_cycle_normal_amplitude=lamp,
            last_cycle_normal_phase_deg=lphase if audit["load"]["normal_amplitude_mpa"] else None,
            last_cycle_registry_amplitude=amp,
            last_cycle_registry_phase_deg=phase if audit["load"]["shear_amplitude_mpa"] else None,
            last_cycle_distribution_L1_change=end["same_phase_L1"],
            mass_residual=audit["max_mass_residual"], well_balance=audit["max_well_balance_residual"],
            registry_balance=audit["max_registry_balance_residual"],
            numerical_repair=audit["numerical_repair"],
            upper_normal_edge_mass=end["upper_normal_edge_mass"], outer_well_mass=end["outer_well_mass"],
            crack_status="not evaluated in reflecting diagnostic",
            physical_plasticity_status="not validated: one-cell thermal normalization and material gates unresolved")
        if control:
            with np.load(run["path"]/"snapshots.npz") as data, np.load(control["path"]/"snapshots.npz") as ref:
                excess = snapshot(data, len(active))-snapshot(ref, len(active))
                prior = snapshot(data, len(active)-1)-snapshot(ref, len(active)-1)
                held_excess = snapshot(data, -1)-snapshot(ref, -1)
                row["end_driven_vs_zero_L1"] = float(np.sum(np.abs(excess)))
                row["hold_driven_vs_zero_L1"] = float(np.sum(np.abs(held_excess)))
                row["last_cycle_excess_distribution_change_L1"] = float(np.sum(abs(excess-prior)))
            for current, base in zip(run["cycles"], control["cycles"]):
                paired_cycles.append(dict(run=run["path"].name, cycle=current["cycle"], segment=current["segment"],
                    excess_outside_mass=current["outside_central_mass"]-base["outside_central_mass"],
                    excess_net_transfer=current["cycle_net_transfer"]-base["cycle_net_transfer"],
                    cycle_net_transfer=current["cycle_net_transfer"],
                    cycle_work_ev=current["external_work_ev_cell"],
                    mean_n=current["mean_well_index"], P0=current["P_+0"],
                    Pminus1=current["P_-1"], Pplus1=current["P_+1"]))
        else:
            row.update(end_driven_vs_zero_L1=None, hold_driven_vs_zero_L1=None,
                       last_cycle_excess_distribution_change_L1=None)
        summary.append(row)
    write_csv(OUT/"cyclic_summary.csv", summary)
    write_csv(OUT/"paired_cycle_summary.csv", paired_cycles)
    save_json(OUT/"summary_execution.json", dict(executed_runs=len(runs),
        purpose="summary of saved actual experiments, not a rerun",
        sources_with_missing_zero_controls=[r["run"] for r in summary if not r["matched_zero_available"]],
        physical_hz=False, production_adoption=False, crack_probability_computed=False))
    for r in summary:
        print(r["run"], "outside", r["end_outside_mass"], "excess", r["excess_outside_mass"],
              "hold_n", r["final_hold_mean_n"], flush=True)


if __name__ == "__main__": main()
