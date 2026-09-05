from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import numpy as np
from .model import ModelParams, TwoRowLJ

@dataclass
class LoadParams:
    force_max: float = 3.2
    force_min: float = 0.0
    period: float = 10.0
    cycles: int = 10
    phase_radians: float = -0.5*np.pi
    value_function: Callable[[float], float] | None = None

    def value(self, t: float) -> float:
        if self.value_function is not None:
            return float(self.value_function(t))
        mid = 0.5*(self.force_max + self.force_min)
        amp = 0.5*(self.force_max - self.force_min)
        return mid + amp*np.sin(2*np.pi*t/self.period + self.phase_radians)

@dataclass
class SolverParams:
    dt: float = 0.01
    n_trajectories: int = 64
    seed: int = 7
    first_passage_stride: int = 10
    record_stride: int = 10

def run_ensemble(
    model_p=ModelParams(),
    load_p=LoadParams(),
    solver_p=SolverParams(),
    model=None,
    *,
    record_callback=None,
    stop_requested=None,
    retain_history: bool = True,
):
    model = TwoRowLJ(model_p) if model is None else model
    rng = np.random.default_rng(solver_p.seed)

    n = model_p.n_cells
    nt = int(round(load_p.cycles*load_p.period/solver_p.dt)) + 1
    a = np.full((solver_p.n_trajectories, n), model.a0, dtype=float)
    s = np.zeros((solver_p.n_trajectories, n), dtype=float)
    s += rng.normal(scale=1e-3, size=s.shape)
    alive = np.ones(solver_p.n_trajectories, dtype=bool)
    invalid = np.zeros(solver_p.n_trajectories, dtype=bool)
    first_passage_time = np.full(solver_p.n_trajectories, np.nan, dtype=float)

    rec_t, rec_f, rec_eps, rec_normal, rec_intrawell, rec_eps_plastic = [], [], [], [], [], []
    rec_alive, rec_plastic, rec_barrier = [], [], []
    rec_opening_eigenvalue, rec_plastic_eigenvalue = [], []
    sqrta = np.sqrt(2*model_p.kT*model_p.mobility_a*solver_p.dt)
    sqrts = np.sqrt(2*model_p.kT*model_p.mobility_s*solver_p.dt)

    last_time = 0.0
    for step in range(nt):
        if stop_requested is not None and stop_requested():
            break
        t = step*solver_p.dt
        last_time = t
        force = load_p.value(t)

        if step % solver_p.record_stride == 0:
            idx = np.where(alive)[0]
            if len(idx):
                total_i, normal_i, intrawell_i, plastic_i = model.strain_components_batch(a[idx], s[idx])
                eps = float(np.mean(total_i))
                normal_strain = float(np.mean(normal_i))
                intrawell_strain = float(np.mean(intrawell_i))
                plastic_strain = float(np.mean(plastic_i))
                plastic = float(np.mean(np.abs(model.well_index(s[idx]))))
                opening_eigenvalue, plastic_eigenvalue, _ = model.modal_stability_batch(a[idx], s[idx])
                min_opening_eigenvalue = (
                    float(np.nanmin(opening_eigenvalue))
                    if np.any(np.isfinite(opening_eigenvalue)) else np.nan
                )
                min_plastic_eigenvalue = (
                    float(np.nanmin(plastic_eigenvalue))
                    if np.any(np.isfinite(plastic_eigenvalue)) else np.nan
                )
                barriers = [model.opening_barrier(float(s[k,0]), force) for k in idx[:min(12,len(idx))]]
                barrier = float(np.mean(barriers)) if barriers else np.nan
            else:
                eps = normal_strain = intrawell_strain = plastic_strain = np.nan
                plastic, barrier = np.nan, 0.0
                min_opening_eigenvalue = min_plastic_eigenvalue = np.nan
            valid_count = int(np.count_nonzero(~invalid))
            survival = float(np.count_nonzero(alive) / valid_count) if valid_count else np.nan
            if retain_history:
                rec_t.append(t); rec_f.append(force); rec_eps.append(eps)
                rec_normal.append(normal_strain); rec_intrawell.append(intrawell_strain)
                rec_eps_plastic.append(plastic_strain)
                rec_alive.append(survival); rec_plastic.append(plastic)
                rec_barrier.append(barrier)
                rec_opening_eigenvalue.append(min_opening_eigenvalue)
                rec_plastic_eigenvalue.append(min_plastic_eigenvalue)
            if record_callback is not None:
                record_callback({
                    "time": t,
                    "force": force,
                    "strain": eps,
                    "normal_strain": normal_strain,
                    "intrawell_strain": intrawell_strain,
                    "plastic_strain": plastic_strain,
                    "survival": survival,
                    "plastic_well_activity": plastic,
                    "opening_barrier": barrier,
                    "min_opening_eigenvalue": min_opening_eigenvalue,
                    "min_plastic_eigenvalue": min_plastic_eigenvalue,
                })

        if step == nt-1:
            break

        idx = np.where(alive)[0]
        if len(idx):
            _, ga, gs = model.energy_gradient_batch(a[idx], s[idx], force)
            good = np.all(np.isfinite(ga), axis=1) & np.all(np.isfinite(gs), axis=1)
            if np.any(~good):
                invalid[idx[~good]] = True
                alive[idx[~good]] = False
            good_idx = idx[good]
            if len(good_idx):
                a[good_idx] += -model_p.mobility_a*ga[good]*solver_p.dt + sqrta*rng.normal(size=(len(good_idx), n))
                s[good_idx] += -model_p.mobility_s*gs[good]*solver_p.dt + sqrts*rng.normal(size=(len(good_idx), n))
                a[good_idx] = np.maximum(a[good_idx], model_p.a_min*1.001)

        if step % solver_p.first_passage_stride == 0:
            idx = np.where(alive)[0]
            if len(idx):
                fail, _, _, _ = model.coupled_opening_escape_batch(a[idx], s[idx], force)
                first_passage_time[idx[fail]] = (step + 1) * solver_p.dt
                alive[idx[fail]] = False
                if not np.any(alive):
                    last_time = (step + 1) * solver_p.dt
                    break

    return {
        "model": model,
        "time": np.asarray(rec_t),
        "force": np.asarray(rec_f),
        "strain": np.asarray(rec_eps),
        "normal_strain": np.asarray(rec_normal),
        "intrawell_strain": np.asarray(rec_intrawell),
        "plastic_strain": np.asarray(rec_eps_plastic),
        "survival": np.asarray(rec_alive),
        "plastic_well_activity": np.asarray(rec_plastic),
        "opening_barrier": np.asarray(rec_barrier),
        "min_opening_eigenvalue": np.asarray(rec_opening_eigenvalue),
        "min_plastic_eigenvalue": np.asarray(rec_plastic_eigenvalue),
        "first_passage_time": first_passage_time,
        "invalid_trajectory": invalid,
        "observation_end_time": last_time,
    }
