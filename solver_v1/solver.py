from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .model import ModelParams, TwoRowLJ

@dataclass
class LoadParams:
    force_max: float = 3.2
    force_min: float = 0.0
    period: float = 10.0
    cycles: int = 10

    def value(self, t: float) -> float:
        mid = 0.5*(self.force_max + self.force_min)
        amp = 0.5*(self.force_max - self.force_min)
        return mid - amp*np.cos(2*np.pi*t/self.period)

@dataclass
class SolverParams:
    dt: float = 0.01
    n_trajectories: int = 64
    seed: int = 7
    first_passage_stride: int = 10
    record_stride: int = 10

def run_ensemble(model_p=ModelParams(), load_p=LoadParams(), solver_p=SolverParams(), model=None):
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

    rec_t, rec_f, rec_eps, rec_alive, rec_plastic, rec_barrier = [], [], [], [], [], []
    sqrta = np.sqrt(2*model_p.kT*model_p.mobility_a*solver_p.dt)
    sqrts = np.sqrt(2*model_p.kT*model_p.mobility_s*solver_p.dt)

    for step in range(nt):
        t = step*solver_p.dt
        force = load_p.value(t)

        if step % solver_p.record_stride == 0:
            idx = np.where(alive)[0]
            if len(idx):
                pp = model.p
                eps = float(np.mean((a[idx]-model.a0)/model.a0 + pp.chi_axial_projection*s[idx]/model.a0))
                plastic = float(np.mean(np.abs(model.well_index(s[idx]))))
                barriers = [model.opening_barrier(float(s[k,0]), force) for k in idx[:min(12,len(idx))]]
                barrier = float(np.mean(barriers)) if barriers else np.nan
            else:
                eps, plastic, barrier = np.nan, np.nan, 0.0
            rec_t.append(t); rec_f.append(force); rec_eps.append(eps)
            valid_count = int(np.count_nonzero(~invalid))
            rec_alive.append(float(np.count_nonzero(alive) / valid_count) if valid_count else np.nan)
            rec_plastic.append(plastic)
            rec_barrier.append(barrier)

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
                amin_b, asad_b, bound_b = model.opening_saddle_batch(s[idx], force)
                lost_spinodal = ~np.all(bound_b, axis=1)
                crossed = np.any(a[idx] >= asad_b, axis=1)
                fail = lost_spinodal | crossed | ~np.all(np.isfinite(amin_b), axis=1)
                first_passage_time[idx[fail]] = (step + 1) * solver_p.dt
                alive[idx[fail]] = False

    return {
        "model": model,
        "time": np.asarray(rec_t),
        "force": np.asarray(rec_f),
        "strain": np.asarray(rec_eps),
        "survival": np.asarray(rec_alive),
        "plastic_well_activity": np.asarray(rec_plastic),
        "opening_barrier": np.asarray(rec_barrier),
        "first_passage_time": first_passage_time,
        "invalid_trajectory": invalid,
        "observation_end_time": (nt - 1) * solver_p.dt,
    }
