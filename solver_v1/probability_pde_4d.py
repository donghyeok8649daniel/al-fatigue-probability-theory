"""Dense N=2 correlated probability-PDE reference solver for Theory Core v1.

State
-----
q = (a1, a2, s1, s2)

This module is the first explicit cell-cell correlation reference.  It solves
exactly the same deterministic Smoluchowski probability equation as the N=1
2D gold-standard solver, but on a small four-dimensional tensor grid.  It is
intentionally limited to modest grids and short convergence cases; the N=3
production solver must use a compressed representation.

No random-number sampling, Monte Carlo resampling, product closure, named life
distribution, or empirical crack-probability law appears here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .model import ModelParams, TwoRowLJ
from .probability_pde_2d import _bernoulli


@dataclass(frozen=True)
class Grid4DParams:
    n_a: int = 13
    n_s: int = 17
    a_lower_factor: float = 0.84
    a_upper: float = 1.65
    s_wells: int = 3


@dataclass(frozen=True)
class PDE4DTimeParams:
    max_dt: float = 5.0e-4
    cfl: float = 0.35
    record_interval: float = 5.0e-3
    negative_tolerance: float = 5.0e-12


@dataclass(frozen=True)
class CyclicLoad4D:
    force_min: float = 0.0
    force_max: float = 3.2
    period: float = 1.0
    cycles: float = 0.1
    phase_radians: float = 0.0
    value_function: Callable[[float], float] | None = None

    @property
    def duration(self) -> float:
        return float(self.period * self.cycles)

    def value(self, time: float) -> float:
        if self.value_function is not None:
            return float(self.value_function(time))
        midpoint = 0.5 * (self.force_max + self.force_min)
        amplitude = 0.5 * (self.force_max - self.force_min)
        return float(
            midpoint
            + amplitude
            * np.sin(2.0 * np.pi * time / self.period + self.phase_radians)
        )


@dataclass(frozen=True)
class Grid4D:
    a: np.ndarray
    s: np.ndarray
    da: float
    ds: float

    @property
    def cell_volume(self) -> float:
        return float(self.da**2 * self.ds**2)

    @property
    def shape(self) -> tuple[int, int, int, int]:
        return (self.a.size, self.a.size, self.s.size, self.s.size)


def build_grid_4d(model: TwoRowLJ, params: Grid4DParams) -> Grid4D:
    if model.p.n_cells != 2:
        raise ValueError("probability_pde_4d requires ModelParams(n_cells=2)")
    if params.n_a < 5 or params.n_s < 5:
        raise ValueError("n_a and n_s must both be at least 5")
    if params.s_wells < 1 or params.s_wells % 2 == 0:
        raise ValueError("s_wells must be a positive odd integer")

    a_low = max(model.p.a_min * 1.001, params.a_lower_factor * model.a0)
    if params.a_upper <= a_low:
        raise ValueError("a_upper must exceed the lower a-domain bound")

    a_edges = np.linspace(a_low, params.a_upper, params.n_a + 1)
    s_half_span = 0.5 * params.s_wells * model.p.b
    s_edges = np.linspace(-s_half_span, s_half_span, params.n_s + 1)
    return Grid4D(
        a=0.5 * (a_edges[:-1] + a_edges[1:]),
        s=0.5 * (s_edges[:-1] + s_edges[1:]),
        da=float(a_edges[1] - a_edges[0]),
        ds=float(s_edges[1] - s_edges[0]),
    )


def _coordinate_mesh(grid: Grid4D) -> tuple[np.ndarray, ...]:
    return np.meshgrid(grid.a, grid.a, grid.s, grid.s, indexing="ij")


def energy_grid_4d(model: TwoRowLJ, grid: Grid4D, force: float) -> np.ndarray:
    a1, a2, s1, s2 = _coordinate_mesh(grid)
    a = np.column_stack((a1.ravel(), a2.ravel()))
    s = np.column_stack((s1.ravel(), s2.ravel()))
    energy, _, _ = model.energy_gradient_batch(a, s, float(force))
    return energy.reshape(grid.shape)


def opening_intact_mask_4d(model: TwoRowLJ, grid: Grid4D, force: float) -> np.ndarray:
    a1, a2, s1, s2 = _coordinate_mesh(grid)
    a = np.column_stack((a1.ravel(), a2.ravel()))
    s = np.column_stack((s1.ravel(), s2.ravel()))
    _, saddle, bound = model.opening_saddle_batch(s, float(force))
    intact = np.all(bound & (a < saddle), axis=1)
    return intact.reshape(grid.shape)


def initial_gibbs_density_4d(
    model: TwoRowLJ,
    grid: Grid4D,
    *,
    preload_force: float = 0.0,
    principal_well_only: bool = True,
) -> np.ndarray:
    """Correlated conditional Gibbs density in the intact N=2 basin."""

    if model.p.kT <= 0.0:
        raise ValueError("finite-temperature Gibbs density requires kT > 0")
    energy = energy_grid_4d(model, grid, preload_force)
    mask = opening_intact_mask_4d(model, grid, preload_force)
    if principal_well_only:
        _, _, s1, s2 = _coordinate_mesh(grid)
        mask &= (np.abs(s1) < 0.5 * model.p.b) & (np.abs(s2) < 0.5 * model.p.b)
    if not np.any(mask):
        raise ValueError("initial N=2 Gibbs basin contains no grid cells")

    g0 = float(np.min(energy[mask]))
    density = np.zeros_like(energy)
    exponent = -(energy - g0) / model.p.kT
    density[mask] = np.exp(np.clip(exponent[mask], -745.0, 0.0))
    z = float(np.sum(density) * grid.cell_volume)
    if not np.isfinite(z) or z <= 0.0:
        raise FloatingPointError("failed to normalize N=2 Gibbs density")
    return density / z


def _sg_rhs_4d(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid4D,
) -> tuple[np.ndarray, float]:
    """Conservative four-dimensional Scharfetter--Gummel semi-discretization."""

    if model.p.kT <= 0.0:
        raise ValueError("Smoluchowski diffusion requires kT > 0")
    beta = 1.0 / model.p.kT
    rhs = np.zeros_like(density)
    outgoing = np.zeros_like(density)

    axis_data = (
        (0, model.p.mobility_a * model.p.kT, grid.da),
        (1, model.p.mobility_a * model.p.kT, grid.da),
        (2, model.p.mobility_s * model.p.kT, grid.ds),
        (3, model.p.mobility_s * model.p.kT, grid.ds),
    )

    for axis, diffusivity, spacing in axis_data:
        left = [slice(None)] * 4
        right = [slice(None)] * 4
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        left_t = tuple(left)
        right_t = tuple(right)

        psi = beta * (energy[right_t] - energy[left_t])
        bp = _bernoulli(psi)
        bm = _bernoulli(-psi)
        flux = (diffusivity / spacing) * (
            bp * density[left_t] - bm * density[right_t]
        )
        rhs[left_t] -= flux / spacing
        rhs[right_t] += flux / spacing
        outgoing[left_t] += (diffusivity / spacing**2) * bp
        outgoing[right_t] += (diffusivity / spacing**2) * bm

    return rhs, float(np.max(outgoing))


def _absorb_4d(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid4D,
    force: float,
) -> tuple[np.ndarray, float]:
    mask = opening_intact_mask_4d(model, grid, force)
    removed = float(np.sum(density[~mask]) * grid.cell_volume)
    if removed == 0.0:
        return density, 0.0
    out = density.copy()
    out[~mask] = 0.0
    return out, removed


def observables_4d(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid4D,
    force: float,
    *,
    first_passage_flux: float = 0.0,
) -> dict[str, float]:
    volume = grid.cell_volume
    survival = float(np.sum(density) * volume)
    a1, a2, s1, s2 = _coordinate_mesh(grid)

    if survival > 0.0:
        p = density / survival
        mean_a1 = float(np.sum(p * a1) * volume)
        mean_a2 = float(np.sum(p * a2) * volume)
        mean_s1 = float(np.sum(p * s1) * volume)
        mean_s2 = float(np.sum(p * s2) * volume)
        cov_a12 = float(np.sum(p * (a1 - mean_a1) * (a2 - mean_a2)) * volume)
        cov_s12 = float(np.sum(p * (s1 - mean_s1) * (s2 - mean_s2)) * volume)
        strain_field = (
            (a1 - model.a0)
            + (a2 - model.a0)
            + model.p.chi_axial_projection * (s1 + s2)
        ) / (2.0 * model.a0)
        strain = float(np.sum(p * strain_field) * volume)

        # Exact one-cell joint marginals from the full correlated N=2 density.
        p1 = np.sum(p, axis=(1, 3)) * grid.da * grid.ds
        p2 = np.sum(p, axis=(0, 2)) * grid.da * grid.ds
        product = p1[:, None, :, None] * p2[None, :, None, :]
        product_l1_error = float(np.sum(np.abs(p - product)) * volume)
    else:
        strain = np.nan
        cov_a12 = np.nan
        cov_s12 = np.nan
        product_l1_error = np.nan

    return {
        "force": float(force),
        "survival": survival,
        "initiation_probability": max(0.0, 1.0 - survival),
        "first_passage_flux": float(first_passage_flux),
        "strain": strain,
        "cov_a12": cov_a12,
        "cov_s12": cov_s12,
        "product_closure_l1_error": product_l1_error,
    }


def run_probability_pde_4d(
    *,
    model_params: ModelParams | None = None,
    grid_params: Grid4DParams = Grid4DParams(),
    time_params: PDE4DTimeParams = PDE4DTimeParams(),
    load: CyclicLoad4D = CyclicLoad4D(),
    preload_force: float = 0.0,
) -> dict[str, object]:
    """Run the dense N=2 correlated probability reference solve."""

    p = model_params or ModelParams(n_cells=2)
    if p.n_cells != 2:
        raise ValueError("N=2 reference solver requires n_cells=2")
    if time_params.max_dt <= 0.0 or not 0.0 < time_params.cfl < 1.0:
        raise ValueError("invalid 4D PDE time controls")
    if time_params.record_interval <= 0.0:
        raise ValueError("record_interval must be positive")

    model = TwoRowLJ(p)
    model._build_opening_table()
    grid = build_grid_4d(model, grid_params)
    density = initial_gibbs_density_4d(model, grid, preload_force=preload_force)

    names = (
        "time",
        "force",
        "survival",
        "initiation_probability",
        "first_passage_flux",
        "strain",
        "cov_a12",
        "cov_s12",
        "product_closure_l1_error",
    )
    records: dict[str, list[float]] = {name: [] for name in names}

    t = 0.0
    next_record = 0.0
    last_flux = 0.0
    duration = load.duration

    def record(now: float, force: float) -> None:
        obs = observables_4d(
            density, model, grid, force, first_passage_flux=last_flux
        )
        records["time"].append(float(now))
        for name, value in obs.items():
            records[name].append(float(value))

    while True:
        force = load.value(t)
        density, instant_loss = _absorb_4d(density, model, grid, force)
        if instant_loss > 0.0:
            last_flux = instant_loss / max(time_params.max_dt, np.finfo(float).eps)

        if t + 1.0e-14 >= next_record:
            record(t, force)
            next_record += time_params.record_interval

        if t >= duration - 1.0e-14 or float(np.sum(density)) == 0.0:
            break

        energy = energy_grid_4d(model, grid, force)
        rhs, max_rate = _sg_rhs_4d(density, energy, model, grid)
        stable_dt = time_params.max_dt
        if max_rate > 0.0:
            stable_dt = min(stable_dt, time_params.cfl / max_rate)
        dt = min(stable_dt, duration - t)
        if dt <= 0.0:
            break

        updated = density + dt * rhs
        minimum = float(np.min(updated))
        if minimum < -time_params.negative_tolerance:
            raise FloatingPointError(
                f"4D probability density lost positivity: min={minimum:.3e}"
            )
        updated = np.maximum(updated, 0.0)

        next_t = t + dt
        next_force = load.value(next_t)
        before = float(np.sum(updated) * grid.cell_volume)
        updated, removed = _absorb_4d(updated, model, grid, next_force)
        after = float(np.sum(updated) * grid.cell_volume)
        last_flux = max(0.0, (before - after) / dt)
        if removed > 0.0:
            last_flux = max(last_flux, removed / dt)
        density = updated
        t = next_t

    if not records["time"] or records["time"][-1] < t - 1.0e-12:
        record(t, load.value(t))

    return {
        **{name: np.asarray(values, dtype=float) for name, values in records.items()},
        "density": density,
        "grid": grid,
        "model": model,
    }
