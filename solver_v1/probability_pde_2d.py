"""Deterministic probability-density reference solver for Theory Core v1.

This module solves the N=1, two-coordinate state q=(a,s) directly as a
Smoluchowski/Fokker--Planck PDE.  It contains no random-number sampling.

The purpose is numerical validation before compressing the full N=3,
six-dimensional probability density with a sparse-grid or tensor method.

Discretization
--------------
A cell-centred conservative finite-volume method is used.  Interior fluxes use
the Scharfetter--Gummel exponential fitting formula,

    J = D/h [ B(beta*DeltaG) P_L - B(-beta*DeltaG) P_R ],

where B(x)=x/(exp(x)-1), D=M*kT.  This flux has the correct discrete Gibbs
stationary ratio on a fixed energy landscape and is much more robust than a
naive centred drift-diffusion stencil in steep LJ regions.

Crack first passage
-------------------
The mechanically defined local opening saddle from ``TwoRowLJ`` defines the
intact normal-opening basin.  Probability outside that basin is absorbed and
is never renormalized into the surviving density.  The cumulative initiated
probability is therefore 1-S, where S is the remaining intact probability
mass.

This N=1 module is a gold-standard/reference calculation, not the production
N=3 solver and not a calibrated pure-Al fatigue-life predictor.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .model import ModelParams, TwoRowLJ


@dataclass(frozen=True)
class Grid2DParams:
    """Cell-centred truncated domain for q=(a,s)."""

    n_a: int = 61
    n_s: int = 121
    a_lower_factor: float = 0.82
    a_upper: float = 1.80
    s_wells: int = 5


@dataclass(frozen=True)
class PDETimeParams:
    """Explicit conservative integration controls.

    ``max_dt`` is only an upper bound.  The solver computes a positivity CFL
    bound from the Scharfetter--Gummel transition rates at every step.
    """

    max_dt: float = 2.0e-3
    cfl: float = 0.45
    record_interval: float = 2.0e-2
    negative_tolerance: float = 2.0e-12


@dataclass(frozen=True)
class CyclicLoad2D:
    force_min: float = 0.0
    force_max: float = 3.2
    period: float = 1.0
    cycles: float = 1.0
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
class Grid2D:
    a: np.ndarray
    s: np.ndarray
    da: float
    ds: float

    @property
    def cell_volume(self) -> float:
        return float(self.da * self.ds)


def _bernoulli(x: np.ndarray) -> np.ndarray:
    """Stable Bernoulli function B(x)=x/(exp(x)-1)."""

    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1.0e-5
    positive = x > 50.0
    negative = x < -50.0
    middle = ~(small | positive | negative)

    xs = x[small]
    out[small] = 1.0 - 0.5 * xs + xs * xs / 12.0 - xs**4 / 720.0

    xp = x[positive]
    out[positive] = xp * np.exp(-xp) / (1.0 - np.exp(-xp))

    xn = x[negative]
    out[negative] = -xn / (1.0 - np.exp(xn))

    xm = x[middle]
    out[middle] = xm / np.expm1(xm)
    return out


def build_grid(model: TwoRowLJ, params: Grid2DParams) -> Grid2D:
    if model.p.n_cells != 1:
        raise ValueError("probability_pde_2d requires ModelParams(n_cells=1)")
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
    a = 0.5 * (a_edges[:-1] + a_edges[1:])
    s = 0.5 * (s_edges[:-1] + s_edges[1:])
    return Grid2D(
        a=a,
        s=s,
        da=float(a_edges[1] - a_edges[0]),
        ds=float(s_edges[1] - s_edges[0]),
    )


def _mesh(grid: Grid2D) -> tuple[np.ndarray, np.ndarray]:
    return np.meshgrid(grid.a, grid.s, indexing="ij")


def energy_grid(model: TwoRowLJ, grid: Grid2D, force: float) -> np.ndarray:
    """Evaluate the full N=1 effective energy G(a,s;force) at cell centres."""

    aa, ss = _mesh(grid)
    energy, _, _ = model.energy_gradient_batch(
        aa.reshape(-1, 1), ss.reshape(-1, 1), float(force)
    )
    return energy.reshape(aa.shape)


def opening_intact_mask(
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> np.ndarray:
    """Return cells that remain inside the mechanically bound opening basin."""

    _, saddle, bound = model.opening_saddle_batch(grid.s[:, None], float(force))
    saddle = saddle[:, 0]
    bound = bound[:, 0]
    return bound[None, :] & (grid.a[:, None] < saddle[None, :])


def initial_gibbs_density(
    model: TwoRowLJ,
    grid: Grid2D,
    *,
    preload_force: float = 0.0,
    principal_well_only: bool = True,
) -> np.ndarray:
    """Conditional Gibbs density in the declared intact initial basin.

    The baseline uses the principal configurational well |s|<b/2.  This is a
    conditional metastable initial ensemble, not an imposed Gaussian spacing
    law and not a product closure.
    """

    if model.p.kT <= 0.0:
        raise ValueError("finite-temperature Gibbs density requires kT > 0")

    g = energy_grid(model, grid, preload_force)
    mask = opening_intact_mask(model, grid, preload_force)
    if principal_well_only:
        mask &= np.abs(grid.s[None, :]) < 0.5 * model.p.b
    if not np.any(mask):
        raise ValueError("initial Gibbs basin contains no grid cells")

    g0 = float(np.min(g[mask]))
    exponent = -(g - g0) / model.p.kT
    density = np.zeros_like(g)
    density[mask] = np.exp(np.clip(exponent[mask], -745.0, 0.0))
    z = float(np.sum(density) * grid.cell_volume)
    if not np.isfinite(z) or z <= 0.0:
        raise FloatingPointError("failed to normalize the initial Gibbs density")
    return density / z


def _sg_rates_and_rhs(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
) -> tuple[np.ndarray, float]:
    """Return conservative SG semi-discrete RHS and maximum outgoing rate."""

    if model.p.kT <= 0.0:
        raise ValueError("Smoluchowski diffusion requires kT > 0")

    beta = 1.0 / model.p.kT
    d_a = model.p.mobility_a * model.p.kT
    d_s = model.p.mobility_s * model.p.kT
    rhs = np.zeros_like(density)
    outgoing = np.zeros_like(density)

    # a-direction faces; external domain faces are reflecting (zero flux).
    psi_a = beta * (energy[1:, :] - energy[:-1, :])
    bp_a = _bernoulli(psi_a)
    bm_a = _bernoulli(-psi_a)
    flux_a = (d_a / grid.da) * (
        bp_a * density[:-1, :] - bm_a * density[1:, :]
    )
    rhs[:-1, :] -= flux_a / grid.da
    rhs[1:, :] += flux_a / grid.da
    outgoing[:-1, :] += (d_a / grid.da**2) * bp_a
    outgoing[1:, :] += (d_a / grid.da**2) * bm_a

    # s-direction faces.  The computational truncation is reflecting; a
    # boundary-mass diagnostic tells us when s_wells must be enlarged.
    psi_s = beta * (energy[:, 1:] - energy[:, :-1])
    bp_s = _bernoulli(psi_s)
    bm_s = _bernoulli(-psi_s)
    flux_s = (d_s / grid.ds) * (
        bp_s * density[:, :-1] - bm_s * density[:, 1:]
    )
    rhs[:, :-1] -= flux_s / grid.ds
    rhs[:, 1:] += flux_s / grid.ds
    outgoing[:, :-1] += (d_s / grid.ds**2) * bp_s
    outgoing[:, 1:] += (d_s / grid.ds**2) * bm_s

    return rhs, float(np.max(outgoing))


def _absorb_outside_opening_basin(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> tuple[np.ndarray, float]:
    mask = opening_intact_mask(model, grid, force)
    removed = float(np.sum(density[~mask]) * grid.cell_volume)
    if removed == 0.0:
        return density, 0.0
    out = density.copy()
    out[~mask] = 0.0
    return out, removed


def observables(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
    *,
    first_passage_flux: float = 0.0,
) -> dict[str, float]:
    """Compute survivor-conditioned macroscopic and diagnostic observables."""

    volume = grid.cell_volume
    survival = float(np.sum(density) * volume)
    aa, ss = _mesh(grid)
    if survival > 0.0:
        conditional = density / survival
        strain_field = ((aa - model.a0) + model.p.chi_axial_projection * ss) / model.a0
        strain = float(np.sum(conditional * strain_field) * volume)
        nwell = np.floor((ss + 0.5 * model.p.b) / model.p.b)
        well_activity = float(np.sum(conditional * np.abs(nwell)) * volume)
    else:
        strain = np.nan
        well_activity = np.nan

    # Probability near artificial truncation boundaries is a convergence flag,
    # not a physical observable.
    edge = max(1, min(2, grid.s.size // 4))
    s_boundary_mass = float(
        (np.sum(density[:, :edge]) + np.sum(density[:, -edge:])) * volume
    )
    a_lower_mass = float(np.sum(density[:edge, :]) * volume)
    a_upper_mass = float(np.sum(density[-edge:, :]) * volume)

    return {
        "force": float(force),
        "survival": survival,
        "initiation_probability": max(0.0, 1.0 - survival),
        "first_passage_flux": float(first_passage_flux),
        "strain": strain,
        "plastic_well_activity": well_activity,
        "s_truncation_boundary_mass": s_boundary_mass,
        "a_lower_boundary_mass": a_lower_mass,
        "a_upper_boundary_mass": a_upper_mass,
    }


def run_probability_pde_2d(
    *,
    model_params: ModelParams | None = None,
    grid_params: Grid2DParams = Grid2DParams(),
    time_params: PDETimeParams = PDETimeParams(),
    load: CyclicLoad2D = CyclicLoad2D(),
    preload_force: float = 0.0,
) -> dict[str, np.ndarray | Grid2D | TwoRowLJ]:
    """Solve the deterministic N=1 probability PDE under a cyclic load."""

    p = model_params or ModelParams(n_cells=1)
    if p.n_cells != 1:
        raise ValueError("N=1 reference solver requires n_cells=1")
    if time_params.max_dt <= 0.0 or not 0.0 < time_params.cfl < 1.0:
        raise ValueError("invalid PDE time controls")
    if time_params.record_interval <= 0.0:
        raise ValueError("record_interval must be positive")
    if load.period <= 0.0 or load.cycles < 0.0:
        raise ValueError("invalid cyclic load duration")

    model = TwoRowLJ(p)
    model._build_opening_table()
    grid = build_grid(model, grid_params)
    density = initial_gibbs_density(
        model, grid, preload_force=preload_force, principal_well_only=True
    )

    records: dict[str, list[float]] = {
        key: []
        for key in (
            "time",
            "force",
            "survival",
            "initiation_probability",
            "first_passage_flux",
            "strain",
            "plastic_well_activity",
            "s_truncation_boundary_mass",
            "a_lower_boundary_mass",
            "a_upper_boundary_mass",
        )
    }

    t = 0.0
    next_record = 0.0
    last_flux = 0.0
    duration = load.duration

    def append_record(now: float, force: float) -> None:
        obs = observables(
            density, model, grid, force, first_passage_flux=last_flux
        )
        records["time"].append(float(now))
        for name, value in obs.items():
            records[name].append(float(value))

    while True:
        force = load.value(t)
        density, instant_loss = _absorb_outside_opening_basin(
            density, model, grid, force
        )
        if instant_loss > 0.0:
            # A moving spinodal/dividing surface can instantaneously remove
            # mass.  Report it over the next numerical interval as a flux-like
            # rate; cumulative initiation remains exactly 1-S.
            last_flux = instant_loss / max(time_params.max_dt, np.finfo(float).eps)

        if t + 1.0e-14 >= next_record:
            append_record(t, force)
            next_record += time_params.record_interval

        if t >= duration - 1.0e-14 or np.sum(density) == 0.0:
            break

        energy = energy_grid(model, grid, force)
        rhs, max_outgoing_rate = _sg_rates_and_rhs(density, energy, model, grid)
        if max_outgoing_rate > 0.0:
            stable_dt = time_params.cfl / max_outgoing_rate
        else:
            stable_dt = time_params.max_dt
        dt = min(time_params.max_dt, stable_dt, duration - t)
        if not np.isfinite(dt) or dt <= 0.0:
            raise FloatingPointError("failed to obtain a positive stable PDE step")

        trial = density + dt * rhs
        min_value = float(np.min(trial))
        if min_value < -time_params.negative_tolerance:
            raise FloatingPointError(
                f"probability density became negative ({min_value:.3e}); refine dt/grid"
            )
        trial = np.maximum(trial, 0.0)

        next_t = t + dt
        next_force = load.value(next_t)
        before_absorb = float(np.sum(trial) * grid.cell_volume)
        trial, removed = _absorb_outside_opening_basin(
            trial, model, grid, next_force
        )
        after_absorb = float(np.sum(trial) * grid.cell_volume)
        # The conservative SG update should preserve mass before absorption.
        numerical_mass_error = abs(before_absorb - float(np.sum(density) * grid.cell_volume))
        if numerical_mass_error > 5.0e-10:
            raise FloatingPointError(
                f"finite-volume mass conservation error {numerical_mass_error:.3e}"
            )

        last_flux = removed / dt
        density = trial
        t = next_t

    # Ensure the terminal state is present even if it does not land exactly on
    # the record interval.
    if not records["time"] or abs(records["time"][-1] - t) > 1.0e-12:
        append_record(t, load.value(t))

    return {
        **{name: np.asarray(values, dtype=float) for name, values in records.items()},
        "density": density,
        "grid": grid,
        "model": model,
    }
