"""N=3 six-dimensional probability tensor / tensor-train prototype.

This module does not yet time-integrate the Smoluchowski equation in tensor
train form.  It constructs the full correlated N=3 Gibbs initial density on a
small verification grid and compresses that density with TT-SVD.  The purpose
is to measure ranks, compression, mass error, positivity error, and correlation
retention before a TT time integrator is trusted.

The physical state is q=(a1,a2,a3,s1,s2,s3).  Rank one is never imposed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import ModelParams, TwoRowLJ
from .tensor_train import TensorTrain, relative_frobenius_error, tt_svd


@dataclass(frozen=True)
class Grid6DParams:
    n_a: int = 7
    n_s: int = 7
    a_lower_factor: float = 0.86
    a_upper: float = 1.50
    s_wells: int = 1


@dataclass(frozen=True)
class Grid6D:
    a: np.ndarray
    s: np.ndarray
    da: float
    ds: float

    @property
    def shape(self) -> tuple[int, ...]:
        return (
            self.a.size,
            self.a.size,
            self.a.size,
            self.s.size,
            self.s.size,
            self.s.size,
        )

    @property
    def cell_volume(self) -> float:
        return float(self.da**3 * self.ds**3)


def build_grid_6d(model: TwoRowLJ, params: Grid6DParams) -> Grid6D:
    if model.p.n_cells != 3:
        raise ValueError("probability_tt_6d requires ModelParams(n_cells=3)")
    if params.n_a < 3 or params.n_s < 3:
        raise ValueError("n_a and n_s must both be at least 3")
    if params.s_wells < 1 or params.s_wells % 2 == 0:
        raise ValueError("s_wells must be a positive odd integer")

    a_low = max(model.p.a_min * 1.001, params.a_lower_factor * model.a0)
    if params.a_upper <= a_low:
        raise ValueError("a_upper must exceed the lower a-domain bound")
    a_edges = np.linspace(a_low, params.a_upper, params.n_a + 1)
    s_half_span = 0.5 * params.s_wells * model.p.b
    s_edges = np.linspace(-s_half_span, s_half_span, params.n_s + 1)
    return Grid6D(
        a=0.5 * (a_edges[:-1] + a_edges[1:]),
        s=0.5 * (s_edges[:-1] + s_edges[1:]),
        da=float(a_edges[1] - a_edges[0]),
        ds=float(s_edges[1] - s_edges[0]),
    )


def _mesh_6d(grid: Grid6D) -> tuple[np.ndarray, ...]:
    return np.meshgrid(
        grid.a,
        grid.a,
        grid.a,
        grid.s,
        grid.s,
        grid.s,
        indexing="ij",
    )


def _state_arrays(grid: Grid6D) -> tuple[np.ndarray, np.ndarray]:
    a1, a2, a3, s1, s2, s3 = _mesh_6d(grid)
    a = np.column_stack((a1.ravel(), a2.ravel(), a3.ravel()))
    s = np.column_stack((s1.ravel(), s2.ravel(), s3.ravel()))
    return a, s


def energy_grid_6d(model: TwoRowLJ, grid: Grid6D, force: float) -> np.ndarray:
    a, s = _state_arrays(grid)
    energy, _, _ = model.energy_gradient_batch(a, s, float(force))
    return energy.reshape(grid.shape)


def opening_intact_mask_6d(model: TwoRowLJ, grid: Grid6D, force: float) -> np.ndarray:
    a, s = _state_arrays(grid)
    _, saddle, bound = model.opening_saddle_batch(s, float(force))
    intact = np.all(bound & (a < saddle), axis=1)
    return intact.reshape(grid.shape)


def initial_gibbs_density_6d(
    model: TwoRowLJ,
    grid: Grid6D,
    *,
    preload_force: float = 0.0,
) -> np.ndarray:
    """Build the full correlated N=3 Gibbs density on a small verification grid."""

    if model.p.kT <= 0.0:
        raise ValueError("finite-temperature Gibbs density requires kT > 0")
    energy = energy_grid_6d(model, grid, preload_force)
    mask = opening_intact_mask_6d(model, grid, preload_force)
    _, _, _, s1, s2, s3 = _mesh_6d(grid)
    mask &= (
        (np.abs(s1) < 0.5 * model.p.b)
        & (np.abs(s2) < 0.5 * model.p.b)
        & (np.abs(s3) < 0.5 * model.p.b)
    )
    finite = np.isfinite(energy)
    mask &= finite
    if not np.any(mask):
        raise ValueError("N=3 Gibbs verification basin contains no valid cells")

    g0 = float(np.min(energy[mask]))
    density = np.zeros_like(energy)
    exponent = -(energy - g0) / model.p.kT
    density[mask] = np.exp(np.clip(exponent[mask], -745.0, 0.0))
    z = float(np.sum(density) * grid.cell_volume)
    if not np.isfinite(z) or z <= 0.0:
        raise FloatingPointError("failed to normalize N=3 Gibbs density")
    return density / z


def compress_initial_gibbs_6d(
    *,
    model_params: ModelParams | None = None,
    grid_params: Grid6DParams = Grid6DParams(),
    preload_force: float = 0.0,
    relative_tolerance: float = 1.0e-7,
    max_rank: int | None = None,
) -> dict[str, object]:
    """Construct and TT-compress the correlated N=3 initial probability field."""

    p = model_params or ModelParams(n_cells=3)
    if p.n_cells != 3:
        raise ValueError("N=3 TT prototype requires n_cells=3")
    model = TwoRowLJ(p)
    model._build_opening_table()
    grid = build_grid_6d(model, grid_params)
    density = initial_gibbs_density_6d(model, grid, preload_force=preload_force)
    tt = tt_svd(
        density,
        relative_tolerance=relative_tolerance,
        max_rank=max_rank,
    )
    reconstructed = tt.reconstruct()

    exact_mass = float(np.sum(density) * grid.cell_volume)
    reconstructed_mass = float(np.sum(reconstructed) * grid.cell_volume)
    negative_mass = float(
        np.sum(np.maximum(-reconstructed, 0.0)) * grid.cell_volume
    )
    positive_mass = float(
        np.sum(np.maximum(reconstructed, 0.0)) * grid.cell_volume
    )

    return {
        "model": model,
        "grid": grid,
        "density": density,
        "tensor_train": tt,
        "reconstructed_density": reconstructed,
        "tt_ranks": tt.ranks,
        "dense_storage": tt.dense_storage,
        "tt_storage": tt.storage,
        "compression_ratio": tt.compression_ratio,
        "relative_frobenius_error": relative_frobenius_error(density, reconstructed),
        "exact_mass": exact_mass,
        "reconstructed_mass": reconstructed_mass,
        "negative_mass": negative_mass,
        "positive_mass": positive_mass,
    }
