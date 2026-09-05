"""Deterministic finite-chain solver for the active normal theory.

The probability measure is supplied as a weighted discrete initial measure and
is pushed forward by the same deterministic ODE for every atom.  No thermal
noise, mobility, fitted life distribution, or Monte Carlo resampling appears
in this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np


@dataclass(frozen=True)
class InitialMeasureAtom:
    """One atom of a declared full-state initial measure."""

    weight: float
    spacing: np.ndarray
    spacing_rate: np.ndarray


@dataclass(frozen=True)
class NormalChainParams:
    n_cells: int = 40
    repulsive_exponent: float = 12.19
    attractive_exponent: float = 6.0

    @property
    def lambda_c(self) -> float:
        m = self.repulsive_exponent
        n = self.attractive_exponent
        return float(((m + 1.0) / (n + 1.0)) ** (1.0 / (m - n)))


@dataclass(frozen=True)
class DeterministicRunParams:
    dt: float = 0.01
    duration: float = 10.0
    record_stride: int = 10


@dataclass(frozen=True)
class LoadParams:
    """Prescribed dimensionless normal traction history."""

    force_max: float = 0.0
    force_min: float = 0.0
    period: float = 10.0
    cycles: int = 10
    phase_radians: float = 0.0
    value_function: Callable[[float], float] | None = None

    def value(self, time: float) -> float:
        if self.value_function is not None:
            return float(self.value_function(time))
        midpoint = 0.5 * (self.force_max + self.force_min)
        amplitude = 0.5 * (self.force_max - self.force_min)
        return float(
            midpoint
            + amplitude * np.sin(2.0 * np.pi * time / self.period + self.phase_radians)
        )


def phi(spacing: np.ndarray, params: NormalChainParams) -> np.ndarray:
    m = params.repulsive_exponent
    n = params.attractive_exponent
    spacing = np.asarray(spacing, dtype=float)
    return (
        spacing**(-m) / (m * (m - n))
        - spacing**(-n) / (n * (m - n))
    )


def phi_prime(spacing: np.ndarray, params: NormalChainParams) -> np.ndarray:
    m = params.repulsive_exponent
    n = params.attractive_exponent
    spacing = np.asarray(spacing, dtype=float)
    return (spacing**(-n - 1.0) - spacing**(-m - 1.0)) / (m - n)


def phi_second(spacing: np.ndarray, params: NormalChainParams) -> np.ndarray:
    m = params.repulsive_exponent
    n = params.attractive_exponent
    spacing = np.asarray(spacing, dtype=float)
    return (
        (m + 1.0) * spacing**(-m - 2.0)
        - (n + 1.0) * spacing**(-n - 2.0)
    ) / (m - n)


def spacing_acceleration(
    spacing: np.ndarray,
    reduced_traction: float,
    params: NormalChainParams,
) -> np.ndarray:
    """Evaluate the exact spacing ODE stated by the active finite chain."""
    spacing = np.asarray(spacing, dtype=float)
    gradient = phi_prime(spacing, params)
    acceleration = np.empty_like(spacing)
    if spacing.shape[-1] == 1:
        acceleration[..., 0] = reduced_traction - gradient[..., 0]
        return acceleration
    acceleration[..., 0] = gradient[..., 1] - gradient[..., 0]
    if spacing.shape[-1] > 2:
        acceleration[..., 1:-1] = (
            gradient[..., 2:] - 2.0 * gradient[..., 1:-1] + gradient[..., :-2]
        )
    acceleration[..., -1] = (
        reduced_traction + gradient[..., -2] - 2.0 * gradient[..., -1]
    )
    return acceleration


def delta_initial_measure(params: NormalChainParams) -> tuple[InitialMeasureAtom, ...]:
    """Return the declared ideal baseline measure mu0=delta_(lambda=1,c=0)."""
    return (
        InitialMeasureAtom(
            weight=1.0,
            spacing=np.ones(params.n_cells, dtype=float),
            spacing_rate=np.zeros(params.n_cells, dtype=float),
        ),
    )


def _measure_arrays(
    params: NormalChainParams,
    initial_measure: Sequence[InitialMeasureAtom] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    atoms = delta_initial_measure(params) if initial_measure is None else tuple(initial_measure)
    if not atoms:
        raise ValueError("initial_measure must contain at least one atom")
    weights = np.asarray([atom.weight for atom in atoms], dtype=float)
    spacing = np.asarray([atom.spacing for atom in atoms], dtype=float)
    rate = np.asarray([atom.spacing_rate for atom in atoms], dtype=float)
    if spacing.shape != (len(atoms), params.n_cells) or rate.shape != spacing.shape:
        raise ValueError("every initial-measure atom must contain n_cells spacings and rates")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("initial-measure weights must be finite and nonnegative")
    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError("initial-measure weights must have positive total mass")
    if not np.all(np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("initial spacings must be finite and positive")
    if not np.all(np.isfinite(rate)):
        raise ValueError("initial spacing rates must be finite")
    return weights / total, spacing.copy(), rate.copy()


def empirical_phase_space_support(
    spacing: np.ndarray,
    rate: np.ndarray,
    measure_weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the exact delta-support representation of ``F(lambda,c,tau)``."""
    spacing = np.asarray(spacing, dtype=float)
    rate = np.asarray(rate, dtype=float)
    weights = np.asarray(measure_weights, dtype=float)
    if spacing.ndim != 2 or rate.shape != spacing.shape or weights.shape != (spacing.shape[0],):
        raise ValueError("incompatible phase-space support arrays")
    site_weights = np.broadcast_to(weights[:, None] / spacing.shape[1], spacing.shape)
    return spacing.reshape(-1).copy(), rate.reshape(-1).copy(), site_weights.reshape(-1).copy()


def run_deterministic_pushforward(
    chain: NormalChainParams,
    run: DeterministicRunParams,
    reduced_traction: Callable[[float], float],
    *,
    initial_measure: Sequence[InitialMeasureAtom] | None = None,
    record_callback=None,
    stop_requested=None,
    retain_history: bool = True,
) -> dict:
    """Push a declared discrete ``mu0`` through the conservative chain flow."""
    if chain.n_cells < 1:
        raise ValueError("n_cells must be positive")
    if run.dt <= 0.0 or run.duration < 0.0 or run.record_stride < 1:
        raise ValueError("invalid deterministic integration controls")

    weights, spacing, rate = _measure_arrays(chain, initial_measure)
    atom_count, site_count = spacing.shape
    site_weight = weights[:, None] / site_count
    crossed = spacing >= chain.lambda_c
    first_passage_time = np.full_like(spacing, np.nan, dtype=float)
    first_passage_time[crossed] = 0.0

    history: dict[str, list[float]] = {
        name: [] for name in (
            "time", "force", "strain", "normal_strain", "intrawell_strain",
            "plastic_strain", "mean_spacing_rate", "global_spacing_rate_variance",
            "survival", "local_survival", "specimen_survival",
            "intrinsic_energy", "mechanical_energy", "min_opening_eigenvalue",
            "min_plastic_eigenvalue",
        )
    }
    support_spacing: list[np.ndarray] = []
    support_rate: list[np.ndarray] = []
    support_intact_atoms: list[np.ndarray] = []
    steps = int(round(run.duration / run.dt))
    last_time = 0.0

    def aggregate(t: float, traction: float) -> dict:
        intact_atoms = ~np.any(crossed, axis=1)
        specimen_survival = float(np.sum(weights * intact_atoms))
        if specimen_survival > 0.0:
            observable_weights = weights * intact_atoms / specimen_survival
        else:
            # This terminal state is recorded immediately at first passage,
            # before any post-failure evolution can contaminate the strain.
            observable_weights = weights
        observable_site_weight = observable_weights[:, None] / site_count
        mean_spacing = float(np.sum(observable_site_weight * spacing))
        mean_rate = float(np.sum(observable_site_weight * rate))
        global_rate_variance = float(
            np.sum(observable_site_weight * (rate - mean_rate) ** 2)
        )
        local_survival = float(np.sum(site_weight * (~crossed)))
        energy = float(np.sum(
            observable_site_weight
            * (phi(spacing, chain) - phi(np.ones_like(spacing), chain))
        ))
        indices = np.arange(site_count)
        spacing_metric = site_count - np.maximum(indices[:, None], indices[None, :])
        kinetic_by_atom = 0.5 * np.einsum("ai,ij,aj->a", rate, spacing_metric, rate)
        mechanical_energy = energy + float(
            np.sum(observable_weights * kinetic_by_atom) / site_count
        )
        return {
            "time": t,
            "force": traction,
            "strain": mean_spacing - 1.0,
            "normal_strain": mean_spacing - 1.0,
            "intrawell_strain": 0.0,
            "plastic_strain": 0.0,
            "mean_spacing_rate": mean_rate,
            "global_spacing_rate_variance": global_rate_variance,
            "survival": specimen_survival,
            "local_survival": local_survival,
            "specimen_survival": specimen_survival,
            "plastic_well_activity": 0.0,
            "opening_barrier": np.nan,
            "intrinsic_energy": energy,
            "mechanical_energy": mechanical_energy,
            "min_opening_eigenvalue": float(np.min(phi_second(spacing, chain))),
            "min_plastic_eigenvalue": np.nan,
            "initial_measure_atom_count": atom_count,
            "spatial_site_count": site_count,
            "probability_resolution": float(np.min(weights[weights > 0.0])),
            "lambda_c": chain.lambda_c,
        }

    for step in range(steps + 1):
        if stop_requested is not None and stop_requested():
            break
        t = step * run.dt
        last_time = t
        traction = float(reduced_traction(t))
        if step % run.record_stride == 0:
            record = aggregate(t, traction)
            if retain_history:
                for name in history:
                    history[name].append(float(record[name]))
                support_spacing.append(spacing.copy())
                support_rate.append(rate.copy())
                support_intact_atoms.append(~np.any(crossed, axis=1))
            if record_callback is not None:
                record_callback(dict(record))
        if step == steps:
            break

        active = ~np.any(crossed, axis=1)
        if not np.any(active):
            break
        acceleration = spacing_acceleration(spacing[active], traction, chain)
        next_spacing = spacing.copy()
        next_spacing[active] = (
            spacing[active] + run.dt * rate[active] + 0.5 * run.dt**2 * acceleration
        )
        if not np.all(np.isfinite(next_spacing[active])) or np.any(next_spacing[active] <= 0.0):
            raise FloatingPointError("deterministic chain left the positive finite spacing domain")
        next_t = (step + 1) * run.dt
        next_traction = float(reduced_traction(next_t))
        next_acceleration = spacing_acceleration(next_spacing[active], next_traction, chain)
        rate[active] += 0.5 * run.dt * (acceleration + next_acceleration)
        spacing = next_spacing
        newly_crossed = (~crossed) & (spacing >= chain.lambda_c)
        first_passage_time[newly_crossed] = next_t
        crossed |= newly_crossed
        if not np.any(~np.any(crossed, axis=1)):
            last_time = next_t
            terminal = aggregate(next_t, next_traction)
            if retain_history:
                for name in history:
                    history[name].append(float(terminal[name]))
                support_spacing.append(spacing.copy())
                support_rate.append(rate.copy())
                support_intact_atoms.append(~np.any(crossed, axis=1))
            if record_callback is not None:
                record_callback(dict(terminal))
            break

    phase_weights = np.broadcast_to(weights[:, None] / site_count, spacing.shape).copy()
    return {
        **{name: np.asarray(values, dtype=float) for name, values in history.items()},
        "spacing_support": np.asarray(support_spacing, dtype=float),
        "spacing_rate_support": np.asarray(support_rate, dtype=float),
        "intact_atom_support": np.asarray(support_intact_atoms, dtype=bool),
        "phase_space_weights": phase_weights,
        "first_passage_time": first_passage_time,
        "measure_weights": weights,
        "observation_end_time": last_time,
        "lambda_c": chain.lambda_c,
    }
