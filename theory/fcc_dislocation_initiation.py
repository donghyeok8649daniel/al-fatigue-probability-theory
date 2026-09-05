"""Mechanistic FCC slip-band crack-initiation bridge for axial HCF.

This module adds an explicitly declared physical law that is absent from the
dimensionless Theory Core v1 demonstration solver: irreversible dislocation
dipole accumulation in a persistent slip band.  It uses the stress-controlled
Tanaka--Mura--Wu (TMW) crack-nucleation relation

    N_c = 2 mu w_s / ((1 - nu) b (Delta tau - 2 k)^2).

The relation is a model assumption, not an identity of the P(a,s,t) theory.
It is kept separate so the application cannot silently present a fitted S--N
curve as an output of the current LJ phase-space dynamics.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .fcc_axial_slip import FCCSlipSystem, fcc_axial_slip_systems


@dataclass(frozen=True)
class AluminumSlipInitiationParameters:
    """Explicit material inputs for the high-purity-Al TMW baseline.

    The defaults are physical-property values rather than S--N fit
    coefficients.  ``friction_stress_mpa=0`` is the ideal high-purity FCC
    lower-bound baseline; a measured cyclic friction stress may be supplied
    without changing the formulation.
    """

    lattice_parameter_nm: float = 0.40495
    surface_energy_j_m2: float = 1.14
    friction_stress_mpa: float = 0.0

    @property
    def burgers_vector_m(self) -> float:
        return self.lattice_parameter_nm * 1.0e-9 / math.sqrt(2.0)

    def validate(self) -> None:
        if not np.isfinite(self.lattice_parameter_nm) or self.lattice_parameter_nm <= 0.0:
            raise ValueError("lattice_parameter_nm must be finite and positive")
        if not np.isfinite(self.surface_energy_j_m2) or self.surface_energy_j_m2 <= 0.0:
            raise ValueError("surface_energy_j_m2 must be finite and positive")
        if not np.isfinite(self.friction_stress_mpa) or self.friction_stress_mpa < 0.0:
            raise ValueError("friction_stress_mpa must be finite and nonnegative")


@dataclass(frozen=True)
class SlipSystemLife:
    system: FCCSlipSystem
    resolved_shear_range_mpa: float
    effective_shear_range_mpa: float
    cycles_to_initiation: float


@dataclass(frozen=True)
class EmpiricalFirstPassageShape:
    """Discrete first-passage law from equally weighted theory trajectories."""

    event_cycle_multipliers: np.ndarray
    probability_mass: np.ndarray
    cumulative_probability: np.ndarray
    censor_cycle_multiplier: float
    censored_probability: float

    def quantile_multiplier(self, probability: float) -> float:
        if not 0.0 < probability < 1.0:
            raise ValueError("probability must lie strictly between zero and one")
        indices = np.flatnonzero(self.cumulative_probability >= probability)
        return float(self.event_cycle_multipliers[indices[0]]) if indices.size else math.nan


def empirical_first_passage_shape(
    first_passage_cycles: np.ndarray,
    observation_end_cycles: float,
) -> EmpiricalFirstPassageShape:
    """Normalize a theory first-passage empirical measure by its median.

    NaN entries are right-censored at ``observation_end_cycles``.  No kernel,
    named life distribution, or tail extrapolation is introduced.  At least
    half of the trajectories must cross so that a finite empirical median is
    available as the bridge to a physical characteristic life.
    """
    samples = np.asarray(first_passage_cycles, dtype=float)
    if samples.ndim != 1 or samples.size < 2:
        raise ValueError("first_passage_cycles must be a one-dimensional ensemble")
    if not np.isfinite(observation_end_cycles) or observation_end_cycles <= 0.0:
        raise ValueError("observation_end_cycles must be finite and positive")
    if np.any(np.isfinite(samples) & ((samples <= 0.0) | (samples > observation_end_cycles))):
        raise ValueError("event cycles must lie inside the observation interval")
    events = np.sort(samples[np.isfinite(samples)])
    if events.size / samples.size < 0.5:
        raise ValueError("at least half the reference trajectories must reach first passage")
    median_rank = int(math.ceil(0.5 * samples.size)) - 1
    median_cycles = float(events[median_rank])
    unique_cycles, counts = np.unique(events, return_counts=True)
    mass = counts.astype(float) / samples.size
    cdf = np.cumsum(mass)
    return EmpiricalFirstPassageShape(
        event_cycle_multipliers=unique_cycles / median_cycles,
        probability_mass=mass,
        cumulative_probability=cdf,
        censor_cycle_multiplier=observation_end_cycles / median_cycles,
        censored_probability=float(1.0 - events.size / samples.size),
    )


def tmw_cycles_to_initiation(
    resolved_shear_range_pa: float,
    shear_modulus_pa: float,
    poisson_ratio: float,
    surface_energy_j_m2: float,
    burgers_vector_m: float,
    friction_stress_pa: float = 0.0,
) -> float:
    """Return deterministic TMW crack-nucleation life in cycles.

    Infinite life here means only that the effective range
    ``Delta tau - 2 k`` is nonpositive within this particular mechanism.  It
    is not a universal fatigue limit and does not exclude other mechanisms.
    """
    values = (
        resolved_shear_range_pa,
        shear_modulus_pa,
        poisson_ratio,
        surface_energy_j_m2,
        burgers_vector_m,
        friction_stress_pa,
    )
    if not all(np.isfinite(value) for value in values):
        raise ValueError("TMW inputs must be finite")
    if resolved_shear_range_pa < 0.0:
        raise ValueError("resolved_shear_range_pa must be nonnegative")
    if shear_modulus_pa <= 0.0:
        raise ValueError("shear_modulus_pa must be positive")
    if not -1.0 < poisson_ratio < 0.5:
        raise ValueError("poisson_ratio must satisfy -1 < nu < 0.5")
    if surface_energy_j_m2 <= 0.0 or burgers_vector_m <= 0.0:
        raise ValueError("surface energy and Burgers vector must be positive")
    if friction_stress_pa < 0.0:
        raise ValueError("friction_stress_pa must be nonnegative")

    effective_range = resolved_shear_range_pa - 2.0 * friction_stress_pa
    if effective_range <= 0.0:
        return math.inf
    numerator = 2.0 * shear_modulus_pa * surface_energy_j_m2
    denominator = (1.0 - poisson_ratio) * burgers_vector_m * effective_range**2
    return numerator / denominator


def fcc_axial_tmw_lives(
    h: int,
    k: int,
    l: int,
    axial_stress_amplitude_mpa: float,
    young_modulus_pa: float,
    poisson_ratio: float,
    material: AluminumSlipInitiationParameters = AluminumSlipInitiationParameters(),
) -> tuple[SlipSystemLife, ...]:
    """Evaluate the TMW life of all FCC systems under sinusoidal axial load.

    The axial stress range is twice the entered amplitude.  Mean-stress and
    multiaxial corrections are intentionally absent; adding either would need
    a separately justified physical law.
    """
    material.validate()
    if not np.isfinite(axial_stress_amplitude_mpa) or axial_stress_amplitude_mpa < 0.0:
        raise ValueError("axial_stress_amplitude_mpa must be finite and nonnegative")
    if not np.isfinite(young_modulus_pa) or young_modulus_pa <= 0.0:
        raise ValueError("young_modulus_pa must be finite and positive")
    shear_modulus_pa = young_modulus_pa / (2.0 * (1.0 + poisson_ratio))
    lives = []
    for system in fcc_axial_slip_systems(h, k, l):
        shear_range_mpa = 2.0 * system.schmid_factor * axial_stress_amplitude_mpa
        effective_mpa = max(0.0, shear_range_mpa - 2.0 * material.friction_stress_mpa)
        cycles = tmw_cycles_to_initiation(
            shear_range_mpa * 1.0e6,
            shear_modulus_pa,
            poisson_ratio,
            material.surface_energy_j_m2,
            material.burgers_vector_m,
            material.friction_stress_mpa * 1.0e6,
        )
        lives.append(SlipSystemLife(system, shear_range_mpa, effective_mpa, cycles))
    return tuple(sorted(lives, key=lambda item: item.cycles_to_initiation))


def axial_tmw_initiation_life(*args, **kwargs) -> SlipSystemLife:
    """Return the earliest FCC slip-system initiation prediction."""
    return fcc_axial_tmw_lives(*args, **kwargs)[0]
