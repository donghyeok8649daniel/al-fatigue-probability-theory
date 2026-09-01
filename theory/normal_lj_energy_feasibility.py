# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D layer-LJ potential의 convexity와 support constraint를 이용해 crack-free energy feasibility bound를 계산한다.
# - 주요 클래스: SafeEnergyInterval
# - 주요 함수/메서드: shifted_lj_energy, safe_energy_interval, safe_energy_margin, safe_distribution_exists
#   first_energy_ceiling_crossing_time, no_compression_bound_counterexample_energy
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Exact energy-feasibility bounds for the active 1D normal-LJ theory.

The active state is P(lambda, t), where lambda=a/a0 is the local normal
spacing.  This module does not prescribe an evolution law for P.

Instead it asks a narrower question:

    Given normalization, mean stretch mu(t), mean configurational energy
    e(t), the LJ normal-opening limit lambda_c, and a mechanically justified
    lower compression bound lambda_L(t), can a crack-free distribution still
    exist?

For the normalized generalized-LJ potential, phi is convex on
(0, lambda_c].  Therefore, for probability measures supported on
[lambda_L, lambda_c] with fixed mean mu, the energy interval is exactly

    phi(mu) <= E[phi(lambda)] <= chord_phi(mu; lambda_L, lambda_c).

The lower bound is Jensen's inequality.  The upper bound follows because a
convex graph lies below its secant chord.  Both bounds are attainable, so the
interval is exact under the stated support constraint.

No Gaussian/Weibull family, damage variable, damping law, or cycle-count
state variable is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from theory.normal_lj_chain import (
    critical_stretch,
    normalized_lj_energy,
)


@dataclass(frozen=True)
class SafeEnergyInterval:
    lower_stretch: float
    mean_stretch: float
    critical_stretch: float
    minimum_energy: float
    maximum_energy: float
    lower_endpoint_weight: float
    critical_endpoint_weight: float


def shifted_lj_energy(
    stretch: np.ndarray | float,
    m: float = 12.19,
    n: float = 6.0,
):
    """LJ energy shifted so the equilibrium energy at lambda=1 is zero."""
    return normalized_lj_energy(stretch, m, n) - float(
        normalized_lj_energy(1.0, m, n)
    )


def safe_energy_interval(
    mean_stretch: float,
    lower_stretch: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
    upper_stretch: Optional[float] = None,
) -> SafeEnergyInterval:
    """Exact energy interval for crack-free measures at a fixed mean.

    Assumptions for this theorem:
      * P >= 0 and integral P = 1;
      * support(P) is contained in [lower_stretch, upper_stretch];
      * integral lambda P(lambda) d lambda = mean_stretch;
      * upper_stretch <= lambda_c, so the LJ potential is convex throughout
        the admissible interval.

    The default upper_stretch is the LJ tangent-instability stretch lambda_c.
    """
    lam_c = critical_stretch(m, n)
    upper = lam_c if upper_stretch is None else float(upper_stretch)
    lower = float(lower_stretch)
    mu = float(mean_stretch)

    if lower <= 0.0:
        raise ValueError("lower_stretch must be positive")
    if upper > lam_c + 1.0e-14:
        raise ValueError("upper_stretch must not exceed the LJ convexity limit lambda_c")
    if not (lower < upper):
        raise ValueError("require lower_stretch < upper_stretch")
    if not (lower <= mu <= upper):
        raise ValueError("mean_stretch must lie inside the admissible support")

    critical_weight = (mu - lower) / (upper - lower)
    lower_weight = 1.0 - critical_weight

    minimum = float(shifted_lj_energy(mu, m, n))
    maximum = (
        lower_weight * float(shifted_lj_energy(lower, m, n))
        + critical_weight * float(shifted_lj_energy(upper, m, n))
    )

    return SafeEnergyInterval(
        lower_stretch=lower,
        mean_stretch=mu,
        critical_stretch=upper,
        minimum_energy=minimum,
        maximum_energy=maximum,
        lower_endpoint_weight=lower_weight,
        critical_endpoint_weight=critical_weight,
    )


def safe_energy_margin(
    mean_energy: float,
    mean_stretch: float,
    lower_stretch: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Return E_safe^max - E for the current continuous-time state."""
    interval = safe_energy_interval(
        mean_stretch,
        lower_stretch,
        m=m,
        n=n,
    )
    return interval.maximum_energy - float(mean_energy)


def safe_distribution_exists(
    mean_energy: float,
    mean_stretch: float,
    lower_stretch: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
    tolerance: float = 1.0e-12,
) -> bool:
    """Whether any crack-free distribution can satisfy the stated moments."""
    interval = safe_energy_interval(
        mean_stretch,
        lower_stretch,
        m=m,
        n=n,
    )
    energy = float(mean_energy)
    return (
        interval.minimum_energy - tolerance
        <= energy
        <= interval.maximum_energy + tolerance
    )


def first_energy_ceiling_crossing_time(
    times: Sequence[float],
    mean_energies: Sequence[float],
    mean_stretches: Sequence[float],
    lower_stretches: Sequence[float],
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> Optional[float]:
    """First sampled time at which no crack-free distribution can exist.

    This is a discrete diagnostic of the continuous definition

        tau_E = inf { t : E(t) > E_safe^max(t) }.
    """
    t = np.asarray(times, dtype=float)
    e = np.asarray(mean_energies, dtype=float)
    mu = np.asarray(mean_stretches, dtype=float)
    lower = np.asarray(lower_stretches, dtype=float)

    if not (t.ndim == e.ndim == mu.ndim == lower.ndim == 1):
        raise ValueError("all inputs must be one-dimensional")
    if not (len(t) == len(e) == len(mu) == len(lower)):
        raise ValueError("all inputs must have equal length")
    if np.any(np.diff(t) < 0.0):
        raise ValueError("times must be nondecreasing")

    for ti, ei, mui, li in zip(t, e, mu, lower):
        interval = safe_energy_interval(mui, li, m=m, n=n)
        if ei > interval.maximum_energy:
            return float(ti)
    return None


def no_compression_bound_counterexample_energy(
    epsilon: float,
    mean_stretch: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Energy of a crack-free two-point measure showing why a third condition is required.

    The measure is supported at epsilon and lambda_c and has the requested
    mean.  As epsilon -> 0+, the generalized-LJ repulsive energy diverges,
    while the support remains entirely at or below lambda_c.
    """
    lam_c = critical_stretch(m, n)
    eps = float(epsilon)
    mu = float(mean_stretch)
    if not (0.0 < eps < mu < lam_c):
        raise ValueError("require 0 < epsilon < mean_stretch < lambda_c")

    weight_critical = (mu - eps) / (lam_c - eps)
    weight_epsilon = 1.0 - weight_critical
    return (
        weight_epsilon * float(shifted_lj_energy(eps, m, n))
        + weight_critical * float(shifted_lj_energy(lam_c, m, n))
    )
