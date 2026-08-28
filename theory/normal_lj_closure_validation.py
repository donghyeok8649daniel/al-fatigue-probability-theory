"""Direct validation utilities for the active 1D layer-LJ distribution closure.

This module compares a deterministic spacing snapshot from the reduced 1D
layer-LJ mechanics against the large-M closure

    p(lambda) = Z^{-1} exp[-alpha lambda - beta psi(lambda)]

at exactly the same empirical mean stretch and mean configurational energy.

The comparison is a NUMERICAL FALSIFICATION DIAGNOSTIC. It does not promote
the saddle-point closure to an exact driven-state law.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from theory.normal_lj_chain import critical_stretch
from theory.normal_lj_distribution import (
    ClosureSolution,
    closure_density,
    shifted_lj_energy,
    solve_distribution_closure,
)


@dataclass(frozen=True)
class ClosureSnapshotComparison:
    represented_spacings: int
    empirical_mean_stretch: float
    empirical_mean_energy: float
    empirical_variance: float
    closure_variance: float
    variance_relative_error: float
    empirical_third_central_moment: float
    closure_third_central_moment: float
    empirical_skewness: float
    closure_skewness: float
    empirical_critical_tail_probability: float
    closure_critical_tail_probability: float
    kolmogorov_distance: float
    alpha: float
    beta: float


def _legendre_interval(order: int, lower: float, upper: float):
    if order < 32:
        raise ValueError("quadrature order must be at least 32")
    if not (0.0 <= lower < upper):
        raise ValueError("require 0 <= lower < upper")
    z, w = np.polynomial.legendre.leggauss(order)
    x = 0.5 * (upper - lower) * z + 0.5 * (upper + lower)
    weights = 0.5 * (upper - lower) * w
    return x, weights


def _legendre_unbounded(order: int):
    if order < 64:
        raise ValueError("quadrature order must be at least 64")
    z, w = np.polynomial.legendre.leggauss(order)
    q = 0.5 * (z + 1.0)
    wq = 0.5 * w
    lam = q / (1.0 - q)
    weights = wq / (1.0 - q) ** 2
    return lam, weights


def closure_third_central_moment(solution: ClosureSolution, *, quadrature_order: int = 640) -> float:
    lam, weights = _legendre_unbounded(int(quadrature_order))
    density = closure_density(
        lam,
        solution.moments.alpha,
        solution.moments.beta,
        quadrature_order=solution.quadrature_order,
    )
    mu = solution.moments.mean_stretch
    return float(np.sum(weights * density * (lam - mu) ** 3))


def closure_cdf(stretch: float, solution: ClosureSolution, *, quadrature_order: int = 128) -> float:
    x = float(stretch)
    if x <= 0.0:
        return 0.0
    lam, weights = _legendre_interval(int(quadrature_order), 0.0, x)
    density = closure_density(
        lam,
        solution.moments.alpha,
        solution.moments.beta,
        quadrature_order=solution.quadrature_order,
    )
    return float(np.sum(weights * density))


def kolmogorov_distance(spacings, solution: ClosureSolution, *, cdf_quadrature_order: int = 128) -> float:
    values = np.sort(np.asarray(spacings, dtype=float))
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spacings must be a one-dimensional sample of length >= 2")
    if np.any(values <= 0.0):
        raise ValueError("all spacings must be positive")
    model_cdf = np.asarray([
        closure_cdf(value, solution, quadrature_order=cdf_quadrature_order)
        for value in values
    ])
    count = len(values)
    empirical_upper = np.arange(1, count + 1, dtype=float) / count
    empirical_lower = np.arange(0, count, dtype=float) / count
    return float(max(
        np.max(np.abs(empirical_upper - model_cdf)),
        np.max(np.abs(empirical_lower - model_cdf)),
    ))


def compare_snapshot_to_closure(
    spacings,
    *,
    closure_quadrature_order: int = 640,
    cdf_quadrature_order: int = 128,
) -> ClosureSnapshotComparison:
    values = np.asarray(spacings, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spacings must be a one-dimensional sample of length >= 2")
    if np.any(values <= 0.0):
        raise ValueError("all spacings must be positive")

    mean = float(np.mean(values))
    mean_energy = float(np.mean(shifted_lj_energy(values)))
    centered = values - mean
    empirical_variance = float(np.mean(centered ** 2))
    empirical_third = float(np.mean(centered ** 3))

    solution = solve_distribution_closure(
        mean,
        mean_energy,
        quadrature_order=int(closure_quadrature_order),
    )
    closure_variance = solution.moments.variance_stretch
    closure_third = closure_third_central_moment(
        solution,
        quadrature_order=int(closure_quadrature_order),
    )
    empirical_skew = empirical_third / empirical_variance ** 1.5 if empirical_variance > 0.0 else 0.0
    closure_skew = closure_third / closure_variance ** 1.5 if closure_variance > 0.0 else 0.0
    lam_c = critical_stretch()

    return ClosureSnapshotComparison(
        represented_spacings=len(values),
        empirical_mean_stretch=mean,
        empirical_mean_energy=mean_energy,
        empirical_variance=empirical_variance,
        closure_variance=closure_variance,
        variance_relative_error=abs(closure_variance - empirical_variance) / max(empirical_variance, 1.0e-300),
        empirical_third_central_moment=empirical_third,
        closure_third_central_moment=closure_third,
        empirical_skewness=empirical_skew,
        closure_skewness=closure_skew,
        empirical_critical_tail_probability=float(np.mean(values >= lam_c)),
        closure_critical_tail_probability=solution.moments.critical_tail_probability,
        kolmogorov_distance=kolmogorov_distance(
            values,
            solution,
            cdf_quadrature_order=int(cdf_quadrature_order),
        ),
        alpha=solution.moments.alpha,
        beta=solution.moments.beta,
    )
