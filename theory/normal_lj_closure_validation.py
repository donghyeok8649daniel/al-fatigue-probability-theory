# === 한국어 파일 안내 시작 ===
# - 파일 역할: deterministic layer-spacing snapshot과 mean/energy 기반 one-point distribution closure를 같은 조건에서 비교한다.
# - 주요 클래스: ClosureSnapshotComparison
# - 주요 함수/메서드: _standard_legendre_rule, _legendre_unbounded, closure_third_central_moment
#   closure_cdf_many, closure_cdf, kolmogorov_distance, compare_snapshot_to_closure
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Direct validation utilities for the active 1D layer-LJ distribution closure.

This module compares deterministic spacing snapshots from the reduced 1D
layer-LJ mechanics against the large-M closure

    p(lambda) = Z^{-1} exp[-alpha lambda - beta psi(lambda)]

at exactly the same empirical mean stretch and mean configurational energy.

The comparison is a NUMERICAL FALSIFICATION DIAGNOSTIC. It does not promote
the saddle-point closure to an exact driven-state law.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

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


@lru_cache(maxsize=16)
def _standard_legendre_rule(order: int):
    if order < 32:
        raise ValueError("quadrature order must be at least 32")
    return np.polynomial.legendre.leggauss(order)


def _legendre_unbounded(order: int):
    if order < 64:
        raise ValueError("quadrature order must be at least 64")
    z, w = _standard_legendre_rule(order)
    q = 0.5 * (z + 1.0)
    wq = 0.5 * w
    lam = q / (1.0 - q)
    weights = wq / (1.0 - q) ** 2
    return lam, weights


def closure_third_central_moment(
    solution: ClosureSolution,
    *,
    quadrature_order: int = 640,
) -> float:
    """Third central stretch moment evaluated by the closure moment rule.

    quadrature_order is retained for API compatibility; the solved closure
    already stores the moment evaluated with its own resolved quadrature.
    """
    del quadrature_order
    return solution.moments.third_central_moment_stretch


def closure_cdf_many(
    stretches,
    solution: ClosureSolution,
    *,
    quadrature_order: int = 128,
):
    """Evaluate closure CDF values in one vectorized Gauss-Legendre pass."""
    x = np.asarray(stretches, dtype=float)
    if x.ndim != 1:
        raise ValueError("stretches must be one-dimensional")

    out = np.zeros_like(x)
    positive = x > 0.0
    if not np.any(positive):
        return out

    z, w = _standard_legendre_rule(int(quadrature_order))
    u = 0.5 * (z + 1.0)
    wu = 0.5 * w

    xp = x[positive]
    lam = xp[:, None] * u[None, :]
    weights = xp[:, None] * wu[None, :]
    density = closure_density(
        lam,
        solution.moments.alpha,
        solution.moments.beta,
        quadrature_order=solution.quadrature_order,
    )
    out[positive] = np.sum(weights * density, axis=1)
    return out


def closure_cdf(
    stretch: float,
    solution: ClosureSolution,
    *,
    quadrature_order: int = 128,
) -> float:
    return float(
        closure_cdf_many(
            np.asarray([stretch], dtype=float),
            solution,
            quadrature_order=quadrature_order,
        )[0]
    )


def kolmogorov_distance(
    spacings,
    solution: ClosureSolution,
    *,
    cdf_quadrature_order: int = 128,
) -> float:
    values = np.sort(np.asarray(spacings, dtype=float))
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spacings must be a one-dimensional sample of length >= 2")
    if np.any(values <= 0.0):
        raise ValueError("all spacings must be positive")

    model_cdf = closure_cdf_many(
        values,
        solution,
        quadrature_order=int(cdf_quadrature_order),
    )
    count = len(values)
    empirical_upper = np.arange(1, count + 1, dtype=float) / count
    empirical_lower = np.arange(0, count, dtype=float) / count

    return float(
        max(
            np.max(np.abs(empirical_upper - model_cdf)),
            np.max(np.abs(empirical_lower - model_cdf)),
        )
    )


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
    empirical_skew = (
        empirical_third / empirical_variance ** 1.5
        if empirical_variance > 0.0
        else 0.0
    )
    closure_skew = (
        closure_third / closure_variance ** 1.5
        if closure_variance > 0.0
        else 0.0
    )
    lam_c = critical_stretch()

    return ClosureSnapshotComparison(
        represented_spacings=len(values),
        empirical_mean_stretch=mean,
        empirical_mean_energy=mean_energy,
        empirical_variance=empirical_variance,
        closure_variance=closure_variance,
        variance_relative_error=abs(closure_variance - empirical_variance)
        / max(empirical_variance, 1.0e-300),
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
