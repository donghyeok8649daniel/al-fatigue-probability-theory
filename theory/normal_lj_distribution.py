# === 한국어 파일 안내 시작 ===
# - 파일 역할: fixed-length/fixed-energy ensemble에서 유도한 과거 large-M one-point spacing closure와 moment 계산을 구현한다.
# - 주요 클래스: ClosureMoments, ClosureSolution
# - 주요 함수/메서드: shifted_lj_energy, _quadrature, _legendre_rule, _logsumexp, _closure_mode, _moment_rule
#   closure_moments, solve_alpha_for_mean, solve_distribution_closure, closure_density
#   energy_derivative_at_fixed_mean
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Large-system distribution closure for the active 1D normal layer-LJ theory.

The closure is derived from an equiprobable fixed-length/fixed-configurational-
energy ensemble for M positive layer spacings. For finite M the one-spacing
marginal is proportional to the density of states of the remaining M-1
spacings. A large-M saddle-point expansion gives

    p(lambda) = Z^{-1} exp[-alpha*lambda - beta*psi(lambda)].

This is a CONTROLLED APPROXIMATION to the driven fatigue state, not an exact
consequence of deterministic cyclic dynamics. alpha and beta are determined
by the prescribed mean stretch and mean configurational energy; they are not
fitted fatigue parameters.

Numerics:
A transformed Gauss-Legendre rule is used on (0, infinity) for ordinary
states. Very sharply concentrated states can be under-resolved by a fixed
global rule. For large alpha and beta, the implementation therefore switches
to a mode-centered finite Gauss-Legendre rule whose width is derived from the
local LJ curvature. This is a numerical resolution strategy, not a new
physical closure.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from theory.normal_lj_chain import (
    critical_stretch,
    normalized_lj_energy,
    normalized_lj_force,
    normalized_lj_stiffness,
)


@dataclass(frozen=True)
class ClosureMoments:
    alpha: float
    beta: float
    log_partition: float
    mean_stretch: float
    mean_energy: float
    variance_stretch: float
    variance_energy: float
    covariance_stretch_energy: float
    third_central_moment_stretch: float
    critical_tail_probability: float


@dataclass(frozen=True)
class ClosureSolution:
    target_mean_stretch: float
    target_mean_energy: float
    moments: ClosureMoments
    quadrature_order: int


def shifted_lj_energy(stretch, m: float = 12.19, n: float = 6.0):
    """Dimensionless LJ energy shifted to zero at the equilibrium stretch."""
    lam = np.asarray(stretch, dtype=float)
    if np.any(lam <= 0.0):
        raise ValueError("stretch must be positive")
    return normalized_lj_energy(lam, m, n) - float(normalized_lj_energy(1.0, m, n))


@lru_cache(maxsize=32)
def _quadrature(order: int, lower: float = 0.0):
    """Gauss-Legendre quadrature on lambda in (lower, infinity)."""
    if order < 64:
        raise ValueError("quadrature order must be at least 64")
    if lower < 0.0:
        raise ValueError("lower must be non-negative")
    z, w = np.polynomial.legendre.leggauss(order)
    x = 0.5 * (z + 1.0)
    wx = 0.5 * w
    lam = float(lower) + x / (1.0 - x)
    jacobian = 1.0 / (1.0 - x) ** 2
    log_measure = np.log(wx) + np.log(jacobian)
    return lam, log_measure


@lru_cache(maxsize=16)
def _legendre_rule(order: int):
    if order < 64:
        raise ValueError("quadrature order must be at least 64")
    return np.polynomial.legendre.leggauss(order)


def _logsumexp(values: np.ndarray) -> float:
    vmax = float(np.max(values))
    return vmax + math.log(float(np.sum(np.exp(values - vmax))))


def _closure_mode(
    alpha: float,
    beta: float,
    *,
    m: float,
    n: float,
) -> float:
    """Mode of exp[-alpha lambda-beta psi(lambda)] on the compressive branch."""
    target_force = -alpha / beta
    lo = 1.0e-8
    hi = 1.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        force_mid = float(normalized_lj_force(mid, m, n))
        if force_mid < target_force:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _moment_rule(
    alpha: float,
    beta: float,
    *,
    m: float,
    n: float,
    quadrature_order: int,
):
    """Return quadrature nodes and log-measure for closure moments.

    A mode-centered rule is used only for sharply concentrated states for
    which a fixed global transformed rule can miss the narrow probability
    peak. The switch thresholds are numerical, not material parameters.
    """
    if beta > 1.0e4 and alpha > 50.0:
        mode = _closure_mode(alpha, beta, m=m, n=n)
        stiffness = float(normalized_lj_stiffness(mode, m, n))
        curvature = beta * stiffness
        if curvature > 0.0:
            sigma_local = 1.0 / math.sqrt(curvature)
            lower = max(1.0e-10, mode - 24.0 * sigma_local)
            upper = mode + 24.0 * sigma_local
            z, w = _legendre_rule(int(quadrature_order))
            lam = 0.5 * (upper - lower) * z + 0.5 * (upper + lower)
            weights = 0.5 * (upper - lower) * w
            return lam, np.log(weights)

    return _quadrature(int(quadrature_order), 0.0)


def closure_moments(
    alpha: float,
    beta: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
    quadrature_order: int = 640,
) -> ClosureMoments:
    """Evaluate normalized moments of exp[-alpha lambda-beta psi(lambda)]."""
    alpha = float(alpha)
    beta = float(beta)
    if alpha <= 0.0:
        raise ValueError("alpha must be positive for tensile-tail integrability")
    if beta <= 0.0:
        raise ValueError("beta must be positive for repulsive-tail integrability")

    lam, log_measure = _moment_rule(
        alpha,
        beta,
        m=m,
        n=n,
        quadrature_order=int(quadrature_order),
    )
    energy = shifted_lj_energy(lam, m, n)
    log_terms = log_measure - alpha * lam - beta * energy
    log_z = _logsumexp(log_terms)
    probabilities = np.exp(log_terms - log_z)

    mean_lam = float(np.sum(probabilities * lam))
    mean_e = float(np.sum(probabilities * energy))
    dl = lam - mean_lam
    de = energy - mean_e
    var_lam = float(np.sum(probabilities * dl * dl))
    var_e = float(np.sum(probabilities * de * de))
    cov = float(np.sum(probabilities * dl * de))
    third_lam = float(np.sum(probabilities * dl ** 3))

    lam_c = critical_stretch(m, n)
    tail_lam, tail_log_measure = _quadrature(int(quadrature_order), float(lam_c))
    tail_energy = shifted_lj_energy(tail_lam, m, n)
    tail_log_terms = tail_log_measure - alpha * tail_lam - beta * tail_energy
    tail = math.exp(_logsumexp(tail_log_terms) - log_z)

    return ClosureMoments(
        alpha=alpha,
        beta=beta,
        log_partition=log_z,
        mean_stretch=mean_lam,
        mean_energy=mean_e,
        variance_stretch=var_lam,
        variance_energy=var_e,
        covariance_stretch_energy=cov,
        third_central_moment_stretch=third_lam,
        critical_tail_probability=tail,
    )


def solve_alpha_for_mean(
    beta: float,
    target_mean_stretch: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
    quadrature_order: int = 640,
    relative_tolerance: float = 1.0e-11,
) -> float:
    """Solve the length multiplier alpha at fixed beta and mean stretch."""
    beta = float(beta)
    mu = float(target_mean_stretch)
    if beta <= 0.0 or mu <= 0.0:
        raise ValueError("beta and target_mean_stretch must be positive")

    lo = 1.0e-10
    hi = 1.0
    while closure_moments(
        hi,
        beta,
        m=m,
        n=n,
        quadrature_order=quadrature_order,
    ).mean_stretch > mu:
        hi *= 2.0
        if hi > 1.0e12:
            raise RuntimeError("failed to bracket alpha")

    for _ in range(120):
        mid = 0.5 * (lo + hi)
        mean_mid = closure_moments(
            mid,
            beta,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        ).mean_stretch
        if mean_mid > mu:
            lo = mid
        else:
            hi = mid
        if abs(hi - lo) <= relative_tolerance * max(1.0, mid):
            break
    return 0.5 * (lo + hi)


def solve_distribution_closure(
    target_mean_stretch: float,
    target_mean_energy: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
    quadrature_order: int = 640,
    relative_tolerance: float = 1.0e-9,
) -> ClosureSolution:
    """Solve alpha and beta from the mean-stretch and energy constraints."""
    mu = float(target_mean_stretch)
    target_e = float(target_mean_energy)
    if mu <= 0.0:
        raise ValueError("target_mean_stretch must be positive")

    minimum_e = float(shifted_lj_energy(mu, m, n))
    if target_e <= minimum_e:
        raise ValueError(
            "target_mean_energy must exceed the homogeneous Jensen minimum; "
            "equality corresponds to the singular delta distribution"
        )

    beta_lo = 1.0e-3
    alpha_lo = solve_alpha_for_mean(
        beta_lo,
        mu,
        m=m,
        n=n,
        quadrature_order=quadrature_order,
    )
    e_lo = closure_moments(
        alpha_lo,
        beta_lo,
        m=m,
        n=n,
        quadrature_order=quadrature_order,
    ).mean_energy
    while e_lo < target_e:
        beta_lo *= 0.1
        if beta_lo < 1.0e-14:
            raise RuntimeError("failed to bracket target energy from above")
        alpha_lo = solve_alpha_for_mean(
            beta_lo,
            mu,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        )
        e_lo = closure_moments(
            alpha_lo,
            beta_lo,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        ).mean_energy

    beta_hi = 1.0
    while True:
        alpha_hi = solve_alpha_for_mean(
            beta_hi,
            mu,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        )
        e_hi = closure_moments(
            alpha_hi,
            beta_hi,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        ).mean_energy
        if e_hi <= target_e:
            break
        beta_hi *= 2.0
        if beta_hi > 1.0e14:
            raise RuntimeError("failed to bracket target energy from below")

    for _ in range(120):
        beta_mid = math.sqrt(beta_lo * beta_hi)
        alpha_mid = solve_alpha_for_mean(
            beta_mid,
            mu,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        )
        e_mid = closure_moments(
            alpha_mid,
            beta_mid,
            m=m,
            n=n,
            quadrature_order=quadrature_order,
        ).mean_energy
        if e_mid > target_e:
            beta_lo = beta_mid
        else:
            beta_hi = beta_mid
        if abs(e_mid - target_e) <= relative_tolerance * max(target_e, 1.0e-14):
            break

    beta = math.sqrt(beta_lo * beta_hi)
    alpha = solve_alpha_for_mean(
        beta,
        mu,
        m=m,
        n=n,
        quadrature_order=quadrature_order,
    )
    moments = closure_moments(
        alpha,
        beta,
        m=m,
        n=n,
        quadrature_order=quadrature_order,
    )
    return ClosureSolution(
        target_mean_stretch=mu,
        target_mean_energy=target_e,
        moments=moments,
        quadrature_order=int(quadrature_order),
    )


def closure_density(
    stretch,
    alpha: float,
    beta: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
    quadrature_order: int = 640,
):
    """Evaluate the normalized p_lambda(lambda) density on positive stretches."""
    lam = np.asarray(stretch, dtype=float)
    if np.any(lam <= 0.0):
        raise ValueError("stretch must be positive")
    moments = closure_moments(
        alpha,
        beta,
        m=m,
        n=n,
        quadrature_order=quadrature_order,
    )
    log_density = (
        -alpha * lam
        - beta * shifted_lj_energy(lam, m, n)
        - moments.log_partition
    )
    return np.exp(log_density)


def energy_derivative_at_fixed_mean(moments: ClosureMoments) -> float:
    """Exact exponential-family identity d<E>/d beta at fixed mean."""
    if moments.variance_stretch <= 0.0:
        return 0.0
    residual = (
        moments.variance_energy
        - moments.covariance_stretch_energy ** 2 / moments.variance_stretch
    )
    return -max(0.0, residual)
