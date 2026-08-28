"""Large-system distribution closure for the active 1D normal layer-LJ theory.

The closure is derived from an equiprobable fixed-length/fixed-configurational-
energy ensemble for M positive layer spacings.  For finite M the one-spacing
marginal is proportional to the density of states of the remaining M-1
spacings.  A large-M saddle-point expansion gives

    p(lambda) = Z^{-1} exp[-alpha*lambda - beta*psi(lambda)].

This is a CONTROLLED APPROXIMATION to the driven fatigue state, not an exact
consequence of deterministic cyclic dynamics.  alpha and beta are determined
by the prescribed mean stretch and mean configurational energy; they are not
fitted fatigue parameters.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from theory.normal_lj_chain import critical_stretch, normalized_lj_energy


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
    """Gauss-Legendre quadrature on lambda in (lower, infinity).

    The map lambda=lower+x/(1-x), x in (0,1), resolves the unbounded domain
    without introducing a finite support cutoff.
    """
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


def _logsumexp(values: np.ndarray) -> float:
    vmax = float(np.max(values))
    return vmax + math.log(float(np.sum(np.exp(values - vmax))))


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

    lam, log_measure = _quadrature(int(quadrature_order), 0.0)
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
    while closure_moments(hi, beta, m=m, n=n, quadrature_order=quadrature_order).mean_stretch > mu:
        hi *= 2.0
        if hi > 1.0e8:
            raise RuntimeError("failed to bracket alpha")

    for _ in range(100):
        mid = 0.5 * (lo + hi)
        mean_mid = closure_moments(
            mid, beta, m=m, n=n, quadrature_order=quadrature_order
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
        beta_lo, mu, m=m, n=n, quadrature_order=quadrature_order
    )
    e_lo = closure_moments(
        alpha_lo, beta_lo, m=m, n=n, quadrature_order=quadrature_order
    ).mean_energy
    while e_lo < target_e:
        beta_lo *= 0.1
        if beta_lo < 1.0e-14:
            raise RuntimeError("failed to bracket target energy from above")
        alpha_lo = solve_alpha_for_mean(
            beta_lo, mu, m=m, n=n, quadrature_order=quadrature_order
        )
        e_lo = closure_moments(
            alpha_lo, beta_lo, m=m, n=n, quadrature_order=quadrature_order
        ).mean_energy

    beta_hi = 1.0
    while True:
        alpha_hi = solve_alpha_for_mean(
            beta_hi, mu, m=m, n=n, quadrature_order=quadrature_order
        )
        e_hi = closure_moments(
            alpha_hi, beta_hi, m=m, n=n, quadrature_order=quadrature_order
        ).mean_energy
        if e_hi <= target_e:
            break
        beta_hi *= 2.0
        if beta_hi > 1.0e12:
            raise RuntimeError("failed to bracket target energy from below")

    for _ in range(100):
        beta_mid = math.sqrt(beta_lo * beta_hi)
        alpha_mid = solve_alpha_for_mean(
            beta_mid, mu, m=m, n=n, quadrature_order=quadrature_order
        )
        e_mid = closure_moments(
            alpha_mid, beta_mid, m=m, n=n, quadrature_order=quadrature_order
        ).mean_energy
        if e_mid > target_e:
            beta_lo = beta_mid
        else:
            beta_hi = beta_mid
        if abs(e_mid - target_e) <= relative_tolerance * max(target_e, 1.0e-14):
            break

    beta = math.sqrt(beta_lo * beta_hi)
    alpha = solve_alpha_for_mean(
        beta, mu, m=m, n=n, quadrature_order=quadrature_order
    )
    moments = closure_moments(
        alpha, beta, m=m, n=n, quadrature_order=quadrature_order
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
        alpha, beta, m=m, n=n, quadrature_order=quadrature_order
    )
    log_density = -alpha * lam - beta * shifted_lj_energy(lam, m, n) - moments.log_partition
    return np.exp(log_density)


def energy_derivative_at_fixed_mean(moments: ClosureMoments) -> float:
    """Exact exponential-family identity d<E>/d beta at fixed mean.

    d<E>/d beta|_mu = -[Var(E)-Cov(lambda,E)^2/Var(lambda)] <= 0.
    """
    if moments.variance_stretch <= 0.0:
        return 0.0
    residual = (
        moments.variance_energy
        - moments.covariance_stretch_energy ** 2 / moments.variance_stretch
    )
    return -max(0.0, residual)
