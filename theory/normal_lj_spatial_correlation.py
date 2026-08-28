"""Spatial-correlation diagnostics for the active 1D normal layer-LJ theory.

A one-point spacing density P(lambda,t) is invariant under any permutation of
the layer-spacing labels. Spatial covariance is not. This module quantifies
that missing ordering information without introducing a fitted constitutive
law.

Definitions for a finite open chain with M represented spacings:
    mu = mean(lambda_i)
    C_k = mean_{i=1..M-k}[(lambda_i-mu)(lambda_{i+k}-mu)]
    rho_k = C_k / C_0

C_0 is exactly the empirical spacing variance.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class CorrelationSummary:
    represented_spacings: int
    variance_c0: float
    rho1: float
    first_zero_crossing_lag: float | None
    first_zero_crossing_scaled_lag: float | None
    positive_correlation_integral_over_m: float


def spatial_covariance(spacings, lag: int) -> float:
    """Return the open-chain lag-k covariance C_k."""
    values = np.asarray(spacings, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spacings must be one-dimensional with length >= 2")
    if np.any(values <= 0.0):
        raise ValueError("all spacings must be positive")
    k = int(lag)
    if not (0 <= k < len(values)):
        raise ValueError("lag must satisfy 0 <= lag < number of spacings")
    centered = values - float(np.mean(values))
    if k == 0:
        return float(np.mean(centered * centered))
    return float(np.mean(centered[:-k] * centered[k:]))


def normalized_spatial_correlation(spacings, lag: int) -> float:
    """Return rho_k=C_k/C_0; zero for a zero-variance sample."""
    c0 = spatial_covariance(spacings, 0)
    if c0 == 0.0:
        return 0.0
    return spatial_covariance(spacings, lag) / c0


def correlation_profile(spacings, max_lag: int | None = None):
    """Return arrays (lag, C_k, rho_k) for an open-chain sample."""
    values = np.asarray(spacings, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spacings must be one-dimensional with length >= 2")
    if max_lag is None:
        max_lag = len(values) - 1
    max_lag = int(max_lag)
    if not (0 <= max_lag < len(values)):
        raise ValueError("max_lag must be in [0, M-1]")
    lags = np.arange(max_lag + 1, dtype=int)
    covariance = np.asarray([spatial_covariance(values, int(k)) for k in lags])
    c0 = covariance[0]
    correlation = np.zeros_like(covariance)
    if c0 != 0.0:
        correlation = covariance / c0
    return lags, covariance, correlation


def first_zero_crossing(spacings) -> tuple[float | None, float | None]:
    """Linearly interpolate the first rho_k=0 crossing and return (k, k/M)."""
    values = np.asarray(spacings, dtype=float)
    _, _, rho = correlation_profile(values)
    for k in range(1, len(rho)):
        if rho[k] <= 0.0:
            r0 = float(rho[k - 1])
            r1 = float(rho[k])
            crossing = float(k) if r1 == r0 else float(k - 1) + r0 / (r0 - r1)
            return crossing, crossing / len(values)
    return None, None


def positive_correlation_integral_over_m(spacings) -> float:
    """Discrete positive-rho area divided by M, up to the first zero crossing."""
    values = np.asarray(spacings, dtype=float)
    _, _, rho = correlation_profile(values)
    total = 0.0
    for k in range(1, len(rho)):
        if rho[k] <= 0.0:
            break
        total += float(rho[k])
    return total / len(values)


def random_permutation_expected_rho(represented_spacings: int) -> float:
    """Exact E[rho_k] for any nonzero lag under a uniform random permutation.

    For a centered finite population d_i with sum_i d_i=0, two distinct
    elements sampled without replacement obey
        E[d_i d_j] = -C_0/(M-1).
    Hence E[rho_k] = -1/(M-1) for every nonzero lag.
    """
    m = int(represented_spacings)
    if m < 2:
        raise ValueError("represented_spacings must be at least 2")
    return -1.0 / (m - 1)


def summarize_spatial_correlation(spacings) -> CorrelationSummary:
    """Return compact correlation diagnostics for a spacing snapshot."""
    values = np.asarray(spacings, dtype=float)
    crossing_k, crossing_eta = first_zero_crossing(values)
    return CorrelationSummary(
        represented_spacings=len(values),
        variance_c0=spatial_covariance(values, 0),
        rho1=normalized_spatial_correlation(values, 1),
        first_zero_crossing_lag=crossing_k,
        first_zero_crossing_scaled_lag=crossing_eta,
        positive_correlation_integral_over_m=positive_correlation_integral_over_m(values),
    )
