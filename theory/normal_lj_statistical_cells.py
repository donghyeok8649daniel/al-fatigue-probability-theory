# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D layer-spacing 상관으로부터 분산 기준 유효 독립개수와 축방향 통계 특성길이를 계산하고, 완전종속/독립 극한을 검증한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: finite_correlation_factor, effective_independent_count, variance_equivalent_axial_length, identical_pair_msd, independent_any_event_probability
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Correlation-based statistical-cell quantities for the active 1D theory.

The module deliberately stays one-dimensional.  It does not define a
transverse statistical area or a three-dimensional characteristic volume.
The representative energy-calibration area A0 is not silently identified
with an independence area.

For a second-order stationary spacing sequence with true lag correlations
rho_k and M represented spacings,

    Var(mean_M) = sigma^2/M * tau_M,

where

    tau_M = 1 + 2 sum_{k=1}^{M-1} (1-k/M) rho_k.

This identity motivates the variance-equivalent independent count
M_eff=M/tau_M and axial statistical length ell_var=a0*tau_M.  M_eff is an
exact variance-equivalent count when true correlations are supplied; it is
not a claim of full probabilistic independence.  Supplying empirical rho_k
turns it into an estimator/diagnostic.
"""
from __future__ import annotations

import math

import numpy as np


def finite_correlation_factor(rho) -> float:
    """Return tau_M from lag correlations rho[0],...,rho[M-1].

    EXACT / IDENTITY when ``rho`` contains the true correlation coefficients
    of a second-order stationary sequence.  ``rho[0]`` must equal one.
    Empirical correlations produce an estimate rather than an exact material
    property.
    """
    values = np.asarray(rho, dtype=float)
    if values.ndim != 1 or values.size < 1:
        raise ValueError("rho must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("rho must be finite")
    if not math.isclose(float(values[0]), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("rho[0] must equal 1")
    if np.any(np.abs(values) > 1.0 + 1.0e-12):
        raise ValueError("correlation coefficients must lie in [-1,1]")

    m = values.size
    if m == 1:
        return 1.0
    k = np.arange(1, m, dtype=float)
    weights = 1.0 - k / float(m)
    tau = 1.0 + 2.0 * float(np.dot(weights, values[1:]))
    if tau <= 0.0:
        raise ValueError(
            "non-positive correlation factor: rho may be inconsistent with a valid covariance sequence"
        )
    return tau


def effective_independent_count(rho) -> float:
    """Return M_eff=M/tau_M, the variance-equivalent independent count.

    This is not a proof that the variables factorize jointly.  Full
    independence requires the appropriate joint density to factorize.
    """
    values = np.asarray(rho, dtype=float)
    tau = finite_correlation_factor(values)
    return float(values.size) / tau


def variance_equivalent_axial_length(
    rho,
    equilibrium_spacing: float = 1.0,
) -> float:
    """Return ell_var=a0*tau_M for the represented 1D sequence.

    ``equilibrium_spacing`` may be physical a0 (metres) or 1 in reduced units.
    The result is the length of one perfectly correlated block that would
    give the same variance of the sample mean as the supplied correlation
    sequence.  It is therefore a variance-equivalent statistical length, not
    automatically a strict independence length for arbitrary events.
    """
    if equilibrium_spacing <= 0.0:
        raise ValueError("equilibrium_spacing must be positive")
    return float(equilibrium_spacing) * finite_correlation_factor(rho)


def identical_pair_msd(first, second) -> float:
    """Return E_sample[(lambda_i-lambda_j)^2] for paired observations.

    In probability theory E[(X-Y)^2]=0 iff X=Y almost surely.  Therefore zero
    is the exact criterion for complete identical dependence, whereas a small
    positive value is only an empirical/numerical closeness diagnostic.
    """
    x = np.asarray(first, dtype=float)
    y = np.asarray(second, dtype=float)
    if x.shape != y.shape or x.size == 0:
        raise ValueError("first and second must be nonempty arrays of equal shape")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("paired observations must be finite")
    return float(np.mean((x - y) ** 2))


def independent_any_event_probability(
    single_cell_probability: float,
    independent_cells: int,
) -> float:
    """Return 1-(1-q)^N for N genuinely independent statistical cells.

    This formula must not be applied to partially dependent or completely
    identical cells merely because their covariance is small.
    """
    q = float(single_cell_probability)
    n = int(independent_cells)
    if not (0.0 <= q <= 1.0):
        raise ValueError("single_cell_probability must lie in [0,1]")
    if n < 0:
        raise ValueError("independent_cells must be non-negative")
    if n == 0:
        return 0.0
    return -math.expm1(n * math.log1p(-q)) if q < 1.0 else 1.0
