# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D layer-spacing 상관으로부터 true-process 분산 기준 유효 독립개수와 finite-snapshot positive-window 통계 특성길이를 구분해 계산하고, 완전종속/독립 및 동일-block event 극한을 검증한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: finite_correlation_factor, effective_independent_count, variance_equivalent_axial_length, positive_window_empirical_correlation_factor, positive_window_empirical_effective_count, positive_window_empirical_axial_length, identical_pair_msd, independent_any_event_probability, identical_block_any_event_probability
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Correlation-based statistical-cell quantities for the active 1D theory.

Two objects are intentionally separated.

1. For a second-order stationary stochastic spacing process with *true* lag
   correlations rho_k,

       Var(mean_M) = sigma^2/M * tau_M,
       tau_M = 1 + 2 sum_{k=1}^{M-1}(1-k/M) rho_k.

   This is an exact second-moment identity.

2. A single deterministic finite-chain snapshot provides only empirical
   open-chain correlation diagnostics.  Because the snapshot is centered by
   its own sample mean, a full all-lag plug-in sum has a zero-sum finite-sample
   artifact and must not be identified with the population tau_M.  For such
   snapshots this module provides a separately labeled positive-window
   estimator that stops at the first non-positive empirical correlation.

The module remains strictly one-dimensional.  The energy-calibration area A0
is not identified with a transverse statistical independence area.
"""
from __future__ import annotations

import math

import numpy as np


def finite_correlation_factor(rho) -> float:
    """Return exact finite-M tau_M from true stationary correlations.

    EXACT / IDENTITY when ``rho`` contains the true correlation coefficients
    rho[0],...,rho[M-1] of a second-order stationary sequence.  This function
    deliberately enforces the population-correlation bounds |rho_k|<=1.

    Do not pass the open-chain finite-snapshot diagnostic rho_k from
    ``normal_lj_spatial_correlation`` here; use the separately labeled
    positive-window empirical estimator below.
    """
    values = np.asarray(rho, dtype=float)
    if values.ndim != 1 or values.size < 1:
        raise ValueError("rho must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("rho must be finite")
    if not math.isclose(float(values[0]), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("rho[0] must equal 1")
    if np.any(np.abs(values) > 1.0 + 1.0e-12):
        raise ValueError("true correlation coefficients must lie in [-1,1]")

    m = values.size
    if m == 1:
        return 1.0
    k = np.arange(1, m, dtype=float)
    weights = 1.0 - k / float(m)
    tau = 1.0 + 2.0 * float(np.dot(weights, values[1:]))
    if tau <= 0.0:
        raise ValueError(
            "non-positive correlation factor: rho is inconsistent with a positive variance of the represented mean"
        )
    return tau


def effective_independent_count(rho) -> float:
    """Return M_eff=M/tau_M for true stationary correlations.

    This is a variance-equivalent independent count, not a proof of joint
    factorization for arbitrary events.
    """
    values = np.asarray(rho, dtype=float)
    tau = finite_correlation_factor(values)
    return float(values.size) / tau


def variance_equivalent_axial_length(
    rho,
    equilibrium_spacing: float = 1.0,
) -> float:
    """Return ell_stat^(2)=a0*tau_M for true stationary correlations."""
    if equilibrium_spacing <= 0.0:
        raise ValueError("equilibrium_spacing must be positive")
    return float(equilibrium_spacing) * finite_correlation_factor(rho)


def positive_window_empirical_correlation_factor(open_chain_rho) -> float:
    """Return a first-positive-lobe estimator of the correlation factor.

    ESTIMATOR / DIAGNOSTIC for a single finite deterministic snapshot.

    The supplied rho values may be the open-chain diagnostic C_k/C_0, whose
    finite-sample values are not required to remain in [-1,1].  The estimator
    uses

        tau_hat_plus = 1 + 2 sum_{k=1}^{K0} (1-k/M) rho_hat_k,

    where K0 is the last strictly positive lag before the first non-positive
    empirical correlation.  If no non-positive lag occurs, all available
    positive lags are included.

    The positive-window rule prevents the exact zero-sum artifact that occurs
    when all lags of a sample-mean-centered finite snapshot are summed.  It is
    not promoted to an exact material constant; convergence with represented
    system size must be checked.
    """
    values = np.asarray(open_chain_rho, dtype=float)
    if values.ndim != 1 or values.size < 1:
        raise ValueError("open_chain_rho must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("open_chain_rho must be finite")
    if not math.isclose(float(values[0]), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("open_chain_rho[0] must equal 1")

    m = values.size
    total = 1.0
    for k in range(1, m):
        rho_k = float(values[k])
        if rho_k <= 0.0:
            break
        total += 2.0 * (1.0 - k / float(m)) * rho_k
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("positive-window empirical correlation factor is invalid")
    return total


def positive_window_empirical_effective_count(open_chain_rho) -> float:
    """Return M/tau_hat_plus for a finite-snapshot positive-window estimate."""
    values = np.asarray(open_chain_rho, dtype=float)
    tau = positive_window_empirical_correlation_factor(values)
    return float(values.size) / tau


def positive_window_empirical_axial_length(
    open_chain_rho,
    equilibrium_spacing: float = 1.0,
) -> float:
    """Return a0*tau_hat_plus for a finite-snapshot positive-window estimate."""
    if equilibrium_spacing <= 0.0:
        raise ValueError("equilibrium_spacing must be positive")
    return float(equilibrium_spacing) * positive_window_empirical_correlation_factor(
        open_chain_rho
    )


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
    """Return 1-(1-q)^N for N genuinely independent statistical cells."""
    q = float(single_cell_probability)
    n = int(independent_cells)
    if not (0.0 <= q <= 1.0):
        raise ValueError("single_cell_probability must lie in [0,1]")
    if n < 0:
        raise ValueError("independent_cells must be non-negative")
    if n == 0:
        return 0.0
    return -math.expm1(n * math.log1p(-q)) if q < 1.0 else 1.0


def identical_block_any_event_probability(
    single_variable_probability: float,
    total_variables: int,
    identical_block_size: int,
) -> float:
    """Exact union probability for independent blocks of identical variables.

    Assume M variables are partitioned into N=M/b mutually independent
    blocks.  All b variables within a block are exactly the same random
    variable.  If its event probability is q, then

        P(any event) = 1-(1-q)^(M/b).

    This is an exact special case and is not asserted for partial dependence.
    """
    q = float(single_variable_probability)
    m = int(total_variables)
    b = int(identical_block_size)
    if not (0.0 <= q <= 1.0):
        raise ValueError("single_variable_probability must lie in [0,1]")
    if m < 0 or b <= 0:
        raise ValueError("total_variables must be non-negative and block size positive")
    if m == 0:
        return 0.0
    if m % b != 0:
        raise ValueError("total_variables must be an integer multiple of identical_block_size")
    return independent_any_event_probability(q, m // b)
