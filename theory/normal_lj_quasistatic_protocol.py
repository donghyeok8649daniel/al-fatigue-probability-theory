# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 이론 계산에 사용하는 Python 모듈이다.
# - 주요 클래스: ProtocolResidualMetrics
# - 주요 함수/메서드: stable_stretch_for_tensile_force, quasistatic_open_chain_spacings, cycle_boundary_force
#   residual_snapshot_metrics
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Quasistatic force-control and protocol diagnostics for the active 1D layer-LJ chain.

This module separates an exact static statement from a dynamical diagnostic.
For an open homogeneous chain under a constant tensile end force f, the
potential in spacing coordinates is

    Pi = sum_i [phi(lambda_i) - f lambda_i].

On the stable branch, phi'' > 0 and phi' is strictly increasing. Therefore
stationarity requires the same unique stable spacing lambda_s(f) in every
represented layer interval. The exact zero-temperature quasistatic empirical
spacing distribution is consequently a delta distribution.

Nonzero spatial variance observed in the conservative cyclic chain is thus a
dynamical/initial-condition effect unless an additional physical ensemble,
thermal fluctuation, disorder, or other source of heterogeneity is supplied.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from theory.normal_lj_chain import (
    NormalLJParameters,
    critical_dimensionless_force,
    critical_stretch,
    normalized_lj_force,
)
from theory.normal_lj_spatial_correlation import normalized_spatial_correlation
from theory.normal_lj_statistical_cells import (
    positive_window_empirical_correlation_factor,
    positive_window_empirical_effective_count,
)


@dataclass(frozen=True)
class ProtocolResidualMetrics:
    represented_spacings: int
    mean_stretch: float
    quasistatic_stretch: float
    mean_offset_from_quasistatic: float
    variance_c0: float
    rms_nonuniformity: float
    rho1: float
    tau_positive_window: float
    m_eff_positive_window: float


def stable_stretch_for_tensile_force(
    force: float,
    m: float = 12.19,
    n: float = 6.0,
    *,
    tolerance: float = 1.0e-13,
    max_iterations: int = 200,
) -> float:
    """Return the unique stable tensile root phi'(lambda)=f.

    Classification: EXACT STATIC CONSTITUTIVE ROOT within the calibrated
    homogeneous 1D layer model.

    The supported interval is 0 <= f <= f_c.  For f<f_c the stable root lies
    in 1 <= lambda < lambda_c, where phi''>0.  At f=f_c the stable and barrier
    roots coalesce at lambda_c.
    """
    f = float(force)
    if tolerance <= 0.0 or max_iterations <= 0:
        raise ValueError("tolerance and max_iterations must be positive")
    f_c = critical_dimensionless_force(m, n)
    if f < 0.0 or f > f_c * (1.0 + 1.0e-12):
        raise ValueError("force must satisfy 0 <= force <= critical force")
    if abs(f) <= tolerance:
        return 1.0
    lam_c = critical_stretch(m, n)
    if abs(f - f_c) <= tolerance:
        return lam_c

    lo = 1.0
    hi = lam_c
    for _ in range(max_iterations):
        mid = 0.5 * (lo + hi)
        value = float(normalized_lj_force(mid, m, n))
        if abs(value - f) <= tolerance:
            return mid
        if value < f:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def quasistatic_open_chain_spacings(
    force: float,
    represented_spacings: int,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Return the exact homogeneous stable force-controlled spacing state."""
    count = int(represented_spacings)
    if count < 1:
        raise ValueError("represented_spacings must be positive")
    lam = stable_stretch_for_tensile_force(force, m, n)
    return np.full(count, lam, dtype=float)


def cycle_boundary_force(parameters: NormalLJParameters, cycle_index: int) -> float:
    """Return the prescribed force at an exact integer cycle boundary.

    This mirrors the loading definition in ``simulate_normal_lj_chain``.
    For the current zero-mean sine protocol, the sinusoidal contribution is
    exactly zero at every integer cycle.  The envelope still multiplies a
    nonzero mean force when one is specified.
    """
    cycle = int(cycle_index)
    if cycle < 0:
        raise ValueError("cycle_index must be non-negative")
    if parameters.omega <= 0.0:
        raise ValueError("omega must be positive")
    period = 2.0 * math.pi / parameters.omega
    t = cycle * period
    if parameters.ramp_cycles <= 0:
        envelope = 1.0
    else:
        ramp_time = parameters.ramp_cycles * period
        if t >= ramp_time:
            envelope = 1.0
        else:
            envelope = 0.5 * (1.0 - math.cos(math.pi * t / ramp_time))
    # Use the exact integer-cycle identity sin(2*pi*cycle)=0 rather than a
    # floating evaluation of sin for this diagnostic.
    return envelope * parameters.mean_force


def residual_snapshot_metrics(
    spacings,
    *,
    quasistatic_stretch: float,
) -> ProtocolResidualMetrics:
    """Measure departure of one deterministic snapshot from its static state.

    ``rho1`` and the positive-window effective count are normalized-shape
    diagnostics.  They can remain finite even when ``variance_c0`` tends to
    zero, so they must not be interpreted without the fluctuation amplitude.
    """
    values = np.asarray(spacings, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spacings must be one-dimensional with length >= 2")
    if np.any(values <= 0.0):
        raise ValueError("all spacings must be positive")
    mean = float(np.mean(values))
    centered = values - mean
    variance = float(np.mean(centered * centered))
    if variance == 0.0:
        rho1 = 0.0
        tau = 1.0
        m_eff = float(len(values))
    else:
        rho = np.asarray(
            [normalized_spatial_correlation(values, k) for k in range(len(values))],
            dtype=float,
        )
        tau = float(positive_window_empirical_correlation_factor(rho))
        m_eff = float(positive_window_empirical_effective_count(rho))
        rho1 = float(rho[1])
    return ProtocolResidualMetrics(
        represented_spacings=len(values),
        mean_stretch=mean,
        quasistatic_stretch=float(quasistatic_stretch),
        mean_offset_from_quasistatic=abs(mean - float(quasistatic_stretch)),
        variance_c0=variance,
        rms_nonuniformity=math.sqrt(variance),
        rho1=rho1,
        tau_positive_window=tau,
        m_eff_positive_window=m_eff,
    )
