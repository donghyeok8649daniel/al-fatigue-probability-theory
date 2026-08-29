# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D normal layer-LJ의 준정적 stable branch와 닫힌 subcritical cycle의 zero-hysteresis/no-accumulation 한계를 계산한다.
# - 주요 클래스: 없음
# - 주요 함수/메서드: stable_stretch_from_force, quasistatic_force_from_phase, quasistatic_cycle, closed_cycle_residual, branch_retrace_error
# - 주의: 실제 피로 누적 법칙을 추가하지 않으며, frequency-independent phase path는 준정적 homogeneous equilibrium 가정 아래의 결과다.
# === 한국어 파일 안내 끝 ===
"""Quasistatic limit of the active one-dimensional normal layer-LJ model.

This module contains no fatigue law and no fitted relaxation time.  It makes
one narrow statement explicit: below the ideal layer-normal instability, the
stable homogeneous equilibrium stretch is a single-valued function of the
instantaneous dimensionless tensile force.  Consequently a quasistatic closed
force cycle retraces the same equilibrium path and cannot by itself accumulate
fatigue damage or hysteresis.

Classification:
- exact for the stated homogeneous equilibrium branch of the normalized LJ
  potential;
- the identification of a laboratory cycle with this quasistatic branch is a
  time-scale approximation and must be checked separately.
"""
from __future__ import annotations

import math

import numpy as np

from theory.normal_lj_chain import (
    critical_dimensionless_force,
    critical_stretch,
    normalized_lj_force,
)


def stable_stretch_from_force(
    force_star: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Return the unique stable tensile root phi'(lambda)=force_star.

    The domain is the tensile stable branch 0 <= force_star <= f_c and
    1 <= lambda <= lambda_c.
    """
    f = float(force_star)
    if f < 0.0:
        raise ValueError("force_star must be non-negative on the tensile branch")
    fc = critical_dimensionless_force(m, n)
    if f > fc:
        raise ValueError("force_star exceeds the stable homogeneous branch")
    if f == 0.0:
        return 1.0
    if f == fc:
        return critical_stretch(m, n)

    lo = 1.0
    hi = critical_stretch(m, n)
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if float(normalized_lj_force(mid, m, n)) < f:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def quasistatic_force_from_phase(
    phase: np.ndarray | float,
    *,
    mean_force: float,
    force_amplitude: float,
):
    """Dimensionless tensile force parameterized by cycle phase theta."""
    theta = np.asarray(phase, dtype=float)
    return mean_force + force_amplitude * np.sin(theta)


def quasistatic_cycle(
    *,
    mean_force: float,
    force_amplitude: float,
    samples: int = 721,
    m: float = 12.19,
    n: float = 6.0,
):
    """Return (phase, force, stable stretch) for one quasistatic cycle.

    Physical frequency does not appear because phase, not laboratory time,
    parameterizes the quasistatic equilibrium path.
    """
    if samples < 3:
        raise ValueError("samples must be at least 3")
    if mean_force - abs(force_amplitude) < 0.0:
        raise ValueError("this helper is restricted to tension-only cycles")
    fc = critical_dimensionless_force(m, n)
    if mean_force + abs(force_amplitude) > fc:
        raise ValueError("cycle exceeds the stable homogeneous branch")

    phase = np.linspace(0.0, 2.0 * math.pi, samples)
    force = quasistatic_force_from_phase(
        phase,
        mean_force=mean_force,
        force_amplitude=force_amplitude,
    )
    stretch = np.asarray(
        [stable_stretch_from_force(float(value), m=m, n=n) for value in force],
        dtype=float,
    )
    return phase, force, stretch


def closed_cycle_residual(
    *,
    mean_force: float,
    force_amplitude: float,
    samples: int = 721,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Absolute stretch mismatch between start and end of a closed cycle."""
    _, _, stretch = quasistatic_cycle(
        mean_force=mean_force,
        force_amplitude=force_amplitude,
        samples=samples,
        m=m,
        n=n,
    )
    return abs(float(stretch[-1] - stretch[0]))


def branch_retrace_error(
    *,
    mean_force: float,
    force_amplitude: float,
    samples: int = 721,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Return max loading/unloading mismatch at matched force levels.

    Since lambda_s is a single-valued function of force on the stable branch,
    the exact mathematical value is zero.  The returned value measures only
    floating-point/root-solve error for the sampled numerical representation.
    """
    _, force, stretch = quasistatic_cycle(
        mean_force=mean_force,
        force_amplitude=force_amplitude,
        samples=samples,
        m=m,
        n=n,
    )
    reconstructed = np.asarray(
        [stable_stretch_from_force(float(value), m=m, n=n) for value in force],
        dtype=float,
    )
    return float(np.max(np.abs(stretch - reconstructed)))
