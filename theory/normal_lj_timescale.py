# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D reduced layer model의 무차원/물리 시간·주파수 변환과 scale separation 진단을 계산한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: lowest_fixed_free_mode_omega_star, lowest_fixed_free_mode_frequency_hz
#   moving_atoms_for_target_frequency, homogeneous_stretch_for_dimensionless_stress
#   local_small_oscillation_frequency_hz, normalized_lj_third_derivative
#   near_critical_distance_for_target_local_frequency
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Time-scale diagnostics for the active normal generalized-LJ mainline.

No damping, fitted relaxation time, or fatigue law is introduced here.

The module asks two falsification questions:
1. Can the lowest conservative fixed-free normal-chain mode naturally lie near
   a target fatigue frequency such as 20 Hz?
2. Can local LJ critical softening reduce an atomic normal vibration to that
   frequency at ordinary stress?

The normal chain is linearized about a stated homogeneous stretch. This is a
CONTROLLED APPROXIMATION to the nonlinear LJ dynamics, but the modal formulas
are exact for that linearized chain.
"""
from __future__ import annotations

import math

from theory.normal_lj_chain import (
    normalized_lj_force,
    normalized_lj_stiffness,
    critical_stretch,
)


def lowest_fixed_free_mode_omega_star(moving_atoms: float, stiffness: float = 1.0) -> float:
    """Lowest dimensionless angular frequency of a fixed-free harmonic chain."""
    if moving_atoms <= 0.0:
        raise ValueError("moving_atoms must be positive")
    if stiffness <= 0.0:
        raise ValueError("stiffness must be positive")
    q1 = math.pi / (2.0 * moving_atoms + 1.0)
    return 2.0 * math.sqrt(stiffness) * math.sin(0.5 * q1)


def lowest_fixed_free_mode_frequency_hz(
    moving_atoms: float,
    atomic_time_scale_s: float,
    stiffness: float = 1.0,
) -> float:
    if atomic_time_scale_s <= 0.0:
        raise ValueError("atomic_time_scale_s must be positive")
    return lowest_fixed_free_mode_omega_star(moving_atoms, stiffness) / (
        2.0 * math.pi * atomic_time_scale_s
    )


def moving_atoms_for_target_frequency(
    target_frequency_hz: float,
    atomic_time_scale_s: float,
    stiffness: float = 1.0,
) -> float:
    """Invert the lowest-mode relation for an effectively continuous atom count."""
    if target_frequency_hz <= 0.0:
        raise ValueError("target_frequency_hz must be positive")
    if atomic_time_scale_s <= 0.0:
        raise ValueError("atomic_time_scale_s must be positive")
    if stiffness <= 0.0:
        raise ValueError("stiffness must be positive")

    omega_star = 2.0 * math.pi * target_frequency_hz * atomic_time_scale_s
    argument = omega_star / (2.0 * math.sqrt(stiffness))
    if not (0.0 < argument < 1.0):
        raise ValueError("target frequency lies outside the linear-chain band")

    q1 = 2.0 * math.asin(argument)
    return 0.5 * (math.pi / q1 - 1.0)


def homogeneous_stretch_for_dimensionless_stress(
    stress_over_E: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Stable homogeneous tensile root phi'(lambda)=stress/E below ideal strength."""
    if stress_over_E < 0.0:
        raise ValueError("stress_over_E must be non-negative")

    lam_c = critical_stretch(m, n)
    f_c = float(normalized_lj_force(lam_c, m, n))
    if stress_over_E > f_c:
        raise ValueError("stress exceeds the stable homogeneous branch")

    lo, hi = 1.0, lam_c
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if float(normalized_lj_force(mid, m, n)) < stress_over_E:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def local_small_oscillation_frequency_hz(
    stretch: float,
    atomic_time_scale_s: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Local harmonic frequency from the LJ tangent stiffness."""
    if atomic_time_scale_s <= 0.0:
        raise ValueError("atomic_time_scale_s must be positive")
    k_star = float(normalized_lj_stiffness(stretch, m, n))
    if k_star <= 0.0:
        return 0.0
    return math.sqrt(k_star) / (2.0 * math.pi * atomic_time_scale_s)


def normalized_lj_third_derivative(
    stretch: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Third derivative of the normalized generalized-LJ potential."""
    return (
        -(m + 1.0) * (m + 2.0) * stretch ** (-m - 3.0)
        + (n + 1.0) * (n + 2.0) * stretch ** (-n - 3.0)
    ) / (m - n)


def near_critical_distance_for_target_local_frequency(
    target_frequency_hz: float,
    atomic_time_scale_s: float,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Linearized estimate of lambda_c-lambda needed for a target local frequency."""
    if target_frequency_hz <= 0.0 or atomic_time_scale_s <= 0.0:
        raise ValueError("frequency and time scale must be positive")
    lam_c = critical_stretch(m, n)
    phi3_c = normalized_lj_third_derivative(lam_c, m=m, n=n)
    required_k = (2.0 * math.pi * target_frequency_hz * atomic_time_scale_s) ** 2
    return required_k / abs(phi3_c)
