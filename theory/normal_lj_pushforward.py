# === 한국어 파일 안내 시작 ===
# - 파일 역할: deterministic spacing field에서 one-point density를 얻는 push-forward 관련 계산 도구를 모은다. harmonic/Taylor 기반 항목은 active 전역 분포 가정으로 사용하지 않는다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: lj_force_taylor_coefficients, linear_mode_frequency, arcsine_density, arcsine_cdf
#   single_mode_moments, two_harmonic_moments, two_harmonic_max_abs_skewness
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Push-forward structure implied by the 1D layer-LJ governing equations.

The empirical one-point spacing density is the push-forward of the uniformly
weighted layer label under the deterministic spacing field.  In a continuum
label xi in [0,1],

    p(lambda,t) = int_0^1 delta(lambda - Lambda(xi,t)) dxi.

For a piecewise monotone spacing profile this becomes a sum of inverse spatial
slopes over all preimages of lambda.  Linearizing the calibrated layer-LJ chain
about lambda=1 gives the discrete wave equation, so a single normal mode has an
arcsine spacing density.  Nonlinear LJ terms generate harmonic distortion and
therefore deform that arcsine law.

These statements provide a governing-equation clue for the FORM of p(lambda,t)
without selecting a named statistical family by fit.
"""
from __future__ import annotations

import math
import numpy as np


def lj_force_taylor_coefficients(m: float = 12.19, n: float = 6.0) -> tuple[float, float, float]:
    """Return coefficients in phi'(1+u)=c1*u+c2*u^2+c3*u^3+O(u^4).

    For the normalized generalized LJ model used by the repository,
    c1 = phi''(1) = 1,
    c2 = phi'''(1)/2,
    c3 = phi''''(1)/6.
    """
    if not (m > n > 1.0):
        raise ValueError("require m > n > 1")
    phi3 = -(m + n + 3.0)
    phi4 = (
        (m + 1.0) * (m + 2.0) * (m + 3.0)
        - (n + 1.0) * (n + 2.0) * (n + 3.0)
    ) / (m - n)
    return 1.0, 0.5 * phi3, phi4 / 6.0


def linear_mode_frequency(wavenumber: float) -> float:
    """Exact dispersion relation of the linearized unit-stiffness spacing chain."""
    q = float(wavenumber)
    return 2.0 * abs(math.sin(0.5 * q))


def arcsine_density(stretch, mean: float, amplitude: float):
    """Density induced by Lambda=mean+amplitude*cos(theta), theta uniform.

    This is not a fitted probability ansatz.  It is the push-forward density of
    a single linear normal mode sampled uniformly in spatial phase.
    """
    a = abs(float(amplitude))
    if a <= 0.0:
        raise ValueError("amplitude must be nonzero")
    lam = np.asarray(stretch, dtype=float)
    z = lam - float(mean)
    out = np.zeros_like(lam, dtype=float)
    inside = np.abs(z) < a
    out[inside] = 1.0 / (math.pi * np.sqrt(a * a - z[inside] * z[inside]))
    return out


def arcsine_cdf(stretch, mean: float, amplitude: float):
    """CDF corresponding to :func:`arcsine_density`."""
    a = abs(float(amplitude))
    if a <= 0.0:
        raise ValueError("amplitude must be nonzero")
    x = np.asarray(stretch, dtype=float)
    z = (x - float(mean)) / a
    out = np.empty_like(z, dtype=float)
    out[z <= -1.0] = 0.0
    out[z >= 1.0] = 1.0
    mid = (z > -1.0) & (z < 1.0)
    out[mid] = 1.0 - np.arccos(z[mid]) / math.pi
    return out


def single_mode_moments(amplitude: float) -> dict[str, float]:
    """Exact central moments of the single-mode push-forward density."""
    a = abs(float(amplitude))
    return {
        "variance": 0.5 * a * a,
        "third_central_moment": 0.0,
        "fourth_central_moment": 3.0 * a ** 4 / 8.0,
        "skewness": 0.0,
        "kurtosis": 1.5,
    }


def two_harmonic_moments(amplitude_1: float, amplitude_2: float) -> dict[str, float]:
    """Exact low moments of Lambda-mu=A*cos(theta)+B*cos(2 theta)."""
    A = float(amplitude_1)
    B = float(amplitude_2)
    variance = 0.5 * (A * A + B * B)
    third = 0.75 * A * A * B
    skew = third / variance ** 1.5 if variance > 0.0 else 0.0
    return {
        "variance": variance,
        "third_central_moment": third,
        "skewness": skew,
    }


def two_harmonic_max_abs_skewness() -> float:
    """Maximum possible |skewness| for A*cos(theta)+B*cos(2 theta)."""
    return math.sqrt(2.0 / 3.0)
