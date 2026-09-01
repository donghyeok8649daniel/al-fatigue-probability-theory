# === 한국어 파일 안내 시작 ===
# - 파일 역할: 원래 nonlinear 1D layer-LJ 식에서 finite-M spacing/velocity empirical measure의 exact transport 및 moment identity를 계산한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: exact_spacing_acceleration, monomial_moment, monomial_moment_rate
#   monomial_moment_second_rate
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Exact transport identities for the nonlinear 1D layer-LJ spacing state.

No Taylor expansion, harmonic ansatz, Gaussian/Weibull family, damping law, or
fatigue-damage variable is used here.

For represented spacings lambda_i and spacing velocities v_i, define the
finite-M phase-space empirical measure

    F_M(lambda, v, t) = M^{-1} sum_i delta(lambda-lambda_i) delta(v-v_i).

The corresponding one-point density is P_M = int F_M dv.  In distributional
form the empirical measures obey exact transport identities.  The dynamics is
second order, so P(lambda,t) alone is generally not a closed state: its flux
requires velocity information and its flux evolution requires neighbor-force
information.
"""
from __future__ import annotations

import numpy as np

from theory.normal_lj_chain import normalized_lj_force


def exact_spacing_acceleration(
    spacing: np.ndarray,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Return exact interior spacing accelerations from the nonlinear LJ law.

    For interior spacing indices i=1,...,M-2,

        ddot(lambda_i)
        = phi'(lambda_{i+1}) - 2 phi'(lambda_i) + phi'(lambda_{i-1}).

    The two boundary entries are returned as NaN because their exact equations
    depend on the chosen boundary loading/fixation and should not be silently
    replaced by a bulk formula.
    """
    lam = np.asarray(spacing, dtype=float)
    if lam.ndim != 1 or lam.size < 3:
        raise ValueError("spacing must be a one-dimensional array with at least 3 entries")
    if np.any(lam <= 0.0):
        raise ValueError("all layer spacings must be positive")

    force = normalized_lj_force(lam, m, n)
    acc = np.full_like(lam, np.nan, dtype=float)
    acc[1:-1] = force[2:] - 2.0 * force[1:-1] + force[:-2]
    return acc


def monomial_moment(spacing: np.ndarray, order: int) -> float:
    """Finite-M empirical raw moment <lambda^order>."""
    lam = np.asarray(spacing, dtype=float)
    if lam.ndim != 1 or lam.size == 0:
        raise ValueError("spacing must be a non-empty one-dimensional array")
    if order < 0:
        raise ValueError("order must be non-negative")
    return float(np.mean(lam ** order))


def monomial_moment_rate(
    spacing: np.ndarray,
    spacing_velocity: np.ndarray,
    order: int,
) -> float:
    """Exact first time derivative of the finite-M moment <lambda^order>.

    Identity:

        d/dt <lambda^r> = r <lambda^(r-1) v>.

    This is purely kinematic and makes no closure approximation.
    """
    lam = np.asarray(spacing, dtype=float)
    vel = np.asarray(spacing_velocity, dtype=float)
    if lam.shape != vel.shape or lam.ndim != 1 or lam.size == 0:
        raise ValueError("spacing and spacing_velocity must be equal-size 1D arrays")
    if order < 0:
        raise ValueError("order must be non-negative")
    if order == 0:
        return 0.0
    return float(order * np.mean((lam ** (order - 1)) * vel))


def monomial_moment_second_rate(
    spacing: np.ndarray,
    spacing_velocity: np.ndarray,
    spacing_acceleration: np.ndarray,
    order: int,
) -> float:
    """Exact second derivative of <lambda^order> for supplied accelerations.

    Identity:

        d2/dt2 <lambda^r>
        = r(r-1)<lambda^(r-2) v^2>
          + r<lambda^(r-1) a>.

    Boundary entries may be excluded by passing arrays that contain only the
    interior spacings/velocities/accelerations for which the acceleration law
    is defined.
    """
    lam = np.asarray(spacing, dtype=float)
    vel = np.asarray(spacing_velocity, dtype=float)
    acc = np.asarray(spacing_acceleration, dtype=float)
    if not (lam.shape == vel.shape == acc.shape) or lam.ndim != 1 or lam.size == 0:
        raise ValueError("spacing, velocity, acceleration must be equal-size 1D arrays")
    if order < 0:
        raise ValueError("order must be non-negative")
    if order == 0:
        return 0.0

    term_acc = order * np.mean((lam ** (order - 1)) * acc)
    if order == 1:
        return float(term_acc)
    term_vel = order * (order - 1) * np.mean((lam ** (order - 2)) * vel * vel)
    return float(term_vel + term_acc)
