# === 한국어 파일 안내 시작 ===
# - 파일 역할: exact phase-space moment balance에서 1D layer-spacing 분포 P(lambda,t)의 순간 함수형/로그기울기 제약을 계산한다.
# - 주요 함수/메서드: density_log_slope, reconstruct_density_from_shape_fields, lj_conditional_acceleration
# - 입력: stretch grid, conditional velocity variance Theta, conditional acceleration, mean-flow material acceleration, neighbor joint densities
# - 출력: d(log P)/d(lambda), 정규화된 P(lambda), full nonlinear LJ 조건부 acceleration
# - 주의: Taylor/harmonic/Gaussian/Weibull 가정을 사용하지 않는다. P의 smooth representation에서 Theta>0인 구간에 대한 exact moment identity다.
# === 한국어 파일 안내 끝 ===
"""Exact instantaneous shape constraint for the 1D layer-spacing density.

Starting from the exact phase-space moment equations

    d_t P + d_lambda(P u) = 0,
    d_t(P u) + d_lambda[P(u^2 + Theta)] = P a_bar,

one obtains, wherever P>0 and Theta>0,

    d_lambda log P
      = (a_bar - D_t u)/Theta - d_lambda log Theta.

Therefore

    P(lambda,t)
      = C(t)/Theta(lambda,t)
        * exp(int (a_bar-D_t u)/Theta d lambda).

No Taylor expansion of the LJ force, finite-harmonic ansatz, named probability
family, damping law, or empirical fatigue variable is introduced.  The formula
is exact at the smooth one-point moment level, but it is not closed because
Theta, D_t u, and a_bar carry phase-space and neighboring-spacing information.
"""
from __future__ import annotations

import numpy as np

from theory.normal_lj_chain import normalized_lj_force


def _as_same_1d(*arrays: np.ndarray) -> tuple[np.ndarray, ...]:
    out = tuple(np.asarray(a, dtype=float) for a in arrays)
    if not out or any(a.ndim != 1 for a in out):
        raise ValueError("all inputs must be one-dimensional arrays")
    if any(a.shape != out[0].shape for a in out[1:]):
        raise ValueError("all input arrays must have the same shape")
    if out[0].size < 2:
        raise ValueError("at least two grid points are required")
    return out


def density_log_slope(
    stretch: np.ndarray,
    conditional_velocity_variance: np.ndarray,
    conditional_acceleration: np.ndarray,
    mean_flow_material_acceleration: np.ndarray,
) -> np.ndarray:
    """Return the exact smooth-moment log-slope d(log P)/d(lambda).

    With u = E[v|lambda], Theta = Var(v|lambda), and
    a_bar = E[ddot(lambda_i)|lambda_i=lambda], the exact first two moment
    balances imply

        d_lambda log P
          = (a_bar - D_t u)/Theta - d_lambda log Theta.

    ``mean_flow_material_acceleration`` is D_t u = d_t u + u d_lambda u.
    This function makes no closure assumption for any supplied field.
    """
    lam, theta, abar, dtu = _as_same_1d(
        stretch,
        conditional_velocity_variance,
        conditional_acceleration,
        mean_flow_material_acceleration,
    )
    if np.any(np.diff(lam) <= 0.0):
        raise ValueError("stretch grid must be strictly increasing")
    if np.any(theta <= 0.0):
        raise ValueError("conditional velocity variance must be strictly positive")

    dlogtheta = np.gradient(np.log(theta), lam, edge_order=2 if lam.size >= 3 else 1)
    return (abar - dtu) / theta - dlogtheta


def reconstruct_density_from_shape_fields(
    stretch: np.ndarray,
    conditional_velocity_variance: np.ndarray,
    conditional_acceleration: np.ndarray,
    mean_flow_material_acceleration: np.ndarray,
) -> np.ndarray:
    """Reconstruct normalized P(lambda) from the exact instantaneous shape law.

    The integration constant is fixed by normalization on the supplied stretch
    grid.  The result is a continuum representation on that grid; it does not
    assert that the finite-M empirical measure itself is smooth.
    """
    lam, theta, abar, dtu = _as_same_1d(
        stretch,
        conditional_velocity_variance,
        conditional_acceleration,
        mean_flow_material_acceleration,
    )
    slope = density_log_slope(lam, theta, abar, dtu)

    # Integrate d(log P)/d lambda with the trapezoidal rule.  The arbitrary
    # additive constant in log P is removed by normalization.
    delta = np.diff(lam)
    increments = 0.5 * (slope[:-1] + slope[1:]) * delta
    logp = np.concatenate(([0.0], np.cumsum(increments)))
    logp -= float(np.max(logp))
    p = np.exp(logp)
    norm = float(np.trapezoid(p, lam))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("density normalization failed")
    return p / norm


def lj_conditional_acceleration(
    stretch: np.ndarray,
    right_neighbor_stretch: np.ndarray,
    right_neighbor_joint_density: np.ndarray,
    left_neighbor_stretch: np.ndarray,
    left_neighbor_joint_density: np.ndarray,
    one_point_density: np.ndarray,
    *,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Return exact interior conditional acceleration from neighbor joint states.

    For central spacing lambda,

        a_bar(lambda)
          = E[phi'(lambda_{i+1})|lambda_i=lambda]
            + E[phi'(lambda_{i-1})|lambda_i=lambda]
            - 2 phi'(lambda).

    ``right_neighbor_joint_density[j,k]`` represents P2+(lambda_j,lambda'_k)
    and similarly for the left neighbor.  The first axis must match ``stretch``.
    No independence or symmetry between left/right neighbors is imposed.
    """
    lam = np.asarray(stretch, dtype=float)
    rp = np.asarray(right_neighbor_stretch, dtype=float)
    lp = np.asarray(left_neighbor_stretch, dtype=float)
    p2r = np.asarray(right_neighbor_joint_density, dtype=float)
    p2l = np.asarray(left_neighbor_joint_density, dtype=float)
    p1 = np.asarray(one_point_density, dtype=float)

    if lam.ndim != 1 or rp.ndim != 1 or lp.ndim != 1 or p1.shape != lam.shape:
        raise ValueError("stretch, neighbor grids, and one-point density have incompatible shapes")
    if p2r.shape != (lam.size, rp.size) or p2l.shape != (lam.size, lp.size):
        raise ValueError("joint-density shapes must be (central_grid, neighbor_grid)")
    if np.any(p1 <= 0.0):
        raise ValueError("one-point density must be positive where conditional acceleration is evaluated")
    if np.any(np.diff(rp) <= 0.0) or np.any(np.diff(lp) <= 0.0):
        raise ValueError("neighbor stretch grids must be strictly increasing")

    force_r = normalized_lj_force(rp, m, n)
    force_l = normalized_lj_force(lp, m, n)
    mean_force_r = np.trapezoid(p2r * force_r[None, :], rp, axis=1) / p1
    mean_force_l = np.trapezoid(p2l * force_l[None, :], lp, axis=1) / p1
    return mean_force_r + mean_force_l - 2.0 * normalized_lj_force(lam, m, n)
