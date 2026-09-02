# === 한국어 파일 안내 시작 ===
# - 파일 역할: reduced state q=(a,s)의 exact velocity-moment hierarchy에서
#   2차 조건부 모멘트 Theta의 진화와 평균 kinetic-energy 항을 계산한다.
# - 주요 함수: theta_material_rate, local_mean_kinetic_energy,
#   conditional_second_velocity_moment, density_shape_residual.
# - 주의: 이 모듈은 closure를 가정하지 않는다. Theta 진화에 필요한 3차
#   central-moment flux divergence는 외부에서 실제 mechanics/ensemble로 제공해야 한다.
# === 한국어 파일 안내 끝 ===
"""Exact local identities for the reduced (a,s) velocity-moment hierarchy.

For q=(a,s), v=dot(q), define

    u_i(q,t)      = E[v_i | q],
    Theta_ij      = Cov(v_i,v_j | q),
    C_ijk         = E[(v_i-u_i)(v_j-u_j)(v_k-u_k) | q].

If the reduced acceleration b(q,t) is deterministic at fixed q,t (for example
b=G^{-1}(Q-grad U0) for a constant declared metric), then the exact first two
central-moment equations include

    D_t u_i = b_i - (1/P) partial_j(P Theta_ij),

and

    D_t Theta_ij
      + (partial_k u_i) Theta_kj
      + Theta_ik (partial_k u_j)
      + (1/P) partial_k(P C_ijk) = 0.

Thus the second-moment system is not closed by U0 and Theta alone: the third
conditional central moment enters exactly.  Setting C=0 is a closure assumption
and is not done here.

The Theta-based density-shape identity

    Theta grad(log P) = A - D_t u - div(Theta)

is also an identity obtained from the same first-moment balance.  It is useful
for reconstruction and consistency checks, but it is not an independent
predictive evolution law if D_t u is itself computed from that balance.
"""

from __future__ import annotations

import numpy as np


def _vector(values: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.ndim != 1 or out.size == 0 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite nonempty vector")
    return out


def _square(values: np.ndarray, name: str, n: int | None = None) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.ndim != 2 or out.shape[0] != out.shape[1]:
        raise ValueError(f"{name} must be square")
    if n is not None and out.shape != (n, n):
        raise ValueError(f"{name} must have shape {(n, n)}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def conditional_second_velocity_moment(
    mean_velocity: np.ndarray,
    covariance: np.ndarray,
) -> np.ndarray:
    """Return E[v_i v_j|q] = u_i u_j + Theta_ij."""
    u = _vector(mean_velocity, "mean_velocity")
    theta = _square(covariance, "covariance", u.size)
    return np.outer(u, u) + theta


def theta_material_rate(
    covariance: np.ndarray,
    velocity_gradient: np.ndarray,
    third_central_moment_flux_divergence_over_density: np.ndarray,
) -> np.ndarray:
    """Return the exact material rate D_t Theta for deterministic b(q,t).

    ``velocity_gradient[i,j]`` means partial u_i / partial q_j.
    ``third_central_moment_flux_divergence_over_density[i,j]`` means

        (1/P) partial_k(P C_ijk).

    No approximation to this third-moment term is made.
    """
    theta = _square(covariance, "covariance")
    n = theta.shape[0]
    grad_u = _square(velocity_gradient, "velocity_gradient", n)
    div_c = _square(
        third_central_moment_flux_divergence_over_density,
        "third_central_moment_flux_divergence_over_density",
        n,
    )
    if not np.allclose(theta, theta.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("covariance must be symmetric")
    return -(grad_u @ theta + theta @ grad_u.T + div_c)


def local_mean_kinetic_energy(
    mean_velocity: np.ndarray,
    covariance: np.ndarray,
    metric: np.ndarray,
) -> float:
    """Return E[T|q] for T=1/2 v^T G v.

    The exact conditional mean is

        1/2 [u^T G u + tr(G Theta)].

    Therefore Theta has a direct mechanical-energy meaning; it is not merely a
    numerical uncertainty parameter.
    """
    u = _vector(mean_velocity, "mean_velocity")
    theta = _square(covariance, "covariance", u.size)
    g = _square(metric, "metric", u.size)
    if not np.allclose(theta, theta.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("covariance must be symmetric")
    if not np.allclose(g, g.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("metric must be symmetric")
    return 0.5 * float(u @ g @ u + np.trace(g @ theta))


def density_shape_residual(
    covariance: np.ndarray,
    grad_log_density: np.ndarray,
    conditional_acceleration: np.ndarray,
    material_acceleration: np.ndarray,
    covariance_divergence: np.ndarray,
) -> np.ndarray:
    """Return residual of Theta grad(log P)=A-D_tu-div(Theta)."""
    grad = _vector(grad_log_density, "grad_log_density")
    theta = _square(covariance, "covariance", grad.size)
    acc = _vector(conditional_acceleration, "conditional_acceleration")
    dtu = _vector(material_acceleration, "material_acceleration")
    div_theta = _vector(covariance_divergence, "covariance_divergence")
    if acc.shape != grad.shape or dtu.shape != grad.shape or div_theta.shape != grad.shape:
        raise ValueError("all vectors must have the same dimension")
    return theta @ grad - (acc - dtu - div_theta)
