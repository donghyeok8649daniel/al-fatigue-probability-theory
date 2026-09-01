# === 한국어 파일 안내 시작 ===
# - 파일 역할: active reduced state q=(a,s)의 일반화좌표 역학을 가정 명시형으로 계산한다.
# - 주요 함수: mass_metric_from_embedding_jacobians, christoffel_first_kind,
#   generalized_acceleration, conditional_acceleration, kinetic_energy,
#   mechanical_energy, external_power.
# - 주의: U0(a,s)만으로 질량 metric, 공간 patch 결합, 소산은 정해지지 않는다.
#   이 모듈은 사용자가 지정한 finite representative embedding에서 유도된 metric을 받는다.
# === 한국어 파일 안내 끝 ===
"""Assumption-explicit mechanics for the reduced spacing--registry state.

Let q=(a,s) be a two-coordinate holonomic description of a *finite*
representative atomic/layer set with microscopic positions R_A(q).  Its
kinetic energy is

    T = 1/2 qdot_i G_ij(q) qdot_j,

where

    G_ij(q) = sum_A m_A d_i R_A . d_j R_A.

With intrinsic reduced potential U0(q) and external generalized force Q,
Euler--Lagrange gives

    G_ij qddot_j + Gamma_{i,jk} qdot_j qdot_k + d_i U0 = Q_i,

    Gamma_{i,jk} = 1/2 (d_j G_ik + d_k G_ij - d_i G_jk).

This is exact *within the declared reduced coordinate embedding*.  It is not a
claim that the repository's infinite-layer U0(a,s) by itself fixes an
inertial mass, a spatial coupling law, or dissipation.  In particular, a
literal embedding in which infinitely many layers all move with the same
collective a and s generally gives a divergent kinetic metric and must not be
used silently.

For the density-shape theory, if Q and grad(U0) depend only on q,t and no
unresolved force is omitted, conditional averaging gives

    A_i = (G^{-1})_iℓ [Q_ℓ - d_ℓ U0
          - Gamma_{ℓ,jk}(u_j u_k + Theta_jk)].

Thus the mechanics supplies the conditional acceleration A used by the exact
smooth-moment identity without introducing a Boltzmann or Fokker--Planck law.
"""

from __future__ import annotations

import numpy as np


def _square_matrix(values: np.ndarray, name: str) -> np.ndarray:
    a = np.asarray(values, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be finite")
    return a


def _positive_definite_metric(metric: np.ndarray) -> np.ndarray:
    g = _square_matrix(metric, "metric")
    if not np.allclose(g, g.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("metric must be symmetric")
    if np.min(np.linalg.eigvalsh(g)) <= 0.0:
        raise ValueError("metric must be positive definite")
    return g


def mass_metric_from_embedding_jacobians(
    jacobians: np.ndarray,
    masses: np.ndarray,
) -> np.ndarray:
    """Return G_ij=sum_A m_A d_i R_A . d_j R_A.

    Parameters
    ----------
    jacobians:
        Array with shape ``(n_objects, spatial_dimension, n_coordinates)``.
        ``jacobians[A, :, i]`` is partial R_A / partial q_i.
    masses:
        Positive masses of the represented finite objects/atoms.

    No diagonal-mass or orthogonality assumption is made.  A nonzero G_as is
    retained if the chosen coordinate embedding generates it.
    """
    j = np.asarray(jacobians, dtype=float)
    m = np.asarray(masses, dtype=float)
    if j.ndim != 3 or j.shape[0] == 0 or j.shape[2] == 0:
        raise ValueError("jacobians must have shape (objects, dimension, coordinates)")
    if m.shape != (j.shape[0],) or np.any(m <= 0.0):
        raise ValueError("masses must be positive and match the number of objects")
    if not np.all(np.isfinite(j)) or not np.all(np.isfinite(m)):
        raise ValueError("jacobians and masses must be finite")
    metric = np.einsum("A,Aki,Akj->ij", m, j, j)
    return _positive_definite_metric(metric)


def christoffel_first_kind(metric_derivatives: np.ndarray) -> np.ndarray:
    """Return Gamma_{i,jk} from d_l G_ij.

    ``metric_derivatives[l,i,j]`` means partial G_ij / partial q_l.
    The output has indices ``Gamma[i,j,k]``.
    """
    dg = np.asarray(metric_derivatives, dtype=float)
    if dg.ndim != 3 or not (dg.shape[0] == dg.shape[1] == dg.shape[2]):
        raise ValueError("metric_derivatives must have shape (n,n,n)")
    if not np.all(np.isfinite(dg)):
        raise ValueError("metric_derivatives must be finite")
    n = dg.shape[0]
    gamma = np.empty((n, n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                gamma[i, j, k] = 0.5 * (
                    dg[j, i, k] + dg[k, i, j] - dg[i, j, k]
                )
    return gamma


def generalized_acceleration(
    velocity: np.ndarray,
    intrinsic_energy_gradient: np.ndarray,
    generalized_force: np.ndarray,
    metric: np.ndarray,
    metric_derivatives: np.ndarray | None = None,
) -> np.ndarray:
    """Return qddot from the reduced Euler--Lagrange equation.

    With a constant metric, the connection term vanishes and

        qddot = G^{-1} (Q - grad U0).

    Supplying ``metric_derivatives`` retains the exact coordinate-metric
    connection term for the declared embedding.
    """
    g = _positive_definite_metric(metric)
    v = np.asarray(velocity, dtype=float)
    grad_u = np.asarray(intrinsic_energy_gradient, dtype=float)
    qforce = np.asarray(generalized_force, dtype=float)
    n = g.shape[0]
    if v.shape != (n,) or grad_u.shape != (n,) or qforce.shape != (n,):
        raise ValueError("velocity, energy gradient, and force must match metric size")
    rhs = qforce - grad_u
    if metric_derivatives is not None:
        gamma = christoffel_first_kind(metric_derivatives)
        if gamma.shape != (n, n, n):
            raise ValueError("metric derivative size must match metric")
        rhs = rhs - np.einsum("ijk,j,k->i", gamma, v, v)
    return np.linalg.solve(g, rhs)


def conditional_acceleration(
    conditional_mean_velocity: np.ndarray,
    conditional_velocity_covariance: np.ndarray,
    intrinsic_energy_gradient: np.ndarray,
    generalized_force: np.ndarray,
    metric: np.ndarray,
    metric_derivatives: np.ndarray | None = None,
) -> np.ndarray:
    """Return A=E[qddot|q] for the declared reduced mechanics.

    This function assumes that, at fixed q and t, the supplied generalized
    force and intrinsic energy gradient do not contain additional unresolved
    random forces.  Under that explicit condition,

        E[qdot_j qdot_k|q] = u_j u_k + Theta_jk

    closes the coordinate-metric connection term exactly.
    """
    g = _positive_definite_metric(metric)
    u = np.asarray(conditional_mean_velocity, dtype=float)
    theta = _square_matrix(
        conditional_velocity_covariance, "conditional_velocity_covariance"
    )
    grad_u = np.asarray(intrinsic_energy_gradient, dtype=float)
    qforce = np.asarray(generalized_force, dtype=float)
    n = g.shape[0]
    if (
        u.shape != (n,)
        or theta.shape != (n, n)
        or grad_u.shape != (n,)
        or qforce.shape != (n,)
    ):
        raise ValueError("conditional fields must match metric size")
    if not np.allclose(theta, theta.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("conditional velocity covariance must be symmetric")
    if np.min(np.linalg.eigvalsh(theta)) < -1.0e-12:
        raise ValueError("conditional velocity covariance must be positive semidefinite")

    rhs = qforce - grad_u
    if metric_derivatives is not None:
        gamma = christoffel_first_kind(metric_derivatives)
        if gamma.shape != (n, n, n):
            raise ValueError("metric derivative size must match metric")
        second_velocity_moment = np.outer(u, u) + theta
        rhs = rhs - np.einsum("ijk,jk->i", gamma, second_velocity_moment)
    return np.linalg.solve(g, rhs)


def kinetic_energy(velocity: np.ndarray, metric: np.ndarray) -> float:
    """Return T=1/2 qdot^T G qdot."""
    g = _positive_definite_metric(metric)
    v = np.asarray(velocity, dtype=float)
    if v.shape != (g.shape[0],):
        raise ValueError("velocity must match metric size")
    return 0.5 * float(v @ g @ v)


def mechanical_energy(
    velocity: np.ndarray,
    metric: np.ndarray,
    intrinsic_energy: float,
) -> float:
    """Return T+U0 for the reduced conservative subsystem."""
    if not np.isfinite(intrinsic_energy):
        raise ValueError("intrinsic_energy must be finite")
    return kinetic_energy(velocity, metric) + float(intrinsic_energy)


def external_power(velocity: np.ndarray, generalized_force: np.ndarray) -> float:
    """Return Q.qdot, the generalized external mechanical power."""
    v = np.asarray(velocity, dtype=float)
    qforce = np.asarray(generalized_force, dtype=float)
    if v.ndim != 1 or qforce.shape != v.shape:
        raise ValueError("velocity and generalized_force must be equal-size vectors")
    return float(qforce @ v)
