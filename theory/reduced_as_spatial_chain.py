# === 한국어 파일 안내 시작 ===
# - 파일 역할: 현재 local U0(a,s)를 반복 cell에 배치한 최소 deterministic spatial chain의
#   힘, 가속도, 에너지, empirical state measure용 진단량을 계산한다.
# - 핵심 가정: V_M=sum_i U0(a_i,s_i)라는 local-state additivity를 새로 선언한다.
#   이는 full many-body/FCC 모델이 아니며, noise/damping/FP/Boltzmann을 넣지 않는다.
# - 주요 함수: intrinsic_gradients, atom_acceleration, spacing_acceleration,
#   registry_acceleration, total_mechanical_energy, external_power,
#   empirical_state_observables.
# === 한국어 파일 안내 끝 ===
"""Minimal deterministic spatial bridge for the active reduced state (a,s).

The purpose of this module is to expose one physically transparent source of
state diversity without prescribing random noise or a named probability law.

For M serial cells, let

    a_i = x_{i+1} - x_i,
    s_i = local unwrapped registry internal coordinate,

and postulate the *local-state additive* intrinsic energy

    V_M = sum_i U0(a_i, s_i).

This is a new reduced-model assumption.  It is not claimed to be an exact
nonuniform many-body lattice energy.  It is the minimal spatial replication of
the already active local state energy U0(a,s), and its uniform state has the
same intrinsic energy per cell U0(a,s).

The left node x_0 is fixed. A generalized normal force Q_a(t) acts at the right
node.  Each serial registry coordinate may receive a declared generalized
force q_s,i(t).  With equal node mass m_x and registry inertia mu_s,

    m_x xddot_j = d_a U0(a_j,s_j) - d_a U0(a_{j-1},s_{j-1})

for interior nodes, while the right boundary receives Q_a. Registry dynamics is

    mu_s sddot_i = q_s,i - d_s U0(a_i,s_i).

No damping or stochastic term is present.  Consequently any spread of
(a_i,s_i) comes from deterministic spatial nonuniformity/boundary propagation
and the mixed dependence of U0 on a and s.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from theory.registry_lattice import (
    MultilayerPotentialParameters,
    dU_da,
    dU_ds,
    u0,
)


@dataclass(frozen=True)
class EmpiricalStateObservables:
    """Finite-cell one-point observables of the deterministic spatial state."""

    mean_a: float
    variance_a: float
    mean_s: float
    variance_s: float
    covariance_as: float
    mean_intrinsic_energy: float
    mean_a_velocity: float
    mean_s_velocity: float
    theta_aa_global: float
    theta_as_global: float
    theta_ss_global: float


def _one_dimensional(name: str, values: np.ndarray, *, positive: bool = False) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.ndim != 1 or out.size == 0 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a non-empty finite one-dimensional array")
    if positive and np.any(out <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return out


def intrinsic_gradients(
    spacing: np.ndarray,
    registry: np.ndarray,
    params: MultilayerPotentialParameters,
) -> tuple[np.ndarray, np.ndarray]:
    """Return arrays (partial_a U0, partial_s U0) cell by cell."""
    a = _one_dimensional("spacing", spacing, positive=True)
    s = _one_dimensional("registry", registry)
    if s.shape != a.shape:
        raise ValueError("spacing and registry must have the same shape")
    grad_a = np.empty_like(a)
    grad_s = np.empty_like(a)
    for i, (ai, si) in enumerate(zip(a, s, strict=True)):
        grad_a[i] = float(dU_da(float(ai), float(si), params))
        grad_s[i] = float(dU_ds(float(ai), float(si), params))
    return grad_a, grad_s


def atom_acceleration(
    spacing: np.ndarray,
    registry: np.ndarray,
    normal_boundary_force: float,
    node_mass: float,
    params: MultilayerPotentialParameters,
) -> np.ndarray:
    """Return accelerations of nodes x_0,...,x_M with x_0 fixed.

    For M cells there are M+1 nodes.  The first returned acceleration is
    exactly zero.  Interior and right-boundary forces follow from

        V_M = sum_i U0(x_{i+1}-x_i, s_i)

    plus the external right-boundary generalized force Q_a.
    """
    if not math.isfinite(normal_boundary_force):
        raise ValueError("normal_boundary_force must be finite")
    if not math.isfinite(node_mass) or node_mass <= 0.0:
        raise ValueError("node_mass must be positive and finite")
    grad_a, _ = intrinsic_gradients(spacing, registry, params)
    m = grad_a.size
    acc = np.zeros(m + 1, dtype=float)
    if m > 1:
        acc[1:m] = (grad_a[1:] - grad_a[:-1]) / node_mass
    acc[m] = (-grad_a[-1] + normal_boundary_force) / node_mass
    return acc


def spacing_acceleration(
    spacing: np.ndarray,
    registry: np.ndarray,
    normal_boundary_force: float,
    node_mass: float,
    params: MultilayerPotentialParameters,
) -> np.ndarray:
    """Return exact ddot(a_i)=xddot_{i+1}-xddot_i for the reduced chain."""
    acc_x = atom_acceleration(
        spacing, registry, normal_boundary_force, node_mass, params
    )
    return np.diff(acc_x)


def registry_acceleration(
    spacing: np.ndarray,
    registry: np.ndarray,
    registry_force: float | np.ndarray,
    registry_inertia: float,
    params: MultilayerPotentialParameters,
) -> np.ndarray:
    """Return sddot_i=(q_s,i-partial_s U0)/mu_s with no damping/noise."""
    if not math.isfinite(registry_inertia) or registry_inertia <= 0.0:
        raise ValueError("registry_inertia must be positive and finite")
    a = _one_dimensional("spacing", spacing, positive=True)
    s = _one_dimensional("registry", registry)
    if s.shape != a.shape:
        raise ValueError("spacing and registry must have the same shape")
    q = np.asarray(registry_force, dtype=float)
    if q.ndim == 0:
        q = np.full_like(a, float(q))
    if q.shape != a.shape or not np.all(np.isfinite(q)):
        raise ValueError("registry_force must be finite scalar or match cell shape")
    _, grad_s = intrinsic_gradients(a, s, params)
    return (q - grad_s) / registry_inertia


def total_intrinsic_energy(
    spacing: np.ndarray,
    registry: np.ndarray,
    params: MultilayerPotentialParameters,
) -> float:
    """Return sum_i U0(a_i,s_i) under the declared local-additivity model."""
    a = _one_dimensional("spacing", spacing, positive=True)
    s = _one_dimensional("registry", registry)
    if s.shape != a.shape:
        raise ValueError("spacing and registry must have the same shape")
    return float(sum(float(u0(float(ai), float(si), params)) for ai, si in zip(a, s, strict=True)))


def total_mechanical_energy(
    node_velocity: np.ndarray,
    registry_velocity: np.ndarray,
    spacing: np.ndarray,
    registry: np.ndarray,
    node_mass: float,
    registry_inertia: float,
    params: MultilayerPotentialParameters,
) -> float:
    """Return T_x+T_s+V_M; x_0 velocity must be supplied explicitly."""
    vx = _one_dimensional("node_velocity", node_velocity)
    vs = _one_dimensional("registry_velocity", registry_velocity)
    a = _one_dimensional("spacing", spacing, positive=True)
    s = _one_dimensional("registry", registry)
    if vx.shape != (a.size + 1,) or vs.shape != a.shape or s.shape != a.shape:
        raise ValueError("velocity/state shapes are inconsistent")
    if abs(vx[0]) > 1.0e-14:
        raise ValueError("left fixed node must have zero velocity")
    if min(node_mass, registry_inertia) <= 0.0:
        raise ValueError("inertias must be positive")
    kinetic_x = 0.5 * node_mass * float(vx @ vx)
    kinetic_s = 0.5 * registry_inertia * float(vs @ vs)
    return kinetic_x + kinetic_s + total_intrinsic_energy(a, s, params)


def external_power(
    right_node_velocity: float,
    registry_velocity: np.ndarray,
    normal_boundary_force: float,
    registry_force: float | np.ndarray,
) -> float:
    """Return Q_a*x_dot_right + sum_i q_s,i*s_dot_i."""
    if not math.isfinite(right_node_velocity) or not math.isfinite(normal_boundary_force):
        raise ValueError("normal boundary quantities must be finite")
    vs = _one_dimensional("registry_velocity", registry_velocity)
    q = np.asarray(registry_force, dtype=float)
    if q.ndim == 0:
        q = np.full_like(vs, float(q))
    if q.shape != vs.shape or not np.all(np.isfinite(q)):
        raise ValueError("registry_force must be finite scalar or match registry velocity")
    return float(normal_boundary_force * right_node_velocity + q @ vs)


def empirical_state_observables(
    spacing: np.ndarray,
    registry: np.ndarray,
    spacing_velocity: np.ndarray,
    registry_velocity: np.ndarray,
    params: MultilayerPotentialParameters,
    *,
    reference_a: float | None = None,
    reference_s: float | None = None,
) -> EmpiricalStateObservables:
    """Return finite-cell G1/G2-like and global velocity-covariance diagnostics.

    The velocity covariance here is the unconditional finite-cell covariance.
    It is **not** the conditional field Theta(a,s,t) required by the divided
    density-shape identity.  Estimating that field requires local conditioning
    or coarse graining in state space.
    """
    a = _one_dimensional("spacing", spacing, positive=True)
    s = _one_dimensional("registry", registry)
    va = _one_dimensional("spacing_velocity", spacing_velocity)
    vs = _one_dimensional("registry_velocity", registry_velocity)
    if not (s.shape == va.shape == vs.shape == a.shape):
        raise ValueError("all cell arrays must have the same shape")

    energies = np.array([
        float(u0(float(ai), float(si), params))
        for ai, si in zip(a, s, strict=True)
    ])
    if (reference_a is None) != (reference_s is None):
        raise ValueError("reference_a and reference_s must be supplied together")
    if reference_a is not None:
        ref = float(u0(float(reference_a), float(reference_s), params))
        energies = energies - ref

    centered_a = a - np.mean(a)
    centered_s = s - np.mean(s)
    centered_va = va - np.mean(va)
    centered_vs = vs - np.mean(vs)
    return EmpiricalStateObservables(
        mean_a=float(np.mean(a)),
        variance_a=float(np.mean(centered_a**2)),
        mean_s=float(np.mean(s)),
        variance_s=float(np.mean(centered_s**2)),
        covariance_as=float(np.mean(centered_a * centered_s)),
        mean_intrinsic_energy=float(np.mean(energies)),
        mean_a_velocity=float(np.mean(va)),
        mean_s_velocity=float(np.mean(vs)),
        theta_aa_global=float(np.mean(centered_va**2)),
        theta_as_global=float(np.mean(centered_va * centered_vs)),
        theta_ss_global=float(np.mean(centered_vs**2)),
    )
