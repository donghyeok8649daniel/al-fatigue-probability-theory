# === 한국어 파일 안내 시작 ===
# - 파일 역할: reduced state x=(a,s)의 exact smooth-moment density shape identity를 계산한다.
# - 주요 함수: density_log_gradient_2d, compatibility_curl_2d, reconstruct_density_path_2d,
#   probability_mass_2d, mean_intrinsic_energy_2d
# - 주의: Boltzmann, Gaussian/Weibull, Fokker-Planck, Markov, 독립성 가정을 사용하지 않는다.
#   다만 divided form은 smooth P>0 및 양의 정부호 conditional velocity covariance를 요구한다.
# === 한국어 파일 안내 끝 ===
"""Exact two-coordinate smooth-moment identity for P(a,s,t).

Let x=(a,s), v=dot(x), P(x,t) the one-point state density,

    u_i(x,t) = E[v_i | x],
    Theta_ij(x,t) = Cov(v_i,v_j | x),
    A_i(x,t) = E[ddot(x_i) | x].

The exact continuity and first velocity-moment balances are

    partial_t P + partial_j(P u_j) = 0,

    partial_t(P u_i)
      + partial_j[P(u_i u_j + Theta_ij)] = P A_i.

Therefore, wherever P>0 and Theta is invertible,

    Theta grad(log P)
      = A - D_t u - div(Theta),

and hence

    grad(log P)
      = Theta^{-1}[A - D_t u - div(Theta)].

This is a moment identity, not a kinetic closure.  In particular it does not
assume Boltzmann equilibrium, a named PDF family, a Markov bath, Fokker--Planck
dynamics, or statistical independence of neighboring microscopic states.

For exact smooth fields the right-hand side must be a gradient field.  In two
coordinates this implies the compatibility condition

    partial_a g_s - partial_s g_a = 0,

where g=grad(log P).  Non-zero numerical curl is therefore a direct diagnostic
of finite-sample error, differentiation error, or an inconsistent closure.
"""

from __future__ import annotations

import numpy as np


def _validate_axes(a: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    aa = np.asarray(a, dtype=float)
    ss = np.asarray(s, dtype=float)
    if aa.ndim != 1 or ss.ndim != 1 or aa.size < 3 or ss.size < 3:
        raise ValueError("a and s must be one-dimensional axes with at least 3 points")
    if np.any(np.diff(aa) <= 0.0) or np.any(np.diff(ss) <= 0.0):
        raise ValueError("a and s axes must be strictly increasing")
    return aa, ss


def _field(name: str, value: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def density_log_gradient_2d(
    a: np.ndarray,
    s: np.ndarray,
    theta_aa: np.ndarray,
    theta_as: np.ndarray,
    theta_ss: np.ndarray,
    conditional_acceleration_a: np.ndarray,
    conditional_acceleration_s: np.ndarray,
    material_acceleration_a: np.ndarray,
    material_acceleration_s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return g=(partial_a log P, partial_s log P) from the exact moment identity.

    Parameters are fields on a rectangular grid with shape ``(len(a), len(s))``.
    ``material_acceleration_*`` denotes D_t u, not the conditional microscopic
    acceleration.  The conditional velocity covariance is

        Theta = [[theta_aa, theta_as],
                 [theta_as, theta_ss]].

    The divided form requires Theta to be positive definite pointwise.  If a
    physical regime has singular Theta, the undivided matrix equation must be
    used instead and this routine intentionally refuses to invert it.
    """
    aa, ss = _validate_axes(a, s)
    shape = (aa.size, ss.size)
    taa = _field("theta_aa", theta_aa, shape)
    tas = _field("theta_as", theta_as, shape)
    tss = _field("theta_ss", theta_ss, shape)
    acc_a = _field("conditional_acceleration_a", conditional_acceleration_a, shape)
    acc_s = _field("conditional_acceleration_s", conditional_acceleration_s, shape)
    dtu_a = _field("material_acceleration_a", material_acceleration_a, shape)
    dtu_s = _field("material_acceleration_s", material_acceleration_s, shape)

    det = taa * tss - tas * tas
    if np.any(taa <= 0.0) or np.any(tss <= 0.0) or np.any(det <= 0.0):
        raise ValueError("conditional velocity covariance must be positive definite")

    edge_a = 2 if aa.size >= 3 else 1
    edge_s = 2 if ss.size >= 3 else 1
    d_taa_da = np.gradient(taa, aa, axis=0, edge_order=edge_a)
    d_tas_ds = np.gradient(tas, ss, axis=1, edge_order=edge_s)
    d_tas_da = np.gradient(tas, aa, axis=0, edge_order=edge_a)
    d_tss_ds = np.gradient(tss, ss, axis=1, edge_order=edge_s)

    div_theta_a = d_taa_da + d_tas_ds
    div_theta_s = d_tas_da + d_tss_ds
    rhs_a = acc_a - dtu_a - div_theta_a
    rhs_s = acc_s - dtu_s - div_theta_s

    grad_a = (tss * rhs_a - tas * rhs_s) / det
    grad_s = (-tas * rhs_a + taa * rhs_s) / det
    return grad_a, grad_s


def compatibility_curl_2d(
    a: np.ndarray,
    s: np.ndarray,
    grad_log_p_a: np.ndarray,
    grad_log_p_s: np.ndarray,
) -> np.ndarray:
    """Return partial_a g_s - partial_s g_a for g=grad(log P).

    Exact smooth fields have zero curl on a simply connected state domain.
    """
    aa, ss = _validate_axes(a, s)
    shape = (aa.size, ss.size)
    ga = _field("grad_log_p_a", grad_log_p_a, shape)
    gs = _field("grad_log_p_s", grad_log_p_s, shape)
    d_gs_da = np.gradient(gs, aa, axis=0, edge_order=2)
    d_ga_ds = np.gradient(ga, ss, axis=1, edge_order=2)
    return d_gs_da - d_ga_ds


def _cumtrapz_axis0(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    increments = 0.5 * (values[:-1, :] + values[1:, :]) * np.diff(x)[:, None]
    return np.vstack((np.zeros((1, values.shape[1])), np.cumsum(increments, axis=0)))


def _cumtrapz_axis1(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    increments = 0.5 * (values[:, :-1] + values[:, 1:]) * np.diff(x)[None, :]
    return np.hstack((np.zeros((values.shape[0], 1)), np.cumsum(increments, axis=1)))


def probability_mass_2d(a: np.ndarray, s: np.ndarray, density: np.ndarray) -> float:
    """Integrate a rectangular-grid density over da ds."""
    aa, ss = _validate_axes(a, s)
    p = _field("density", density, (aa.size, ss.size))
    return float(np.trapezoid(np.trapezoid(p, ss, axis=1), aa, axis=0))


def reconstruct_density_path_2d(
    a: np.ndarray,
    s: np.ndarray,
    grad_log_p_a: np.ndarray,
    grad_log_p_s: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Reconstruct normalized P from an integrable log-gradient on a grid.

    Two rectangular paths from the lower-left reference state are evaluated:
    (a then s) and (s then a).  Their maximum log-density difference is
    returned as ``path_mismatch``.  For exact compatible fields it vanishes up
    to quadrature/differentiation error.  The returned density uses the mean of
    the two log-potentials before normalization.
    """
    aa, ss = _validate_axes(a, s)
    shape = (aa.size, ss.size)
    ga = _field("grad_log_p_a", grad_log_p_a, shape)
    gs = _field("grad_log_p_s", grad_log_p_s, shape)

    int_a = _cumtrapz_axis0(ga, aa)
    int_s = _cumtrapz_axis1(gs, ss)

    # Path 1: lower-left -> (a, s_min) -> (a, s)
    logp_as = int_a[:, [0]] + int_s
    # Path 2: lower-left -> (a_min, s) -> (a, s)
    logp_sa = int_s[[0], :] + int_a
    mismatch = float(np.max(np.abs(logp_as - logp_sa)))

    logp = 0.5 * (logp_as + logp_sa)
    logp -= float(np.max(logp))
    p = np.exp(logp)
    mass = probability_mass_2d(aa, ss, p)
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("density normalization failed")
    return p / mass, mismatch


def mean_intrinsic_energy_2d(
    a: np.ndarray,
    s: np.ndarray,
    density: np.ndarray,
    delta_u0: np.ndarray,
    *,
    normalization_tolerance: float = 5.0e-3,
) -> float:
    """Return G2, integral Delta U0(a,s) P(a,s,t) da ds, on a grid."""
    aa, ss = _validate_axes(a, s)
    shape = (aa.size, ss.size)
    p = _field("density", density, shape)
    energy = _field("delta_u0", delta_u0, shape)
    mass = probability_mass_2d(aa, ss, p)
    if abs(mass - 1.0) > normalization_tolerance:
        raise ValueError("density must be normalized before computing mean energy")
    return float(np.trapezoid(np.trapezoid(energy * p, ss, axis=1), aa, axis=0))
