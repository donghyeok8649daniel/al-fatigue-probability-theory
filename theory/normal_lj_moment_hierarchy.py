# === 한국어 파일 안내 시작 ===
# - 파일 역할: active 1D normal LJ probability theory의 exact P-u-Theta moment
#   identities에서 빠지면 안 되는 acceleration-covariance source와 spacing-coordinate
#   kinetic metric을 계산한다.
# - Boltzmann, Gaussian PDF, Fokker-Planck, damping, empirical damage closure를 쓰지 않는다.
# === 한국어 파일 안내 끝 ===
"""Exact helper identities for the active 1D P-u-Theta formulation.

The microscopic chain is expressed in atomic positions x_j with unit masses,
while the probability projection uses spacing coordinates lambda_i and spacing
rates c_i = dot(lambda_i).  These are not independent unit-mass coordinates.

For the projected phase-space density F(lambda,c,t), define

    u      = E[c | lambda]
    Theta  = Var(c | lambda)
    C3     = E[(c-u)^3 | lambda]
    Psi    = Cov(c, ddot(lambda) | lambda)

Then the exact general second-central-moment equation is

    D_t Theta
      + 2 Theta * d_lambda u
      + (1/P) d_lambda(P C3)
      = 2 Psi.

The right-hand side vanishes only under an additional conditional-acceleration
assumption; it is not generally zero for the spatial LJ chain because the
spacing acceleration depends on neighbouring spacings.
"""

from __future__ import annotations

import math
import numpy as np


def spacing_mass_metric(number_of_spacings: int) -> np.ndarray:
    """Return the exact dimensionless mass metric in spacing coordinates.

    With the left atom fixed and

        x_j = sum_{k=1}^j lambda_k,

    the unit-mass atomic kinetic energy is

        T = 1/2 c^T G c,

    where

        G_kl = M - max(k,l) + 1

    for M spacings, using one-based mathematical indices.
    """
    m = int(number_of_spacings)
    if m < 1:
        raise ValueError("number_of_spacings must be positive")
    idx = np.arange(m)
    # zero-based equivalent of M-max(k,l)+1 for one-based indices
    return (m - np.maximum(idx[:, None], idx[None, :])).astype(float)


def spacing_rate_kinetic_energy(spacing_rates: np.ndarray) -> float:
    """Return exact unit-mass chain kinetic energy from spacing rates."""
    c = np.asarray(spacing_rates, dtype=float)
    if c.ndim != 1 or c.size < 1 or not np.all(np.isfinite(c)):
        raise ValueError("spacing_rates must be a finite nonempty vector")
    g = spacing_mass_metric(c.size)
    return 0.5 * float(c @ g @ c)


def conditional_second_spacing_rate_moment(
    mean_spacing_rate: np.ndarray | float,
    theta: np.ndarray | float,
):
    """Return E[c^2|lambda] = u^2 + Theta pointwise."""
    u = np.asarray(mean_spacing_rate, dtype=float)
    th = np.asarray(theta, dtype=float)
    if u.shape != th.shape:
        raise ValueError("mean_spacing_rate and theta must have the same shape")
    if np.any(~np.isfinite(u)) or np.any(~np.isfinite(th)) or np.any(th < 0.0):
        raise ValueError("inputs must be finite and theta non-negative")
    return u * u + th


def theta_material_rate(
    theta: np.ndarray | float,
    mean_spacing_rate_gradient: np.ndarray | float,
    third_central_flux_divergence_over_density: np.ndarray | float,
    spacing_rate_acceleration_covariance: np.ndarray | float,
):
    """Return the exact material rate D_t Theta.

    Implements

        D_t Theta
          = 2 Psi
            - 2 Theta d_lambda u
            - (1/P) d_lambda(P C3).

    No assumption Psi=0 or C3=0 is made.
    """
    th = np.asarray(theta, dtype=float)
    grad_u = np.asarray(mean_spacing_rate_gradient, dtype=float)
    div_c3 = np.asarray(third_central_flux_divergence_over_density, dtype=float)
    psi = np.asarray(spacing_rate_acceleration_covariance, dtype=float)
    if not (th.shape == grad_u.shape == div_c3.shape == psi.shape):
        raise ValueError("all fields must have the same shape")
    if any(np.any(~np.isfinite(x)) for x in (th, grad_u, div_c3, psi)):
        raise ValueError("all fields must be finite")
    if np.any(th < 0.0):
        raise ValueError("theta must be non-negative")
    return 2.0 * psi - 2.0 * th * grad_u - div_c3


def theta_balance_residual(
    material_rate_theta: np.ndarray | float,
    theta: np.ndarray | float,
    mean_spacing_rate_gradient: np.ndarray | float,
    third_central_flux_divergence_over_density: np.ndarray | float,
    spacing_rate_acceleration_covariance: np.ndarray | float,
):
    """Return residual of the exact general Theta balance."""
    dth = np.asarray(material_rate_theta, dtype=float)
    rhs_rate = theta_material_rate(
        theta,
        mean_spacing_rate_gradient,
        third_central_flux_divergence_over_density,
        spacing_rate_acceleration_covariance,
    )
    if dth.shape != rhs_rate.shape or np.any(~np.isfinite(dth)):
        raise ValueError("material_rate_theta has incompatible shape or values")
    return dth - rhs_rate


def irreversible_cycle_work(
    external_work: float,
    mechanical_energy_change: float,
) -> float:
    """Return D_irr = W_ext - Delta E_mech from the cycle first law.

    This is an accounting identity once a true irreversible mechanism exists.
    For the present conservative baseline the physically expected value is zero
    up to numerical error; a positive result must not be imposed by clipping.
    """
    w = float(external_work)
    de = float(mechanical_energy_change)
    if not (math.isfinite(w) and math.isfinite(de)):
        raise ValueError("work and energy change must be finite")
    return w - de
