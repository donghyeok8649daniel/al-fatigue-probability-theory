from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.special import gamma, kv, kvp


@dataclass(frozen=True)
class FourierLatticeConfig:
    """Numerical controls for the Poisson-summed infinite lower lattice.

    ``tol`` controls truncation of the exponentially convergent reciprocal-
    lattice series. ``max_modes`` is only a safety ceiling; convergence is
    tolerance-led rather than based on a fixed real-space image cutoff.
    """

    tol: float = 1.0e-13
    max_modes: int = 64
    consecutive_small: int = 2


def _as_arrays(a, s) -> tuple[np.ndarray, np.ndarray]:
    aa, ss = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(s, dtype=float))
    if np.any(aa <= 0.0):
        raise ValueError("interlayer spacing a must be positive")
    return aa, ss


def _zero_mode(a: np.ndarray, p: int, b: float) -> tuple[np.ndarray, np.ndarray]:
    """m=0 Poisson term A_p(a) and exact dA_p/da."""

    value = (
        math.sqrt(math.pi)
        * gamma(p - 0.5)
        / (b * gamma(p))
        * a ** (1 - 2 * p)
    )
    deda = (1 - 2 * p) * value / a
    return value, deda


def _zero_mode_second_deda(a: np.ndarray, p: int, b: float) -> np.ndarray:
    """Exact d^2 A_p / da^2 for the m=0 Poisson term."""

    value, _ = _zero_mode(a, p, b)
    return (1 - 2 * p) * (-2 * p) * value / (a * a)


def _mode_coefficient_and_deda(
    a: np.ndarray,
    p: int,
    m: int,
    b: float,
) -> tuple[np.ndarray, np.ndarray]:
    """B_{p,m}(a) and its exact derivative with respect to a."""

    nu = p - 0.5
    c = 2.0 * math.pi * m / b
    z = c * a
    prefactor = (
        4.0
        * math.sqrt(math.pi)
        / (b * gamma(p))
        * (math.pi * m / b) ** nu
    )
    a_power = a ** (-nu)
    kval = kv(nu, z)
    coefficient = prefactor * a_power * kval
    dcoefficient = prefactor * (
        -nu * a ** (-nu - 1.0) * kval
        + a_power * c * kvp(nu, z, 1)
    )
    return coefficient, dcoefficient


def _mode_coefficient_second_deda(
    a: np.ndarray,
    p: int,
    m: int,
    b: float,
) -> np.ndarray:
    """Exact second a-derivative of B_{p,m}(a)."""

    nu = p - 0.5
    c = 2.0 * math.pi * m / b
    z = c * a
    prefactor = (
        4.0
        * math.sqrt(math.pi)
        / (b * gamma(p))
        * (math.pi * m / b) ** nu
    )
    kval = kv(nu, z)
    kprime = kvp(nu, z, 1)
    ksecond = kvp(nu, z, 2)
    return prefactor * (
        nu * (nu + 1.0) * a ** (-nu - 2.0) * kval
        - 2.0 * nu * a ** (-nu - 1.0) * c * kprime
        + a ** (-nu) * c * c * ksecond
    )


def power_lattice_sum_and_gradient(
    a,
    s,
    *,
    p: int,
    b: float = 1.0,
    config: FourierLatticeConfig = FourierLatticeConfig(),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    r"""Evaluate the staggered infinite-row power lattice sum.

    The real-space sum is

        S_p(a,s) = sum_n [a^2 + ((n+1/2)b - s)^2]^{-p}.

    Poisson summation gives

        S_p = A_p(a)
              + sum_{m>=1} B_{p,m}(a)
                cos[2*pi*m*(1/2 - s/b)].

    The half-cell phase is essential: ``s=0`` is the staggered registry used
    by ``TwoRowLJ`` where the upper row is located at x=(i+1/2)b+s_i.

    Returns ``(S_p, dS_p/da, dS_p/ds, modes_used)``.
    """

    if p <= 0:
        raise ValueError("p must be positive")
    if b <= 0.0:
        raise ValueError("b must be positive")
    if (
        config.tol <= 0.0
        or config.max_modes < 1
        or config.consecutive_small < 1
    ):
        raise ValueError("invalid Fourier lattice controls")

    aa, ss = _as_arrays(a, s)
    total, deda = _zero_mode(aa, p, b)
    deds = np.zeros_like(total)

    small_count = 0
    modes_used = config.max_modes
    for m in range(1, config.max_modes + 1):
        coefficient, dcoefficient = _mode_coefficient_and_deda(aa, p, m, b)
        wave_number = 2.0 * math.pi * m / b
        phase = 2.0 * math.pi * m * (0.5 - ss / b)
        cphase = np.cos(phase)
        sphase = np.sin(phase)

        energy_term = coefficient * cphase
        a_term = dcoefficient * cphase
        s_term = wave_number * coefficient * sphase

        total = total + energy_term
        deda = deda + a_term
        deds = deds + s_term

        scale = max(
            1.0,
            float(np.max(np.abs(total))),
            float(np.max(np.abs(deda))),
            float(np.max(np.abs(deds))),
        )
        term_size = max(
            float(np.max(np.abs(energy_term))),
            float(np.max(np.abs(a_term))),
            float(np.max(np.abs(s_term))),
        )
        if term_size <= config.tol * scale:
            small_count += 1
            if small_count >= config.consecutive_small:
                modes_used = m
                break
        else:
            small_count = 0

    return total, deda, deds, modes_used


def power_lattice_normal_second_derivative(
    a,
    s,
    *,
    p: int,
    b: float = 1.0,
    config: FourierLatticeConfig = FourierLatticeConfig(),
) -> tuple[np.ndarray, int]:
    """Evaluate the exact Poisson-summed d^2 S_p / da^2.

    This is used only to identify the tangent normal stiffness of the canonical
    lattice potential. It does not linearize the subsequent nonlinear PDE.
    """

    if p <= 0:
        raise ValueError("p must be positive")
    if b <= 0.0:
        raise ValueError("b must be positive")
    if (
        config.tol <= 0.0
        or config.max_modes < 1
        or config.consecutive_small < 1
    ):
        raise ValueError("invalid Fourier lattice controls")

    aa, ss = _as_arrays(a, s)
    second = _zero_mode_second_deda(aa, p, b)
    small_count = 0
    modes_used = config.max_modes

    for m in range(1, config.max_modes + 1):
        phase = 2.0 * math.pi * m * (0.5 - ss / b)
        term = _mode_coefficient_second_deda(aa, p, m, b) * np.cos(phase)
        second = second + term
        scale = max(1.0, float(np.max(np.abs(second))))
        term_size = float(np.max(np.abs(term)))
        if term_size <= config.tol * scale:
            small_count += 1
            if small_count >= config.consecutive_small:
                modes_used = m
                break
        else:
            small_count = 0

    return second, modes_used


def two_row_lj_infinite_energy_gradient(
    a,
    s,
    *,
    epsilon: float,
    sigma_lj: float,
    b: float = 1.0,
    config: FourierLatticeConfig = FourierLatticeConfig(),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Infinite lower-row LJ energy W(a,s) and its exact first derivatives.

    W(a,s) is the Poisson-summed form of the full real-space interaction
    between one upper cell and the infinite staggered lower row.  The Bessel
    functions are therefore part of the energy kernel itself, and its
    derivatives drive the Smoluchowski drift.
    """

    if epsilon <= 0.0 or sigma_lj <= 0.0:
        raise ValueError("epsilon and sigma_lj must be positive")

    s6, s6a, s6s, modes6 = power_lattice_sum_and_gradient(
        a, s, p=6, b=b, config=config
    )
    s3, s3a, s3s, modes3 = power_lattice_sum_and_gradient(
        a, s, p=3, b=b, config=config
    )

    factor = 4.0 * epsilon
    energy = factor * (sigma_lj**12 * s6 - sigma_lj**6 * s3)
    deda = factor * (sigma_lj**12 * s6a - sigma_lj**6 * s3a)
    deds = factor * (sigma_lj**12 * s6s - sigma_lj**6 * s3s)
    return energy, deda, deds, max(modes6, modes3)


def two_row_lj_infinite_normal_stiffness(
    a,
    s,
    *,
    epsilon: float,
    sigma_lj: float,
    b: float = 1.0,
    config: FourierLatticeConfig = FourierLatticeConfig(),
) -> tuple[np.ndarray, int]:
    """Exact normal tangent stiffness d^2 W/da^2 of the infinite LJ row."""

    if epsilon <= 0.0 or sigma_lj <= 0.0:
        raise ValueError("epsilon and sigma_lj must be positive")
    s6aa, modes6 = power_lattice_normal_second_derivative(
        a, s, p=6, b=b, config=config
    )
    s3aa, modes3 = power_lattice_normal_second_derivative(
        a, s, p=3, b=b, config=config
    )
    stiffness = 4.0 * epsilon * (
        sigma_lj**12 * s6aa - sigma_lj**6 * s3aa
    )
    return stiffness, max(modes6, modes3)


def two_row_lj_direct_reference(
    a,
    s,
    *,
    epsilon: float,
    sigma_lj: float,
    b: float = 1.0,
    images: int = 2048,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Large real-space sum retained only to validate the analytic kernel."""

    if images < 4:
        raise ValueError("images must be at least 4")
    aa, ss = _as_arrays(a, s)
    n = np.arange(-images, images + 1, dtype=float)
    dx = (n + 0.5) * b - ss[..., None]
    rr = np.sqrt(aa[..., None] ** 2 + dx**2)

    x = sigma_lj / rr
    phi = 4.0 * epsilon * (x**12 - x**6)
    dphi = 4.0 * epsilon * (
        -12.0 * sigma_lj**12 * rr**(-13)
        + 6.0 * sigma_lj**6 * rr**(-7)
    )
    energy = np.sum(phi, axis=-1)
    deda = np.sum(dphi * aa[..., None] / rr, axis=-1)
    deds = np.sum(-dphi * dx / rr, axis=-1)
    return energy, deda, deds
