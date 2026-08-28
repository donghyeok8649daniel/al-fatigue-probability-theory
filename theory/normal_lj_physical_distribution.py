# === 한국어 파일 안내 시작 ===
# - 파일 역할: full nonlinear layer-LJ에 비선형 탄성 안정성과 통계역학을 적용해 zero-T, fixed-length canonical, tensile metastable P의 물리 함수형을 계산한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: physical_energy_scale, inverse_reduced_temperature, force_biased_potential, _bisect_root
#   quasistatic_stable_spacing, metastable_stationary_points, metastable_barrier_height
#   metastable_gibbs_density, metastable_tail_probability, fixed_length_two_spacing_density
#   fixed_length_canonical_density
# - 주의: thermal/canonical 및 metastable 분포는 명시된 ensemble·시간척도 가정 아래에서만 물리적으로 해석한다. 피로수명 법칙이나 escape rate를 임의로 넣지 않는다.
# === 한국어 파일 안내 끝 ===
"""Physical-statistical forms of the full nonlinear 1D layer-LJ spacing state.

This module does not fit a named distribution and does not Taylor-expand the
LJ force. It separates three physically distinct statements:

1. zero-temperature quasistatic homogeneous force balance;
2. equilibrium statistical mechanics at fixed length / temperature;
3. a metastable intact-basin Gibbs distribution under tensile force.

For the fixed-length canonical ensemble, the exact finite-M one-spacing
marginal can be evaluated from the convolution recursion

    Z_M(L) = integral z(lambda) Z_{M-1}(L-lambda) d lambda,
    z(lambda) = exp[-chi phi(lambda)].

The code uses the shifted potential phi(lambda)-phi(1), which changes each
Z_M only by an M-dependent constant and therefore leaves the normalized
fixed-M spacing marginal unchanged while improving numerical conditioning.
"""
from __future__ import annotations

import math

import numpy as np

from theory.normal_lj_chain import (
    critical_dimensionless_force,
    critical_stretch,
    normalized_lj_energy,
    normalized_lj_force,
    normalized_lj_stiffness,
)

BOLTZMANN_CONSTANT_J_PER_K = 1.380649e-23


def physical_energy_scale(
    youngs_modulus_pa: float,
    reference_area_m2: float,
    equilibrium_spacing_m: float,
) -> float:
    """Return E0 = E A0 a0 in joules."""
    if min(youngs_modulus_pa, reference_area_m2, equilibrium_spacing_m) <= 0.0:
        raise ValueError("E, A0, and a0 must all be positive")
    return youngs_modulus_pa * reference_area_m2 * equilibrium_spacing_m


def inverse_reduced_temperature(
    youngs_modulus_pa: float,
    reference_area_m2: float,
    equilibrium_spacing_m: float,
    temperature_k: float,
) -> float:
    """Return chi = E0/(k_B T)."""
    if temperature_k <= 0.0:
        raise ValueError("temperature must be positive")
    e0 = physical_energy_scale(
        youngs_modulus_pa, reference_area_m2, equilibrium_spacing_m
    )
    return e0 / (BOLTZMANN_CONSTANT_J_PER_K * temperature_k)


def force_biased_potential(
    stretch: np.ndarray | float,
    dimensionless_force: float,
    m: float = 12.19,
    n: float = 6.0,
):
    """Return w_f(lambda)=phi(lambda)-f*lambda using the full nonlinear LJ law."""
    lam = np.asarray(stretch, dtype=float)
    if np.any(lam <= 0.0):
        raise ValueError("stretch must be positive")
    return normalized_lj_energy(lam, m, n) - dimensionless_force * lam


def _bisect_root(func, lo: float, hi: float, *, iterations: int = 120) -> float:
    flo = float(func(lo))
    fhi = float(func(hi))
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        raise ValueError("root is not bracketed")
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        fm = float(func(mid))
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def quasistatic_stable_spacing(
    dimensionless_force: float,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Return the stable homogeneous spacing satisfying phi'(lambda)=f."""
    fc = critical_dimensionless_force(m, n)
    if dimensionless_force < 0.0 or dimensionless_force >= fc:
        raise ValueError("require 0 <= dimensionless_force < f_c")
    if dimensionless_force == 0.0:
        return 1.0
    lam_c = critical_stretch(m, n)
    return _bisect_root(
        lambda lam: normalized_lj_force(lam, m, n) - dimensionless_force,
        1.0,
        lam_c,
    )


def metastable_stationary_points(
    dimensionless_force: float,
    m: float = 12.19,
    n: float = 6.0,
) -> tuple[float, float]:
    """Return (stable_spacing, barrier_spacing) for 0 < f < f_c."""
    fc = critical_dimensionless_force(m, n)
    if not (0.0 < dimensionless_force < fc):
        raise ValueError("metastable tensile basin requires 0 < f < f_c")
    stable = quasistatic_stable_spacing(dimensionless_force, m, n)
    lam_c = critical_stretch(m, n)

    def residual(lam: float) -> float:
        return float(normalized_lj_force(lam, m, n) - dimensionless_force)

    lo = lam_c
    hi = max(2.0 * lam_c, 2.0)
    while residual(hi) > 0.0:
        hi *= 1.5
        if hi > 1.0e8:
            raise RuntimeError("failed to bracket the metastable barrier root")
    barrier = _bisect_root(residual, lo, hi)
    return stable, barrier


def metastable_barrier_height(
    dimensionless_force: float,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Return Delta w = w_f(lambda_b)-w_f(lambda_s), dimensionless."""
    stable, barrier = metastable_stationary_points(dimensionless_force, m, n)
    ws = float(force_biased_potential(stable, dimensionless_force, m, n))
    wb = float(force_biased_potential(barrier, dimensionless_force, m, n))
    return wb - ws


def metastable_gibbs_density(
    stretch: np.ndarray,
    dimensionless_force: float,
    inverse_temperature: float,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Return the normalized intact-basin metastable Gibbs density."""
    lam = np.asarray(stretch, dtype=float)
    if lam.ndim != 1 or lam.size < 3 or np.any(np.diff(lam) <= 0.0):
        raise ValueError("stretch must be a strictly increasing 1D grid")
    if inverse_temperature <= 0.0:
        raise ValueError("inverse_temperature chi must be positive")

    stable, barrier = metastable_stationary_points(dimensionless_force, m, n)
    mask = (lam > 0.0) & (lam < barrier)
    if np.count_nonzero(mask) < 2:
        raise ValueError("grid must resolve the intact metastable basin")

    density = np.zeros_like(lam)
    ws = float(force_biased_potential(stable, dimensionless_force, m, n))
    wf = force_biased_potential(lam[mask], dimensionless_force, m, n)
    logw = -inverse_temperature * (wf - ws)
    logw -= float(np.max(logw))
    density[mask] = np.exp(logw)
    norm = float(np.trapezoid(density, lam))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError("metastable density normalization failed")
    return density / norm


def metastable_tail_probability(
    stretch: np.ndarray,
    density: np.ndarray,
    m: float = 12.19,
    n: float = 6.0,
) -> float:
    """Integrate the instantaneous metastable mass from lambda_c upward."""
    lam = np.asarray(stretch, dtype=float)
    p = np.asarray(density, dtype=float)
    if lam.shape != p.shape or lam.ndim != 1 or np.any(np.diff(lam) <= 0.0):
        raise ValueError("stretch and density must be matching increasing 1D arrays")
    lam_c = critical_stretch(m, n)
    if lam_c <= lam[0]:
        return float(np.trapezoid(p, lam))
    if lam_c >= lam[-1]:
        return 0.0
    idx = int(np.searchsorted(lam, lam_c))
    tail_lam = np.concatenate(([lam_c], lam[idx:]))
    p_at_c = float(np.interp(lam_c, lam, p))
    tail_p = np.concatenate(([p_at_c], p[idx:]))
    return float(np.trapezoid(tail_p, tail_lam))


def _uniform_grid_step(grid: np.ndarray) -> float:
    diffs = np.diff(grid)
    if grid.ndim != 1 or grid.size < 4 or np.any(diffs <= 0.0):
        raise ValueError("grid must be a strictly increasing 1D array")
    dx = float(diffs[0])
    if not np.allclose(diffs, dx, rtol=2.0e-10, atol=2.0e-13):
        raise ValueError("fixed-length convolution requires a uniform grid")
    if abs(float(grid[0])) > max(1.0e-13, 1.0e-10 * dx):
        raise ValueError("fixed-length convolution grid must start at zero")
    return dx


def _shifted_boltzmann_weight(
    grid: np.ndarray,
    inverse_temperature: float,
    m: float,
    n: float,
) -> np.ndarray:
    weight = np.zeros_like(grid, dtype=float)
    positive = grid > 0.0
    phi_eq = float(normalized_lj_energy(1.0, m, n))
    psi = normalized_lj_energy(grid[positive], m, n) - phi_eq
    weight[positive] = np.exp(-inverse_temperature * psi)
    return weight


def _fft_convolution_prefix(a: np.ndarray, b: np.ndarray, dx: float) -> np.ndarray:
    n = a.size
    nfft = 1 << (2 * n - 1).bit_length()
    conv = np.fft.irfft(np.fft.rfft(a, nfft) * np.fft.rfft(b, nfft), nfft)
    conv = conv[:n] * dx
    conv[conv < 0.0] = 0.0
    return conv


def fixed_length_canonical_density(
    grid: np.ndarray,
    total_stretch: float,
    spacing_count: int,
    inverse_temperature: float,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Numerically evaluate the exact finite-M fixed-length canonical marginal.

    The physical formula is

        P_M(lambda|L,chi)
          = z(lambda) Z_{M-1}(L-lambda) / Z_M(L),

    where Z_M is the M-fold convolution of z(lambda)=exp[-chi phi(lambda)].
    The additive shift phi -> phi-phi(1) used numerically cancels exactly from
    the normalized fixed-M marginal.

    The convolution is a numerical quadrature of the exact ensemble formula;
    grid refinement is therefore required for quantitative use.
    """
    lam = np.asarray(grid, dtype=float)
    dx = _uniform_grid_step(lam)
    if spacing_count < 2:
        raise ValueError("spacing_count must be at least 2")
    if inverse_temperature <= 0.0 or total_stretch <= 0.0:
        raise ValueError("inverse_temperature and total_stretch must be positive")
    if total_stretch > float(lam[-1]):
        raise ValueError("grid must extend at least to total_stretch")

    z = _shifted_boltzmann_weight(lam, inverse_temperature, m, n)
    prev = z.copy()
    peak = float(np.max(prev))
    if peak <= 0.0:
        raise ValueError("Boltzmann weight vanished on the supplied grid")
    prev /= peak

    for _count in range(2, spacing_count):
        prev = _fft_convolution_prefix(z, prev, dx)
        peak = float(np.max(prev))
        if not math.isfinite(peak) or peak <= 0.0:
            raise ValueError("partition recursion lost numerical support")
        prev /= peak

    partner = total_stretch - lam
    reservoir = np.interp(partner, lam, prev, left=0.0, right=0.0)
    one_spacing_weight = z * reservoir
    one_spacing_weight[(lam <= 0.0) | (lam >= total_stretch)] = 0.0
    norm = float(np.trapezoid(one_spacing_weight, lam))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError("fixed-length canonical marginal normalization failed")
    return one_spacing_weight / norm


def fixed_length_two_spacing_density(
    stretch: np.ndarray,
    total_stretch: float,
    inverse_temperature: float,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Exact canonical M=2 one-spacing density at fixed total stretch."""
    lam = np.asarray(stretch, dtype=float)
    if lam.ndim != 1 or lam.size < 3 or np.any(np.diff(lam) <= 0.0):
        raise ValueError("stretch must be a strictly increasing 1D grid")
    if total_stretch <= 0.0 or inverse_temperature <= 0.0:
        raise ValueError("total_stretch and inverse_temperature must be positive")

    partner = total_stretch - lam
    mask = (lam > 0.0) & (partner > 0.0)
    density = np.zeros_like(lam)
    if np.count_nonzero(mask) < 2:
        raise ValueError("grid does not resolve the positive fixed-length domain")

    logw = -inverse_temperature * (
        normalized_lj_energy(lam[mask], m, n)
        + normalized_lj_energy(partner[mask], m, n)
    )
    logw -= float(np.max(logw))
    density[mask] = np.exp(logw)
    norm = float(np.trapezoid(density, lam))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError("fixed-length density normalization failed")
    return density / norm
