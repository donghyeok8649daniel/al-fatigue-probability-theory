# === 한국어 파일 안내 시작 ===
# - 파일 역할: full nonlinear layer-LJ에 비선형 탄성 안정성과 통계역학을 적용해 zero-T, fixed-length canonical, tensile metastable P의 물리 함수형을 계산한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: physical_energy_scale, inverse_reduced_temperature, force_biased_potential, _bisect_root
#   quasistatic_stable_spacing, metastable_stationary_points, metastable_barrier_height
#   metastable_gibbs_density, metastable_tail_probability, fixed_length_two_spacing_density
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Physical-statistical forms of the full nonlinear 1D layer-LJ spacing state.

This module does not fit a named distribution and does not Taylor-expand the
LJ force.  It separates three physically distinct statements:

1. zero-temperature quasistatic homogeneous force balance;
2. equilibrium statistical mechanics at fixed length / temperature;
3. a metastable intact-basin Gibbs distribution under tensile force.

The normalized layer potential is phi(lambda).  If the physical effective
layer-patch energy is U = E0 * phi + const and phi''(1)=1, the calibration
E = (a0/A0) U''(a0) gives E0 = E A0 a0 exactly.

A tensile force f = sigma/E biases the reduced potential to

    w_f(lambda) = phi(lambda) - f lambda.

For 0 < f < f_c this function has a stable stationary point below lambda_c
and an unstable barrier point above lambda_c.  The full-domain tensile Gibbs
integral diverges because w_f -> -infinity as lambda -> infinity; therefore a
normalizable tensile distribution can only be an intact-basin/metastable
object unless length is constrained.
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
    """Return E0 = E A0 a0 in joules.

    This follows exactly from the model normalization phi''(1)=1 and the
    calibration E=(a0/A0) U''(a0), provided U=E0*phi(a/a0)+const.
    """
    if min(youngs_modulus_pa, reference_area_m2, equilibrium_spacing_m) <= 0.0:
        raise ValueError("E, A0, and a0 must all be positive")
    return youngs_modulus_pa * reference_area_m2 * equilibrium_spacing_m


def inverse_reduced_temperature(
    youngs_modulus_pa: float,
    reference_area_m2: float,
    equilibrium_spacing_m: float,
    temperature_k: float,
) -> float:
    """Return chi = E0/(k_B T).

    The representative layer area A0 is a physical input, not a fitting
    parameter.  No aluminum temperature prediction is meaningful until A0 is
    defined consistently with the coarse-grained layer potential.
    """
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
    """Return the stable homogeneous spacing satisfying phi'(lambda)=f.

    For this active tensile implementation 0 <= f < f_c is required.  At
    f=0 the stable state is exactly lambda=1.  The result is the T=0
    quasistatic homogeneous state, for which P is a delta measure at this
    spacing; applying it to a finite-rate driven state is a quasistatic
    approximation rather than an exact dynamic claim.
    """
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
    """Return (stable_spacing, barrier_spacing) for 0 < f < f_c.

    The stable point lies in (1, lambda_c), while the barrier point lies above
    lambda_c.  No polynomial/Taylor approximation is used.
    """
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
    """Return the normalized intact-basin metastable Gibbs density.

    P_ms(lambda) is proportional to exp[-chi(w_f(lambda)-w_f(lambda_s))]
    on 0 < lambda < lambda_b and is zero outside that basin.  This is a
    controlled local-equilibrium/metastable approximation, not a global
    equilibrium distribution and not a fatigue escape-rate model.
    """
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


def fixed_length_two_spacing_density(
    stretch: np.ndarray,
    total_stretch: float,
    inverse_temperature: float,
    m: float = 12.19,
    n: float = 6.0,
) -> np.ndarray:
    """Exact canonical M=2 one-spacing density at fixed total stretch.

    For two positive spacings with lambda_1+lambda_2=L,

        P(lambda_1|L,chi) ~ exp[-chi(phi(lambda_1)+phi(L-lambda_1))].

    This is an exact configurational canonical result for the stated M=2
    fixed-length ensemble.  It is not a driven nonequilibrium fatigue law.
    """
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
