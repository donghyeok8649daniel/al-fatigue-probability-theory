# === 한국어 파일 안내 시작 ===
# - 파일 역할: 현재 reduced U0(a,s)에서 registry symmetry manifold와 normal/registry
#   curvature instability를 direct (k,p) sum으로 계산한다.
# - FCC, Boltzmann, Fokker-Planck, damping, random noise를 사용하지 않는다.
# - 핵심 함수: curvature_aa_direct, curvature_ss_direct, curvature_ass_direct,
#   find_curvature_zero.
# === 한국어 파일 안내 끝 ===
"""Registry symmetry and local stability of the active reduced U0(a,s).

For

    U0(a,s)=C_mn eps [sigma^m H_m(a,s)-sigma^n H_n(a,s)]

with

    H_q=sum_{k>=1,p in Z} [(k a)^2+(p b+s)^2]^(-q/2),

this module evaluates local normal and registry curvatures directly from the
termwise derivatives. Integer truncations are numerical convergence controls,
not physical cutoffs.

At a registry-symmetry point s0 for which dU/ds(a,s0)=0 for every a, the
zero-drive reduced equation mu_s sddot=-dU/ds has the exact invariant solution
s(t)=s0, sdot(t)=0.  Small perturbations xi=s-s0 obey

    mu_s xi_ddot + U_ss(a(t),s0) xi = 0

at linear order.  Thus negative U_ss is a static registry instability, while a
periodic positive U_ss(a(t),s0) can in principle produce Hill/parametric
instability for a nonzero physical seed.  Exact zero perturbation remains zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from theory.registry_lattice import MultilayerPotentialParameters, generalized_lj_coefficient


@dataclass(frozen=True)
class CurvatureRoots:
    normal_zero_a: float
    registry_zero_a: float


def _validate(a: float, s: float, params: MultilayerPotentialParameters, kmax: int, pmax: int) -> None:
    params.validate()
    if not (math.isfinite(a) and a > 0.0 and math.isfinite(s)):
        raise ValueError("a must be positive and a,s finite")
    if kmax < 1 or pmax < 1:
        raise ValueError("kmax and pmax must be positive")


def _grids(a: float, s: float, params: MultilayerPotentialParameters, kmax: int, pmax: int):
    _validate(a, s, params, kmax, pmax)
    k = np.arange(1, kmax + 1, dtype=float)[:, None]
    p = np.arange(-pmax, pmax + 1, dtype=float)[None, :]
    y = p * params.b + s
    r2 = (k * a) ** 2 + y**2
    return k, y, r2


def _inverse_power_aa(q: float, a: float, k: np.ndarray, r2: np.ndarray) -> float:
    return float(np.sum(
        -q * k**2 * r2 ** (-(q + 2.0) / 2.0)
        + q * (q + 2.0) * k**4 * a**2 * r2 ** (-(q + 4.0) / 2.0)
    ))


def _inverse_power_ss(q: float, y: np.ndarray, r2: np.ndarray) -> float:
    return float(np.sum(
        -q * r2 ** (-(q + 2.0) / 2.0)
        + q * (q + 2.0) * y**2 * r2 ** (-(q + 4.0) / 2.0)
    ))


def _inverse_power_ass(q: float, a: float, k: np.ndarray, y: np.ndarray, r2: np.ndarray) -> float:
    return float(np.sum(
        q * (q + 2.0) * k**2 * a * r2 ** (-(q + 4.0) / 2.0)
        - q * (q + 2.0) * (q + 4.0) * k**2 * a * y**2
        * r2 ** (-(q + 6.0) / 2.0)
    ))


def _combine(dm: float, dn: float, params: MultilayerPotentialParameters) -> float:
    c = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj
    return float(c * (params.sigma_lj**params.m * dm - params.sigma_lj**params.n * dn))


def curvature_aa_direct(
    a: float,
    s: float,
    params: MultilayerPotentialParameters,
    *,
    kmax: int = 120,
    pmax: int = 300,
) -> float:
    """Return d^2 U0 / da^2 from a direct truncated double sum."""
    k, y, r2 = _grids(a, s, params, kmax, pmax)
    return _combine(
        _inverse_power_aa(params.m, a, k, r2),
        _inverse_power_aa(params.n, a, k, r2),
        params,
    )


def curvature_ss_direct(
    a: float,
    s: float,
    params: MultilayerPotentialParameters,
    *,
    kmax: int = 120,
    pmax: int = 300,
) -> float:
    """Return d^2 U0 / ds^2 from a direct truncated double sum."""
    k, y, r2 = _grids(a, s, params, kmax, pmax)
    del k
    return _combine(
        _inverse_power_ss(params.m, y, r2),
        _inverse_power_ss(params.n, y, r2),
        params,
    )


def curvature_ass_direct(
    a: float,
    s: float,
    params: MultilayerPotentialParameters,
    *,
    kmax: int = 120,
    pmax: int = 300,
) -> float:
    """Return d/da(d^2 U0/ds^2), the normal modulation of registry stiffness."""
    k, y, r2 = _grids(a, s, params, kmax, pmax)
    return _combine(
        _inverse_power_ass(params.m, a, k, y, r2),
        _inverse_power_ass(params.n, a, k, y, r2),
        params,
    )


def find_curvature_zero(function, lower: float, upper: float, *, iterations: int = 80) -> float:
    """Bisection root for a sign-changing scalar curvature function."""
    lo, hi = float(lower), float(upper)
    flo, fhi = float(function(lo)), float(function(hi))
    if not (math.isfinite(flo) and math.isfinite(fhi) and flo * fhi < 0.0):
        raise ValueError("root interval must contain a finite strict sign change")
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        fm = float(function(mid))
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def local_registry_frequency_squared(curvature_ss: float, registry_inertia: float) -> float:
    """Return omega_s^2=U_ss/mu_s for the linearized registry coordinate."""
    if not math.isfinite(registry_inertia) or registry_inertia <= 0.0:
        raise ValueError("registry_inertia must be positive and finite")
    return float(curvature_ss / registry_inertia)
