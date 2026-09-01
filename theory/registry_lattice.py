# === 한국어 파일 안내 시작 ===
# - 파일 역할: multiplicity 없는 다층 generalized-LJ 위치에너지 U0(a,s)를 계산한다.
# - 주요 클래스: RegistryLattice, MultilayerPotentialParameters
# - 주요 함수/메서드: h_q_direct, h_q_bessel, h_q_polylog, u0, dU_da, dU_ds
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Exact multiplicity-free multilayer lattice energy.

``W(d,s)`` is a row--row interaction kernel.  The intrinsic local-fatigue
potential is the interaction of one reference layer with layers at
``a, 2a, 3a, ...``::

    U0(a,s) = sum_{k>=1} W(k*a,s).

There is no multiplicity ``k`` in front of ``W`` and the same unwrapped
collective registry ``s`` is used at every normal layer.  The alternative
``sum k W(k*a,s)`` counts all cross-interface layer pairs and is not the
counting convention implemented here.

The Mellin--Poisson--Bessel and Bessel--Lambert formulae below are exact
identities.  Integer cutoffs are numerical evaluation parameters, not
physical interaction cutoffs.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Iterable

import mpmath
import numpy as np


def generalized_lj_coefficient(m: float, n: float) -> float:
    """Return C_mn when epsilon_LJ is the isolated-pair well depth."""
    if not (m > n > 0.0):
        raise ValueError("require m > n > 0")
    return float(m / (m - n) * (m / n) ** (n / (m - n)))


def v_mn(
    r: float | np.ndarray,
    m: float,
    n: float,
    epsilon_lj: float,
    sigma_lj: float,
) -> float | np.ndarray:
    """Generalized Lennard--Jones pair potential in energy units."""
    values = np.asarray(r, dtype=float)
    if np.any(values <= 0.0) or epsilon_lj <= 0.0 or sigma_lj <= 0.0:
        raise ValueError("r, epsilon_lj and sigma_lj must be positive")
    result = generalized_lj_coefficient(m, n) * epsilon_lj * (
        (sigma_lj / values) ** m - (sigma_lj / values) ** n
    )
    return float(result) if values.ndim == 0 else result


@dataclass(frozen=True)
class RegistryLattice:
    """Dimensionless parameters for U0/(C_mn*epsilon_LJ)."""
    normal_ratio: float = 1.0
    sigma_ratio: float = 1.0
    m: float = 12.19
    n: float = 6.0
    bessel_modes: int = 24
    layer_modes: int = 48

    def validate(self) -> None:
        if not (self.m > self.n > 2.0):
            raise ValueError("full multilayer energy requires m > n > 2")
        if self.normal_ratio <= 0.0 or self.sigma_ratio <= 0.0:
            raise ValueError("normal_ratio and sigma_ratio must be positive")
        if self.bessel_modes < 1 or self.layer_modes < 1:
            raise ValueError("Bessel and layer mode counts must be positive")


@dataclass(frozen=True)
class MultilayerPotentialParameters:
    """Physical parameters of the intrinsic multilayer potential."""
    b: float = 1.0
    epsilon_lj: float = 1.0
    sigma_lj: float = 1.0
    m: float = 12.0
    n: float = 6.0
    bessel_modes: int = 24
    layer_modes: int = 48

    def validate(self) -> None:
        if self.b <= 0.0 or self.epsilon_lj <= 0.0 or self.sigma_lj <= 0.0:
            raise ValueError("b, epsilon_lj and sigma_lj must be positive")
        RegistryLattice(
            1.0, self.sigma_lj / self.b, self.m, self.n,
            self.bessel_modes, self.layer_modes,
        ).validate()

    def reduced(self, a: float) -> RegistryLattice:
        self.validate()
        if a <= 0.0:
            raise ValueError("a must be positive")
        return RegistryLattice(
            a / self.b, self.sigma_lj / self.b, self.m, self.n,
            self.bessel_modes, self.layer_modes,
        )


def b_q_direct(q: float, delta: float, eta: float, pmax: int = 10000) -> float:
    """Symmetrically truncated direct definition of the single-row B_q."""
    if not (q > 1.0 and eta > 0.0 and pmax >= 1):
        raise ValueError("single-row B_q requires q>1, eta>0 and pmax>=1")
    p = np.arange(-pmax, pmax + 1, dtype=float)
    return float(np.sum(((p + delta) ** 2 + eta**2) ** (-0.5 * q)))


@lru_cache(maxsize=256)
def inverse_power_bessel_coefficients(
    q: float, eta: float, modes: int
) -> tuple[float, np.ndarray]:
    """Fourier coefficients of exact B_q, truncated only in ell."""
    if not (q > 1.0 and eta > 0.0 and modes >= 1):
        raise ValueError("single-row B_q requires q>1, eta>0 and modes>=1")
    half_q = 0.5 * q
    nu = 0.5 * (q - 1.0)
    prefactor = math.sqrt(math.pi) / float(mpmath.gamma(half_q))
    zero = prefactor * float(mpmath.gamma(nu)) * eta ** (1.0 - q)
    coefficients = np.empty(modes, dtype=float)
    for index in range(modes):
        ell = index + 1
        coefficients[index] = (
            4.0 * prefactor * (math.pi * ell / eta) ** nu
            * float(mpmath.besselk(nu, 2.0 * math.pi * ell * eta))
        )
    return zero, coefficients


def shifted_inverse_power_bessel(
    q: float, delta: float | np.ndarray, eta: float, modes: int = 24,
) -> float | np.ndarray:
    """Evaluate exact B_q by its reciprocal Fourier--Bessel series."""
    zero, coefficients = inverse_power_bessel_coefficients(q, eta, modes)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def h_q_direct(
    q: float, delta: float, eta: float, kmax: int = 200, pmax: int = 400,
) -> float:
    """Direct H_q=sum_k sum_p [(p+delta)^2+k^2 eta^2]^-q/2."""
    if not (q > 2.0 and eta > 0.0 and kmax >= 1 and pmax >= 1):
        raise ValueError("absolute multilayer H_q requires q>2 and positive cutoffs")
    p = np.arange(-pmax, pmax + 1, dtype=float)
    total = 0.0
    for k in range(1, kmax + 1):
        total += float(np.sum(((p + delta) ** 2 + (k * eta) ** 2) ** (-0.5 * q)))
    return total


def h_q_eta_derivatives_direct(
    q: float, delta: float, eta: float, kmax: int = 120, pmax: int = 240,
) -> tuple[float, float]:
    """Termwise direct first and second eta derivatives of H_q."""
    if not (q > 2.0 and eta > 0.0 and kmax >= 1 and pmax >= 1):
        raise ValueError("H_q derivatives require q>2 and positive cutoffs")
    p = np.arange(-pmax, pmax + 1, dtype=float)[None, :]
    k = np.arange(1, kmax + 1, dtype=float)[:, None]
    r2 = (p + delta) ** 2 + (k * eta) ** 2
    first = -q * eta * np.sum(k**2 * r2 ** (-0.5 * q - 1.0))
    second = (
        -q * np.sum(k**2 * r2 ** (-0.5 * q - 1.0))
        + q * (q + 2.0) * eta**2
        * np.sum(k**4 * r2 ** (-0.5 * q - 2.0))
    )
    return float(first), float(second)


def h_q_delta_derivative_direct(
    q: float, delta: float, eta: float, kmax: int = 120, pmax: int = 240,
) -> float:
    """Termwise direct delta derivative of H_q."""
    if not (q > 2.0 and eta > 0.0 and kmax >= 1 and pmax >= 1):
        raise ValueError("H_q derivative requires q>2 and positive cutoffs")
    p = np.arange(-pmax, pmax + 1, dtype=float)[None, :]
    k = np.arange(1, kmax + 1, dtype=float)[:, None]
    shifted = p + delta
    r2 = shifted**2 + (k * eta) ** 2
    return float(-q * np.sum(shifted * r2 ** (-0.5 * q - 1.0)))


@lru_cache(maxsize=8192)
def bessel_lambert(nu: float, x: float, kmax: int = 48) -> float:
    """Numerically evaluate Kcal_nu(x)=sum k^-nu K_nu(k x)."""
    if nu <= 0.0 or x <= 0.0 or kmax < 1:
        raise ValueError("require nu>0, x>0 and kmax>=1")
    return float(mpmath.fsum(
        mpmath.power(k, -nu) * mpmath.besselk(nu, k * x)
        for k in range(1, kmax + 1)
    ))


@lru_cache(maxsize=2048)
def bessel_lambert_polylog(q: int, x: float) -> float:
    """Exact half-integer K-Lambert closure for q=6 or q=12.

    K_5/2 has coefficients (1,3,3); K_11/2 has
    (1,15,105,420,945,945).  The k powers are summed as polylogarithms.
    """
    if x <= 0.0 or q not in (6, 12):
        raise ValueError("polylog closure is implemented for q=6 or q=12")
    coefficients = (1, 3, 3) if q == 6 else (1, 15, 105, 420, 945, 945)
    base_order = q // 2
    value = mpmath.fsum(
        c * mpmath.polylog(base_order + j, mpmath.e ** (-x)) / x**j
        for j, c in enumerate(coefficients)
    )
    return float(mpmath.sqrt(mpmath.pi / (2.0 * x)) * value)


@lru_cache(maxsize=256)
def h_q_bessel_coefficients(
    q: float, eta: float, modes: int = 24, layer_modes: int = 48,
) -> tuple[float, np.ndarray]:
    """Fourier coefficients of exact H_q using the Bessel--Lambert series."""
    if not (q > 2.0 and eta > 0.0 and modes >= 1 and layer_modes >= 1):
        raise ValueError("H_q requires q>2, eta>0 and positive mode counts")
    half_q = 0.5 * q
    nu = 0.5 * (q - 1.0)
    prefactor = math.sqrt(math.pi) / float(mpmath.gamma(half_q))
    zero = (
        prefactor * float(mpmath.gamma(nu)) * eta ** (1.0 - q)
        * float(mpmath.zeta(q - 1.0))
    )
    coefficients = np.empty(modes, dtype=float)
    for index in range(modes):
        ell = index + 1
        coefficients[index] = (
            4.0 * prefactor * (math.pi * ell / eta) ** nu
            * bessel_lambert(nu, 2.0 * math.pi * ell * eta, layer_modes)
        )
    return zero, coefficients


def h_q_bessel(
    q: float, delta: float | np.ndarray, eta: float,
    modes: int = 24, layer_modes: int = 48,
) -> float | np.ndarray:
    """Evaluate H_q from its exact Bessel--Lambert representation."""
    zero, coefficients = h_q_bessel_coefficients(q, eta, modes, layer_modes)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def h_q_polylog(
    q: int, delta: float | np.ndarray, eta: float, modes: int = 24,
) -> float | np.ndarray:
    """Evaluate H_6 or H_12 with the exact polylog K-Lambert closure."""
    if q not in (6, 12) or eta <= 0.0 or modes < 1:
        raise ValueError("require q in {6,12}, eta>0 and modes>=1")
    nu = 0.5 * (q - 1.0)
    prefactor = math.sqrt(math.pi) / float(mpmath.gamma(0.5 * q))
    zero = (
        prefactor * float(mpmath.gamma(nu)) * eta ** (1.0 - q)
        * float(mpmath.zeta(q - 1))
    )
    coefficients = np.array([
        4.0 * prefactor * (math.pi * ell / eta) ** nu
        * bessel_lambert_polylog(q, 2.0 * math.pi * ell * eta)
        for ell in range(1, modes + 1)
    ])
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def delta_h_q(
    q: float, delta: float | np.ndarray, delta0: float, eta: float,
    modes: int = 24, layer_modes: int = 48,
) -> float | np.ndarray:
    """Slip-excess Delta H_q; its registry-independent zero mode cancels."""
    _, coefficients = h_q_bessel_coefficients(q, eta, modes, layer_modes)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = np.sum(coefficients * (
        np.cos(2.0 * math.pi * values[..., None] * ell)
        - np.cos(2.0 * math.pi * delta0 * ell)
    ), axis=-1)
    return float(result) if values.ndim == 0 else result


@lru_cache(maxsize=256)
def _h_q_eta_derivative_coefficients(
    q: float, eta: float, modes: int, layer_modes: int,
) -> tuple[float, np.ndarray]:
    """Analytic eta derivative of the Bessel--Lambert H_q coefficients."""
    if not (q > 2.0 and eta > 0.0):
        raise ValueError("require q>2 and eta>0")
    nu = 0.5 * (q - 1.0)
    prefactor = math.sqrt(math.pi) / float(mpmath.gamma(0.5 * q))
    zero = (
        prefactor * float(mpmath.gamma(nu)) * float(mpmath.zeta(q - 1.0))
        * (1.0 - q) * eta ** (-q)
    )
    coefficients = np.empty(modes)
    for index in range(modes):
        ell = index + 1
        total = mpmath.mpf("0")
        for k in range(1, layer_modes + 1):
            x = 2.0 * math.pi * ell * k * eta
            kval = mpmath.besselk(nu, x)
            kprime = -0.5 * (
                mpmath.besselk(nu - 1.0, x) + mpmath.besselk(nu + 1.0, x)
            )
            total += k ** (-nu) * (
                -nu * eta ** (-nu - 1.0) * kval
                + eta ** (-nu) * (2.0 * math.pi * ell * k) * kprime
            )
        coefficients[index] = float(4.0 * prefactor * (math.pi * ell) ** nu * total)
    return zero, coefficients


def h_q_eta_derivative_bessel(
    q: float, delta: float | np.ndarray, eta: float,
    modes: int = 24, layer_modes: int = 48,
) -> float | np.ndarray:
    """Analytic partial H_q / partial eta from the exact series."""
    zero, coefficients = _h_q_eta_derivative_coefficients(q, eta, modes, layer_modes)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def h_q_delta_derivative_bessel(
    q: float, delta: float | np.ndarray, eta: float,
    modes: int = 24, layer_modes: int = 48,
) -> float | np.ndarray:
    """Analytic partial H_q / partial delta."""
    _, coefficients = h_q_bessel_coefficients(q, eta, modes, layer_modes)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = np.sum(
        -2.0 * math.pi * ell * coefficients
        * np.sin(2.0 * math.pi * values[..., None] * ell), axis=-1,
    )
    return float(result) if values.ndim == 0 else result


@lru_cache(maxsize=128)
def registry_energy_coefficients(lattice: RegistryLattice) -> tuple[float, np.ndarray]:
    """Return U0/(C_mn epsilon_LJ) Fourier coefficients."""
    lattice.validate()
    zero_m, coeff_m = h_q_bessel_coefficients(
        lattice.m, lattice.normal_ratio, lattice.bessel_modes, lattice.layer_modes
    )
    zero_n, coeff_n = h_q_bessel_coefficients(
        lattice.n, lattice.normal_ratio, lattice.bessel_modes, lattice.layer_modes
    )
    sm, sn = lattice.sigma_ratio**lattice.m, lattice.sigma_ratio**lattice.n
    return sm * zero_m - sn * zero_n, sm * coeff_m - sn * coeff_n


def registry_energy(
    delta: float | np.ndarray, lattice: RegistryLattice = RegistryLattice(),
) -> float | np.ndarray:
    """Return intrinsic multilayer U0/(C_mn epsilon_LJ)."""
    zero, coefficients = registry_energy_coefficients(lattice)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, coefficients.size + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def registry_energy_derivative(
    delta: float | np.ndarray, lattice: RegistryLattice = RegistryLattice(),
) -> float | np.ndarray:
    """Return partial[U0/(C_mn epsilon_LJ)]/partial delta."""
    _, coefficients = registry_energy_coefficients(lattice)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, coefficients.size + 1, dtype=float)
    result = np.sum(
        -2.0 * math.pi * ell * coefficients
        * np.sin(2.0 * math.pi * values[..., None] * ell), axis=-1,
    )
    return float(result) if values.ndim == 0 else result


def registry_normal_derivative(
    delta: float | np.ndarray, lattice: RegistryLattice = RegistryLattice(),
) -> float | np.ndarray:
    """Return partial[U0/(C_mn epsilon_LJ)]/partial eta."""
    lattice.validate()
    dm = h_q_eta_derivative_bessel(
        lattice.m, delta, lattice.normal_ratio,
        lattice.bessel_modes, lattice.layer_modes,
    )
    dn = h_q_eta_derivative_bessel(
        lattice.n, delta, lattice.normal_ratio,
        lattice.bessel_modes, lattice.layer_modes,
    )
    return lattice.sigma_ratio**lattice.m * dm - lattice.sigma_ratio**lattice.n * dn


def u0(
    a: float, s: float | np.ndarray, params: MultilayerPotentialParameters,
) -> float | np.ndarray:
    """Physical intrinsic U0(a,s), excluding all external work."""
    lattice = params.reduced(a)
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj
    return scale * registry_energy(np.asarray(s, dtype=float) / params.b, lattice)


def u0_direct(
    a: float, s: float, params: MultilayerPotentialParameters,
    kmax: int = 200, pmax: int = 400,
) -> float:
    """Physical direct double-sum U0 with declared numerical cutoffs."""
    params.validate()
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj
    delta, eta = s / params.b, a / params.b
    return scale * (
        (params.sigma_lj / params.b) ** params.m
        * h_q_direct(params.m, delta, eta, kmax, pmax)
        - (params.sigma_lj / params.b) ** params.n
        * h_q_direct(params.n, delta, eta, kmax, pmax)
    )


def v_slip(
    a: float, s: float | np.ndarray, params: MultilayerPotentialParameters,
    s0: float = 0.0,
) -> float | np.ndarray:
    """Exact slip excess U0(a,s)-U0(a,s0) from the common U0."""
    return u0(a, s, params) - u0(a, s0, params)


def dU_ds(
    a: float, s: float | np.ndarray, params: MultilayerPotentialParameters,
) -> float | np.ndarray:
    """Physical intrinsic registry force partial U0/partial s."""
    lattice = params.reduced(a)
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj
    return scale / params.b * registry_energy_derivative(
        np.asarray(s, dtype=float) / params.b, lattice
    )


def dU_da(
    a: float, s: float | np.ndarray, params: MultilayerPotentialParameters,
) -> float | np.ndarray:
    """Physical intrinsic normal force partial U0/partial a."""
    lattice = params.reduced(a)
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj
    return scale / params.b * registry_normal_derivative(
        np.asarray(s, dtype=float) / params.b, lattice
    )


def dU_da_direct(
    a: float, s: float, params: MultilayerPotentialParameters,
    kmax: int = 120, pmax: int = 240,
) -> float:
    """Direct termwise intrinsic normal force for independent verification."""
    params.validate()
    delta, eta = s / params.b, a / params.b
    dm, _ = h_q_eta_derivatives_direct(params.m, delta, eta, kmax, pmax)
    dn, _ = h_q_eta_derivatives_direct(params.n, delta, eta, kmax, pmax)
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj / params.b
    return float(scale * (
        (params.sigma_lj / params.b) ** params.m * dm
        - (params.sigma_lj / params.b) ** params.n * dn
    ))


def d2U_da2_direct(
    a: float, s: float, params: MultilayerPotentialParameters,
    kmax: int = 120, pmax: int = 240,
) -> float:
    """Direct termwise normal curvature used to classify stationary roots."""
    params.validate()
    delta, eta = s / params.b, a / params.b
    _, d2m = h_q_eta_derivatives_direct(params.m, delta, eta, kmax, pmax)
    _, d2n = h_q_eta_derivatives_direct(params.n, delta, eta, kmax, pmax)
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj / params.b**2
    return float(scale * (
        (params.sigma_lj / params.b) ** params.m * d2m
        - (params.sigma_lj / params.b) ** params.n * d2n
    ))


def slip_barrier(
    a: float, params: MultilayerPotentialParameters,
    s0: float = 0.0, samples: int = 4096,
) -> float:
    """Full-energy barrier height in one registry period at fixed a."""
    if samples < 64:
        raise ValueError("samples must be at least 64")
    s = s0 + np.linspace(0.0, params.b, samples, endpoint=False)
    energy = np.asarray(u0(a, s, params))
    return float(np.max(energy) - np.min(energy))


def normal_stationary_points(
    generalized_force: float,
    s: float,
    params: MultilayerPotentialParameters,
    a_min: float,
    a_max: float,
    samples: int = 400,
    kmax: int = 120,
    pmax: int = 240,
) -> list[tuple[float, float]]:
    """Roots of partial_a U0=Q_a classified by numerical curvature.

    The returned tuples are ``(a_root, curvature)``.  Positive curvature is
    the bonded stable root and negative curvature is the outer metastable
    barrier.  No fitted critical distance is used.
    """
    if not (a_min > 0.0 and a_max > a_min and samples >= 32):
        raise ValueError("invalid stationary-point search interval")
    grid = np.linspace(a_min, a_max, samples)
    residual = np.array([
        dU_da_direct(a, s, params, kmax, pmax) - generalized_force for a in grid
    ])
    roots: list[tuple[float, float]] = []
    for left, right, f_left, f_right in zip(grid[:-1], grid[1:], residual[:-1], residual[1:]):
        if f_left == 0.0:
            root = float(left)
        elif f_left * f_right > 0.0:
            continue
        else:
            lo, hi, flo = float(left), float(right), float(f_left)
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                fm = dU_da_direct(mid, s, params, kmax, pmax) - generalized_force
                if flo * fm <= 0.0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            root = 0.5 * (lo + hi)
        if roots and abs(root - roots[-1][0]) < 1.0e-9 * params.b:
            continue
        curvature = d2U_da2_direct(root, s, params, kmax, pmax)
        roots.append((root, curvature))
    return roots


def preferred_registry(lattice: RegistryLattice = RegistryLattice()) -> float:
    """Return the numerically verified minimum delta in one period."""
    delta = np.linspace(0.0, 1.0, 4096, endpoint=False)
    return float(delta[int(np.argmin(registry_energy(delta, lattice)))])


def plastic_well_index(
    s: float | np.ndarray, b: float, s0: float = 0.0,
) -> int | np.ndarray:
    """Integer z in s=s0+z*b+tilde{s}, tilde{s} in [-b/2,b/2)."""
    if b <= 0.0:
        raise ValueError("b must be positive")
    values = np.asarray(s, dtype=float)
    result = np.floor((values - s0) / b + 0.5).astype(int)
    return int(result) if values.ndim == 0 else result


def schmid_factor(
    loading_axis: Iterable[float],
    plane_normal: Iterable[float],
    slip_direction: Iterable[float],
) -> float:
    """Return signed M=(l.n)(l.d) for one declared slip system."""
    load = np.asarray(tuple(loading_axis), dtype=float)
    normal = np.asarray(tuple(plane_normal), dtype=float)
    direction = np.asarray(tuple(slip_direction), dtype=float)
    if load.shape != (3,) or normal.shape != (3,) or direction.shape != (3,):
        raise ValueError("all crystallographic vectors must have three components")
    norms = np.linalg.norm(load), np.linalg.norm(normal), np.linalg.norm(direction)
    if min(norms) <= 0.0:
        raise ValueError("crystallographic vectors must be nonzero")
    load, normal, direction = (
        load / norms[0], normal / norms[1], direction / norms[2]
    )
    if abs(float(np.dot(normal, direction))) > 1.0e-10:
        raise ValueError("slip direction must lie in the slip plane")
    return float(np.dot(load, normal) * np.dot(load, direction))
