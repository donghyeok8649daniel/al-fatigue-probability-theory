# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 이론 계산에 사용하는 Python 모듈이다.
# - 주요 클래스: RegistryLattice
# - 주요 함수/메서드: RegistryLattice.validate, inverse_power_bessel_coefficients, shifted_inverse_power_bessel
#   registry_energy_coefficients, registry_energy, registry_energy_derivative, preferred_registry
#   schmid_factor
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Exact one-registry lattice energy for the active ideal-slip extension.

The geometry is one-dimensional: each row is indexed by one integer and the
only slip coordinate is the scalar registry ``s``.  If normal opening ``a`` is
also evolved, the probability state space is two-dimensional even though the
lattice and slip mechanism remain reduced one-dimensional objects.

This module deliberately does *not* add the collinear-chain energy U_infinity
to the two-row cross energy W.  They describe different pair geometries.  W is
the complete configuration-dependent interaction energy per commensurate row
repeat for this ideal-interface branch; same-row terms are constants when the
row repeat b is fixed.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import mpmath
import numpy as np


@dataclass(frozen=True)
class RegistryLattice:
    """Parameters of the generalized-m,n two-row registry landscape."""

    normal_ratio: float = 1.0  # eta=a/b
    sigma_ratio: float = 1.0  # sigma_LJ/b
    m: float = 12.19
    n: float = 6.0
    bessel_modes: int = 24

    def validate(self) -> None:
        if not (self.m > self.n > 1.0):
            raise ValueError("require m > n > 1")
        if self.normal_ratio <= 0.0 or self.sigma_ratio <= 0.0:
            raise ValueError("normal_ratio and sigma_ratio must be positive")
        if self.bessel_modes < 1:
            raise ValueError("bessel_modes must be positive")


def inverse_power_bessel_coefficients(
    q: float, eta: float, modes: int
) -> tuple[float, np.ndarray]:
    """Return the exact Fourier--Bessel coefficients truncated numerically.

    The returned values represent

        Z_q(delta,eta) = zero + sum_l coeff[l-1] cos(2*pi*l*delta).

    The identity itself is exact.  ``modes`` controls exponentially convergent
    reciprocal-space evaluation and is not a physical interaction cutoff.
    """

    if not (q > 1.0 and eta > 0.0 and modes >= 1):
        raise ValueError("require q>1, eta>0 and modes>=1")
    nu = 0.5 * q
    prefactor = math.sqrt(math.pi) / float(mpmath.gamma(nu))
    zero = (
        prefactor
        * float(mpmath.gamma(nu - 0.5))
        * eta ** (1.0 - q)
    )
    order = 0.5 * (q - 1.0)
    coefficients = np.empty(modes, dtype=float)
    for index in range(modes):
        ell = index + 1
        coefficients[index] = (
            4.0
            * prefactor
            * (math.pi * ell / eta) ** order
            * float(mpmath.besselk(order, 2.0 * math.pi * ell * eta))
        )
    return zero, coefficients


def shifted_inverse_power_bessel(
    q: float,
    delta: float | np.ndarray,
    eta: float,
    modes: int = 24,
) -> float | np.ndarray:
    """Evaluate the shifted Epstein--Hurwitz sum by its Bessel-K identity."""

    zero, coefficients = inverse_power_bessel_coefficients(q, eta, modes)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, modes + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def registry_energy_coefficients(
    lattice: RegistryLattice,
) -> tuple[float, np.ndarray]:
    """Return W/epsilon_c Fourier coefficients for one row repeat."""

    lattice.validate()
    zero_m, coeff_m = inverse_power_bessel_coefficients(
        lattice.m, lattice.normal_ratio, lattice.bessel_modes
    )
    zero_n, coeff_n = inverse_power_bessel_coefficients(
        lattice.n, lattice.normal_ratio, lattice.bessel_modes
    )
    sm = lattice.sigma_ratio**lattice.m
    sn = lattice.sigma_ratio**lattice.n
    return sm * zero_m - sn * zero_n, sm * coeff_m - sn * coeff_n


def registry_energy(
    delta: float | np.ndarray, lattice: RegistryLattice = RegistryLattice()
) -> float | np.ndarray:
    """Return W(a,s)/epsilon_c per row repeat, with delta=s/b."""

    zero, coefficients = registry_energy_coefficients(lattice)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, coefficients.size + 1, dtype=float)
    result = zero + np.sum(
        coefficients * np.cos(2.0 * math.pi * values[..., None] * ell), axis=-1
    )
    return float(result) if values.ndim == 0 else result


def registry_energy_derivative(
    delta: float | np.ndarray, lattice: RegistryLattice = RegistryLattice()
) -> float | np.ndarray:
    """Return d(W/epsilon_c)/d(delta) from the same Bessel series."""

    _, coefficients = registry_energy_coefficients(lattice)
    values = np.asarray(delta, dtype=float)
    ell = np.arange(1, coefficients.size + 1, dtype=float)
    result = np.sum(
        -2.0
        * math.pi
        * ell
        * coefficients
        * np.sin(2.0 * math.pi * values[..., None] * ell),
        axis=-1,
    )
    return float(result) if values.ndim == 0 else result


def preferred_registry(lattice: RegistryLattice = RegistryLattice()) -> float:
    """Return the global minimum in one period to grid-verification accuracy."""

    # The central-force two-row energy is even and one-periodic, but retaining
    # a dense search avoids silently assuming whether 0 or 1/2 is preferred.
    # An even count samples both symmetry candidates delta=0 and delta=1/2.
    delta = np.linspace(0.0, 1.0, 4096, endpoint=False)
    return float(delta[int(np.argmin(registry_energy(delta, lattice)))])


def schmid_factor(
    loading_axis: Iterable[float],
    plane_normal: Iterable[float],
    slip_direction: Iterable[float],
) -> float:
    """Return signed resolved-shear projection (l.n)(l.d).

    All vectors are supplied in one common crystal frame.  The plane normal
    and slip direction must be orthogonal.  The signed value is retained so
    that forward and reverse slip are distinguishable.
    """

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
