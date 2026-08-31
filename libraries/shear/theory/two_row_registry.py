# 이 파일은 두 평행 원자열의 반복단위당 registry 에너지를 수치 검산한다.
# 주요 함수는 well_depth_scale, inverse_power_q2_closed, two_row_cross_energy_per_repeat_direct이다.
# 이 모듈은 비활성 archive이며 과학적 분류와 한계는 docs/SLIP_LATTICE_ENERGY_REVIEW.md를 따른다.

"""Diagnostic evaluators for the archived two-row registry derivation.

The mathematical object is the cross-row interaction energy per upper atom,
equivalently per commensurate row repeat.  It is not the total energy of two
infinite rows and it is not an active aluminum plasticity law.

The direct sums below use an explicit symmetric numerical half-width.  That
half-width is a convergence-control parameter for theorem checks, not a
physical interaction cutoff.  The q=2 helper supplies an independent exact
closed form for checking the counting and periodicity of the shifted sum.
"""

from __future__ import annotations

import math


def _check_exponents(m: float, n: float) -> None:
    if not (m > n > 1.0):
        raise ValueError("require m > n > 1")


def well_depth_scale(m: float, n: float) -> float:
    """Return C_mn mapping pair-well depth to coefficient epsilon.

    If ``epsilon_well`` is the positive pair-well depth and the zero crossing
    is ``sigma``, the repository coefficient convention is recovered with
    ``epsilon_coefficient = C_mn * epsilon_well``.
    """
    _check_exponents(m, n)
    return m / (m - n) * (m / n) ** (n / (m - n))


def pair_equilibrium_ratio(m: float, n: float) -> float:
    """Return r_e/sigma for one generalized-m,n pair."""
    _check_exponents(m, n)
    return (m / n) ** (1.0 / (m - n))


def pair_potential_coefficient(
    r: float,
    epsilon_coefficient: float,
    sigma: float,
    m: float,
    n: float,
) -> float:
    """Generalized pair energy in the active coefficient convention."""
    _check_exponents(m, n)
    if min(r, epsilon_coefficient, sigma) <= 0.0:
        raise ValueError("r, epsilon_coefficient and sigma must be positive")
    return epsilon_coefficient * ((sigma / r) ** m - (sigma / r) ** n)


def inverse_power_direct(
    q: float,
    delta: float,
    eta: float,
    half_width: int,
) -> float:
    """Symmetrically truncate sum_p [((p+delta)^2+eta^2)]^(-q/2).

    This is a diagnostic numerical sum.  Convergence must be demonstrated by
    increasing ``half_width``; it is not an exact finite-range interaction.
    """
    if not (q > 1.0 and eta > 0.0 and half_width >= 1):
        raise ValueError("require q>1, eta>0 and half_width>=1")
    return math.fsum(
        ((p + delta) ** 2 + eta**2) ** (-0.5 * q)
        for p in range(-half_width, half_width + 1)
    )


def inverse_power_q2_closed(delta: float, eta: float) -> float:
    """Exact closed form of the shifted q=2 sum.

    sum_p 1/[(p+delta)^2+eta^2]
      = pi/eta * sinh(2*pi*eta)
        / [cosh(2*pi*eta)-cos(2*pi*delta)].
    """
    if eta <= 0.0:
        raise ValueError("eta must be positive")
    x = 2.0 * math.pi * eta
    return (
        math.pi
        / eta
        * math.sinh(x)
        / (math.cosh(x) - math.cos(2.0 * math.pi * delta))
    )


def two_row_cross_energy_per_repeat_direct(
    a: float,
    s: float,
    b: float,
    epsilon_coefficient: float,
    sigma: float,
    m: float = 12.19,
    n: float = 6.0,
    half_width: int = 4096,
) -> float:
    """Diagnostic direct sum for cross-row energy per upper atom/repeat."""
    _check_exponents(m, n)
    if min(a, b, epsilon_coefficient, sigma) <= 0.0:
        raise ValueError("a, b, epsilon_coefficient and sigma must be positive")
    delta = s / b
    eta = a / b
    return epsilon_coefficient * (
        (sigma / b) ** m * inverse_power_direct(m, delta, eta, half_width)
        - (sigma / b) ** n * inverse_power_direct(n, delta, eta, half_width)
    )


def registry_force_per_repeat_direct(
    a: float,
    s: float,
    b: float,
    epsilon_coefficient: float,
    sigma: float,
    m: float = 12.19,
    n: float = 6.0,
    half_width: int = 4096,
) -> float:
    """Return dW/ds for the same explicit diagnostic direct sum."""
    _check_exponents(m, n)
    if min(a, b, epsilon_coefficient, sigma) <= 0.0:
        raise ValueError("a, b, epsilon_coefficient and sigma must be positive")

    def derivative(q: float) -> float:
        return math.fsum(
            -q
            * sigma**q
            * (p * b + s)
            * (a * a + (p * b + s) ** 2) ** (-0.5 * q - 1.0)
            for p in range(-half_width, half_width + 1)
        )

    return epsilon_coefficient * (derivative(m) - derivative(n))

