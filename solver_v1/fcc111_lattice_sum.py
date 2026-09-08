"""Two-dimensional triangular-lattice Poisson sums for FCC(111) planes.

The reciprocal formulas use the Fourier convention

    f_hat(G) = integral_R2 f(x) exp(-i G.x) dx

and ``sum_R f(R+delta) = A_cell^-1 sum_G f_hat(G) exp(i G.delta)``.
Direct real-space routines are independent validation references only.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.special import gamma, kv, kvp, zeta

from .fcc111_geometry import FCC111Geometry, fcc111_geometry_from_b


@dataclass(frozen=True)
class ReciprocalSumConfig:
    tol: float = 1.0e-13
    max_shell_index: int = 24
    consecutive_small_shells: int = 3

    def validate(self) -> None:
        if not np.isfinite(self.tol) or self.tol <= 0.0:
            raise ValueError("tol must be finite and positive")
        if self.max_shell_index < 2:
            raise ValueError("max_shell_index must be at least 2")
        if self.consecutive_small_shells < 1:
            raise ValueError("consecutive_small_shells must be positive")


@dataclass(frozen=True)
class ReciprocalShell:
    quadratic_index: int
    magnitude: float
    vectors: np.ndarray

    @property
    def degeneracy(self) -> int:
        return int(self.vectors.shape[0])


@dataclass(frozen=True)
class PlaneKernelResult:
    value: float
    d_d: float
    grad_delta: np.ndarray
    d2_dd: float
    mixed_d_delta: np.ndarray
    hess_delta: np.ndarray
    shells_used: int
    reciprocal_vectors_used: int
    estimated_tail_absolute: float


@dataclass(frozen=True)
class DirectPlaneResult:
    value: float
    d_d: float
    grad_delta: np.ndarray
    d2_dd: float
    mixed_d_delta: np.ndarray
    hess_delta: np.ndarray
    points_used: int
    continuum_tail_estimate: float


def reciprocal_quadratic_index(h: int, k: int) -> int:
    """Triangular reciprocal metric index ``h^2-h*k+k^2``."""

    return int(h) ** 2 - int(h) * int(k) + int(k) ** 2


@lru_cache(maxsize=32)
def _reciprocal_shells_cached(b: float, max_shell_index: int) -> tuple[ReciprocalShell, ...]:
    geometry = fcc111_geometry_from_b(float(b))
    q_limit = int(max_shell_index) ** 2
    span = 2 * int(max_shell_index)
    grouped: dict[int, list[np.ndarray]] = {}
    for h in range(-span, span + 1):
        for k in range(-span, span + 1):
            q = reciprocal_quadratic_index(h, k)
            if q == 0 or q > q_limit:
                continue
            grouped.setdefault(q, []).append(geometry.reciprocal_vector(h, k))
    shells = []
    for q in sorted(grouped):
        vectors = np.asarray(grouped[q], dtype=float)
        magnitudes = np.linalg.norm(vectors, axis=1)
        if np.ptp(magnitudes) > 5.0e-13 * max(1.0, float(magnitudes[0])):
            raise FloatingPointError("reciprocal shell contains unequal magnitudes")
        shells.append(
            ReciprocalShell(
                quadratic_index=q,
                magnitude=float(np.mean(magnitudes)),
                vectors=vectors,
            )
        )
    return tuple(shells)


def triangular_reciprocal_shells(
    geometry: FCC111Geometry,
    max_shell_index: int,
) -> tuple[ReciprocalShell, ...]:
    """Return complete equal-|G| shells through the declared metric radius."""

    return _reciprocal_shells_cached(float(geometry.b), int(max_shell_index))


def triangular_epstein_zeta(p: float, b: float = 1.0) -> float:
    r"""Return ``sum_(m,n)!=(0,0) |m*a1+n*a2|^(-2p)`` exactly.

    For the Eisenstein triangular form,
    ``Z(p)=6*zeta(p)*L_{-3}(p)`` with
    ``L_{-3}=3^-p [zeta(p,1/3)-zeta(p,2/3)]``.
    """

    exponent = float(p)
    spacing = float(b)
    if exponent <= 1.0:
        raise ValueError("triangular Epstein sum requires p > 1")
    if spacing <= 0.0:
        raise ValueError("b must be positive")
    l_minus_three = 3.0 ** (-exponent) * (
        zeta(exponent, 1.0 / 3.0) - zeta(exponent, 2.0 / 3.0)
    )
    return float(6.0 * zeta(exponent, 1.0) * l_minus_three / spacing ** (2.0 * exponent))


def power_kernel_fourier_2d(d: float, g: float, p: float) -> float:
    r"""Fourier transform of ``(d^2+r^2)^(-p)`` in two dimensions."""

    distance = float(d)
    wave = float(g)
    exponent = float(p)
    if distance <= 0.0 or exponent <= 1.0 or wave < 0.0:
        raise ValueError("require d>0, p>1, and g>=0")
    if wave == 0.0:
        return float(math.pi / (exponent - 1.0) * distance ** (2.0 - 2.0 * exponent))
    return float(
        2.0
        * math.pi
        / gamma(exponent)
        * (wave / (2.0 * distance)) ** (exponent - 1.0)
        * kv(exponent - 1.0, distance * wave)
    )


def exponential_kernel_fourier_2d(d: float, g: float, kappa: float) -> float:
    r"""Fourier transform of ``exp(-kappa*sqrt(d^2+r^2))`` in 2D."""

    distance = float(d)
    wave = float(g)
    decay = float(kappa)
    if distance < 0.0 or wave < 0.0 or decay <= 0.0:
        raise ValueError("require d>=0, g>=0, and kappa>0")
    q = math.sqrt(decay * decay + wave * wave)
    return float(2.0 * math.pi * decay * math.exp(-distance * q) * (1.0 + distance * q) / q**3)


def _phase_sums(vectors: np.ndarray, delta: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    phase = vectors @ delta
    cosine = np.cos(phase)
    sine = np.sin(phase)
    return (
        float(np.sum(cosine)),
        -np.sum(sine[:, None] * vectors, axis=0),
        -np.einsum("i,ij,ik->jk", cosine, vectors, vectors),
    )


def _power_coefficient_derivatives(
    d: float, g: float, p: float, area: float
) -> tuple[float, float, float]:
    nu = float(p) - 1.0
    prefactor = 2.0 * math.pi / (float(area) * gamma(float(p))) * (g / 2.0) ** nu
    kval = kv(nu, d * g)
    first_k = kvp(nu, d * g, 1)
    second_k = kvp(nu, d * g, 2)
    coefficient = prefactor * d ** (-nu) * kval
    first = prefactor * (
        -nu * d ** (-nu - 1.0) * kval
        + d ** (-nu) * g * first_k
    )
    second = prefactor * (
        nu * (nu + 1.0) * d ** (-nu - 2.0) * kval
        - 2.0 * nu * d ** (-nu - 1.0) * g * first_k
        + d ** (-nu) * g * g * second_k
    )
    return float(coefficient), float(first), float(second)


def plane_power_sum_reciprocal(
    d: float,
    delta,
    *,
    p: float,
    geometry: FCC111Geometry,
    config: ReciprocalSumConfig = ReciprocalSumConfig(),
) -> PlaneKernelResult:
    """Poisson-summed triangular plane power kernel and exact derivatives."""

    config.validate()
    distance = float(d)
    exponent = float(p)
    shift = np.asarray(delta, dtype=float)
    if distance <= 0.0 or exponent <= 1.0 or shift.shape != (2,):
        raise ValueError("require d>0, p>1, and a two-component delta")
    area = float(geometry.atomic_cell_area)
    alpha = 2.0 - 2.0 * exponent
    value = math.pi / (area * (exponent - 1.0)) * distance**alpha
    d_d = alpha * value / distance
    d2_dd = alpha * (alpha - 1.0) * value / distance**2
    grad = np.zeros(2)
    mixed = np.zeros(2)
    hess = np.zeros((2, 2))
    small = 0
    shells_used = 0
    vectors_used = 0
    recent_bounds: list[float] = []
    for shell in triangular_reciprocal_shells(geometry, config.max_shell_index):
        coefficient, first, second = _power_coefficient_derivatives(
            distance, shell.magnitude, exponent, area
        )
        cosine_sum, phase_gradient, phase_hessian = _phase_sums(shell.vectors, shift)
        value_term = coefficient * cosine_sum
        d_term = first * cosine_sum
        dd_term = second * cosine_sum
        grad_term = coefficient * phase_gradient
        mixed_term = first * phase_gradient
        hess_term = coefficient * phase_hessian
        value += value_term
        d_d += d_term
        d2_dd += dd_term
        grad += grad_term
        mixed += mixed_term
        hess += hess_term
        shells_used += 1
        vectors_used += shell.degeneracy

        g = shell.magnitude
        bound = shell.degeneracy * max(
            abs(coefficient), abs(first), abs(second),
            abs(coefficient) * g, abs(first) * g, abs(coefficient) * g * g,
        )
        recent_bounds.append(float(bound))
        recent_bounds = recent_bounds[-config.consecutive_small_shells :]
        scale = max(
            1.0, abs(value), abs(d_d), abs(d2_dd),
            float(np.max(np.abs(grad))), float(np.max(np.abs(mixed))),
            float(np.max(np.abs(hess))),
        )
        if bound <= config.tol * scale:
            small += 1
            if small >= config.consecutive_small_shells:
                break
        else:
            small = 0
    else:
        raise RuntimeError("reciprocal power sum did not meet tolerance")

    return PlaneKernelResult(
        value=float(value),
        d_d=float(d_d),
        grad_delta=grad,
        d2_dd=float(d2_dd),
        mixed_d_delta=mixed,
        hess_delta=hess,
        shells_used=shells_used,
        reciprocal_vectors_used=vectors_used,
        estimated_tail_absolute=float(max(recent_bounds, default=0.0)),
    )


def _direct_lattice_points(geometry: FCC111Geometry, radius_index: int) -> np.ndarray:
    radius = int(radius_index)
    if radius < 2:
        raise ValueError("radius_index must be at least 2")
    span = 2 * radius
    points = []
    limit = radius * radius
    for m in range(-span, span + 1):
        for n in range(-span, span + 1):
            if m * m + m * n + n * n <= limit:
                points.append(geometry.lattice_vector(m, n))
    return np.asarray(points, dtype=float)


def plane_power_sum_direct(
    d: float,
    delta,
    *,
    p: float,
    geometry: FCC111Geometry,
    radius_index: int = 128,
) -> DirectPlaneResult:
    """Independent complete-disk real-space reference with analytic derivatives."""

    distance = float(d)
    exponent = float(p)
    shift = np.asarray(delta, dtype=float)
    if distance <= 0.0 or exponent <= 1.0 or shift.shape != (2,):
        raise ValueError("require d>0, p>1, and a two-component delta")
    lateral = _direct_lattice_points(geometry, radius_index) + shift
    radius2 = np.sum(lateral * lateral, axis=1)
    squared = distance * distance + radius2
    base = squared ** (-exponent)
    next_power = squared ** (-exponent - 1.0)
    next_two = squared ** (-exponent - 2.0)
    value = float(np.sum(base))
    d_d = float(np.sum(-2.0 * exponent * distance * next_power))
    grad = np.sum(-2.0 * exponent * next_power[:, None] * lateral, axis=0)
    d2_dd = float(
        np.sum(
            -2.0 * exponent * next_power
            + 4.0 * exponent * (exponent + 1.0) * distance**2 * next_two
        )
    )
    mixed = np.sum(
        4.0 * exponent * (exponent + 1.0)
        * distance * next_two[:, None] * lateral,
        axis=0,
    )
    hess = (
        -2.0 * exponent * np.sum(next_power) * np.eye(2)
        + 4.0 * exponent * (exponent + 1.0)
        * np.einsum("i,ij,ik->jk", next_two, lateral, lateral)
    )
    radial_cut = max(0.5 * geometry.b, radius_index * geometry.b - np.linalg.norm(shift))
    tail = math.pi / (
        geometry.atomic_cell_area * (exponent - 1.0)
    ) * (distance * distance + radial_cut * radial_cut) ** (1.0 - exponent)
    return DirectPlaneResult(
        value=value,
        d_d=d_d,
        grad_delta=grad,
        d2_dd=d2_dd,
        mixed_d_delta=mixed,
        hess_delta=hess,
        points_used=int(lateral.shape[0]),
        continuum_tail_estimate=float(tail),
    )


def _exponential_coefficient_derivatives(
    d: float, g: float, kappa: float, area: float
) -> tuple[float, float, float]:
    q = math.sqrt(kappa * kappa + g * g)
    exponential = math.exp(-d * q)
    common = 2.0 * math.pi * kappa / area
    coefficient = common * exponential * (1.0 + d * q) / q**3
    first = -common * d * exponential / q
    second = common * exponential * (d * q - 1.0) / q
    return float(coefficient), float(first), float(second)


def plane_exponential_sum_reciprocal(
    d: float,
    delta,
    *,
    kappa: float,
    amplitude: float,
    geometry: FCC111Geometry,
    config: ReciprocalSumConfig = ReciprocalSumConfig(),
) -> PlaneKernelResult:
    """Poisson-summed ``amplitude*exp(-kappa*r)`` plane density."""

    config.validate()
    distance = float(d)
    decay = float(kappa)
    scale_amplitude = float(amplitude)
    shift = np.asarray(delta, dtype=float)
    if distance <= 0.0 or decay <= 0.0 or scale_amplitude < 0.0 or shift.shape != (2,):
        raise ValueError("require d>0, kappa>0, amplitude>=0, and 2D delta")
    c0, c0d, c0dd = _exponential_coefficient_derivatives(
        distance, 0.0, decay, geometry.atomic_cell_area
    )
    value = scale_amplitude * c0
    d_d = scale_amplitude * c0d
    d2_dd = scale_amplitude * c0dd
    grad = np.zeros(2)
    mixed = np.zeros(2)
    hess = np.zeros((2, 2))
    small = 0
    shells_used = 0
    vectors_used = 0
    recent_bounds: list[float] = []
    for shell in triangular_reciprocal_shells(geometry, config.max_shell_index):
        coefficient, first, second = _exponential_coefficient_derivatives(
            distance, shell.magnitude, decay, geometry.atomic_cell_area
        )
        coefficient *= scale_amplitude
        first *= scale_amplitude
        second *= scale_amplitude
        cosine_sum, phase_gradient, phase_hessian = _phase_sums(shell.vectors, shift)
        value += coefficient * cosine_sum
        d_d += first * cosine_sum
        d2_dd += second * cosine_sum
        grad += coefficient * phase_gradient
        mixed += first * phase_gradient
        hess += coefficient * phase_hessian
        shells_used += 1
        vectors_used += shell.degeneracy
        g = shell.magnitude
        bound = shell.degeneracy * max(
            abs(coefficient), abs(first), abs(second),
            abs(coefficient) * g, abs(first) * g, abs(coefficient) * g * g,
        )
        recent_bounds.append(float(bound))
        recent_bounds = recent_bounds[-config.consecutive_small_shells :]
        actual_scale = max(
            1.0, abs(value), abs(d_d), abs(d2_dd),
            float(np.max(np.abs(grad))), float(np.max(np.abs(mixed))),
            float(np.max(np.abs(hess))),
        )
        if bound <= config.tol * actual_scale:
            small += 1
            if small >= config.consecutive_small_shells:
                break
        else:
            small = 0
    else:
        raise RuntimeError("reciprocal exponential sum did not meet tolerance")
    return PlaneKernelResult(
        value=float(value), d_d=float(d_d), grad_delta=grad,
        d2_dd=float(d2_dd), mixed_d_delta=mixed, hess_delta=hess,
        shells_used=shells_used, reciprocal_vectors_used=vectors_used,
        estimated_tail_absolute=float(max(recent_bounds, default=0.0)),
    )


def plane_exponential_sum_direct(
    d: float,
    delta,
    *,
    kappa: float,
    amplitude: float,
    geometry: FCC111Geometry,
    radius_index: int = 32,
) -> DirectPlaneResult:
    """Independent real-space exponential plane sum and derivatives."""

    distance = float(d)
    decay = float(kappa)
    prefactor = float(amplitude)
    shift = np.asarray(delta, dtype=float)
    if distance <= 0.0 or decay <= 0.0 or prefactor < 0.0 or shift.shape != (2,):
        raise ValueError("require d>0, kappa>0, amplitude>=0, and 2D delta")
    lateral = _direct_lattice_points(geometry, radius_index) + shift
    squared = distance**2 + np.sum(lateral * lateral, axis=1)
    radius = np.sqrt(squared)
    values = prefactor * np.exp(-decay * radius)
    radial_first_over_r = -decay * values / radius
    curvature_factor = values * (decay**2 + decay / radius) / squared
    value = float(np.sum(values))
    d_d = float(np.sum(radial_first_over_r * distance))
    grad = np.sum(radial_first_over_r[:, None] * lateral, axis=0)
    d2_dd = float(
        np.sum(curvature_factor * distance**2 + radial_first_over_r)
    )
    mixed = np.sum(
        curvature_factor[:, None] * distance * lateral, axis=0
    )
    hess = (
        np.einsum("i,ij,ik->jk", curvature_factor, lateral, lateral)
        + np.sum(radial_first_over_r) * np.eye(2)
    )
    radial_cut = max(0.5 * geometry.b, radius_index * geometry.b - np.linalg.norm(shift))
    # Continuum density tail: 2*pi/A * int_R^inf r exp(-k*r) dr.
    tail = (
        2.0 * math.pi * prefactor / geometry.atomic_cell_area
        * math.exp(-decay * radial_cut)
        * (radial_cut / decay + 1.0 / decay**2)
    )
    return DirectPlaneResult(
        value=value, d_d=d_d, grad_delta=grad, d2_dd=d2_dd,
        mixed_d_delta=mixed, hess_delta=hess,
        points_used=int(lateral.shape[0]), continuum_tail_estimate=float(tail),
    )


def same_plane_exponential_density_direct(
    *,
    kappa: float,
    amplitude: float,
    geometry: FCC111Geometry,
    radius_index: int = 32,
) -> tuple[float, int, float]:
    """Self-excluded same-plane exponential density validation sum."""

    points = _direct_lattice_points(geometry, radius_index)
    radii = np.linalg.norm(points, axis=1)
    nonself = radii > 0.0
    value = float(np.sum(float(amplitude) * np.exp(-float(kappa) * radii[nonself])))
    radial_cut = radius_index * geometry.b
    tail = (
        2.0 * math.pi * float(amplitude) / geometry.atomic_cell_area
        * math.exp(-float(kappa) * radial_cut)
        * (radial_cut / float(kappa) + 1.0 / float(kappa) ** 2)
    )
    return value, int(np.sum(nonself)), float(tail)


def same_plane_power_sum_direct(
    *,
    p: float,
    geometry: FCC111Geometry,
    radius_index: int = 128,
) -> tuple[float, int, float]:
    """Self-excluded complete-disk reference for the in-plane power sum."""

    exponent = float(p)
    if exponent <= 1.0:
        raise ValueError("same-plane power sum requires p > 1")
    points = _direct_lattice_points(geometry, radius_index)
    radius2 = np.sum(points * points, axis=1)
    nonself = radius2 > 0.0
    value = float(np.sum(radius2[nonself] ** (-exponent)))
    radial_cut = radius_index * geometry.b
    tail = math.pi / (
        geometry.atomic_cell_area * (exponent - 1.0)
    ) * radial_cut ** (2.0 - 2.0 * exponent)
    return value, int(np.sum(nonself)), float(tail)
