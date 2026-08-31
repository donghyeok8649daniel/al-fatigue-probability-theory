# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 이론 계산에 사용하는 Python 모듈이다.
# - 주요 클래스: TransportConfig, TransportHistory
# - 주요 함수/메서드: _validate, domain_max, _grid, conditional_equilibrium, _face_coefficients
#   _absorbing_coefficient, step, _barrier_dimensionless, solve
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Finite-volume Smoluchowski transport with reflecting or absorbing escape.

All equations are dimensionless in lambda=a/a0, energy/E0 and time/t_r.
The density is conditional-normalized only for reflecting boundaries.  With an
absorbing upper boundary it is the unnormalised intact density rho and its
finite-volume mass is the survival probability.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from theory.normal_lj_chain import normalized_lj_energy, normalized_lj_force, critical_stretch
from theory.normal_lj_probability_dynamics import _chang_cooper_weight, _solve_tridiagonal


@dataclass(frozen=True)
class TransportConfig:
    inverse_temperature: float = 80.0
    lambda_min: float = 0.78
    lambda_max: float = 1.34
    cells: int = 220
    boundary: str = "reflecting"
    initiation_definition: str = "fixed_coordinate"
    m: float = 12.19
    n: float = 6.0


@dataclass(frozen=True)
class TransportHistory:
    time: np.ndarray
    force: np.ndarray
    stretch: np.ndarray
    density: np.ndarray
    survival: np.ndarray
    initiation: np.ndarray
    outflux: np.ndarray
    hazard: np.ndarray
    mean: np.ndarray
    variance: np.ndarray
    skewness: np.ndarray
    mean_shifted_energy: np.ndarray
    tail_conditional: np.ndarray
    work: np.ndarray
    entropy_production: np.ndarray


def _validate(c: TransportConfig) -> None:
    if c.boundary not in {"reflecting", "absorbing"}:
        raise ValueError("boundary must be reflecting or absorbing")
    if c.initiation_definition not in {"fixed_coordinate", "tangent_instability"}:
        raise ValueError("unknown initiation definition")
    if c.initiation_definition == "tangent_instability" and c.boundary != "absorbing":
        raise ValueError("tangent-instability first passage requires an absorbing boundary")
    if not (c.m > c.n > 1 and c.inverse_temperature > 0):
        raise ValueError("invalid exponents or inverse temperature")
    if not (0 < c.lambda_min < 1 < c.lambda_max) or c.cells < 40:
        raise ValueError("invalid domain or grid")


def domain_max(c: TransportConfig) -> float:
    """Return the actual upper boundary implied by the declared definition."""
    if c.initiation_definition == "tangent_instability":
        return critical_stretch(c.m, c.n)
    return c.lambda_max


def _grid(c: TransportConfig) -> tuple[np.ndarray, float]:
    upper = domain_max(c)
    dx = (upper - c.lambda_min) / c.cells
    return c.lambda_min + (np.arange(c.cells) + 0.5) * dx, dx


def conditional_equilibrium(x: np.ndarray, dx: float, force: float,
                            c: TransportConfig) -> np.ndarray:
    potential = normalized_lj_energy(x, c.m, c.n) - force * (x - 1.0)
    logp = -c.inverse_temperature * potential
    logp -= np.max(logp)
    p = np.exp(logp)
    return p / (np.sum(p) * dx)


def _face_coefficients(x: np.ndarray, dx: float, force: float,
                       c: TransportConfig) -> tuple[np.ndarray, np.ndarray]:
    phi = normalized_lj_energy(x, c.m, c.n)
    drift = force - np.diff(phi) / dx
    diffusion = 1.0 / c.inverse_temperature
    delta = _chang_cooper_weight(drift * dx / diffusion)
    return (drift * (1.0 - delta) + diffusion / dx,
            drift * delta - diffusion / dx)


def _absorbing_coefficient(x_last: float, dx: float, force: float,
                           c: TransportConfig) -> float:
    """Coefficient J_out=c p_last to a zero-density boundary face."""
    h = 0.5 * dx
    phi_last = float(normalized_lj_energy(x_last, c.m, c.n))
    phi_face = float(normalized_lj_energy(domain_max(c), c.m, c.n))
    drift = force - (phi_face - phi_last) / h
    diffusion = 1.0 / c.inverse_temperature
    delta = float(_chang_cooper_weight(np.asarray([drift * h / diffusion]))[0])
    return drift * (1.0 - delta) + diffusion / h


def step(density: np.ndarray, x: np.ndarray, dx: float, dt: float, force: float,
         c: TransportConfig) -> tuple[np.ndarray, float]:
    """Backward-Euler conservative step; return rho and end-step outflux."""
    left, right = _face_coefficients(x, dx, force, c)
    ratio = dt / dx
    n = x.size
    lower = np.empty(n - 1)
    diagonal = np.ones(n)
    upper = np.empty(n - 1)
    diagonal[0] += ratio * left[0]
    upper[0] = ratio * right[0]
    for i in range(1, n - 1):
        lower[i - 1] = -ratio * left[i - 1]
        diagonal[i] += ratio * (left[i] - right[i - 1])
        upper[i] = ratio * right[i]
    lower[-1] = -ratio * left[-1]
    diagonal[-1] -= ratio * right[-1]
    out_coefficient = 0.0
    if c.boundary == "absorbing":
        out_coefficient = _absorbing_coefficient(float(x[-1]), dx, force, c)
        diagonal[-1] += ratio * out_coefficient
    updated = _solve_tridiagonal(lower, diagonal, upper, density)
    if np.min(updated) < -1e-11:
        raise RuntimeError("negative finite-volume density")
    updated = np.maximum(updated, 0.0)
    return updated, out_coefficient * float(updated[-1])


def _barrier_dimensionless(force: float, m: float, n: float) -> float | None:
    if force <= 0:
        return None
    lc = critical_stretch(m, n)
    fc = float(normalized_lj_force(lc, m, n))
    if force >= fc:
        return lc
    lo, hi = lc, 2 * lc
    while float(normalized_lj_force(hi, m, n)) > force:
        hi *= 2
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if float(normalized_lj_force(mid, m, n)) > force:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve(time: np.ndarray, force: np.ndarray, c: TransportConfig = TransportConfig(),
          max_dt: float = 0.02) -> TransportHistory:
    _validate(c)
    t = np.asarray(time, dtype=float)
    f = np.asarray(force, dtype=float)
    if t.ndim != 1 or f.shape != t.shape or t.size < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("time and force must be matching, increasing arrays")
    x, dx = _grid(c)
    rho = conditional_equilibrium(x, dx, float(f[0]), c)
    densities = np.empty((t.size, c.cells)); densities[0] = rho
    out = np.zeros(t.size)
    for k in range(1, t.size):
        count = max(1, math.ceil((t[k] - t[k - 1]) / max_dt))
        local_dt = (t[k] - t[k - 1]) / count
        flux_integral = 0.0
        for j in range(count):
            fj = f[k - 1] + (j + 0.5) / count * (f[k] - f[k - 1])
            rho, jout = step(rho, x, dx, local_dt, float(fj), c)
            flux_integral += local_dt * jout
        out[k] = flux_integral / (t[k] - t[k - 1])
        densities[k] = rho

    survival = np.sum(densities, axis=1) * dx
    if c.boundary == "reflecting":
        densities /= survival[:, None]
        survival[:] = 1.0
    mean = np.empty(t.size); variance = np.empty(t.size); skew = np.empty(t.size)
    energy = np.empty(t.size); tail = np.empty(t.size); entropy_prod = np.empty(t.size)
    phi = normalized_lj_energy(x, c.m, c.n)
    psi = phi - float(normalized_lj_energy(1.0, c.m, c.n))
    lc = critical_stretch(c.m, c.n)
    for k, row in enumerate(densities):
        mass = survival[k]
        cond = row / mass if mass > 0 else row
        mean[k] = np.sum(x * cond) * dx
        variance[k] = np.sum((x - mean[k]) ** 2 * cond) * dx
        third = np.sum((x - mean[k]) ** 3 * cond) * dx
        skew[k] = third / variance[k] ** 1.5 if variance[k] > 0 else 0.0
        energy[k] = np.sum(psi * cond) * dx
        barrier = _barrier_dimensionless(float(f[k]), c.m, c.n)
        upper = barrier if barrier is not None else c.lambda_max
        tail[k] = np.sum(cond[(x >= lc) & (x <= upper)]) * dx
        left, right = _face_coefficients(x, dx, float(f[k]), c)
        currents = left * row[:-1] + right * row[1:]
        face_p = np.maximum(0.5 * (row[:-1] + row[1:]), 1e-300)
        entropy_prod[k] = np.sum(currents ** 2 / face_p) * dx
    work = np.zeros(t.size)
    work[1:] = np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(mean))
    hazard = np.divide(out, survival, out=np.zeros_like(out), where=survival > 0)
    return TransportHistory(t, f, x, densities, survival, 1.0 - survival, out,
                            hazard, mean, variance, skew, energy, tail, work,
                            entropy_prod)
