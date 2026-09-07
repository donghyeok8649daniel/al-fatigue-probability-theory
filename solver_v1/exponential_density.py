"""Poisson-summed exponential environment density for two infinite rows.

The canonical evaluation is the reciprocal-space Bessel series.  Direct
real-space summation is provided only as an independent validation utility.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.special import kv


@dataclass(frozen=True)
class ExponentialDensityParams:
    """Analytic EAM-style density-kernel parameters.

    ``rho0`` and ``rho_ref`` share the chosen electron-density unit. ``r_e``,
    ``b``, ``a``, and ``s`` share the model length unit. ``beta_rho`` is
    dimensionless, while ``kappa=beta_rho/r_e`` has inverse-length units.
    """

    rho0: float = 1.0
    beta_rho: float = 3.0
    r_e: float = 1.0
    tol: float = 1.0e-13
    max_modes: int = 128
    consecutive_small: int = 3

    def validate(self) -> None:
        if not np.isfinite(self.rho0) or self.rho0 < 0.0:
            raise ValueError("rho0 must be finite and nonnegative")
        if not np.isfinite(self.beta_rho) or self.beta_rho <= 0.0:
            raise ValueError("beta_rho must be finite and positive")
        if not np.isfinite(self.r_e) or self.r_e <= 0.0:
            raise ValueError("r_e must be finite and positive")
        if not np.isfinite(self.tol) or self.tol <= 0.0:
            raise ValueError("tol must be finite and positive")
        if self.max_modes < 1 or self.consecutive_small < 1:
            raise ValueError("reciprocal-mode controls must be positive")

    @property
    def kappa(self) -> float:
        return float(self.beta_rho / self.r_e)

    @property
    def C_rho(self) -> float:
        return float(self.rho0 * math.exp(self.beta_rho))


@dataclass(frozen=True)
class ExponentialDensityResult:
    value: np.ndarray
    d_da: np.ndarray
    d_ds: np.ndarray
    d2_daa: np.ndarray
    d2_das: np.ndarray
    d2_dss: np.ndarray
    modes_used: int
    estimated_tail_absolute: float


def exponential_fourier_transform(a, wave_number, kappa: float):
    r"""Exact Fourier transform of ``exp(-kappa*sqrt(a^2+x^2))``.

    With the convention ``g_hat(k)=integral g(x) exp(-i*k*x) dx``, the result
    is ``2*a*kappa*K_1(a*q)/q`` where ``q=sqrt(kappa^2+k^2)``.
    """

    aa, kk = np.broadcast_arrays(
        np.asarray(a, dtype=float), np.asarray(wave_number, dtype=float)
    )
    decay = float(kappa)
    if np.any(aa <= 0.0) or not np.isfinite(decay) or decay <= 0.0:
        raise ValueError("a and kappa must be positive")
    q = np.sqrt(decay * decay + kk * kk)
    return 2.0 * aa * decay * kv(1, aa * q) / q


def same_row_exponential_density(
    *, b: float, params: ExponentialDensityParams
) -> float:
    """Exact density from all same-row sites except the reference atom."""

    params.validate()
    spacing = float(b)
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("b must be finite and positive")
    if params.rho0 == 0.0:
        return 0.0
    return float(2.0 * params.C_rho / np.expm1(params.kappa * spacing))


def _transform_a_derivatives(
    a: np.ndarray, q: float, kappa: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return transform H_q and its first two ``a`` derivatives."""

    z = q * a
    h = 2.0 * a * kappa * kv(1, z) / q
    # d[a*K1(q*a)]/da = -q*a*K0(q*a)
    h_a = -2.0 * kappa * a * kv(0, z)
    h_aa = 2.0 * kappa * (a * q * kv(1, z) - kv(0, z))
    return h, h_a, h_aa


def exponential_cross_density(
    a,
    s,
    *,
    b: float,
    params: ExponentialDensityParams = ExponentialDensityParams(),
) -> ExponentialDensityResult:
    r"""Evaluate the infinite cross-row density and analytic derivatives.

    For ``x_n=(n+1/2)b-s`` and ``q_m=sqrt(kappa^2+(2*pi*m/b)^2)``, Poisson
    summation gives

    ``rho_cross = C/b H_0 + 2*C/b sum_{m>=1} H_m cos(theta_m)``,

    where ``H_m=2*a*kappa*K1(a*q_m)/q_m`` and
    ``theta_m=2*pi*m*(1/2-s/b)``.
    """

    params.validate()
    spacing = float(b)
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("b must be finite and positive")
    aa, ss = np.broadcast_arrays(
        np.asarray(a, dtype=float), np.asarray(s, dtype=float)
    )
    if np.any(~np.isfinite(aa)) or np.any(aa <= 0.0) or np.any(~np.isfinite(ss)):
        raise ValueError("a must be positive and a,s must be finite")
    if params.rho0 == 0.0:
        zero = np.zeros_like(aa)
        return ExponentialDensityResult(
            zero, zero.copy(), zero.copy(), zero.copy(), zero.copy(),
            zero.copy(), 0, 0.0,
        )

    decay = params.kappa
    amplitude = params.C_rho
    h0, h0_a, h0_aa = _transform_a_derivatives(aa, decay, decay)
    value = amplitude * h0 / spacing
    d_da = amplitude * h0_a / spacing
    d_ds = np.zeros_like(value)
    d2_daa = amplitude * h0_aa / spacing
    d2_das = np.zeros_like(value)
    d2_dss = np.zeros_like(value)

    last_measure: float | None = None
    small_count = 0
    modes_used = params.max_modes
    tail_estimate = math.inf
    for mode in range(1, params.max_modes + 1):
        wave = 2.0 * math.pi * mode / spacing
        q = math.sqrt(decay * decay + wave * wave)
        h, h_a, h_aa = _transform_a_derivatives(aa, q, decay)
        phase = wave * (0.5 * spacing - ss)
        cosine = np.cos(phase)
        sine = np.sin(phase)
        prefactor = 2.0 * amplitude / spacing

        term = prefactor * h * cosine
        term_a = prefactor * h_a * cosine
        term_s = prefactor * h * wave * sine
        term_aa = prefactor * h_aa * cosine
        term_as = prefactor * h_a * wave * sine
        term_ss = -prefactor * h * wave * wave * cosine
        value += term
        d_da += term_a
        d_ds += term_s
        d2_daa += term_aa
        d2_das += term_as
        d2_dss += term_ss

        measure = max(
            float(np.max(np.abs(term))),
            float(np.max(np.abs(term_a))),
            float(np.max(np.abs(term_s))),
            float(np.max(np.abs(term_aa))),
            float(np.max(np.abs(term_as))),
            float(np.max(np.abs(term_ss))),
        )
        scale = max(
            1.0,
            float(np.max(np.abs(value))),
            float(np.max(np.abs(d_da))),
            float(np.max(np.abs(d_ds))),
            float(np.max(np.abs(d2_daa))),
            float(np.max(np.abs(d2_das))),
            float(np.max(np.abs(d2_dss))),
        )
        if last_measure is not None and last_measure > 0.0:
            ratio = min(0.999999, measure / last_measure)
            tail_estimate = measure * ratio / max(1.0 - ratio, 1.0e-15)
        else:
            tail_estimate = measure
        if measure <= params.tol * scale and tail_estimate <= params.tol * scale:
            small_count += 1
            if small_count >= params.consecutive_small:
                modes_used = mode
                break
        else:
            small_count = 0
        last_measure = measure

    return ExponentialDensityResult(
        value=np.asarray(value),
        d_da=np.asarray(d_da),
        d_ds=np.asarray(d_ds),
        d2_daa=np.asarray(d2_daa),
        d2_das=np.asarray(d2_das),
        d2_dss=np.asarray(d2_dss),
        modes_used=int(modes_used),
        estimated_tail_absolute=float(tail_estimate),
    )


def direct_exponential_cross_density(
    a,
    s,
    *,
    b: float,
    params: ExponentialDensityParams = ExponentialDensityParams(),
    images: int = 1000,
):
    """Large finite real-space sum for validation, never canonical use."""

    params.validate()
    if images < 1:
        raise ValueError("images must be positive")
    aa, ss = np.broadcast_arrays(
        np.asarray(a, dtype=float), np.asarray(s, dtype=float)
    )
    result = np.zeros_like(aa)
    for index in range(-int(images), int(images) + 1):
        x = (index + 0.5) * float(b) - ss
        radius = np.sqrt(aa * aa + x * x)
        result += params.C_rho * np.exp(-params.kappa * radius)
    return result
