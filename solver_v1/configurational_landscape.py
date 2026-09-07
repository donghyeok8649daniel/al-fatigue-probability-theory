"""Bound-basin configurational free-energy branch diagnostics.

The fixed normal-domain potential of mean force can be contaminated by the
tilted large-``a`` tail under tensile load.  This module therefore audits the
metastable configurational branch using the truncated bound basin below the
actual opening saddle.  Truncated Gibbs is a diagnostic, not the absorbing
fast-process QSD and not a replacement for the full 2D PDE.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from .model import TwoRowLJ
from .reduced_fast_a import bound_basin_conditional


@dataclass(frozen=True)
class ConfigurationalBarrier:
    force: float
    minimum_s: float
    saddle_s: float
    barrier: float
    minimum_curvature: float
    saddle_curvature: float
    opening_barrier_at_minimum: float
    opening_barrier_at_saddle: float
    s: np.ndarray
    F_eff: np.ndarray
    dF_eff_ds: np.ndarray
    d2F_eff_ds2: np.ndarray


def bound_configurational_barrier(
    model: TwoRowLJ,
    force: float,
    *,
    s_lower: float = -0.10,
    s_upper: float = 0.50,
    n_s: int = 121,
    quadrature_order: int = 48,
) -> ConfigurationalBarrier | None:
    """Track the principal metastable minimum and its forward ``s`` saddle.

    ``None`` means that a distinct positive-curvature minimum and subsequent
    negative-curvature saddle cannot both be resolved before the bound normal
    basin terminates.  It does not, by itself, prove finite-time well transfer.
    """

    if n_s < 9:
        raise ValueError("n_s must be at least 9")
    s = np.linspace(float(s_lower), float(s_upper), int(n_s))
    F = np.full(s.shape, np.nan)
    derivative = np.full(s.shape, np.nan)
    for index, coordinate in enumerate(s):
        try:
            conditional = bound_basin_conditional(
                model,
                float(coordinate),
                float(force),
                quadrature_order=int(quadrature_order),
            )
        except ValueError:
            continue
        F[index] = conditional.F_eff
        derivative[index] = conditional.dF_eff_ds

    valid_indices = np.flatnonzero(np.isfinite(F) & np.isfinite(derivative))
    if valid_indices.size < 9:
        return None
    # Use the contiguous valid component containing s=0.  Opening instability
    # can truncate other parts of the nominal registry period.
    splits = np.split(valid_indices, np.where(np.diff(valid_indices) != 1)[0] + 1)
    component = min(splits, key=lambda part: np.min(np.abs(s[part])))
    sv = s[component]
    Fv = F[component]
    dv = derivative[component]
    if sv.size < 9:
        return None
    d_spline = CubicSpline(sv, dv)
    F_spline = CubicSpline(sv, Fv)
    roots = np.asarray(d_spline.roots(extrapolate=False), dtype=float)
    roots = roots[(roots > sv[0]) & (roots < sv[-1])]
    curvature = d_spline(roots, 1)
    minima = roots[curvature > 0.0]
    if minima.size == 0:
        return None
    minimum_s = float(minima[np.argmin(np.abs(minima))])
    forward = roots[(roots > minimum_s) & (d_spline(roots, 1) < 0.0)]
    if forward.size == 0:
        return None
    saddle_s = float(np.min(forward))
    minimum_curvature = float(d_spline(minimum_s, 1))
    saddle_curvature = float(d_spline(saddle_s, 1))
    barrier = float(F_spline(saddle_s) - F_spline(minimum_s))
    if barrier <= 0.0:
        return None
    second = np.asarray(d_spline(sv, 1), dtype=float)
    return ConfigurationalBarrier(
        force=float(force),
        minimum_s=minimum_s,
        saddle_s=saddle_s,
        barrier=barrier,
        minimum_curvature=minimum_curvature,
        saddle_curvature=saddle_curvature,
        opening_barrier_at_minimum=float(
            model.opening_barrier(minimum_s, float(force))
        ),
        opening_barrier_at_saddle=float(
            model.opening_barrier(saddle_s, float(force))
        ),
        s=sv,
        F_eff=Fv,
        dF_eff_ds=dv,
        d2F_eff_ds2=second,
    )
