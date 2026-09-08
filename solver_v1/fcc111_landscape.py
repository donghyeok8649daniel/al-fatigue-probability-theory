"""Static barrier and stability diagnostics for the full FCC(111) reference."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar, root, root_scalar

from .fcc111_full_energy import FullFCC111StackEnergy
from .fcc111_geometry import registry_path


@dataclass(frozen=True)
class NormalBarrier:
    force: float
    s: float
    minimum_a: float
    saddle_a: float
    barrier: float
    critical_force: float
    critical_a: float
    status: str


@dataclass(frozen=True)
class RegistryBarrier:
    force: float
    minimum_s: float
    saddle_s: float
    minimum_a: float
    saddle_a: float
    barrier: float
    minimum_relaxed_curvature: float
    saddle_relaxed_curvature: float
    status: str


@dataclass(frozen=True)
class RegistrySpinodal:
    force: float
    a: float
    s: float
    relaxed_curvature: float
    residual_norm: float


def loaded_energy(
    model: FullFCC111StackEnergy, a: float, s: float, force: float
) -> float:
    if not np.isfinite(model.a0):
        raise ValueError("loaded energy requires a stable zero-load reference")
    extension = (float(a) - model.a0) + model.p.chi_axial_projection * float(s)
    return float(model.local_energy(float(a), float(s)) - float(force) * extension)


def normal_stationary_points(
    model: FullFCC111StackEnergy,
    s: float,
    force: float,
    *,
    samples: int = 96,
    upper_over_b: float = 8.0,
) -> list[tuple[float, float]]:
    """Return ``(a,W_aa)`` roots of ``W_a=f`` on the declared branch range."""

    lower = max(float(model.p.a_min), 0.20 * model.p.b)
    upper = max(float(model.p.a_max), float(upper_over_b) * model.p.b)
    grid = np.geomspace(lower, upper, int(samples))
    residual = np.array(
        [model.local_deda(float(value), float(s)) - float(force) for value in grid]
    )
    roots: list[tuple[float, float]] = []
    for left, right, fleft, fright in zip(
        grid[:-1], grid[1:], residual[:-1], residual[1:]
    ):
        if not np.isfinite(fleft) or not np.isfinite(fright):
            continue
        if fleft == 0.0:
            root = float(left)
        elif fleft * fright > 0.0:
            continue
        else:
            solved = root_scalar(
                lambda value: model.local_deda(float(value), float(s)) - force,
                bracket=(float(left), float(right)),
                method="brentq",
                xtol=1.0e-12,
                rtol=1.0e-12,
            )
            if not solved.converged:
                continue
            root = float(solved.root)
        if not roots or abs(root - roots[-1][0]) > 1.0e-8:
            roots.append((root, float(model.local_hessian(root, float(s))[0, 0])))
    return roots


def opening_spinodal(
    model: FullFCC111StackEnergy,
    s: float,
    *,
    upper_over_b: float = 8.0,
) -> tuple[float, float]:
    """Fixed-registry ideal normal force maximum ``max_a W_a(a,s)``."""

    lower = max(float(model.p.a_min), 0.20 * model.p.b)
    upper = max(float(model.p.a_max), float(upper_over_b) * model.p.b)
    result = minimize_scalar(
        lambda value: -model.local_deda(float(value), float(s)),
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": 2.0e-10},
    )
    if not result.success:
        raise RuntimeError("opening-spinodal optimization failed")
    critical_a = float(result.x)
    critical_force = float(model.local_deda(critical_a, float(s)))
    return critical_force, critical_a


def opening_barrier(
    model: FullFCC111StackEnergy, s: float, force: float
) -> NormalBarrier | None:
    """Uniform-plane normal barrier; distinct from a single cleavage interface."""

    if not np.isfinite(model.a0):
        return None
    critical_force, critical_a = opening_spinodal(model, float(s))
    if force >= critical_force:
        return NormalBarrier(
            force=float(force), s=float(s), minimum_a=float("nan"),
            saddle_a=float("nan"), barrier=0.0,
            critical_force=critical_force, critical_a=critical_a,
            status="normal_bound_basin_absent",
        )
    roots = normal_stationary_points(model, float(s), float(force))
    minima = [(a, curvature) for a, curvature in roots if curvature > 0.0]
    if not minima:
        return None
    minimum_a = min(minima, key=lambda item: abs(item[0] - model.a0))[0]
    saddles = [a for a, curvature in roots if curvature < 0.0 and a > minimum_a]
    if force <= 1.0e-14:
        barrier = model.isolated_plane_energy() - model.local_energy(minimum_a, float(s))
        saddle_a = float("inf")
        status = "zero_load_uniform_separation_work"
    elif saddles:
        saddle_a = min(saddles)
        barrier = loaded_energy(model, saddle_a, float(s), force) - loaded_energy(
            model, minimum_a, float(s), force
        )
        status = "finite_uniform_opening_saddle"
    else:
        return None
    return NormalBarrier(
        force=float(force), s=float(s), minimum_a=float(minimum_a),
        saddle_a=float(saddle_a), barrier=float(max(0.0, barrier)),
        critical_force=critical_force, critical_a=critical_a, status=status,
    )


def _stable_normal_root(
    model: FullFCC111StackEnergy,
    s: float,
    force: float,
    previous: float | None = None,
) -> float | None:
    # A stable branch is a negative-to-positive zero of W_a-f.  Searching only
    # the bound-basin neighborhood is both safer and much cheaper than finding
    # the distant opening saddle for every registry sample.
    target = model.a0 if previous is None else float(previous)
    lower = max(float(model.p.a_min), 0.20 * model.p.b)
    upper = max(2.0 * model.a0, 1.75 * model.p.b)
    local_left = max(lower, 0.70 * target)
    local_right = min(upper, 1.35 * target)
    pieces = (
        np.linspace(local_left, local_right, 22),
        np.linspace(lower, upper, 42),
    )
    for grid in pieces:
        residual = np.array(
            [model.local_deda(float(value), float(s)) - force for value in grid]
        )
        candidates: list[float] = []
        for left, right, fleft, fright in zip(
            grid[:-1], grid[1:], residual[:-1], residual[1:]
        ):
            if not np.isfinite(fleft) or not np.isfinite(fright):
                continue
            if fleft <= 0.0 <= fright:
                solved = root_scalar(
                    lambda value: model.local_deda(float(value), float(s)) - force,
                    bracket=(float(left), float(right)), method="brentq",
                    xtol=1.0e-11, rtol=1.0e-11,
                )
                if solved.converged:
                    root = float(solved.root)
                    if model.local_hessian(root, float(s))[0, 0] > 0.0:
                        candidates.append(root)
        if candidates:
            return float(min(candidates, key=lambda value: abs(value - target)))
    return None


def relaxed_registry_barrier(
    model: FullFCC111StackEnergy,
    force: float,
    *,
    samples: int = 41,
) -> RegistryBarrier | None:
    """Barrier along the normal-relaxed homogeneous scalar registry path."""

    if not np.isfinite(model.a0):
        return None
    period = registry_path(model.path_id).period_over_b * model.p.b
    coordinates = np.linspace(-0.15 * period, 0.85 * period, int(samples))
    normal = np.full_like(coordinates, np.nan)
    energy = np.full_like(coordinates, np.nan)
    derivative = np.full_like(coordinates, np.nan)
    curvature = np.full_like(coordinates, np.nan)
    previous = None
    for index, coordinate in enumerate(coordinates):
        root = _stable_normal_root(model, float(coordinate), float(force), previous)
        if root is None:
            previous = None
            continue
        previous = root
        hessian = model.local_hessian(root, float(coordinate))
        normal[index] = root
        energy[index] = loaded_energy(model, root, float(coordinate), force)
        derivative[index] = (
            model.local_deds(root, float(coordinate))
            - force * model.p.chi_axial_projection
        )
        curvature[index] = hessian[1, 1] - hessian[0, 1] ** 2 / hessian[0, 0]
    valid = np.flatnonzero(
        np.isfinite(normal) & np.isfinite(energy) & np.isfinite(derivative)
    )
    if valid.size < 9:
        return None
    splits = np.split(valid, np.where(np.diff(valid) != 1)[0] + 1)
    component = min(splits, key=lambda part: np.min(np.abs(coordinates[part])))
    x = coordinates[component]
    d = derivative[component]
    e = energy[component]
    a_values = normal[component]
    if x.size < 9:
        return None
    derivative_spline = CubicSpline(x, d)
    energy_spline = CubicSpline(x, e)
    normal_spline = CubicSpline(x, a_values)
    roots = derivative_spline.roots(extrapolate=False)
    roots = roots[(roots > x[0]) & (roots < x[-1])]
    second = derivative_spline(roots, 1)
    minima = roots[second > 0.0]
    if minima.size == 0:
        return None
    minimum_s = float(minima[np.argmin(np.abs(minima))])
    saddles = roots[(roots > minimum_s) & (second < 0.0)]
    if saddles.size == 0:
        return None
    saddle_s = float(np.min(saddles))
    barrier = float(energy_spline(saddle_s) - energy_spline(minimum_s))
    if barrier <= 0.0:
        return None
    return RegistryBarrier(
        force=float(force), minimum_s=minimum_s, saddle_s=saddle_s,
        minimum_a=float(normal_spline(minimum_s)),
        saddle_a=float(normal_spline(saddle_s)), barrier=barrier,
        minimum_relaxed_curvature=float(derivative_spline(minimum_s, 1)),
        saddle_relaxed_curvature=float(derivative_spline(saddle_s, 1)),
        status="normal_relaxed_homogeneous_registry_barrier",
    )


def configurational_spinodal(
    model: FullFCC111StackEnergy,
    *,
    initial: tuple[float, float, float] | None = None,
) -> RegistrySpinodal | None:
    """Solve branch stationarity plus vanishing relaxed ``s`` curvature.

    This is the static homogeneous-registry spinodal.  It does not by itself
    imply finite-time probability transfer.
    """

    if not np.isfinite(model.a0):
        return None
    critical_force, _ = opening_spinodal(model, 0.0)
    if initial is None:
        initial = (1.12 * model.a0, 0.14 * model.p.b, 0.78 * critical_force)

    def equations(values):
        a, s, force = map(float, values)
        hessian = model.local_hessian(a, s)
        relaxed = hessian[1, 1] - hessian[0, 1] ** 2 / hessian[0, 0]
        return np.array(
            [
                model.local_deda(a, s) - force,
                model.local_deds(a, s) - force * model.p.chi_axial_projection,
                relaxed,
            ]
        )

    solved = root(equations, np.asarray(initial, dtype=float), method="hybr")
    if not solved.success:
        return None
    a, s, force = map(float, solved.x)
    residual = equations(solved.x)
    if (
        a <= 0.0
        or force <= 0.0
        or force >= 1.05 * critical_force
        or np.linalg.norm(residual) > 1.0e-7
    ):
        return None
    hessian = model.local_hessian(a, s)
    relaxed = hessian[1, 1] - hessian[0, 1] ** 2 / hessian[0, 0]
    return RegistrySpinodal(
        force=force,
        a=a,
        s=s,
        relaxed_curvature=float(relaxed),
        residual_norm=float(np.linalg.norm(residual)),
    )
