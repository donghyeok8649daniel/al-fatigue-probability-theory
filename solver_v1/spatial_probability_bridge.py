"""Research-only surface quadrature bridge; NOT a spatial mechanics solver.

Patch index is not an independent-region count. Correlation area enters only
survival aggregation. No energy, mobility, or production PDE is modified.
"""
from dataclasses import dataclass
import numpy as np

from .interface_static_scenarios import resolved_tensor_tractions


def patch_tractions_mpa(stress_mpa, normals, slip_directions):
    """Return (patch, normal/shear1/shear2); never discard second shear."""
    stress = np.asarray(stress_mpa, dtype=float)
    normals = np.asarray(normals, dtype=float)
    slips = np.asarray(slip_directions, dtype=float)
    if stress.ndim != 3 or stress.shape[1:] != (3, 3) or len(stress) == 0:
        raise ValueError('stress must have shape (patch,3,3)')
    if normals.shape != (len(stress), 3) or slips.shape != normals.shape:
        raise ValueError('one explicit orientation per patch required')
    return np.array([resolved_tensor_tractions(s, n, m)
                     for s, n, m in zip(stress, normals, slips)])


@dataclass(frozen=True)
class SurfaceAggregation:
    mathematical_extrapolation: np.ndarray
    numerically_resolved_extrapolation: np.ndarray
    all_active_patches_resolved: np.ndarray
    effective_region_count: float


def aggregate_surface_history(probability, patch_area_mm2, *,
                              correlation_area_mm2, numerical_floor,
                              convergence_certified):
    """Inputs (time,patch); independent-equivalent-region hypothesis only.

    Numerical certification is NOT physical/material/independence calibration.
    Every positive-area patch must exceed its floor and have a supplied actual
    convergence certificate. Missing/zero local signals conservatively withhold
    numerical certification, but the mathematical extrapolation remains visible.
    """
    p = np.asarray(probability, dtype=float)
    area = np.asarray(patch_area_mm2, dtype=float)
    ac = float(correlation_area_mm2)
    if p.ndim != 2 or 0 in p.shape or area.shape != (p.shape[1],):
        raise ValueError('nonempty (time,patch) probability and patch areas required')
    if not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError('probability must be finite in [0,1]')
    if np.any(np.diff(p, axis=0) < 0):
        raise ValueError('cumulative opening probability must be nondecreasing')
    if not np.all(np.isfinite(area)) or np.any(area < 0) or not np.any(area > 0):
        raise ValueError('finite nonnegative areas with positive total required')
    if not np.isfinite(ac) or ac <= 0:
        raise ValueError('positive finite correlation area required')
    floor = np.broadcast_to(np.asarray(numerical_floor, dtype=float), p.shape)
    cert = np.asarray(convergence_certified)
    if cert.dtype != np.dtype(bool):
        raise ValueError('certification must be boolean, not numerical signal size')
    cert = np.broadcast_to(cert, p.shape)
    if not np.all(np.isfinite(floor)) or np.any(floor < 0):
        raise ValueError('finite nonnegative numerical floors required')
    active = area > 0
    weights = area[active] / ac
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0) or not np.isfinite(weights.sum()):
        raise ValueError('effective region count overflow or underflow')
    with np.errstate(divide='ignore'):
        log_survival = np.log1p(-p[:, active]) @ weights
    extrapolation = -np.expm1(log_survival)
    resolved = np.all((p[:, active] > floor[:, active]) & cert[:, active], axis=1)
    return SurfaceAggregation(extrapolation, np.where(resolved, extrapolation, np.nan),
                              resolved, float(weights.sum()))
