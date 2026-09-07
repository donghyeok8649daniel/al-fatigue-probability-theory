"""Exact registry-well identities and optional crystallographic geometry.

This module adds no empirical plasticity law.  It only evaluates identities
implied by the existing configurational coordinate ``s`` and its probability
flux.  Geometry helpers are diagnostic and do not alter default solver
parameters.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def registry_decomposition(s, b: float) -> tuple[np.ndarray, np.ndarray]:
    """Return ``n, xi`` with ``s=b*n+xi`` and ``xi in [-b/2,b/2)``."""

    spacing = float(b)
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("b must be finite and positive")
    values = np.asarray(s, dtype=float)
    n = np.floor((values + 0.5 * spacing) / spacing).astype(np.int64)
    xi = values - spacing * n
    # Protect the declared half-open interval against one-ulp arithmetic at
    # translated boundaries without changing the well convention.
    upper = xi >= 0.5 * spacing
    lower = xi < -0.5 * spacing
    if np.any(upper):
        n = n.copy()
        xi = xi.copy()
        n[upper] += 1
        xi[upper] -= spacing
    if np.any(lower):
        n = n.copy()
        xi = xi.copy()
        n[lower] -= 1
        xi[lower] += spacing
    return n, xi


@dataclass(frozen=True)
class RegistryRateTerms:
    unnormalized_registry_moment: float
    conditional_mean_well_index: float
    net_registry_flow_rate: float
    gross_configurational_activity_rate: float
    opening_registry_loss_rate: float
    conditional_net_registry_flow_rate: float
    conditional_selective_opening_rate: float


def registry_rate_terms(
    *,
    well_indices,
    well_populations,
    interwell_net_flux,
    interwell_gross_flux,
    opening_absorption_rate_by_well,
) -> RegistryRateTerms:
    """Evaluate exact unnormalized and conditional registry-moment terms.

    Positive interwell flux points toward increasing ``s``.  With reflecting
    outer ``s`` boundaries,

        d sum(n P_n)/dt = sum(J_{n+1/2}) - sum(n A_dot_n).

    The conditional mean additionally contains the normalization contribution
    ``mean_n * sum(A_dot_n) / S``.
    """

    wells = np.asarray(well_indices, dtype=float).reshape(-1)
    populations = np.asarray(well_populations, dtype=float).reshape(-1)
    net = np.asarray(interwell_net_flux, dtype=float).reshape(-1)
    gross = np.asarray(interwell_gross_flux, dtype=float).reshape(-1)
    absorbed = np.asarray(opening_absorption_rate_by_well, dtype=float).reshape(-1)
    if wells.size != populations.size or wells.size != absorbed.size:
        raise ValueError("well arrays must have identical lengths")
    if net.size != max(0, wells.size - 1) or gross.size != net.size:
        raise ValueError("one interwell flux is required per adjacent well boundary")
    survival = float(np.sum(populations))
    moment = float(wells @ populations)
    net_rate = float(np.sum(net))
    gross_rate = float(np.sum(gross))
    opening_registry_loss = float(wells @ absorbed)
    if survival > 0.0:
        mean = moment / survival
        conditional_flow = net_rate / survival
        conditional_opening = (
            -opening_registry_loss + mean * float(np.sum(absorbed))
        ) / survival
    else:
        mean = np.nan
        conditional_flow = np.nan
        conditional_opening = np.nan
    return RegistryRateTerms(
        unnormalized_registry_moment=moment,
        conditional_mean_well_index=float(mean),
        net_registry_flow_rate=net_rate,
        gross_configurational_activity_rate=gross_rate,
        opening_registry_loss_rate=opening_registry_loss,
        conditional_net_registry_flow_rate=float(conditional_flow),
        conditional_selective_opening_rate=float(conditional_opening),
    )


@dataclass(frozen=True)
class SlipProjection:
    axial_shear_projection: float
    chi_candidate: float


def crystallographic_slip_projection(
    loading_direction,
    slip_direction,
    slip_plane_normal,
    *,
    a0: float,
    slip_plane_spacing: float,
    tolerance: float = 1.0e-10,
) -> SlipProjection:
    """Evaluate the candidate single-slip axial projection for supplied axes.

    Input vectors must already be unit vectors and the slip direction must lie
    in the slip plane.  The returned candidate is

        chi = (a0/h) (e dot m) (e dot n_slip).

    It is not installed into ``ModelParams`` by this function.
    """

    vectors = [
        np.asarray(value, dtype=float).reshape(-1)
        for value in (loading_direction, slip_direction, slip_plane_normal)
    ]
    if any(vector.shape != (3,) or np.any(~np.isfinite(vector)) for vector in vectors):
        raise ValueError("crystallographic directions must be finite 3-vectors")
    e, m, normal = vectors
    if any(abs(float(np.linalg.norm(vector)) - 1.0) > tolerance for vector in vectors):
        raise ValueError("crystallographic directions must be normalized")
    if abs(float(m @ normal)) > tolerance:
        raise ValueError("slip direction must lie in the slip plane")
    a0_value = float(a0)
    spacing = float(slip_plane_spacing)
    if not np.isfinite(a0_value) or a0_value <= 0.0:
        raise ValueError("a0 must be finite and positive")
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("slip_plane_spacing must be finite and positive")
    projection = float((e @ m) * (e @ normal))
    return SlipProjection(
        axial_shear_projection=projection,
        chi_candidate=float((a0_value / spacing) * projection),
    )
