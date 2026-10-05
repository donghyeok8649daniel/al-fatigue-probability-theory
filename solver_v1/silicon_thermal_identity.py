"""Boundary-safe classical canonical temperature diagnostic for fixed grips.

For pi(r) proportional to exp(-U/kBT) in a finite Cartesian domain,
<g.grad U> = kBT <div g>, if g has zero normal flow at ALL domain boundaries.
Pair-distance walls require a taper as well as the box-face factors. Forces
must be derivatives of the TOTAL energy in the free Cartesian coordinates.
This identity tests sampling/normalization within a declared model; it cannot
certify real silicon, quantum statistics, Markov closure, or a physical clock.
"""
from __future__ import annotations

import numpy as np

from solver_v1.silicon_conditional_research import KB_EV_K


def thermal_test_field(positions, free, domain_center, reference_center, *,
                       halfwidth_A, minimum_pair_A, taper_width_A):
    """Return g (Angstrom) and its dimensionless free-coordinate divergence.

    g_i=(r_i-reference_i)*(1-((r_i-domain_center_i)/R)^2)*chi(r).
    chi is a product of C1 smoothstep pair factors: zero at the pair wall and
    one beyond wall+taper. No displacement or negative curvature is clipped.
    The computational domain is not a physical intact/crack basin.
    """
    r, center, ref = (np.asarray(x, float) for x in
                      (positions, domain_center, reference_center))
    mask = np.asarray(free)
    radius, wall, width = map(float, (halfwidth_A, minimum_pair_A, taper_width_A))
    if (r.ndim != 2 or r.shape[1] != 3 or center.shape != r.shape
            or ref.shape != r.shape or mask.shape != (len(r),)
            or mask.dtype != np.bool_ or not mask.any()
            or not np.isfinite(r).all() or not np.isfinite(center).all()
            or not np.isfinite(ref).all() or not np.isfinite([radius, wall, width]).all()
            or min(radius, wall, width) <= 0
            or not np.array_equal(center[~mask], r[~mask])):
        raise ValueError('finite fixed-grip Cartesian geometry and positive domain controls required')
    x = r-center
    if np.any(abs(x[mask]) >= radius):
        raise ValueError('position outside declared Cartesian box')
    i, j = np.triu_indices(len(r), 1)
    delta = r[i]-r[j]
    distances = np.linalg.norm(delta, axis=1)
    if not len(distances) or np.any(distances < wall):
        raise ValueError('position outside declared pair-distance domain')
    active = distances < wall+width
    chi, grad_chi = 1., np.zeros_like(r)
    if active.any():
        ai, aj = i[active], j[active]
        ad, av = distances[active], delta[active]
        z = (ad-wall)/width
        factors = z*z*(3-2*z)
        chi = float(np.prod(factors))
        if np.any(factors == 0):
            # smoothstep and its first derivative both vanish at the wall.
            chi = 0.
        else:
            coefficients = chi*(6*z*(1-z)/width)/factors
            contributions = coefficients[:, None]*av/ad[:, None]
            np.add.at(grad_chi, ai, contributions)
            np.add.at(grad_chi, aj, -contributions)
    y = r-ref
    box_factor = 1-(x[mask]/radius)**2
    base = y[mask]*box_factor
    diagonal = box_factor-2*y[mask]*x[mask]/radius**2
    field = np.zeros_like(r)
    field[mask] = chi*base
    divergence_components = chi*diagonal+grad_chi[mask]*base
    divergence = float(divergence_components.sum())
    return dict(field_A=field, divergence=divergence, pair_taper=chi,
                divergence_components=divergence_components,
                active_pair_tapers=int(active.sum()))


def canonical_temperature_terms(positions, forces, free, domain_center,
                                reference_center, *, temperature_K,
                                halfwidth_A, minimum_pair_A, taper_width_A):
    """Return numerator/denominator/residual; average BEFORE taking the ratio.

    A nonzero sample residual is not a temperature change or physical heat.
    Finite-chain bias, autocorrelation, bad forces and domain errors can cause it.
    The iid standard error is invalid for correlated chain draws.
    """
    r, f = np.asarray(positions, float), np.asarray(forces, float)
    t = float(temperature_K)
    if f.shape != r.shape or not np.isfinite(f).all() or not np.isfinite(t) or t <= 0:
        raise ValueError('finite total-energy forces and positive temperature required')
    result = thermal_test_field(r, free, domain_center, reference_center,
                               halfwidth_A=halfwidth_A, minimum_pair_A=minimum_pair_A,
                               taper_width_A=taper_width_A)
    numerator = float(-np.sum(result['field_A']*f))
    denominator = result['divergence']
    result.update(numerator_eV=numerator, denominator=denominator,
                  numerator_components_eV=-result['field_A'][np.asarray(free)]*f[np.asarray(free)],
                  residual_eV=numerator-KB_EV_K*t*denominator,
                  thermal_energy_eV=KB_EV_K*t)
    return result
