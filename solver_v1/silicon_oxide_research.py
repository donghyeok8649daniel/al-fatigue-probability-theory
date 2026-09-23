"""Axial composite audit for oxide-coated Si; not fracture/FE predictions."""
from __future__ import annotations

import numpy as np


def axial_composite(*, areas, young_moduli, eigenstrains, nominal_stress):
    """Common axial strain, linear phases, zero prescribed transverse stress.

    U/L = sum A_i E_i (eps-eps_i*)²/2. Force equilibrium gives eps.
    All areas/moduli use consistent units; output stress has the modulus unit.
    This one-dimensional reference omits 3D Poisson, corner and shear fields.
    """
    a, e, e0 = (np.asarray(x, float) for x in (areas, young_moduli, eigenstrains))
    sigma = float(nominal_stress)
    if (a.ndim != 1 or not len(a) or e.shape != a.shape or e0.shape != a.shape
            or not all(np.all(np.isfinite(x)) for x in (a, e, e0))
            or np.any(a < 0) or np.sum(a) <= 0 or np.any(e <= 0)
            or not np.isfinite(sigma)):
        raise ValueError('finite phase arrays, nonnegative areas and positive moduli required')
    total, stiffness = float(np.sum(a)), float(a@e)
    strain = (sigma*total+float((a*e)@e0))/stiffness
    phase_stress = e*(strain-e0)
    force = float(a@phase_stress)
    energy = float(.5*np.sum(a*e*(strain-e0)**2))
    return dict(strain=float(strain), phase_stress=phase_stress,
        effective_modulus=stiffness/total, total_force=force,
        force_residual=force-sigma*total, elastic_energy_per_length=energy)


def rectangular_shell(*, outer_width, outer_height, thickness_per_face):
    w, h, t = map(float, (outer_width, outer_height, thickness_per_face))
    if not np.all(np.isfinite([w, h, t])) or w <= 0 or h <= 0 or t < 0 or 2*t >= min(w,h):
        raise ValueError('positive outer dimensions and a resolved positive core required')
    core = (w-2*t)*(h-2*t)
    return np.array([core, w*h-core])


def inverse_equivalent_flaw_nm(stress_GPa, *, toughness_MPa_sqrt_m=1., geometry_factor=1.):
    """Reconstruction of paper Eq4, NOT an independently observed initial flaw."""
    sigma, k, y = map(float, (stress_GPa, toughness_MPa_sqrt_m, geometry_factor))
    if not np.all(np.isfinite([sigma, k, y])) or min(sigma,k,y) <= 0:
        raise ValueError('positive finite stress/toughness/geometry required')
    return (k/(y*sigma*1e3))**2/np.pi*1e9


def precipitate_radius_nm(time_minutes, *, diffusion_cm2_s=8e-11,
                          concentration_initial_cm3=7e17, concentration_boundary_cm3=3e17,
                          concentration_precipitate_cm3=4.6e22):
    """Paper Eq3 with explicit high-temperature process inputs (1100 C)."""
    t,d,c0,cb,cp = map(float, (time_minutes,diffusion_cm2_s,concentration_initial_cm3,
                             concentration_boundary_cm3,concentration_precipitate_cm3))
    if (not np.all(np.isfinite([t,d,c0,cb,cp])) or t < 0 or d < 0
            or cb < 0 or c0 < cb or cp <= cb):
        raise ValueError('nonnegative time/diffusion and physical concentration ordering required')
    return float(np.sqrt(2*d*(c0-cb)/(cp-cb)*t*60)*1e7)
