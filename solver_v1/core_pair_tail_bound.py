"""Conservative INFINITE transverse-row LJ force remainder, no energy change.

Atomic rows remain infinite. For the screw exterior used here Uy=Uz=0, all
nonzero transverse displacements lie in the represented free disk. This
bound covers the pair force only, not nonlinear environmental terms, domain
relaxation, or physical material accuracy. It is a diagnostic, not a repair.
"""
import numpy as np
from scipy.special import gammaln,zeta

from .current_core_coefficient_basis import power_row_value_gradient


def pair_force_tail_bound(*,ring,b,h,u,v,maximum_transverse_displacement):
    """Per-free-row norm bound for omitted max(|j|,|l|)>ring offsets.

    For phi=c*r^-m, ||Hess phi||<=|c|m(m+1)r^(-m-2).
    Unimodal row quadrature is bounded by 2*max(f)+integral(f)/b at ANY phase.
    There are 8*k offsets on the integer square shell k, and transverse
    distance is at least alpha*k-2*max|U_perp|, alpha=sigma_min(row map).
    Exact Hurwitz-zeta sums complete all shells; no finite cutoff is fitted.
    """
    numbers=[b,h,u,v,maximum_transverse_displacement]
    if (int(ring)!=ring or ring<1 or np.any(~np.isfinite(numbers))
            or min(b,h)<=0 or min(u,v,maximum_transverse_displacement)<0):
        raise ValueError('integer ring and finite physical geometry/pair/displacement bounds required')
    d=np.sqrt(3)*b/2
    alpha=float(np.linalg.svd([[d,d/3],[0.,h]],compute_uv=False).min())
    transverse=2*maximum_transverse_displacement
    shift=transverse/alpha;start=ring+1-shift
    if start<=0:raise ValueError('the omitted neighborhood does not enclose transverse displacements')
    # For each pair the infinite row is b-periodic in x; choose a displacement
    # representative |delta_x|<=b/2 before the mean-value bound.
    displacement=np.hypot(b/2,transverse)
    def shell_power(exponent):
        return 8*alpha**(-exponent)*(zeta(exponent-1,start)+shift*zeta(exponent,start))
    terms=[]
    for coefficient,power in ((u,12),(v,6)):
        nu=(power+2)/2
        integral=np.sqrt(np.pi)*np.exp(gammaln(nu-.5)-gammaln(nu))/b
        bound=displacement*coefficient*power*(power+1)*(2*shell_power(power+2)+integral*shell_power(power+1))
        terms.append(float(bound))
    return dict(pair_gradient_norm_bound_eV_L0=sum(terms),repulsive_bound=terms[0],
        attractive_bound=terms[1],row_map_minimum_singular_value=alpha,
        maximum_pair_transverse_displacement=transverse,
        pair_phase_representative_displacement_bound=displacement,
        exact_infinite_shell_majorant=True,environmental_force_tail_certified=False,
        infinite_free_domain_certified=False,force_modified=False)


def explicit_pair_gradient(core,field,*,u,v):
    """Independent one-center pair gradient, with no embedding/site assembler."""
    full=core.full_field(field)
    lookup={tuple(i):k for k,i in enumerate(core.indices)}
    destination=np.array([[lookup[tuple(i+offset)] for offset in core.offsets]
                          for i in core.indices[core.free_ids]])
    vectors=core.reference[None]+full[destination]-full[core.free_ids,None]
    flat=vectors.reshape(-1,3)
    gradient=np.zeros_like(flat)
    for coefficient,power in ((u,6),(-v,3)):
        _,term=power_row_value_gradient(flat,power,core.rows.b)
        gradient-=coefficient*term
    return gradient.reshape(len(field),len(core.offsets),3).sum(axis=1)
