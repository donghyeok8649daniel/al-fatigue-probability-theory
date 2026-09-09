"""Declared finite-part partitions; do not fit an atomic core radius.

Smooth radial quadrature reduces lattice-shell jumps. Its continuum constant
is explicitly subtracted. This cannot cure an unstable state, an incorrect
potential, or free-domain/transverse-tail error.
"""
import numpy as np


def radial_window(t,inner_fraction=.5):
    if not np.isfinite(inner_fraction) or not 0<inner_fraction<1:
        raise ValueError('numerical window fraction must lie strictly between zero and one')
    t=np.asarray(t,float)
    if np.any(~np.isfinite(t)) or np.any(t<0):
        raise ValueError('finite nonnegative normalized radii required')
    weights=np.zeros_like(t); weights[t<=inner_fraction]=1.
    middle=(t>inner_fraction)&(t<1.)
    x=(t[middle]-inner_fraction)/(1-inner_fraction)
    weights[middle]=1-3*x*x+2*x*x*x
    return weights


def log_window_constant(inner_fraction=.5,*,tolerance=2e-15):
    """Integral_0^1 (w(t)-1)/t dt, from an absolutely convergent series.

    With delta=1-alpha, c_w=-6 delta sum_n delta^n/[(n+1)(n+3)(n+4)].
    Consecutive absolute term ratios are smaller than delta, giving a
    geometric bound on the omitted tail. This avoids near-alpha=1 cancellation.
    """
    if not np.isfinite(inner_fraction) or not 0<inner_fraction<1:
        raise ValueError('window fraction in (0,1) required')
    if not np.isfinite(tolerance) or tolerance<=0:
        raise ValueError('positive finite quadrature tolerance required')
    delta=1-inner_fraction; total=0.; power=delta
    for n in range(100000):
        term=-6*power/((n+1)*(n+3)*(n+4)); total+=term
        power*=delta
        next_term=6*power/((n+2)*(n+4)*(n+5))
        tail=next_term/inner_fraction
        if tail<=tolerance:
            return dict(value=total,estimated_series_tail=tail,terms_used=n+1)
    raise RuntimeError('log-window series failed the specified tolerance')


def finite_part_samples(distances,remainder_energy,radii,*,log_coefficient,b,
                        free_radius,inner_fractions=(.5,.75)):
    distance=np.asarray(distances,float); energy=np.asarray(remainder_energy,float)
    if (distance.ndim!=1 or energy.shape!=distance.shape or np.any(distance<0)
            or not np.all(np.isfinite(distance)) or not np.all(np.isfinite(energy))):
        raise ValueError('matching finite site distances and unaltered remainder energies required')
    if not np.all(np.isfinite([log_coefficient,b,free_radius])) or min(log_coefficient,b,free_radius)<=0:
        raise ValueError('positive same-potential logarithmic coefficient and geometric scales required')
    rows=[]
    for radius in radii:
        if not np.isfinite(radius) or radius<=0 or radius>free_radius:
            raise ValueError('integration radius must lie within the actual free domain')
        for alpha in inner_fractions:
            if alpha*radius<b:
                continue
            constant=log_window_constant(alpha)
            weighted=float(radial_window(distance/radius,alpha)@energy)
            subtraction=log_coefficient*(np.log(radius/b)+constant['value'])
            rows.append(dict(radius_over_L0=float(radius),inner_fraction=float(alpha),
                weighted_remainder_eV=weighted,continuum_subtraction_eV=float(subtraction),
                finite_part_at_reference_b_eV=weighted-subtraction,
                window_constant=constant['value'],window_series_tail=constant['estimated_series_tail'],
                window_series_terms=constant['terms_used'],log_coefficient_eV=log_coefficient,
                core_radius_fitted=False,convergence_certified=False))
    return rows
