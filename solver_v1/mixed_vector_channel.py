"""Coherent radial vector moment, analytic static research (not production).

Energy per site: lambda*x**z*|Q1(k0)+t Q1(k1)|², lambda>=0.
Moments, not energies, mix first. Signed t is a radial angular weight, NOT a
negative electron density. The scalar environment stays positive/unchanged.
"""
import numpy as np
from .tail_constrained_material import normalized_amplitude


def mix_site_moments(first, second, weight):
    a, b = np.asarray(first, float), np.asarray(second, float)
    if (a.ndim != 3 or b.ndim != 3 or a.shape[1:] != (10, 3) or b.shape[1:] != (10, 3)
            or not len(a) or not len(b) or np.any(~np.isfinite(a))
            or np.any(~np.isfinite(b)) or not np.isfinite(weight)):
        raise ValueError('finite three-component per-site jets and radial weight required')
    result = np.zeros((max(len(a), len(b)), 10, 3))
    result[:len(a)] += a
    result[:len(b)] += weight*b
    return result


def mixed_vector_bulk_column(first, second, q_cubic, weight):
    """Harmonic cubic column with analytic omitted-neighbor norm envelope."""
    if (not np.isfinite(weight) or first.stretch != 1 or second.stretch != 1
            or not np.array_equal(first.R, second.R)):
        raise ValueError('same pristine cubic neighbor geometry required')
    q = first.wavevector(q_cubic)
    phase_weight = -2*np.sin(first.R@q/2)**2
    variations, bounds = [], []
    for model in (first, second):
        idx = next(i for i, (rank, _) in enumerate(model.moments) if rank == 1)
        variation = np.einsum('n,nij->ij', phase_weight, model.gradients[idx])
        variations.append(variation)
        k, tail = model.rank1_decay, model.tail
        def envelope(power):
            return normalized_amplitude(k)*(tail.exponential(k, power)+k*tail.exponential(k, 1+power))
        bounds.append(min(2*envelope(0), .5*np.linalg.norm(q)**2*envelope(2)))
    L = variations[0]+weight*variations[1]
    dL = bounds[0]+abs(weight)*bounds[1]
    return 2*L@L.T, float(2*(2*np.linalg.norm(L, 2)*dL+dL*dL))
