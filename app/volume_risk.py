"""Explicit, uncalibrated volume aggregation of local opening probability.

V_c is an EXTERNAL statistical correlation volume, not FEM cell volume,
atomic volume, or a conversion of A_c. This is a mathematical independent
equivalent-region hypothesis. It never changes mechanics or the local PDE.
"""
import numpy as np


def volume_risk(local_probability, cell_volumes_mm3, correlation_volume_mm3):
    """Return cumulative hazard/mm^3, cell P, and standardized 1 mm^3 P.

    h = -log(1-p_local)/V_c; P(V) = 1-exp(-h*V). Using the survival hazard
    instead of P(V)/V preserves additivity and subdivision invariance even
    for non-rare probabilities. p=1 has infinite hazard and P(V)=1.
    """
    p, volume = np.broadcast_arrays(np.asarray(local_probability, dtype=float),
                                   np.asarray(cell_volumes_mm3, dtype=float))
    vc = float(correlation_volume_mm3)
    if not np.isfinite(vc) or vc <= 0:
        raise ValueError('risk.need_volume')
    if (not p.size or not np.isfinite(p).all() or np.any((p < 0) | (p > 1))
            or not np.isfinite(volume).all() or np.any(volume <= 0)):
        raise ValueError('risk.invalid')
    with np.errstate(divide='ignore', over='ignore', under='ignore', invalid='raise'):
        log_h = np.log(-np.log1p(-p))-np.log(vc)
        density = np.exp(log_h)
        cell_hazard = np.exp(log_h+np.log(volume))
        cell_probability = -np.expm1(-cell_hazard)
        unit_probability = -np.expm1(-density)  # reference volume exactly 1 mm^3
    return density, cell_probability, unit_probability
