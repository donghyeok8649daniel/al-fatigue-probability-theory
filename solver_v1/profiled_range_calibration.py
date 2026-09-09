"""Exact active-face coefficient profiling for the UNCHANGED v12 family.

This removes inner-optimizer noise before a deterministic radial least-squares
continuation. It changes neither targets, discrepancy scales nor energy terms.
Mathematical constrained feasibility does not certify a physical material.
"""
from itertools import combinations

import numpy as np
from scipy.linalg import null_space

from .vector_material_calibration import COEFFICIENTS


def fit_profiled_coefficients(matrix,observations):
    matrix=np.asarray(matrix,float)
    if matrix.shape!=(len(observations),8) or not np.all(np.isfinite(matrix)):
        raise ValueError('finite matrix for the existing eight coefficients required')
    if [i for i,o in enumerate(observations) if o.role=='exact']!=[0,1]:
        raise ValueError('only the existing reference force/cohesion are exact')
    target=np.array([o.target for o in observations]); scale=np.array([o.scale for o in observations])
    equal=matrix[:2]/scale[:2,None]; rhs=target[:2]/scale[:2]
    if np.linalg.matrix_rank(equal)!=2:
        raise ValueError('two independent reference constraints required')
    chosen=[i for i,o in enumerate(observations) if o.role=='fit']
    constrained=(0,1,2,4,5)  # u,v,A,C,D3 >=0; B,D1,D2 unchanged signed sector
    candidates=[]
    for count in range(len(constrained)+1):
        for active in combinations(constrained,count):
            free=[i for i in range(8) if i not in active]
            E=equal[:,free]
            if np.linalg.matrix_rank(E)!=2:
                continue
            offset=np.linalg.lstsq(E,rhs,rcond=None)[0]; tangent=null_space(E)
            design=matrix[np.ix_(chosen,free)]@tangent/scale[chosen,None]
            right=(target[chosen]-matrix[np.ix_(chosen,free)]@offset)/scale[chosen]
            z=np.linalg.lstsq(design,right,rcond=None)[0]
            c=np.zeros(8); c[free]=offset+tangent@z
            # Active coordinates are exactly zero by parameterization. No
            # clipping is applied to the evaluated potential or its derivatives.
            floor=np.zeros(8)
            floor[free]=64*np.finfo(float).eps*np.linalg.cond(E)*(
                abs(offset)+abs(tangent)@abs(z))
            residual=(matrix@c-target)/scale
            if np.any(c[list(constrained)] < -floor[list(constrained)]) or np.max(abs(residual[:2]))>1e-8:
                continue
            candidates.append((float(residual[chosen]@residual[chosen]),c,floor,active,free))
    if not candidates:
        raise ArithmeticError('no feasible active-face profile; do not substitute a fake residual')
    loss,c,floor,active,free=min(candidates,key=lambda row:row[0])
    residual=(matrix@c-target)/scale
    gradient=2*(matrix[chosen]/scale[chosen,None]).T@residual[chosen]
    multiplier=np.linalg.lstsq(equal[:,free].T,-gradient[free],rcond=None)[0]
    normal=gradient+equal.T@multiplier
    stationarity=float(np.max(abs(normal[free])))
    bound_multiplier_min=float(min(normal[list(active)])) if active else None
    return dict(coefficients=c,coefficient_order=COEFFICIENTS,predictions=matrix@c,
        residuals=residual,squared_loss=loss,admissible=True,
        strictly_positive_LJ_resolved=bool(np.all(c[:2]>floor[:2])),
        selected_rows=chosen,active_zero_coefficients=[COEFFICIENTS[i] for i in active],
        coefficient_roundoff_floor=floor,feasible_active_faces=len(candidates),
        equality_residual=float(np.max(abs(residual[:2]))),
        kkt_free_stationarity=stationarity,kkt_min_active_multiplier=bound_multiplier_min,
        success=True,method='convex active-face enumeration at fixed ranges',
        physical_candidate_accepted=False,total_stability_checked=False)
