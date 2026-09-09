"""Independent elastic-constant metric for material-to-yield research.

This does not fit yield. It exposes how errors in three measured cubic
constants propagate through the exact strain-mode transformation. The legacy
5%-per-mode metric is retained for comparison, not silently overwritten.
"""
from dataclasses import replace
from itertools import combinations

import numpy as np
from scipy.linalg import null_space

from .aluminum_calibration import EV_J, GPA_PA
from .full_fcc_calibration_audit import atomic_volume_m3
from .vector_material_calibration import COEFFICIENTS


def cubic_to_mode_matrix():
    """[C11,C12,C44] GPa -> hydro/normal/simple-shear Hessians, eV."""
    return atomic_volume_m3()*GPA_PA/EV_J*np.array(
        [[3., 6., 0.], [1/3, 2/3, 4/3], [1/3, -1/3, 1/3]])


def cubic_metric_problem(matrix, observations, *, relative_scale=.05):
    """Change both predictions AND targets into independent cubic constants.

    The diagonal 5% C metric is an explicitly DIFFERENT discrepancy model
    from independent 5% H modes. Neither is experimental error covariance.
    All held-out observations and exact constraints remain unchanged.
    """
    if not np.isfinite(relative_scale) or relative_scale <= 0:
        raise ValueError('positive finite discrepancy scale required')
    m = np.asarray(matrix, float).copy()
    if m.shape != (len(observations), 8):
        raise ValueError('eight energy coefficients and matching observations required')
    if [o.bulk_index for o in observations[:5]] != list(range(5)):
        raise ValueError('original five independent bulk rows must be first')
    inverse = np.linalg.inv(cubic_to_mode_matrix())
    m[2:5] = inverse@m[2:5]
    target = inverse@np.array([o.target for o in observations[2:5]])
    obs = list(observations)
    for i, name in enumerate(('C11_GPa', 'C12_GPa', 'C44_GPa')):
        obs[i+2] = replace(obs[i+2], name=name, target=float(target[i]),
            scale=relative_scale*abs(float(target[i])), units='GPa', bulk_index=None)
    return m, obs


def transformed_mode_covariance(mode_scales):
    """Preserve the OLD objective if changing coordinates without new weights."""
    scale = np.asarray(mode_scales, float)
    if scale.shape != (3,) or np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError('three positive mode scales required')
    inverse = np.linalg.inv(cubic_to_mode_matrix())
    return inverse@np.diag(scale**2)@inverse.T


def fit_exact_bulk(matrix, observations):
    """Five exact bulk targets, interface least squares, same coefficient sector.

    A fixed-decay convex subproblem. All null directions, equality residuals
    and bound failures are exposed. An infeasible case returns no candidate.
    This is a tradeoff diagnostic, not a claim targets have zero uncertainty.
    """
    matrix = np.asarray(matrix, float)
    if matrix.shape != (len(observations), 8) or np.any(~np.isfinite(matrix)):
        raise ValueError('finite matched eight-column matrix required')
    target = np.array([o.target for o in observations])
    scales = np.array([o.scale for o in observations])
    equal = matrix[:5]/scales[:5, None]
    rhs_equal = target[:5]/scales[:5]
    singular = np.linalg.svd(equal, compute_uv=False)
    numerical_floor = np.finfo(float).eps*max(equal.shape)*singular[0]
    rank = int(np.sum(singular > numerical_floor))
    if rank != 5:
        raise ValueError('five independently resolved bulk constraints required')
    constrained = np.array([0, 1, 2, 4, 5])  # u,v,A,C,D3; signed B,D1,D2
    chosen = [i for i,o in enumerate(observations) if i>=5 and o.role=='fit']
    # Only three null coordinates remain. Enumerating independent active faces
    # (<=3 of5 bounds) solves this convex QP with linear algebra. Explicitly
    # setting an ACTIVE coefficient to0 is its equality parameterization, not
    # post-fit force clipping. This avoids SLSQP's near-zero-bound failures.
    candidates=[]
    for number in range(4):
        for active in combinations(constrained,number):
            free=[i for i in range(8) if i not in active]
            E=equal[:,free]
            if np.linalg.matrix_rank(E)!=5: continue
            offset=np.linalg.lstsq(E,rhs_equal,rcond=None)[0]
            tangent=null_space(E)
            design=matrix[np.ix_(chosen,free)]@tangent/scales[chosen,None]
            rhs=(target[chosen]-matrix[np.ix_(chosen,free)]@offset)/scales[chosen]
            z=np.linalg.lstsq(design,rhs,rcond=None)[0]
            values=offset+tangent@z
            c=np.zeros(8); c[free]=values
            floor=np.zeros(8)
            floor[free]=32*np.finfo(float).eps*np.linalg.cond(E)*(
                abs(offset)+abs(tangent)@abs(z))
            residual=(matrix@c-target)/scales
            if np.any(c[constrained]<-floor[constrained]) or np.max(abs(residual[:5]))>=1e-8:
                continue
            candidates.append((float(residual[chosen]@residual[chosen]),c,floor,active))
    if not candidates:
        return dict(admissible=False, coefficients=None,
            message='no feasible active-face least-squares solution',
            exact_bulk_rank=rank,equality_singular_values=singular)
    loss,c,floor,active=min(candidates,key=lambda x:x[0])
    residual = (matrix@c-target)/scales
    error = float(np.max(abs(residual[:5])))
    full_tangent=null_space(equal)
    raw_design=matrix[chosen]@full_tangent/scales[chosen,None]
    design=raw_design/np.linalg.norm(raw_design,axis=0)
    return dict(admissible=True, coefficients=c, coefficient_order=COEFFICIENTS,
        success=True, message='convex active-face enumeration; no radial global claim', predictions=matrix@c,
        squared_loss=loss, residuals=residual,active_zero_coefficients=[COEFFICIENTS[i] for i in active],
        feasible_active_faces=len(candidates),coefficient_roundoff_floor=floor,
        strictly_positive_LJ_resolved=bool(np.all(c[:2]>floor[:2])),
        exact_bulk_rank=rank, equality_singular_values=singular,
        exact_bulk_residual_max=error, reduced_jacobian=design,
        singular_values=np.linalg.svd(design,compute_uv=False),
        selected_rows=chosen, physical_candidate_accepted=False)
