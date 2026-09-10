"""Fixed-shape least-squares LOWER BOUND, not an adopted material fit.

Allow every sign in the null space of the declared exact bulk/interface
constraints. If even this relaxed problem cannot reproduce independent core
forces, imposing physical inequalities cannot improve its minimum residual.
No fatigue, yield, core-barrier target or parameter search is inserted.
"""
import numpy as np
from scipy.linalg import null_space


def force_compatibility_lower_bound(design,target,exact,exact_target,baseline,*,scales):
    design,target,exact,exact_target,baseline,scales=map(np.asarray,
        (design,target,exact,exact_target,baseline,scales))
    n=len(baseline)
    if (design.ndim!=2 or not len(design) or design.shape[1]!=n or target.shape!=(len(design),)
            or exact.ndim!=2 or exact.shape[1]!=n or exact_target.shape!=(len(exact),)
            or scales.shape!=(n,) or np.any(scales<=0)
            or any(np.any(~np.isfinite(a)) for a in (design,target,exact,exact_target,baseline,scales))):
        raise ValueError('finite consistently dimensioned force/constraint arrays and positive coefficient scales required')
    scaled_exact=exact*scales
    row_scale=np.linalg.norm(scaled_exact,axis=1)
    if np.any(row_scale==0):
        raise ValueError('zero exact-constraint rows require an explicit gauge audit')
    A=scaled_exact/row_scale[:,None]
    correction=np.linalg.lstsq(A,(exact_target-exact@baseline)/row_scale,rcond=None)[0]
    particular=baseline+scales*correction
    if np.max(abs(exact@particular-exact_target)/row_scale,initial=0.) > 1e-10:
        raise ValueError('inconsistent exact constraints cannot define a force lower bound')
    Z=null_space(A)
    projected=design@(scales[:,None]*Z)
    shift=np.linalg.lstsq(projected,target-design@particular,rcond=None)[0]
    best=particular+scales*(Z@shift)
    before=design@baseline-target; after=design@best-target
    return dict(coefficients=best,baseline_coefficients=baseline,coefficient_scales=scales,
        exact_rank=n-Z.shape[1],null_dimension=Z.shape[1],
        exact_singular_values=np.linalg.svd(A,compute_uv=False),
        projected_force_singular_values=np.linalg.svd(projected,compute_uv=False),
        initial_force_rms=float(np.sqrt(np.mean(before**2))),
        unconstrained_lower_bound_force_rms=float(np.sqrt(np.mean(after**2))),
        maximum_initial_force_error=float(np.max(abs(before))),
        maximum_lower_bound_force_error=float(np.max(abs(after))),
        normalized_exact_residual=float(np.max(abs(exact@best-exact_target)/row_scale)),
        force_predictions=design@best,baseline_force_predictions=design@baseline,
        least_squares_normal_residual=float(np.max(abs(projected.T@after),initial=0.)),
        physical_inequalities_enforced=False,relaxed_core_recomputed=False,material_accepted=False,
        fixed_shape_only=True,global_family_impossibility_claimed=False)


def constrained_force_profile(design,target,exact,exact_target,*,nonnegative,
                              operators=None,tails=None):
    """Exact coefficient QP for frozen source forces, not a relaxed material fit.

    Every force component has the same eV/L0 units and unit weight; no yield,
    lifetime or desired barrier is a target. Other source observables must be
    audited separately. Sampled bulk spectral inequalities are optional.
    """
    from .tail_constrained_material import convex_profile,spectral_profile
    from .vector_material_calibration import MaterialObservation
    F,t,E,rhs=map(np.asarray,(design,target,exact,exact_target))
    if (F.ndim!=2 or not len(F) or E.ndim!=2 or F.shape[1]!=E.shape[1]
            or t.shape!=(len(F),) or rhs.shape!=(len(E),)):
        raise ValueError('matching nonempty force and exact constraint arrays required')
    matrix=np.vstack([E,F])
    observations=[MaterialObservation(f'exact_{i}',float(value),1.,'declared source unit','exact')
                  for i,value in enumerate(rhs)]
    observations += [MaterialObservation(f'core_gradient_{i}',float(value),1.,'eV/L0','fit')
                     for i,value in enumerate(t)]
    if (operators is None)!=(tails is None):
        raise ValueError('spectral operators and tail bounds must be supplied together')
    if operators is None:
        fit=convex_profile(matrix,observations,nonnegative=nonnegative)
    else:
        fit=spectral_profile(matrix,observations,operators,tails,nonnegative=nonnegative)
    fit.update(force_rms=float(np.sqrt(fit['squared_loss']/len(F))),
               relaxed_core_recomputed=False,material_accepted=False,
               fixed_shape_only=True,physical_yield_used=False)
    return fit
