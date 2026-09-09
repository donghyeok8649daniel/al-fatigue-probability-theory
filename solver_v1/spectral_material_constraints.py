"""Necessary finite-q stability constraints on the same linear energy family.

At fixed radial ranges H(q;c)=sum_j c_j H_j(q). A negative eigenvector v
supplies the linear separating halfspace sum_j c_j v*H_j(q)v >= 0. This is
not an empirical strength criterion. The direct-radius Hessians are independent
validation; they do not replace canonical LJ/Bessel energy or its derivatives.
"""
from itertools import combinations

import numpy as np
from scipy.linalg import null_space

from .range_resolved_material import build_range_surface,RangeBulkValidation
from .vector_material_calibration import COEFFICIENTS


class BulkHessianCoefficientBasis:
    def __init__(self,decays,*,radius=8.):
        self.decays=tuple(decays); self.radius=float(radius)
        self.origin=np.ones(8)
        self.models=[RangeBulkValidation(build_range_surface(*decays,c),cutoff=radius)
            for c in [self.origin,*(self.origin+v for v in np.eye(8))]]

    def matrices(self,fractional_cubic):
        q=self.models[0].crystallographic_wavevector(fractional_cubic)
        base=self.models[0].evaluate(q)['matrix']
        columns=np.array([model.evaluate(q)['matrix']-base for model in self.models[1:]])
        # This is a coefficient difference of an exactly linear operator,
        # not finite differentiation of coordinates or a fitted surrogate.
        if np.max(abs(columns.sum(axis=0)-base))>2e-10*max(1.,np.max(abs(base))):
            raise ArithmeticError('finite-q coefficient linearity failed')
        return columns


def stability_halfspace(columns,coefficients):
    columns=np.asarray(columns,float); c=np.asarray(coefficients,float)
    if columns.shape!=(8,3,3) or c.shape!=(8,):
        raise ValueError('eight symmetric three-coordinate Hessian columns required')
    matrix=np.einsum('c,cij->ij',c,columns)
    values,vectors=np.linalg.eigh(matrix); v=vectors[:,0]
    row=np.einsum('i,cij,j->c',v,columns,v)
    return dict(minimum_eigenvalue=float(values[0]),polarization=v,halfspace=row)


def fit_spectral_faces(matrix,observations,halfspaces):
    """Convex coefficient profile with explicit spectral separating planes.

    Enumerating independent active faces avoids interpreting a failed QP
    optimizer flag as a physical result. No fitted force or coefficient is
    clipped. Pair coefficients on their zero boundary are not admissible LJ.
    """
    matrix=np.asarray(matrix,float); extra=np.asarray(halfspaces,float).reshape((-1,8))
    if matrix.shape!=(len(observations),8) or not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(extra)):
        raise ValueError('finite observation and halfspace matrices required')
    if [i for i,o in enumerate(observations) if o.role=='exact']!=[0,1]:
        raise ValueError('the same two reference constraints are required')
    target=np.array([o.target for o in observations]); scales=np.array([o.scale for o in observations])
    equal=matrix[:2]/scales[:2,None]; rhs=target[:2]/scales[:2]
    if np.linalg.matrix_rank(equal)!=2:
        raise ValueError('independent reference constraints required')
    bounded=(0,1,2,4,5)
    G=np.vstack([np.eye(8)[list(bounded)],extra]); norms=np.linalg.norm(G,axis=1)
    if np.any(norms==0): raise ValueError('nonzero spectral halfspace required')
    G=G/norms[:,None]; selected=[i for i,o in enumerate(observations) if o.role=='fit']
    candidates=[]
    for count in range(min(6,len(G))+1):
        for active in combinations(range(len(G)),count):
            zero=[bounded[i] for i in active if i<len(bounded)]
            free=[i for i in range(8) if i not in zero]
            spectral=[i for i in active if i>=len(bounded)]
            E=np.vstack([equal[:,free],G[np.ix_(spectral,free)]])
            right=np.r_[rhs,np.zeros(len(spectral))]
            if np.linalg.matrix_rank(E)!=len(E): continue
            offset=np.linalg.lstsq(E,right,rcond=None)[0]; T=null_space(E)
            design=matrix[np.ix_(selected,free)]@T/scales[selected,None]
            z=np.linalg.lstsq(design,(target[selected]-matrix[np.ix_(selected,free)]@offset)/scales[selected],rcond=None)[0]
            c=np.zeros(8); c[free]=offset+T@z
            floor=np.zeros(8)
            floor[free]=64*np.finfo(float).eps*np.linalg.cond(E)*(abs(offset)+abs(T)@abs(z))
            residual=(matrix@c-target)/scales
            if np.min(G@c+abs(G)@floor)<0 or np.max(abs(residual[:2]))>1e-8: continue
            candidates.append((float(residual[selected]@residual[selected]),c,floor,active))
    if not candidates:
        raise ArithmeticError('no feasible spectral constrained coefficient profile')
    loss,c,floor,active=min(candidates,key=lambda row:row[0])
    residual=(matrix@c-target)/scales
    grad=2*(matrix[selected]/scales[selected,None]).T@residual[selected]
    operator=np.column_stack([equal.T,-G[list(active)].T])
    multipliers=np.linalg.lstsq(operator,-grad,rcond=None)[0]
    return dict(coefficients=c,coefficient_order=COEFFICIENTS,squared_loss=loss,
        residuals=residual,predictions=matrix@c,selected_rows=selected,
        active_constraints=list(active),feasible_active_faces=len(candidates),
        coefficient_roundoff_floor=floor,minimum_scaled_halfspace=float(np.min(G@c)),
        kkt_stationarity=float(np.max(abs(grad+operator@multipliers))),
        kkt_min_inequality_multiplier=float(np.min(multipliers[2:])) if active else None,
        strictly_positive_LJ_resolved=bool(np.all(c[:2]>floor[:2])),
        admissible=True,physical_candidate_accepted=False)
