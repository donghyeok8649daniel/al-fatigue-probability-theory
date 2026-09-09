"""One-shape-parameter exponential-density ablation, static research only.

After the v15 same-family joint fit fails interface-shape validation, test a
positive mixture of TWO already-defined infinite exponential environments:

 x_i = (1-w) rho_i(k_scalar)/rho_bulk(k_scalar)
             + w rho_i(k_odd)/rho_bulk(k_odd), 0 <= w <= 1.

Only AFTER this site-wise sum is F(x_i) applied. This is NOT a sum of two
embedding energies. Both density normalizations are fixed at reference FCC,
not reset during deformation. LJ and angular terms are unchanged. The second
decay reuses k_odd as an explicit minimal-family assumption, not a new fitted
physical length. w=0 gives the unchanged v14/v15 quartic family.
"""
from functools import lru_cache

import numpy as np

from .full_fcc_calibration_audit import bulk_basis
from .joint_fcc_interface_calibration import MinimalConvexEmbedding
from .quartic_angular_material import (
    QuarticSymmetryInterface, QuarticSymmetryObservationCache,
    QuarticSymmetryTailBulkBasis,
)
from .range_resolved_material import build_range_surface
from .tail_constrained_material import normalized_amplitude
from .vector_interface_reference import _embedding_jet, VectorInterfaceEvaluation


def mixing_weight(value):
    value=float(value)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('finite density mixing weight in [0,1] required')
    return value


def combine_density_jets(first, second, weight):
    """Pad missing exponentially converged depths, never duplicate site energy."""
    weight=mixing_weight(weight)
    first,second=np.asarray(first),np.asarray(second)
    if (first.ndim!=2 or second.ndim!=2 or first.shape[1]!=10 or second.shape[1]!=10
            or np.any(~np.isfinite(first)) or np.any(~np.isfinite(second))):
        raise ValueError('two finite per-site ten-component density jets required')
    mixed=np.zeros((max(len(first),len(second)),10))
    mixed[:len(first)]+=(1-weight)*first
    mixed[:len(second)]+=weight*second
    return mixed


class MixedDensityQuarticInterface(QuarticSymmetryInterface):
    def __init__(self,decays,coefficients,*,tolerance=2e-11):
        if len(decays)!=4:raise ValueError('three decays and one mixing weight required')
        self.weight=mixing_weight(decays[3])
        super().__init__(decays[:3],coefficients,tolerance=tolerance)
        self.shapes=tuple(map(float,decays))
        self.second_density=build_range_surface(decays[1],decays[1],decays[2],
                                                np.ones(8),tolerance=tolerance)
        self.normalized_embedding=MinimalConvexEmbedding(*self.coefficients[2:5],1.)

    def density_changes(self,state):
        a,*diag0=self.base._density(np.asarray(state,float))
        b,*diag1=self.second_density._density(np.asarray(state,float))
        return combine_density_jets(a/self.face.rho_bulk,
            b/self.second_density.face.rho_bulk,self.weight),dict(first=diag0,second=diag1)

    def evaluate(self,state):
        value=super().evaluate(state)
        if self.weight==0:return value
        changes,diag=self.density_changes(state)
        components=dict(value.components)
        # REPLACE the old scalar contribution, not add another embedding.
        components['embedding']=_embedding_jet(changes,1.,self.normalized_embedding)
        diagnostics=dict(value.diagnostics);diagnostics['mixed_density']=diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()),components,diagnostics)


def mixed_bulk_embedding_columns(first_decay, second_decay, weight, *, strain_step=8e-4):
    """Three scalar columns for force/cohesion/hydro/normal/shear.

    The x affine derivatives are combined BEFORE the F chain rule. Cubic
    symmetry implies x_hydro'=3 x_normal', x_shear'=0. The existing infinite-
    sum five-point hydrostatic x'' is retained and independently refined;
    normal and shear derivatives are analytic Poisson derivatives.
    F(0)= -B+C, so cohesion columns remain (1,-1,1).
    """
    weight=mixing_weight(weight)
    first=bulk_basis(float(first_decay),strain_step)
    second=bulk_basis(float(second_decay),strain_step)
    linear=(1-weight)*first[:,3]+weight*second[:,3]
    normal=linear[0]
    gradient=np.array([3*normal,normal,0.])
    sqrt=np.r_[-.5*normal,1.,.25*gradient**2-.5*linear[2:5]]
    quadratic=np.r_[0.,1.,2*gradient**2]
    if weight==0:
        # Preserve legacy endpoint to roundoff, including its refined FD
        # hydrostatic sqrt evaluation convention.
        sqrt=first[:,2].copy()
    elif weight==1:
        sqrt=second[:,2].copy()
    return np.column_stack([sqrt,linear,quadratic])


class MixedDensityObservationCache(QuarticSymmetryObservationCache):
    @lru_cache(maxsize=128)
    def density_data(self,decay):
        probe=build_range_surface(decay,4.,5.,np.ones(8))
        states={o.state for o in self.observations if o.state is not None}
        return {state:probe._density(np.asarray(state))[0]/probe.face.rho_bulk for state in states}

    def matrix(self,shapes):
        if len(shapes)!=4:raise ValueError('three decays and one mixing weight required')
        k0,k1,k2,w=map(float,shapes);w=mixing_weight(w)
        matrix=super().matrix((k0,k1,k2)).copy()
        if w==0:return matrix
        bulk=mixed_bulk_embedding_columns(k0,k1,w)
        first,second=self.density_data(k0),self.density_data(k1)
        jets={}
        for state in first:
            changes=combine_density_jets(first[state],second[state],w)
            jets[state]=np.column_stack([_embedding_jet(changes,1.,MinimalConvexEmbedding(*c,1.))
                                        for c in np.eye(3)])
        for i,o in enumerate(self.observations):
            matrix[i,2:5]=(bulk[o.bulk_index] if o.bulk_index is not None else
                          np.asarray(o.jet_weights)@jets[o.state])
        return matrix


class MixedDensityTailBulkBasis(QuarticSymmetryTailBulkBasis):
    def __init__(self,shapes,*,radius=12.):
        if len(shapes)!=4:raise ValueError('three decays and one mixing weight required')
        self.weight=mixing_weight(shapes[3])
        super().__init__(shapes[:3],radius=radius)
        k=float(shapes[1]);f=normalized_amplitude(k)*np.exp(-k*self.r)
        self.second_gradient=-k*f[:,None]*self.e
        self.second_hessian=f[:,None,None]*(k*k*self.rr
            -k/self.r[:,None,None]*(np.eye(3)-self.rr))

    def evaluate(self,q_cubic):
        columns,errors=super().evaluate(q_cubic)
        if self.weight==0:return columns,errors
        q=self.wavevector(q_cubic);qn=np.linalg.norm(q);phase=self.R@q
        cosine=2*np.sin(phase/2)**2;sine=np.sin(phase);w=self.weight
        rh=np.einsum('n,nij->ij',cosine,(1-w)*self.rho_hessian+w*self.second_hessian)
        rg=sine@((1-w)*self.rho_gradient+w*self.second_gradient)
        columns[2:5]=np.array([-rh+.25*np.outer(rg,rg),2*rh,2*np.outer(rg,rg)])
        def envelopes(k):
            t=self.tail;amplitude=normalized_amplitude(k)
            dh=amplitude*(k*k+k/t.lower)*min(2*t.exponential(k,0),.5*qn*qn*t.exponential(k,2))
            dg=amplitude*k*min(t.exponential(k,0),qn*t.exponential(k,1))
            return np.array([dh,dg])
        dh,dg=(1-w)*envelopes(self.decays[0])+w*envelopes(self.decays[1])
        square=2*np.linalg.norm(rg)*dg+dg*dg
        errors[2:5]=[dh+.25*square,2*dh,2*square]
        return columns,errors
