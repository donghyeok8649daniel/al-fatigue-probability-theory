"""Nested one-amplitude nonlinear angular response, static research only.

J3=sum_i ||Q3_i||^4 is summed PER ATOM. Q3_i is the same infinite exponential
STF moment, with all surrounding atoms summed before taking its norm. At
centrosymmetric bulk Q3=0 for any affine strain. At a pristine cut Q3=O(delta q)
so J3 is fourth order: it changes neither bulk nor perfect-interface Hessians.
It can test whether the rigid quadratic relation between harmonic interface
stiffness and finite fault energy is responsible for the measured fit conflict.
LJ, exponential kernels, Poisson/Bessel derivatives and all old defaults stay.
The single K3>=0 is a nested research coefficient, not an empirical yield law.
"""
from functools import lru_cache

import numpy as np

from .symmetry_resolved_material import SymmetryResolvedInterface, SymmetryObservationCache, SymmetryTailBulkBasis
from .range_resolved_material import build_range_surface
from .vector_interface_reference import HESSIAN_INDICES, VectorInterfaceEvaluation


def per_site_quartic_jet(depths,saturation=0.):
    """H(I)=I^2/(1+alpha I); alpha=0 is the independently tested quartic.

    This rational shape is a separate one-parameter nested ablation. H'=
    I(2+alpha I)/(1+alpha I)^2 and H''=2/(1+alpha I)^3. It is nonnegative,
    increasing/convex for I>=0,alpha>=0; no clipping or target-dependent switch.
    """
    if not np.isfinite(saturation) or saturation<0: raise ValueError('alpha>=0 required')
    t=depths[:,0];g=depths[:,1:4];h=depths[:,HESSIAN_INDICES]
    norm=np.einsum('ri,ri->r',t,t)
    ng=2*np.einsum('rc,ric->ri',t,g)
    nh=2*(np.einsum('ric,rjc->rij',g,g)+np.einsum('rc,rijc->rij',t,h))
    den=1+saturation*norm
    first=norm*(2+saturation*norm)/den**2;second=2/den**3
    grad=np.einsum('r,ri->i',first,ng)
    hess=np.einsum('r,ri,rj->ij',second,ng,ng)+np.einsum('r,rij->ij',first,nh)
    return np.r_[np.sum(norm**2/den),grad,hess[0],hess[1,1:],hess[2,2]]


class QuarticMoment:
    def __init__(self,base,*,saturation=0.):
        if not np.isfinite(saturation) or saturation<0:
            raise ValueError('finite alpha>=0 required')
        self.base=base;self.moment=base.moments[0][1];self.baseline={}
        self.saturation=float(saturation)
        if self.moment.invariant.rank!=3: raise ValueError('derived rank3 environment required')

    @lru_cache(maxsize=128)
    def jet(self,state):
        a,x,y=state;obj=self.moment.invariant;terms=[];small=0;maxshell=0
        for layer in range(1,obj.max_layers+1):
            delta=self.base.geometry.abc_shift(layer)
            if layer not in self.baseline:
                self.baseline[layer]=self.moment.evaluate(layer*self.base.h,delta)[0][0].copy()
            value,shells,_=self.moment.evaluate(a+(layer-1)*self.base.h,delta+np.array([x,y]))
            value[0]-=self.baseline[layer];terms.append(value);maxshell=max(maxshell,shells)
            last=float(np.max(abs(value)));small=small+1 if last<obj.tolerance else 0
            if small>=4:
                depths=np.cumsum(np.array(terms)[::-1],axis=0)[::-1]
                return 2*per_site_quartic_jet(depths,self.saturation),dict(layers=layer,shells=maxshell,last=last)
        raise ArithmeticError('quartic angular neighborhood did not converge')


class QuarticSymmetryInterface(SymmetryResolvedInterface):
    def __init__(self,decays,coefficients,*,tolerance=2e-11,saturation=0.):
        c=np.asarray(coefficients,float)
        if c.shape!=(10,) or np.any(~np.isfinite(c)) or c[9]<0:
            raise ValueError('ten finite coefficients and K3>=0 required')
        super().__init__(decays,c[:9],tolerance=tolerance)
        self.coefficients=c.copy()
        self.quartic_amplitude=float(c[9]);self.quartic=QuarticMoment(self.base,saturation=saturation)

    def evaluate(self,state):
        value=super().evaluate(state)
        if self.quartic_amplitude==0: return value
        jet,diag=self.quartic.jet(tuple(state))
        components=dict(value.components);components['angular_rank3_quartic']=self.quartic_amplitude*jet
        diagnostics=dict(value.diagnostics);diagnostics['angular_rank3_quartic']=diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()),components,diagnostics)


class QuarticSymmetryObservationCache(SymmetryObservationCache):
    @lru_cache(maxsize=128)
    def quartic(self,odd,saturation=0.):
        model=build_range_surface(3.,odd,5.,np.ones(8));moment=QuarticMoment(model,saturation=saturation)
        return np.array([0. if o.bulk_index is not None else
            np.asarray(o.jet_weights)@moment.jet(o.state)[0] for o in self.observations])

    def matrix(self,decays):
        if len(decays) not in (3,4):
            raise ValueError('three radial decays and optional one angular shape required')
        radial=decays[:3];alpha=float(decays[3]) if len(decays)==4 else 0.
        return np.column_stack([super().matrix(radial),self.quartic(float(decays[1]),alpha)])


class QuarticSymmetryTailBulkBasis(SymmetryTailBulkBasis):
    def evaluate(self,q_cubic):
        columns,tail=super().evaluate(q_cubic)
        return np.concatenate([columns,np.zeros((1,3,3))]),np.r_[tail,0.]
