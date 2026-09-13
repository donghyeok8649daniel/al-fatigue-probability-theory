"""True spatial gradient of the SAME exponential site density: research only.

g_i=grad_(central position) sum_j f(|r_j-r_i|), NOT sum_j r_ij f(r_ij).
LJ, density normalization, per-site counting and production models unchanged.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from .polynomial_exponential_density import transform_terms
from .fcc111_geometry import fcc111_geometry_from_b
from .fcc111_lattice_sum import triangular_reciprocal_shells
from .tail_constrained_material import normalized_amplitude


ORDERS=((0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2))


@dataclass(frozen=True)
class ExponentialSpatialKernel:
    k: float
    C: float

    def __post_init__(self):
        if not np.isfinite(self.k+self.C) or min(self.k,self.C)<=0:
            raise ValueError('positive finite decay/amplitude required')

    def radial(self,r,order=0):
        r=np.asarray(r,float)
        if order not in (0,1,2,3) or np.any(r<0) or np.any(~np.isfinite(r)):
            raise ValueError('finite radius and derivative order0..3 required')
        return (-self.k)**order*self.C*np.exp(-self.k*r)

    def plane_transform(self,d,G,d_order=0):
        if d_order not in (0,1,2,3) or not np.isfinite(d) or d<=0 or not np.isfinite(G) or G<0:
            raise ValueError('positive separation/nonnegative wavenumber and order0..3 required')
        Q=np.hypot(self.k,G)
        return 2*np.pi*self.C*np.exp(-d*Q)*sum(v*self.k**h*d**j*Q**(-p)
            for h,j,p,v in transform_terms(0,d_order))


def gradient_plane(kernel,d,delta,geometry,*,tolerance=2e-12):
    """30 analytic entries: three gradient components, each with a full 10-jet.

Central-position gradient is negative relative-neighbor gradient. Three
small reciprocal envelopes are empirical stopping diagnostics, not a bound.
"""
    delta=np.asarray(delta,float)
    if delta.shape!=(2,) or np.any(~np.isfinite(delta)) or not np.isfinite(tolerance) or tolerance<=0:
        raise ValueError('finite two-dimensional registry and positive tolerance required')
    result=np.zeros((10,3));small=0
    shells=[(0.,np.zeros((1,2)))]+[(s.magnitude,s.vectors) for s in triangular_reciprocal_shells(geometry,32)]
    for count,(magnitude,vectors) in enumerate(shells):
        normal=[kernel.plane_transform(d,magnitude,n)/geometry.atomic_cell_area for n in range(4)]
        coefficients=np.zeros((10,3,len(vectors)),complex)
        for row,alpha in enumerate(ORDERS):
            for component in range(3):
                order=np.array(alpha)+np.eye(3,dtype=int)[component]
                coefficients[row,component]=-normal[order[0]]*(1j*vectors[:,0])**order[1]*(1j*vectors[:,1])**order[2]
        result+=np.real(coefficients@np.exp(1j*(vectors@delta)))
        last=float(np.max(np.sum(abs(coefficients),axis=2)))
        small=small+1 if count and last<tolerance*.1 else 0
        if small>=3:
            return result,dict(shells=count,last_shell_envelope=last)
    raise ArithmeticError('density-gradient reciprocal series did not converge')


class SpatialDensityGradient:
    def __init__(self,decay,*,tolerance=2e-12):
        self.geometry=fcc111_geometry_from_b(1.)
        self.h=self.geometry.h111;self.tolerance=tolerance
        self.kernel=ExponentialSpatialKernel(decay,normalized_amplitude(decay))

    @lru_cache(maxsize=192)
    def site_jets(self,state):
        a,x,y=map(float,state)
        if a<=0 or not np.all(np.isfinite(state)):
            raise ValueError('finite state with positive normal spacing required')
        changes=[];small=0
        for layer in range(1,81):
            delta=self.geometry.abc_shift(layer)
            old,_=gradient_plane(self.kernel,layer*self.h,delta,self.geometry,tolerance=self.tolerance)
            new,_=gradient_plane(self.kernel,a+(layer-1)*self.h,delta+[x,y],self.geometry,tolerance=self.tolerance)
            # Bulk total gradient is zero. Intra-half-crystal contributions
            # do not move; cancel them analytically with the old cross terms.
            new[0]-=old[0]
            changes.append(new)
            last=float(np.max(abs(new)))
            small=small+1 if last<self.tolerance*.1 else 0
            if small>=3:
                return np.cumsum(np.asarray(changes)[::-1],axis=0)[::-1],dict(layers=layer,last_layer=last)
        raise ArithmeticError('density-gradient interface layers did not converge')


def gradient_bulk_column(bulk,q_cubic):
    """D |grad x|²/x at cubic x=1,g=0: H_D=2 r_h.T r_h = O(q^4).

This is different from the existing embedding F'' r_g r_g.T term.
"""
    if bulk.stretch!=1:
        raise ValueError('pristine cubic reference required')
    q=bulk.wavevector(q_cubic);qn=np.linalg.norm(q)
    rh=np.einsum('n,nij->ij',2*np.sin(bulk.R@q/2)**2,bulk.rho_hessian)
    k=bulk.decays[0];tail=bulk.tail
    dh=normalized_amplitude(k)*(k*k+k/tail.lower)*min(2*tail.exponential(k,0),.5*qn*qn*tail.exponential(k,2))
    return 2*rh.T@rh,float(2*(2*np.linalg.norm(rh,2)*dh+dh*dh))
