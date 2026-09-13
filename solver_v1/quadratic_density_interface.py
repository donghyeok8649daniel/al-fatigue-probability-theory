"""Research interface with positive quadratic-envelope exponential density.

LJ and angular kernels/gauges stay unchanged. The new full scalar density is
used BOTH in the embedding and in rank1 coordination screening. No production
registration; no existing calibrated parameters are overwritten.
"""
from functools import lru_cache
import numpy as np
from scipy.special import gamma, gammaincc
from .polynomial_exponential_density import QuadraticEnvelopeDensity, plane_sum
from .fcc111_geometry import fcc111_geometry_from_b
from .coordination_screening import CoordinationScreenedInterface, screened_site_jet
from .joint_fcc_interface_calibration import MinimalConvexEmbedding
from .vector_interface_reference import VectorInterfaceEvaluation, _embedding_jet


def pack(value):
    return np.r_[value.value, value.d_d, value.grad_delta, value.d2_dd,
                 value.mixed_d_delta, value.hess_delta[0], value.hess_delta[1, 1]]


class QuadraticScalarEnvironment:
    def __init__(self, kernel, *, tolerance=2e-14):
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('positive tolerance required')
        self.geometry = fcc111_geometry_from_b(1.)
        self.h = self.geometry.h111
        self.tolerance = tolerance
        self.kernel = kernel
        raw = self.bulk()[0]
        if raw <= 0:
            raise ArithmeticError('positive full FCC environment required')
        # A fixed amplitude gauge, used at every later deformation.
        self.kernel = QuadraticEnvelopeDensity(kernel.C/raw, kernel.k, kernel.center, kernel.width)
        self.reference_density = self.bulk()[0]

    def same_plane(self, kernel):
        """Infinite triangular shell series, with an analytic Voronoi tail.

        This state-independent normalization excludes self; no arbitrary
        fixed neighbor cutoff replaces the cross-plane Poisson identity.
        """
        g = self.geometry; cover = g.b/np.sqrt(3)
        c, w, k = abs(kernel.center), kernel.width, kernel.k
        # (r+cover)*[(r+c)^2+w²] upper-envelope coefficients.
        poly = [cover*(c*c+w*w), c*c+w*w+2*cover*c, 2*c+cover, 1.]
        for radius in range(4, 66, 2):
            low = radius-2*cover
            if k*low <= 2:
                continue
            tail = 2*np.pi/g.atomic_cell_area*kernel.C*sum(
                a*gamma(n+1)*gammaincc(n+1,k*low)/k**(n+1) for n,a in enumerate(poly))
            if tail < self.tolerance*.1:
                extent = int(np.ceil(2*radius/g.b))+1
                m,n = np.meshgrid(np.arange(-extent,extent+1),np.arange(-extent,extent+1),indexing='ij')
                r = g.b*np.sqrt(m*m+m*n+n*n)
                keep = (r > 0)&(r <= radius)
                return float(np.sum(kernel.radial(r[keep]))), float(tail)
        raise ArithmeticError('same-plane density tail not resolved')

    def bulk(self, stretch=1.):
        k = self.kernel
        if stretch != 1:
            k = QuadraticEnvelopeDensity(k.C*stretch**2,k.k*stretch,k.center/stretch,k.width/stretch)
        same, tail = self.same_plane(k)
        # density, normal strain gradient/Hessian, engineering shear gradient/Hessian.
        result = np.array([same,0.,0.,0.,0.]); small=0
        for layer in range(1,81):
            d=layer*self.h
            v=plane_sum(k,d,self.geometry.abc_shift(layer),self.geometry,tolerance=self.tolerance)
            term=2*np.array([v.value,d*v.d_d,d*d*v.d2_dd,d*v.grad_delta[0],d*d*v.hess_delta[0,0]])
            result+=term
            last=float(np.max(abs(term)))
            small=small+1 if last < self.tolerance*.1 else 0
            if small>=3:
                return (*result,dict(layers=layer,last_layer=last,same_plane_tail=tail))
        raise ArithmeticError('bulk layer series not converged')

    def bulk_columns(self, step=4e-4):
        rho, first, second, shear_first, shear_second, _ = self.bulk()
        first/=rho; second/=rho; shear_first/=rho; shear_second/=rho
        values={i:self.bulk(1+i*step)[0]/rho for i in (-2,-1,1,2)}
        hydro_second=(-values[2]+16*values[1]-30+16*values[-1]-values[-2])/(12*step**2)
        gradients=np.array([3*first,first,shear_first])
        hessians=np.array([hydro_second,second,shear_second])
        return np.column_stack([np.r_[-.5*first,1.,.25*gradients**2-.5*hessians],
                                np.r_[first,-1.,hessians],np.r_[0.,1.,2*gradients**2]])

    @lru_cache(maxsize=256)
    def site_changes(self,state):
        a,x,y=map(float,state)
        if a<=0 or not np.all(np.isfinite(state)):
            raise ValueError('finite positive opening coordinate required')
        changes=[];small=0
        for layer in range(1,81):
            delta=self.geometry.abc_shift(layer)
            old=plane_sum(self.kernel,layer*self.h,delta,self.geometry,tolerance=self.tolerance).value
            jet=pack(plane_sum(self.kernel,a+(layer-1)*self.h,delta+[x,y],self.geometry,tolerance=self.tolerance))
            jet[0]-=old
            changes.append(jet/self.reference_density)
            last=float(np.max(abs(changes[-1])))
            small=small+1 if last<self.tolerance*.1 else 0
            if small>=3:
                return np.cumsum(np.asarray(changes)[::-1],axis=0)[::-1],dict(layers=layer,last_layer=last)
        raise ArithmeticError('interface density layers did not converge')


class QuadraticDensityInterface(CoordinationScreenedInterface):
    def __init__(self,shape,coefficients,kernel,*,tolerance=2e-14):
        if 2*shape[4]+min(shape[5],0.)*kernel.k<=0:
            raise ValueError('screened isolated-neighbor energy must decay with the NEW scalar kernel')
        super().__init__(shape,coefficients,tolerance=tolerance,law='power')
        self.scalar_environment=QuadraticScalarEnvironment(kernel,tolerance=tolerance)

    def evaluate(self,state):
        old=super().evaluate(state)
        density,diag=self.scalar_environment.site_changes(tuple(state))
        _,moments,_=super().site_inputs(tuple(state))
        components=dict(old.components)
        components['embedding']=_embedding_jet(density,1.,MinimalConvexEmbedding(*self.coefficients[2:5],1.))
        components['angular_1']=2*self.coefficients[6]*screened_site_jet(density,moments,self.screening,law='power')
        diagnostics=dict(old.diagnostics)
        diagnostics['quadratic_density_research']=dict(**diag,production_admissible=False)
        return VectorInterfaceEvaluation.from_jet(sum(components.values()),components,diagnostics)


def bulk_scalar_columns(bulk, environment, q_cubic):
    """Independent Bloch columns of the SAME normalized radial kernel.

    The infinite interface uses Poisson summation; this direct harmonic
    operator has analytical omitted-neighbor bounds for admissibility checks.
    """
    if bulk.stretch != 1:
        raise ValueError('this research operator currently covers pristine cubic bulk only')
    kernel=environment.kernel
    r=bulk.r; e=bulk.e; k=kernel.k; c=abs(kernel.center); w=kernel.width
    fp=kernel.radial(r,1)/environment.reference_density
    fpp=kernel.radial(r,2)/environment.reference_density
    gradient=fp[:,None]*e
    hessian=(fpp-fp/r)[:,None,None]*bulk.rr+(fp/r)[:,None,None]*np.eye(3)
    q=bulk.wavevector(q_cubic); phase=bulk.R@q; qn=np.linalg.norm(q)
    rh=np.einsum('n,nij->ij',2*np.sin(phase/2)**2,hessian)
    rg=np.sin(phase)@gradient
    columns=np.array([-rh+.25*np.outer(rg,rg),2*rh,2*np.outer(rg,rg)])
    tail=bulk.tail
    first=[2*c+k*(c*c+w*w),2+2*k*c,k]
    second=[2+4*k*c+k*k*(c*c+w*w),4*k+2*k*k*c,k*k]
    def integral(poly,power):
        return kernel.C/environment.reference_density*sum(v*tail.exponential(k,n+power) for n,v in enumerate(poly))
    dg=min(integral(first,0),qn*integral(first,1))
    def Hbound(power):
        return integral(second,power)+integral(first,power)/tail.lower
    dh=min(2*Hbound(0),.5*qn*qn*Hbound(2))
    square=2*np.linalg.norm(rg)*dg+dg*dg
    return columns,np.array([dh+.25*square,2*dh,2*square])
