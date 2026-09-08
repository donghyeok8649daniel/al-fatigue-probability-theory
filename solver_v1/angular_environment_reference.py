"""Minimal analytic angular-environment RESEARCH extension.

LJ and the scalar EAM energy are unchanged. The first nested test adds a
nonnegative third-moment invariant. Explicit first/second-moment companions
are separate research extensions; signed coefficients require total stability
checks. Infinite moments follow from derivatives of the SAME exponential
Poisson transform, never from a production neighbor cutoff.

This is not standard/calibrated MEAM and is not registered with the PDE.
"""
from __future__ import annotations

from functools import lru_cache
import itertools
import math
import numpy as np

from .fcc111_lattice_sum import triangular_reciprocal_shells


INDICES=tuple(itertools.product(range(3),repeat=3))


@lru_cache(maxsize=64)
def _radial_terms(t_order,normal_order,z_power):
    """Exact polynomial recurrence for d^j q^-k exp(-d q).

    D_t=(2q)^-1 D_q with t=|G|^2 and q=sqrt(kappa^2+t).
    Fhat/(2 pi kappa)=exp(-d q)(q^-3+d q^-2).
    """
    terms={(0,3):1.,(1,2):1.}
    for _ in range(t_order):
        new={}
        for (j,k),v in terms.items():
            new[j+1,k+1]=new.get((j+1,k+1),0.)-.5*v
            new[j,k+2]=new.get((j,k+2),0.)-.5*k*v
        terms=new
    terms={(j+z_power,k):v for (j,k),v in terms.items()}
    for _ in range(normal_order):
        new={}
        for (j,k),v in terms.items():
            if j:
                new[j-1,k]=new.get((j-1,k),0.)+j*v
            new[j,k-1]=new.get((j,k-1),0.)-v
        terms=new
    return tuple((j,k,v) for (j,k),v in terms.items())


def _radial(d,q,kappa,t_order,normal_order,z_power):
    return 2*math.pi*kappa*np.exp(-d*q)*sum(
        v*d**j*q**(-k) for j,k,v in _radial_terms(t_order,normal_order,z_power))


def _tensor_fourier(d,vectors,kappa,normal_order):
    """i^order D_G^order Fhat times d^z; exact 27 Cartesian components."""
    q=np.sqrt(kappa*kappa+np.sum(vectors*vectors,axis=1))
    radial={(t,z):_radial(d,q,kappa,t,normal_order,z)
            for z in range(4) for t in range(4-z)}
    values=[]
    for indices in INDICES:
        xy=[i for i in indices if i!=2]; z=3-len(xy)
        if not xy:
            value=radial[0,z]
        elif len(xy)==1:
            value=2j*vectors[:,xy[0]]*radial[1,z]
        elif len(xy)==2:
            i,j=xy
            value=-(2*(i==j)*radial[1,z]+4*vectors[:,i]*vectors[:,j]*radial[2,z])
        else:
            i,j,k=xy
            value=-1j*(4*((i==j)*vectors[:,k]+(i==k)*vectors[:,j]+(j==k)*vectors[:,i])*radial[2,0]
                       +8*vectors[:,i]*vectors[:,j]*vectors[:,k]*radial[3,0])
        values.append(value)
    return np.asarray(values).T


def traceless_third(tensor):
    """Symmetric traceless projection, rotationally covariant in 3D."""
    raw=np.asarray(tensor,dtype=float)
    trace=np.einsum("...ijj->...i",raw)
    result=raw.copy()
    for i,j,k in INDICES:
        result[...,i,j,k]-=((j==k)*trace[...,i]+(i==k)*trace[...,j]+(i==j)*trace[...,k])/5
    return result


def _vector_fourier(d,vectors,kappa,normal_order):
    q=np.sqrt(kappa*kappa+np.sum(vectors*vectors,axis=1))
    first=_radial(d,q,kappa,1,normal_order,0)
    normal=_radial(d,q,kappa,0,normal_order,1)
    return np.column_stack((2j*vectors[:,0]*first,2j*vectors[:,1]*first,normal))


def traceless_second(tensor):
    raw=np.asarray(tensor,dtype=float)
    return raw-np.trace(raw,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3


def _second_fourier(d,vectors,kappa,normal_order):
    q=np.sqrt(kappa*kappa+np.sum(vectors*vectors,axis=1))
    h1=_radial(d,q,kappa,1,normal_order,0)
    h2=_radial(d,q,kappa,2,normal_order,0)
    dh1=_radial(d,q,kappa,1,normal_order,1)
    d2h0=_radial(d,q,kappa,0,normal_order,2)
    out=np.empty((len(vectors),3,3),dtype=complex)
    for i in range(2):
        for j in range(2):
            out[:,i,j]=-(2*(i==j)*h1+4*vectors[:,i]*vectors[:,j]*h2)
        out[:,i,2]=out[:,2,i]=2j*vectors[:,i]*dh1
    out[:,2,2]=d2h0
    return out.reshape(len(vectors),9)


def plane_third_moment(d,delta,*,geometry,kappa,amplitude,direction,
                       tolerance=2e-11,max_shell_index=32,rank=3):
    """Six derivatives of the STF moment; explicit rank=1,2,3 research APIs.

    The historical function name is kept for the already tested rank-3 API.
    Rank 2 has cubic, not general affine, zero-bulk symmetry.
    """
    if d<=0 or kappa<=0 or tolerance<=0:
        raise ValueError("positive separation, decay, tolerance required")
    if rank not in (1,2,3):
        raise ValueError("only derived moments rank 1, 2 or 3 are supported")
    transform={1:_vector_fourier,2:_second_fourier,3:_tensor_fourier}[rank]
    total=np.zeros((6,3**rank)); small=0; last=math.inf
    shells=[np.zeros((1,2)),*(shell.vectors for shell in
              triangular_reciprocal_shells(geometry,max_shell_index))]
    for count,vectors in enumerate(shells):
        phase=np.exp(1j*(vectors@delta)); g=vectors@direction
        c0=transform(d,vectors,kappa,0)
        c1=transform(d,vectors,kappa,1)
        c2=transform(d,vectors,kappa,2)
        values=np.array([np.real(np.sum(c*factor[:,None],axis=0)) for c,factor in
                         ((c0,phase),(c1,phase),(c0,1j*g*phase),
                          (c2,phase),(c1,1j*g*phase),(c0,-g*g*phase))])
        total+=values
        # Absolute mode-coefficient bound, NOT cancellation of an odd shell.
        last=float(sum(np.sum(np.abs(c)*factor[:,None],axis=0).max() for c,factor in
                       ((c0,np.maximum(1,g*g)),(c1,np.maximum(1,np.abs(g))),
                        (c2,np.ones(len(g))))))*abs(amplitude)/geometry.atomic_cell_area
        if count>0 and last<tolerance*max(1.,np.max(np.abs(total))*abs(amplitude)/geometry.atomic_cell_area):
            small+=1
            if small>=3:
                scaled=total*amplitude/geometry.atomic_cell_area
                if rank==3:
                    tensor=traceless_third(scaled.reshape(6,3,3,3))
                elif rank==2:
                    tensor=traceless_second(scaled.reshape(6,3,3))
                else:
                    tensor=scaled
                return tensor,count,last
        else:
            small=0
    raise RuntimeError("environment-moment reciprocal shells did not converge")


class AngularInterfaceInvariant:
    """Il=2 sum_depth ||STF(moment_depth)||^2 for one active interface.

    The normalized radial kernel is exp[-k(r/L0-1)]/rho_ref, with rho_ref fixed
    at target bulk. The Cartesian moment uses reduced r/L0. Its odd parity
    makes rank 1/3 EXACTLY zero in affine centrosymmetric FCC bulk. Rank 2 is
    zero in isotropically scaled CUBIC FCC only; this interface implementation
    rejects a noncubic bulk reference rather than dropping its background Q.
    A fault/surface changes each local environment, not one global tensor.
    """
    def __init__(self,interface,*,tolerance=2e-10,max_layers=128,angular_decay=None,rank=3):
        self.interface=interface
        self.geometry=interface.geometry
        self.h=interface.h; self.direction=interface.direction
        if rank not in (1,2,3):
            raise ValueError("only moment ranks 1, 2 and 3 derived")
        if rank==2 and not np.isclose(self.h/self.geometry.b,math.sqrt(2/3),rtol=1e-10,atol=0.):
            raise ValueError("rank-2 interface difference requires cubic FCC bulk; nonzero background Q not implemented")
        self.rank=rank
        self.kappa=interface.bulk.density_params.kappa
        self.amplitude=interface.bulk.density_params.C_rho/interface.bulk.embedding.rho_ref
        if angular_decay is not None:
            # A separate radial range is an additional RESEARCH parameter,
            # normalized once at target b=1, never at each deformed state.
            from .aluminum_full_fcc_calibration import FullFCCParameters, build_full_fcc_calibration_model
            probe=build_full_fcc_calibration_model(FullFCCParameters(.1,1.,float(angular_decay),1.))
            self.kappa=float(angular_decay)
            self.amplitude=math.exp(self.kappa)/probe.embedding.rho_ref
        self.tolerance=tolerance; self.max_layers=max_layers
        self._base={}

    def _plane(self,d,delta):
        return plane_third_moment(d,delta,geometry=self.geometry,kappa=self.kappa,
            amplitude=self.amplitude,direction=self.direction,tolerance=self.tolerance*.1,rank=self.rank)

    def evaluate(self,a,s):
        terms=[]; small=0; maximum_shells=0
        for k in range(1,self.max_layers+1):
            if k not in self._base:
                self._base[k]=self._plane(k*self.h,self.interface._baseline_delta(k))[0][0]
            active,shells,tail=self._plane(a+(k-1)*self.h,self.interface._active_delta(k,s))
            active[0]-=self._base[k]
            terms.append(active); maximum_shells=max(maximum_shells,shells)
            bound=float(np.max(np.abs(active)))
            if bound<self.tolerance:
                small+=1
                if small>=4:
                    break
            else:
                small=0
        else:
            raise RuntimeError("angular interface neighborhood did not converge")
        depths=np.cumsum(np.asarray(terms)[::-1],axis=0)[::-1].reshape(k,6,-1)
        t,ta,ts,taa,tas,tss=(depths[:,i] for i in range(6))
        dot=lambda x,y:float(np.sum(x*y))
        result=2*np.array([dot(t,t),2*dot(t,ta),2*dot(t,ts),
                          2*(dot(ta,ta)+dot(t,taa)),2*(dot(ta,ts)+dot(t,tas)),
                          2*(dot(ts,ts)+dot(t,tss))])
        return result,dict(layers_used=k,shells_used=maximum_shells,last_layer_proxy=bound)

    def direct_value(self,a,s,*,radius=24,layers=24):
        """Independent real-space invariant; never used in a model generator."""
        m,n=np.meshgrid(np.arange(-radius,radius+1),np.arange(-radius,radius+1),indexing="ij")
        rxy=m.ravel()[:,None]*self.geometry.a1+n.ravel()[:,None]*self.geometry.a2
        def plane(d,delta):
            xy=rxy+delta
            r=np.column_stack((xy,np.full(len(xy),d)))
            weights=self.amplitude*np.exp(-self.kappa*np.linalg.norm(r,axis=1))
            if self.rank==1:
                return np.einsum("n,ni->i",weights,r)
            if self.rank==2:
                return traceless_second(np.einsum("n,ni,nj->ij",weights,r,r))
            raw=np.einsum("n,ni,nj,nk->ijk",weights,r,r,r)
            return traceless_third(raw)
        changes=[plane(a+(k-1)*self.h,self.interface._active_delta(k,s))-
                 plane(k*self.h,self.interface._baseline_delta(k)) for k in range(1,layers+1)]
        depths=np.cumsum(np.asarray(changes)[::-1],axis=0)[::-1]
        return float(2*np.sum(depths*depths))


class AngularInterfaceResearchSurface:
    """W_scalar+D3 I3+D1 I1+D2 I2; absent from production registry.

    Companion coefficients default to zero. A signed coefficient is not a
    material-stability certification or a change of probability theory.
    """
    def __init__(self,interface,amplitude_ev,*,tolerance=2e-10,angular_decay=None,
                 vector_amplitude_ev=0.,quadrupole_amplitude_ev=0.):
        if not np.isfinite(amplitude_ev) or amplitude_ev<0:
            raise ValueError("nonnegative finite angular energy amplitude required")
        self.interface=interface; self.amplitude_ev=float(amplitude_ev)
        if not np.isfinite(vector_amplitude_ev):
            raise ValueError("vector moment amplitude must be finite; total stability is checked separately")
        self.vector_amplitude_ev=float(vector_amplitude_ev)
        self.angular=AngularInterfaceInvariant(interface,tolerance=tolerance,angular_decay=angular_decay)
        self.vector=(AngularInterfaceInvariant(interface,tolerance=tolerance,angular_decay=angular_decay,rank=1)
                     if vector_amplitude_ev else None)
        if not np.isfinite(quadrupole_amplitude_ev):
            raise ValueError("quadrupole coefficient must be finite; total stability checked separately")
        self.quadrupole_amplitude_ev=float(quadrupole_amplitude_ev)
        self.quadrupole=(AngularInterfaceInvariant(interface,tolerance=tolerance,angular_decay=angular_decay,rank=2)
                         if quadrupole_amplitude_ev else None)
        self.h=interface.h; self.period=interface.period

    def packed(self,a,s):
        v=self.interface.evaluate(a,s)
        scalar=np.array([v.energy,v.d_da,v.d_ds,v.d2_daa,v.d2_das,v.d2_dss])
        total=scalar+self.amplitude_ev*self.angular.evaluate(a,s)[0]
        if self.vector is not None:
            total+=self.vector_amplitude_ev*self.vector.evaluate(a,s)[0]
        if self.quadrupole is not None:
            total+=self.quadrupole_amplitude_ev*self.quadrupole.evaluate(a,s)[0]
        return total
