"""Independent finite-wavevector STATIC Hessian validation, without mass/time.

Direct neighbor radii here are a refinement-controlled validation of an
infinite analytic potential, NOT a new cutoff energy or PDE generator.
Positive uniform cubic moduli alone do not guarantee finite-q stability.
"""
from __future__ import annotations

import math
import numpy as np

from .angular_environment_reference import traceless_third, traceless_second
from .aluminum_full_fcc_calibration import FullFCCParameters, build_full_fcc_calibration_model


def neighbors_in_plane_frame(bulk,cutoff):
    geometry=bulk.geometry; h=bulk.a0
    extent=int(np.ceil(2*cutoff/geometry.b))+2
    m,n=np.meshgrid(np.arange(-extent,extent+1),np.arange(-extent,extent+1),indexing="ij")
    xy=m.ravel()[:,None]*geometry.a1+n.ravel()[:,None]*geometry.a2
    arrays=[]
    for layer in range(-int(np.ceil(cutoff/h)),int(np.ceil(cutoff/h))+1):
        planar=xy+geometry.abc_shift(layer)
        points=np.column_stack((planar,np.full(len(planar),layer*h)))
        length=np.linalg.norm(points,axis=1)
        selected=(length>1e-10)&(length<cutoff*(1+1e-12))
        arrays.append(points[selected])
    return np.vstack(arrays)


class StaticBulkHessian:
    """Unweighted force-constant matrix K(q), units eV/L0^2, NOT Hz."""
    def __init__(self,bulk,*,cutoff=8.,D3=0.,D1=0.,D2=0.,angular_decay=None):
        self.bulk=bulk; self.cutoff=float(cutoff)
        self.R=neighbors_in_plane_frame(bulk,cutoff)
        self.r=np.linalg.norm(self.R,axis=1); direction=self.R/self.r[:,None]
        rr=direction[:,:,None]*direction[:,None,:]; eye=np.eye(3)
        p=bulk.p; r=self.r
        pair_prime=4*p.epsilon*(-12*p.sigma_lj**12/r**13+6*p.sigma_lj**6/r**7)
        pair_second=4*p.epsilon*(156*p.sigma_lj**12/r**14-42*p.sigma_lj**6/r**8)
        density=bulk.density_params
        rho=density.C_rho*np.exp(-density.kappa*r)
        rho_prime=-density.kappa*rho; rho_second=density.kappa**2*rho
        rho_bulk=bulk.full_environment_density(bulk.a0,0).value
        first=float(bulk.embedding.first_derivative(rho_bulk))
        self.second=float(bulk.embedding.second_derivative(rho_bulk))
        effective_prime=pair_prime+2*first*rho_prime
        effective_second=pair_second+2*first*rho_second
        self.pair_hessian=effective_second[:,None,None]*rr+(effective_prime/r)[:,None,None]*(eye-rr)
        self.density_gradient=rho_prime[:,None]*direction
        self.D3=float(D3); self.D1=float(D1); self.D2=float(D2)
        self.grad3=None; self.grad1=None; self.grad2=None
        if D3 or D1 or D2:
            k=float(angular_decay or density.kappa)
            probe=build_full_fcc_calibration_model(FullFCCParameters(.1,1.,k,1.))
            amplitude=math.exp(k)/probe.embedding.rho_ref
            weight=amplitude*np.exp(-k*r)
            self.grad1=weight[:,None,None]*(eye-k*self.R[:,:,None]*direction[:,None,:])
            if D2:
                raw=np.empty((len(r),3,3,3))
                for b in range(3):
                    for i in range(3):
                        for j in range(3):
                            raw[:,b,i,j]=weight*((i==b)*self.R[:,j]+(j==b)*self.R[:,i]
                                -k*direction[:,b]*self.R[:,i]*self.R[:,j])
                self.grad2=traceless_second(raw).reshape(len(r),3,9)
            if D3:
                # Leading axes (atom, differentiated coordinate), trailing
                # three tensor axes receive the traceless projection.
                raw=np.empty((len(r),3,3,3,3))
                for b in range(3):
                    for i in range(3):
                        for j in range(3):
                            for l in range(3):
                                raw[:,b,i,j,l]=weight*((i==b)*self.R[:,j]*self.R[:,l]+
                                    (j==b)*self.R[:,i]*self.R[:,l]+(l==b)*self.R[:,i]*self.R[:,j]
                                    -k*direction[:,b]*self.R[:,i]*self.R[:,j]*self.R[:,l])
                self.grad3=traceless_third(raw).reshape(len(r),3,27)

    def evaluate(self,q_plane_frame):
        phase=self.R@np.asarray(q_plane_frame,dtype=float)
        cosine=np.cos(phase); sine=np.sin(phase)
        # Pair redistribution F' is exact; the remaining F'' contribution is
        # the square of the collective scalar-density variation.
        pair=np.einsum("n,nij->ij",1-cosine,self.pair_hessian)
        g=np.einsum("n,ni->i",sine,self.density_gradient)
        scalar=self.second*np.outer(g,g)
        angular=np.zeros((3,3))
        if self.grad3 is not None:
            L3=np.einsum("n,nij->ij",cosine-1,self.grad3)
            angular+=2*self.D3*(L3@L3.T)
        if self.D1:
            L1=np.einsum("n,nij->ij",cosine-1,self.grad1)
            angular+=2*self.D1*(L1.T@L1)
        if self.grad2 is not None:
            # Even environment moment -> odd bond derivative. The Fourier
            # variation is imaginary (sin), not the odd moment's cos-1.
            L2=np.einsum("n,nij->ij",sine,self.grad2)
            angular+=2*self.D2*(L2@L2.T)
        total=pair+scalar+angular
        return dict(matrix=total,eigenvalues=np.linalg.eigvalsh(total),
                    pair=pair,scalar_density=scalar,angular=angular)

    def crystallographic_wavevector(self,fractional_cubic):
        q=2*math.pi/self.bulk.geometry.lattice_constant*np.asarray(fractional_cubic)
        geometry=self.bulk.geometry
        return np.array([q@geometry.e1,q@geometry.e2,q@geometry.e3])
