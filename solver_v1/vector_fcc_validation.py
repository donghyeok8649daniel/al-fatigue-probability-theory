"""Independent finite real-atom references; NEVER the canonical energy path."""
from itertools import product

import numpy as np

from .nonlinear_fcc_screw import stf_basis


def direct_row_channels(rows,vector,*,images=384,include_odd_quadrupole=False):
    """Value/gradient/Hessian from direct atom distances and product rules.

    LJ tail needs image-count refinement. Electronic-density tails decay
    exponentially. This shares neither reciprocal coefficients nor their
    derivatives with VectorRowKernel; only the unchanged physical parameters
    and orthonormal STF coordinate convention are shared.
    """
    x,y,z=np.asarray(vector,float)
    rvec=np.column_stack([rows.b*np.arange(-images,images+1)+x,
                         np.full(2*images+1,y),np.full(2*images+1,z)])
    radius=np.linalg.norm(rvec,axis=1);unit=rvec/radius[:,None]
    outer=unit[:,:,None]*unit[:,None,:];identity=np.eye(3)[None,:,:]
    p=rows.bulk.p;t12=(p.sigma_lj/radius)**12;t6=(p.sigma_lj/radius)**6
    pair=4*p.epsilon*(t12-t6)
    first=4*p.epsilon*(-12*t12+6*t6)/radius
    second=4*p.epsilon*(156*t12-42*t6)/radius**2
    grad=first[:,None]*unit
    hess=second[:,None,None]*outer+(first/radius)[:,None,None]*(identity-outer)
    def exponential(amplitude,kappa):
        value=amplitude*np.exp(-kappa*radius)
        derivative=-kappa*value[:,None]*unit
        curvature=value[:,None,None]*((kappa*kappa+kappa/radius)[:,None,None]*outer
                                     -(kappa/radius)[:,None,None]*identity)
        return value,derivative,curvature
    den=rows.bulk.density_params;density,dgrad,dhess=exponential(den.C_rho,den.kappa)
    values=[pair.sum(),density.sum()];gradients=[grad.sum(axis=0),dgrad.sum(axis=0)]
    hessians=[hess.sum(axis=0),dhess.sum(axis=0)]
    ranks=[(rank, getattr(rows.surface, name, None) or rows.surface.angular)
           for rank,name in ((1,'vector'),(2,'quadrupole'),(3,'angular'))]
    if include_odd_quadrupole:
        ranks.append((2,rows.surface.angular))
    for rank,angular in ranks:
        weight,wgrad,whess=exponential(angular.amplitude,angular.kappa)
        raw=[];rawgrad=[];rawhess=[]
        for indices in product(range(3),repeat=rank):
            powers=np.array([indices.count(axis) for axis in range(3)])
            mon=np.prod(rvec**powers,axis=1)
            mg=np.zeros((len(radius),3));mh=np.zeros((len(radius),3,3))
            for i in range(3):
                if powers[i]:
                    lower=powers-np.eye(3,dtype=int)[i]
                    mg[:,i]=powers[i]*np.prod(rvec**lower,axis=1)
                    for j in range(3):
                        if lower[j]:
                            mh[:,i,j]=powers[i]*lower[j]*np.prod(rvec**(lower-np.eye(3,dtype=int)[j]),axis=1)
            raw.append(np.sum(weight*mon))
            rawgrad.append(np.sum(wgrad*mon[:,None]+weight[:,None]*mg,axis=0))
            rawhess.append(np.sum(whess*mon[:,None,None]+wgrad[:,:,None]*mg[:,None,:]
                                 +mg[:,:,None]*wgrad[:,None,:]+weight[:,None,None]*mh,axis=0))
        basis=stf_basis(rank)
        values.extend(np.array(raw)@basis)
        gradients.extend(np.einsum('ti,tc->ci',rawgrad,basis))
        hessians.extend(np.einsum('tij,tc->cij',rawhess,basis))
    return dict(value=np.array(values),gradient=np.array(gradients),hessian=np.array(hessians))
