"""Three identifiable SW amplitudes, with exact conservative force features.

The original radial shapes and angular preference are fixed. epsilon is fixed
to remove its scale gauge; independent repulsion, attraction and angle weights
map back to A, B and lambda. This is a research fit, never a production model.
Lengths are Angstrom; energies eV; forces eV/Angstrom.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import lsq_linear
from .silicon_environment_research import _parameters


def amplitude_parameters(parameters, amplitudes):
    p = _parameters(parameters)
    a = np.asarray(amplitudes, float)
    if a.shape != (3,) or not np.all(np.isfinite(a)) or np.any(a <= 0):
        raise ValueError('three strictly positive SW amplitudes required')
    p['A'] *= a[1]
    p['B'] *= a[0]/a[1]
    p['lambda'] *= a[2]
    return p


def force_features(positions, cell, periodic, parameters):
    """Return E[3] and F[atom,Cartesian,3] for original SW components.

    Periodic images are explicit, including nonzero images of the same atom.
    Sum the three channels to recover original SW; no minimum-image shortcut.
    The forces differentiate the same energy, including all centered angles.
    """
    from ase.neighborlist import primitive_neighbor_list
    p = _parameters(parameters)
    x, cell, pbc = np.asarray(positions,float), np.asarray(cell,float), np.asarray(periodic)
    if (x.ndim != 2 or x.shape[1] != 3 or not len(x)
            or cell.shape != (3,3) or pbc.shape != (3,) or pbc.dtype != np.dtype(bool)
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(cell))
            or abs(np.linalg.det(cell)) < 1e-10):
        raise ValueError('finite nonempty positions, nonsingular cell and bool PBC required')
    rc = p['a']*p['sigma']
    i,j,r,v = primitive_neighbor_list('ijdD',pbc,cell,x,rc,self_interaction=False)
    if np.any(r <= 1e-12):
        raise ValueError('coincident atoms in SW fit geometry')
    order=np.argsort(i,kind='stable');i,j,r,v=[a[order] for a in (i,j,r,v)]
    n=v/r[:,None]
    z=p['sigma']/(r-rc)
    pair=np.stack((p['A']*p['epsilon']*p['B']*(p['sigma']/r)**p['p']*np.exp(z),
                   -p['A']*p['epsilon']*(p['sigma']/r)**p['q']*np.exp(z)),axis=1)
    derivative=pair*(-np.array([p['p'],p['q']])[None,:]/r[:,None]
                     -p['sigma']/(r-rc)[:,None]**2)
    energy=np.r_[.5*np.sum(pair,axis=0),0.]
    force=np.zeros((len(x),3,3))
    for channel in range(2):
        gradient=.5*derivative[:,channel,None]*n
        np.add.at(force[:,:,channel],i,gradient)
        np.add.at(force[:,:,channel],j,-gradient)
    bounds=np.searchsorted(i,np.arange(len(x)+1))
    for center in range(len(x)):
        a,b=bounds[center:center+2]
        if b-a < 2:
            continue
        left,right=np.triu_indices(b-a,k=1);left+=a;right+=a
        cosine=np.einsum('ij,ij->i',n[left],n[right])
        delta=cosine-p['costheta0']
        weight=p['lambda']*p['epsilon']*np.exp(p['gamma']*(z[left]+z[right]))
        energy[2]+=np.sum(weight*delta**2)
        gl=weight[:,None]*(-p['gamma']*p['sigma']/(r[left]-rc)[:,None]**2
              *delta[:,None]**2*n[left]
              +2*delta[:,None]*(n[right]-cosine[:,None]*n[left])/r[left,None])
        gr=weight[:,None]*(-p['gamma']*p['sigma']/(r[right]-rc)[:,None]**2
              *delta[:,None]**2*n[right]
              +2*delta[:,None]*(n[left]-cosine[:,None]*n[right])/r[right,None])
        force[center,:,2]+=np.sum(gl+gr,axis=0)
        np.add.at(force[:,:,2],j[left],-gl)
        np.add.at(force[:,:,2],j[right],-gr)
    if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(force)):
        raise FloatingPointError('nonfinite SW features')
    return energy,force


def fit_amplitudes(design,target):
    """Nonnegative convex least squares, with rank and KKT diagnostics.

    Boundary solutions are returned as diagnostics; they cannot be mapped to
    an admissible positive SW candidate. No ridge term or clipped coefficients.
    """
    x,y=np.asarray(design,float),np.asarray(target,float)
    if x.ndim != 2 or x.shape[1] != 3 or y.shape != (len(x),):
        raise ValueError('three-channel design and matching target required')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('finite fit arrays required')
    singular=np.linalg.svd(x,compute_uv=False)
    if len(singular)!=3 or singular[-1] <= singular[0]*1e-12:
        raise ValueError('rank-deficient SW amplitudes; no regularized fit')
    result=lsq_linear(x,y,bounds=(0,np.inf),method='bvls',tol=1e-12)
    if not result.success:
        raise RuntimeError('amplitude fit did not converge')
    residual=x@result.x-y
    gradient=x.T@residual
    kkt=np.where(result.active_mask == -1,np.minimum(gradient,0),gradient)
    return result.x,dict(amplitudes=result.x.tolist(),rank=3,
        singular_values=singular.tolist(),condition_number=float(singular[0]/singular[-1]),
        active_mask=result.active_mask.tolist(),objective_squared=float(residual@residual),
        original_objective_squared=float(np.sum((x@np.ones(3)-y)**2)),
        kkt_max_abs=float(np.max(abs(kkt))),iterations=int(result.nit),
        regularization=False,epsilon_gauge_fixed=True)
