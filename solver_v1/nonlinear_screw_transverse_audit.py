"""Forces omitted by the fixed-transverse nonlinear screw subspace.

This is a diagnostic at the computed anti-plane states, NOT a vector-core
solver. Keep the zero reciprocal mode: unlike x slip, y/z derivatives do not
annihilate it. Infinite atomic rows are still summed analytically by Bessel.
"""
from itertools import product

import numpy as np
from scipy.special import kv

from .discrete_fcc_screw import exponential_row_transform_derivatives
from .lattice_bessel import _mode_coefficient_and_deda, _zero_mode
from .nonlinear_fcc_screw import stf_basis


def exponential_transform_radial_derivatives(radius,g,kappa):
    """d/dA of (H,H_G,H_GG,H_GGG), using d(A^nu K_nu)/dA."""
    q=np.sqrt(kappa*kappa+g*g)
    d=[-2*radius*kappa*(-radius/2)**j*q**(-j)*kv(j,radius*q) for j in range(4)]
    return np.stack([d[0],2*g*d[1],2*d[1]+4*g*g*d[2],12*g*d[2]+8*g**3*d[3]],axis=-1)


def transverse_row_coefficients(rows,j,l,mode,*,y=None,z=None):
    """Channel value and y/z derivative at fixed x-phase (mode 0 included).

    Coefficients include the reference ABC x phase. Their last derivative
    axis has length 2. This does NOT differentiate affine x shear when moving
    an atom in physical y/z: its x position is held fixed for that force.
    """
    j,l=np.broadcast_arrays(np.atleast_1d(j),np.atleast_1d(l))
    y=rows.d*(j+l/3) if y is None else np.broadcast_to(y,j.shape)
    z=rows.h*l if z is None else np.broadcast_to(z,j.shape)
    radius=np.hypot(y,z)
    if np.any(radius<=0) or mode<0:
        raise ValueError('nonzero transverse row and nonnegative mode required')
    g=2*np.pi*mode/rows.b
    power=lambda p: _zero_mode(radius,p,rows.b) if mode==0 else _mode_coefficient_and_deda(radius,p,mode,rows.b)
    c6,dc6=power(6);c3,dc3=power(3);p=rows.bulk.p
    pair=4*p.epsilon*(p.sigma_lj**12*c6-p.sigma_lj**6*c3)
    dpair=4*p.epsilon*(p.sigma_lj**12*dc6-p.sigma_lj**6*dc3)
    factor=(1 if mode==0 else 2)/rows.b
    density=rows.bulk.density_params
    den=factor*density.C_rho*exponential_row_transform_derivatives(radius,g,density.kappa)[0]
    dden=factor*density.C_rho*exponential_transform_radial_derivatives(radius,g,density.kappa)[...,0]
    radial=np.column_stack([y/radius,z/radius])
    values=[pair[:,None],den[:,None]]
    derivatives=[dpair[:,None,None]*radial[:,None,:],dden[:,None,None]*radial[:,None,:]]
    angular=rows.surface.angular
    h=np.stack(exponential_row_transform_derivatives(radius,g,angular.kappa),axis=-1)
    dh=exponential_transform_radial_derivatives(radius,g,angular.kappa)
    for rank in (1,2,3):
        raw=[];dy=[];dz=[]
        for indices in product(range(3),repeat=rank):
            nx,ny,nz=(indices.count(axis) for axis in range(3))
            common=factor*angular.amplitude*1j**nx
            raw.append(common*h[...,nx]*y**ny*z**nz)
            dy.append(common*(dh[...,nx]*(y/radius)*y**ny*z**nz+
                      (ny*h[...,nx]*y**(ny-1)*z**nz if ny else 0)))
            dz.append(common*(dh[...,nx]*(z/radius)*y**ny*z**nz+
                      (nz*h[...,nx]*y**ny*z**(nz-1) if nz else 0)))
        values.append(np.array(raw).T@stf_basis(rank))
        derivatives.append(np.stack([np.array(dy).T@stf_basis(rank),np.array(dz).T@stf_basis(rank)],axis=-1))
    phase=((-1.)**(mode*(j+l)))
    return np.concatenate(values,axis=1)*phase[:,None],np.concatenate(derivatives,axis=1)*phase[:,None,None]


def omitted_transverse_forces(cell,field,gamma=0.,*,ring=14):
    """Gradient in the two frozen physical coordinates, eV/L0 per line repeat.

    A nonzero result rejects full vector force equilibrium. Refinement of the
    transverse validation ring is required; it is NOT an acceptance threshold.
    """
    rows=cell.rows
    j,l=np.meshgrid(np.arange(-ring,ring+1),np.arange(-ring,ring+1),indexing='ij')
    keep=(j!=0)|(l!=0);j,l=j[keep],l[keep]
    current=cell.evaluate(field,gamma,fields=True)
    weights=np.empty_like(current['channels']);weights[...,0]=.5
    weights[...,1]=rows.embedding.first_derivative(current['density'])
    weights[...,2:]=2*rows.angular_weights*current['channels'][...,2:]
    gradient=np.zeros(cell.shape+(2,))
    reverse=np.ix_((-np.arange(cell.shape[0]))%cell.shape[0],
                   (-np.arange(cell.shape[1]))%cell.shape[1])
    for mode in range(len(rows.g)+1):
        g=2*np.pi*mode/rows.b
        _,derivatives=transverse_row_coefficients(rows,j,l,mode)
        derivatives*=np.exp(1j*g*gamma*rows.h*l)[:,None,None]
        stencil=np.zeros(cell.shape+(17,2),complex)
        np.add.at(stencil,((-j)%cell.shape[0],(-l)%cell.shape[1]),derivatives)
        operator=np.fft.fft2(stencil,axes=(0,1))
        phase=np.exp(1j*g*field)[...,None,None]
        w=weights[...,None]
        direct=cell._apply(operator,phase)
        adjoint=cell._apply(operator[reverse],w*phase.conj())
        gradient+=(phase*adjoint-w*phase.conj()*direct).real.sum(axis=-2)
    return dict(gradient=gradient,maximum_y_gradient=float(np.max(abs(gradient[...,0]))),
                maximum_z_gradient=float(np.max(abs(gradient[...,1]))),
                net_gradient=np.sum(gradient,axis=(0,1)),ring=ring,
                full_vector_equilibrium_certified=False)
