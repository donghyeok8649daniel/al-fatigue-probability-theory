"""Exact fixed-shape coefficient decomposition at a DECLARED atomic state.

Energy/force columns are not a surrogate and no state is relaxed here. The
infinite LJ/Bessel pair remains u/r^12-v/r^6. Source core forces can constrain
these coefficients without using experimental yield or lowering a barrier.
"""
import numpy as np

from .coordination_screening import screening_factor
from .vector_fcc_rows import power_mode_radial


NAMES = ('u','v','A','B','C','D3','D1','D2','D_Eg','K3')


def power_row_value_gradient(vectors, power, b, *, tolerance=2e-13, max_modes=40):
    vectors = np.asarray(vectors,float)
    if vectors.ndim!=2 or vectors.shape[1]!=3 or np.any(~np.isfinite(vectors)) or power not in (3,6):
        raise ValueError('finite row vectors and LJ power3/6 required')
    if not np.isfinite(tolerance) or tolerance<=0 or not np.isfinite(b) or b<=0:
        raise ValueError('positive numerical tolerance and row period required')
    x,y,z = vectors.T; radius=np.hypot(y,z)
    if np.any(radius<=0):
        raise ValueError('separate coincident-row limit required')
    x=x-b*np.floor(x/b+.5)
    value=np.zeros(len(x)); gradient=np.zeros_like(vectors); small=0
    for mode in range(max_modes+1):
        C,Ca,_ = power_mode_radial(radius,power,mode,b)
        g=2*np.pi*mode/b; phase=g*x
        value+=C*np.cos(phase)
        gradient[:,0]-=g*C*np.sin(phase)
        gradient[:,1]+=Ca*y/radius*np.cos(phase)
        gradient[:,2]+=Ca*z/radius*np.cos(phase)
        bound=float(np.max(np.maximum.reduce([abs(C),abs(g*C),abs(Ca)])))
        small=small+1 if mode and bound<tolerance else 0
        if small>=2:
            return value,gradient
    raise ArithmeticError('coefficient pair Bessel sum did not converge')


def site_nonpair_basis(law,channels):
    z=np.asarray(channels,float)
    rho=law.rho_bulk+z[:,1]; x=rho/law.rho_ref
    if np.any(x<=0):
        raise ValueError('positive actual site density required')
    E=np.zeros((len(z),8)); W=np.zeros((len(z),22,8))
    E[:,:3]=np.column_stack([1-np.sqrt(x),x-1,(x-1)**2])
    W[:,1,:3]=np.column_stack([-.5/np.sqrt(x),np.ones(len(x)),2*(x-1)])/law.rho_ref
    def invariant(name):
        B=law.maps[name]; Q=z@B.T
        return np.sum(Q*Q,axis=1),2*Q@B
    I3,d3=invariant('Q3'); I1,d1=invariant('Q1')
    E[:,3]=I3; W[:,:,3]=d3
    g,gx,_=screening_factor(x,law.model.screening,law=law.model.law)
    E[:,4]=g*I1; W[:,:,4]=g[:,None]*d1
    W[:,1,4]+=gx*I1/law.rho_ref
    for index,name,norm in ((5,'Q2',law.model.even_reference_curvature),(6,'QE',1.)):
        I,dI=invariant(name); den=1+law.model.alpha*I/norm
        E[:,index]=I/den; W[:,:,index]=dI/den[:,None]**2
    E[:,7]=I3*I3; W[:,:,7]=2*I3[:,None]*d3
    return E,W


def current_core_coefficients(core,field):
    """Return energy[10], free-row gradient[n,3,10], and site environments."""
    full=core.full_field(field); ns,nr=core.destination.shape
    vectors=core.reference[None]+full[core.destination]-full[core.energy_ids,None]
    flat=vectors.reshape(-1,3)
    values=np.empty((len(flat),22)); gradients=np.empty((len(flat),22,3))
    for begin in range(0,len(flat),core.chunk_size):
        end=min(begin+core.chunk_size,len(flat)); out=core.kernel.evaluate(flat[begin:end])
        values[begin:end]=out['value']; gradients[begin:end]=out['gradient']
    channels=(values.reshape(ns,nr,22)-core.reference_channels).sum(axis=1)
    E,W=site_nonpair_basis(core.site_law,channels)
    local=np.einsum('ncp,nrci->nrip',W,gradients.reshape(ns,nr,22,3))
    all_gradient=core._scatter(local).reshape(len(core.indices),3,8)
    all_gradient[core.energy_ids]-=local.sum(axis=1)
    energy=np.zeros(10); energy[2:]=E.sum(axis=0)
    force=np.zeros((len(core.free_ids),3,10)); force[:,:,2:]=all_gradient[core.free_ids]
    for index,power,sign in ((0,6,1.),(1,3,-1.)):
        base,_=power_row_value_gradient(core.reference,power,core.rows.b)
        actual,grad=power_row_value_gradient(flat,power,core.rows.b)
        energy[index]=sign*.5*np.sum(actual.reshape(ns,nr)-base)
        bonds=sign*.5*grad.reshape(ns,nr,3)
        assembled=core._scatter(bonds)
        assembled[core.energy_ids]-=bonds.sum(axis=1)
        force[:,:,index]=assembled[core.free_ids]
    return dict(energy=energy,gradient=force,channels=channels,
                coefficient_names=NAMES,states_relaxed=False)
