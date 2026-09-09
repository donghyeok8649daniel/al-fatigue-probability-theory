"""First density/angular mixed invariant; one-coefficient research ablation.

With x_i the normalized scalar environment and I_i=||Q3_i||^2, the next mixed
Taylor term in a differentiable per-atom F(x,I) is E_xI (x_i-1) I_i.
This is NOT a replacement LJ pair. Full infinite neighbors are summed before
the site invariant/product. Both x-1 and Q3 vanish in pristine cubic bulk;
the term has no harmonic bulk/finite-q/pristine-interface contribution.

Require [[C,E_xI/2],[E_xI/2,K3]] >=0 for the existing C(x-1)^2+K3 I^2
sector. Thus the added cross term cannot make THIS nonlinear sector negative.
This is not a proof of global stability of all other signed angular sectors.
"""
from functools import lru_cache

import numpy as np

from .quartic_angular_material import (
    QuarticSymmetryInterface,QuarticSymmetryObservationCache,QuarticSymmetryTailBulkBasis,
    QuarticMoment,
)
from .range_resolved_material import build_range_surface
from .vector_interface_reference import VectorInterfaceEvaluation,HESSIAN_INDICES
from .tail_constrained_material import spectral_profile


def cross_site_jet(density_changes, tensor_jets):
    """Sum (x_i-1)||Q3_i||^2 and its full three-coordinate analytic jet."""
    density_changes=np.asarray(density_changes,float);tensor_jets=np.asarray(tensor_jets,float)
    if (density_changes.ndim!=2 or density_changes.shape[1]!=10 or tensor_jets.ndim!=3
            or tensor_jets.shape[1]!=10 or np.any(~np.isfinite(density_changes))
            or np.any(~np.isfinite(tensor_jets))):
        raise ValueError('finite scalar/tensor per-site jets required')
    n=max(len(density_changes),len(tensor_jets))
    y=np.zeros((n,10));y[:len(density_changes)]=density_changes
    t=np.zeros((n,10,tensor_jets.shape[2]));t[:len(tensor_jets)]=tensor_jets
    Q=t[:,0];Qg=t[:,1:4];Qh=t[:,HESSIAN_INDICES]
    I=np.einsum('nc,nc->n',Q,Q)
    Ig=2*np.einsum('nc,nic->ni',Q,Qg)
    Ih=2*(np.einsum('nic,njc->nij',Qg,Qg)+np.einsum('nc,nijc->nij',Q,Qh))
    gradient=np.einsum('ni,n->i',y[:,1:4],I)+np.einsum('n,ni->i',y[:,0],Ig)
    H=(np.einsum('nij,n->ij',y[:,HESSIAN_INDICES],I)
       +np.einsum('ni,nj->ij',y[:,1:4],Ig)+np.einsum('ni,nj->ij',Ig,y[:,1:4])
       +np.einsum('n,nij->ij',y[:,0],Ih))
    return np.r_[y[:,0]@I,gradient,H[0],H[1,1:],H[2,2]]


def nonlinear_sector_matrix(coefficients):
    c=np.asarray(coefficients,float)
    if c.shape!=(11,) or np.any(~np.isfinite(c)):raise ValueError('eleven finite coefficients required')
    return np.array([[c[4],c[10]/2],[c[10]/2,c[9]]])


class CrossDensityAngularInterface(QuarticSymmetryInterface):
    def __init__(self,decays,coefficients,*,tolerance=2e-11):
        c=np.asarray(coefficients,float);block=nonlinear_sector_matrix(c)
        scaling=1/np.sqrt(np.maximum(abs(np.diag(block)),1.))
        if np.linalg.eigvalsh(scaling[:,None]*block*scaling[None,:])[0]<-1e-9:
            raise ValueError('density/angular nonlinear sector is not positive semidefinite')
        super().__init__(decays,c[:10],tolerance=tolerance)
        self.coefficients=c.copy()

    def cross_jet(self,state):
        density,*_=self.base._density(np.asarray(state,float))
        tensors,diag=self.quartic.site_jets(tuple(state))
        return 2*cross_site_jet(density/self.face.rho_bulk,tensors),diag

    def evaluate(self,state):
        result=super().evaluate(state)
        if self.coefficients[10]==0:return result
        jet,diag=self.cross_jet(state)
        components=dict(result.components);components['density_angular_cross']=self.coefficients[10]*jet
        diagnostics=dict(result.diagnostics);diagnostics['density_angular_cross']=diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()),components,diagnostics)


class CrossDensityAngularCache(QuarticSymmetryObservationCache):
    @lru_cache(maxsize=128)
    def cross(self,scalar,odd):
        probe=build_range_surface(scalar,odd,5.,np.ones(8));moment=QuarticMoment(probe)
        jets={}
        for observation in self.observations:
            if observation.state is not None and observation.state not in jets:
                state=observation.state
                density=probe._density(np.asarray(state))[0]/probe.face.rho_bulk
                tensors,_=moment.site_jets(state)
                jets[state]=2*cross_site_jet(density,tensors)
        return np.array([0. if o.bulk_index is not None else np.asarray(o.jet_weights)@jets[o.state]
                         for o in self.observations])

    def matrix(self,decays):
        if len(decays)!=3:raise ValueError('unchanged three radial ranges required')
        return np.column_stack([super().matrix(decays),self.cross(float(decays[0]),float(decays[1]))])


class CrossDensityAngularBulk(QuarticSymmetryTailBulkBasis):
    def evaluate(self,q):
        columns,tails=super().evaluate(q)
        return np.concatenate([columns,np.zeros((1,3,3))]),np.r_[tails,0.]


def cross_sector_profile(matrix,observations,columns,tails,*,nonnegative,max_sector_cuts=24):
    """Spectral constraints AND an independent two-invariant convex PSD cone."""
    cuts=[];history=[]
    for iteration in range(max_sector_cuts+1):
        fit=spectral_profile(matrix,observations,columns,tails,nonnegative=nonnegative,inequalities=cuts)
        block=nonlinear_sector_matrix(fit['coefficients'])
        scale=1/np.sqrt(np.maximum(abs(np.diag(block)),1.))
        normalized=scale[:,None]*block*scale[None,:]
        eigen,vectors=np.linalg.eigh(normalized)
        history.append(dict(loss=fit['squared_loss'],normalized_sector_eigenvalues=eigen))
        if eigen[0]>=-2e-10:
            fit.update(environment_sector_verified=True,environment_sector_cut_iterations=iteration,
                environment_sector_history=history,environment_sector_normalized_margin=float(eigen[0]))
            return fit
        v=scale*vectors[:,0]
        row=np.zeros(11);row[4]=v[0]**2;row[9]=v[1]**2;row[10]=v[0]*v[1]
        cuts.append(row)
    raise ArithmeticError('mixed nonlinear PSD sector did not converge; candidate rejected')
