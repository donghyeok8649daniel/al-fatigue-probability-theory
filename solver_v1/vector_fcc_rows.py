"""Three-component infinite-row FCC STATIC research, not production fatigue.

Actual row radii and per-site environments are recomputed at every state.
The infinite x sum is still Poisson/Bessel, including its transverse zero mode.
No continuum/core splice, empirical yield law, kinetic time or parameter fit.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
import math

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, cg, eigsh
from scipy.special import k0, k1, zeta

from .lattice_bessel import _zero_mode
from .nonlinear_fcc_screw import NonlinearScrewRows, stf_basis


def power_mode_radial(radius, p, mode, b):
    """Exact half-integer Bessel form; returns C,C_A,C_AA for p=3,6."""
    if mode==0:
        value,first=_zero_mode(radius,p,b)
        return value,first,-2*p*first/radius
    g=2*np.pi*mode/b;z=g*radius
    def polynomial(n):
        coefficients=[math.factorial(n+k)/(math.factorial(k)*math.factorial(n-k)*2**k)
                      for k in range(n+1)]
        return np.polynomial.polynomial.polyval(1/z,coefficients)
    common=2**(2-p)*np.pi/(b*math.gamma(p))*g**(p-1)*radius**(-p)*np.exp(-z)
    value=common*polynomial(p-1)
    first=-g*common*polynomial(p)
    return value,first,g*g*value-2*p*first/radius


def exponential_radial(radius,g,kappa):
    """H,H_g,H_gg,H_ggg and their first TWO radial derivatives.

    Integer K recurrence uses only K0,K1; no numerical differentiation.
    """
    q=np.sqrt(kappa*kappa+g*g);z=radius*q
    kval=[k0(z),k1(z)]
    for n in (1,2,3):
        kval.append(kval[n-1]+2*n/z*kval[n])
    values=[];first=[];second=[]
    for j in range(4):
        common=2*radius*kappa*(-radius/2)**j*q**(-j-1)
        values.append(common*kval[j+1])
        first.append(-common*q*kval[j])
        second.append(first[-1]/radius+common*q*q*kval[abs(j-1)])
    def ordinary(array):
        return np.stack([array[0],2*g*array[1],2*array[1]+4*g*g*array[2],
                         12*g*array[2]+8*g**3*array[3]],axis=-1)
    return ordinary(values),ordinary(first),ordinary(second)


@lru_cache(maxsize=3)
def symmetric_monomials(rank):
    """Combine equal raw tensor monomials before the orthonormal STF map."""
    exponents=[];vectors=[];basis=stf_basis(rank)
    for index,indices in enumerate(product(range(3),repeat=rank)):
        exponent=tuple(indices.count(axis) for axis in range(3))
        if exponent not in exponents:
            exponents.append(exponent);vectors.append(np.zeros(2*rank+1))
        vectors[exponents.index(exponent)]+=basis[index]
    return tuple(exponents),np.array(vectors)


def _radial_polynomial(h,ha,haa,y,z,ny,nz,order):
    """h(A)*y^ny*z^nz and its transverse derivatives (all arrays flat)."""
    radius=np.hypot(y,z);uy=y/radius;uz=z/radius
    mon=y**ny*z**nz
    my=ny*y**(ny-1)*z**nz if ny else np.zeros_like(y)
    mz=nz*y**ny*z**(nz-1) if nz else np.zeros_like(y)
    value=h*mon
    gradient=np.column_stack([ha*uy*mon+h*my,ha*uz*mon+h*mz])
    if order<2:
        return value,gradient,None
    myy=ny*(ny-1)*y**(ny-2)*z**nz if ny>1 else np.zeros_like(y)
    mzz=nz*(nz-1)*y**ny*z**(nz-2) if nz>1 else np.zeros_like(y)
    myz=ny*nz*y**(ny-1)*z**(nz-1) if ny and nz else np.zeros_like(y)
    radial=haa-ha/radius
    hessian=np.stack([radial*uy*uy*mon+ha/radius*mon+2*ha*uy*my+h*myy,
                     radial*uy*uz*mon+ha*(uy*mz+uz*my)+h*myz,
                     radial*uz*uz*mon+ha/radius*mon+2*ha*uz*mz+h*mzz],axis=-1)
    return value,gradient,hessian


class VectorRowKernel:
    """Infinite-row channels [pair,rho,STF1,STF2,STF3] at actual (x,y,z).

    Analytic reciprocal series, adaptive per-row positive modes. Tolerance is
    an absolute per-channel derivative envelope, not a material error bar.
    Transverse row/domain omission is a separate convergence study.
    """
    def __init__(self,rows,*,tolerance=2e-12,max_modes=40):
        if not np.isfinite(tolerance) or tolerance<=0 or max_modes<3:
            raise ValueError('positive finite tolerance and at least 3 modes required')
        self.rows=rows;self.tolerance=tolerance;self.max_modes=max_modes
        self.parity=np.r_[1.,1.,-np.ones(3),np.ones(5),-np.ones(7)]

    def _mode(self,vectors,mode,order):
        x,y,z=vectors.T;radius=np.hypot(y,z);g=2*np.pi*mode/self.rows.b
        c6,d6,h6=power_mode_radial(radius,6,mode,self.rows.b)
        c3,d3,h3=power_mode_radial(radius,3,mode,self.rows.b)
        p=self.rows.bulk.p
        pair=[4*p.epsilon*(p.sigma_lj**12*a-p.sigma_lj**6*b) for a,b in zip((c6,d6,h6),(c3,d3,h3))]
        factor=(1 if mode==0 else 2)/self.rows.b
        den=self.rows.bulk.density_params
        density=[factor*den.C_rho*v[:,0] for v in exponential_radial(radius,g,den.kappa)]
        values=[];gradients=[];hessians=[]
        def pack(h,ha,haa,ny=0,nz=0,coefficient=1.):
            val,grad,hess=_radial_polynomial(h,ha,haa,y,z,ny,nz,order)
            val=val*coefficient;grad=grad*coefficient
            gradient=np.column_stack([1j*g*val,grad])
            tensor=None
            if order==2:
                hess=hess*coefficient
                tensor=np.zeros((len(x),3,3),complex)
                tensor[:,0,0]=-g*g*val
                tensor[:,0,1:]=1j*g*grad;tensor[:,1:,0]=1j*g*grad
                tensor[:,1,1]=hess[:,0];tensor[:,1,2]=hess[:,1]
                tensor[:,2,1]=hess[:,1];tensor[:,2,2]=hess[:,2]
            return val,gradient,tensor
        for scalar in (pair,density):
            val,grad,hess=pack(*scalar)
            values.append(val[:,None]);gradients.append(grad[:,None,:])
            if order==2:
                hessians.append(hess[:,None,:,:])
        angular=self.rows.surface.angular
        H,Ha,Haa=exponential_radial(radius,g,angular.kappa)
        for rank in (1,2,3):
            exponents,basis=symmetric_monomials(rank)
            raw=[];rawgrad=[];rawhess=[]
            for nx,ny,nz in exponents:
                val,grad,hess=pack(H[:,nx],Ha[:,nx],Haa[:,nx],ny,nz,
                                   factor*angular.amplitude*1j**nx)
                raw.append(val);rawgrad.append(grad)
                if order==2:
                    rawhess.append(hess)
            values.append(np.stack(raw,axis=1)@basis)
            gradients.append(np.einsum('nti,tc->nci',np.stack(rawgrad,axis=1),basis,optimize=True))
            if order==2:
                hessians.append(np.einsum('ntij,tc->ncij',np.stack(rawhess,axis=1),basis,optimize=True))
        values=np.concatenate(values,axis=1)
        gradients=np.concatenate(gradients,axis=1)
        phase=np.exp(1j*g*x)
        envelope=np.maximum(np.max(abs(values),axis=1),np.max(abs(gradients),axis=(1,2)))
        tensor=None
        if order==2:
            tensor=np.concatenate(hessians,axis=1)
            envelope=np.maximum(envelope,np.max(abs(tensor),axis=(1,2,3)))
            tensor=np.real(tensor*phase[:,None,None,None])
        return np.real(values*phase[:,None]),np.real(gradients*phase[:,None,None]),tensor,envelope

    def evaluate(self,vectors,*,order=1):
        vectors=np.asarray(vectors,float)
        if vectors.ndim!=2 or vectors.shape[1]!=3 or not np.all(np.isfinite(vectors)) or order not in (1,2):
            raise ValueError('finite (n,3) vectors and derivative order 1 or 2 required')
        if np.any(np.hypot(vectors[:,1],vectors[:,2])<=0):
            raise ValueError('coincident transverse rows need a separate limiting representation')
        vectors=vectors.copy()
        vectors[:,0]-=self.rows.b*np.floor(vectors[:,0]/self.rows.b+.5)
        value,gradient,hessian,_=self._mode(vectors,0,order)
        active=np.arange(len(vectors));small=np.zeros(len(vectors),int)
        used=np.zeros(len(vectors),int);last=np.zeros(len(vectors))
        for mode in range(1,self.max_modes+1):
            v,g,h,envelope=self._mode(vectors[active],mode,order)
            value[active]+=v;gradient[active]+=g
            if order==2:
                hessian[active]+=h
            used[active]=mode;last[active]=envelope
            small[active]=np.where(envelope<self.tolerance,small[active]+1,0)
            active=active[small[active]<2]
            if len(active)==0:
                break
        if len(active):
            raise ArithmeticError('vector row reciprocal series exhausted; no repaired force returned')
        return dict(value=value,gradient=gradient,hessian=hessian,
                    modes_used=int(used.max(initial=0)),maximum_last_mode_envelope=float(last.max(initial=0.)))


class VectorPeriodicCell:
    """Three local row displacements, fixed transverse periodic cell vectors.

    A homogeneous xz shear can be strain- or stress-controlled. The other cell
    strains remain fixed: local vector force balance is NOT zero normal-cell
    stress. Row-radius and actual per-site environments are never frozen.
    ``ring`` controls transverse difference neighborhoods, not the infinite
    atomic-row sum. Its algebraic G=0 tail requires independent refinement.
    """
    def __init__(self,rows,shape,*,ring=6,tolerance=2e-12,chunk_size=4096):
        self.rows=rows;self.scalar=rows.periodic_cell(shape)
        self.shape=self.scalar.shape;self.size=self.scalar.size
        if int(ring)!=ring or ring<1 or int(chunk_size)!=chunk_size or chunk_size<1:
            raise ValueError('positive integer ring/chunk controls required')
        self.ring=int(ring);self.chunk_size=int(chunk_size)
        self.kernel=VectorRowKernel(rows,tolerance=tolerance)
        j,l=np.meshgrid(np.arange(-ring,ring+1),np.arange(-ring,ring+1),indexing='ij')
        half=(l>0)|((l==0)&(j>0));self.j=j[half];self.l=l[half]
        self.reference=np.column_stack([rows.b*(self.j+self.l)/2,
                                       rows.d*(self.j+self.l/3),rows.h*self.l])
        self.neighbors=len(self.j)
        sites=np.arange(self.size).reshape(self.shape)
        self.destination=np.stack([np.roll(sites,(-int(j),-int(l)),axis=(0,1)).ravel()
                                   for j,l in zip(self.j,self.l)],axis=1)
        self._flat_destination=self.destination.ravel()
        self.reference_channels=self.kernel.evaluate(self.reference)['value']
        self.volume=self.size*rows.b*rows.d*rows.h

    def _scatter(self,values):
        """Sum oriented bond contributions into their destination sites."""
        values=np.asarray(values);tail=values.shape[2:]
        flat=values.reshape((self.size*self.neighbors,-1))
        out=np.stack([np.bincount(self._flat_destination,weights=flat[:,c],minlength=self.size)
                      for c in range(flat.shape[1])],axis=-1)
        return out.reshape((self.size,)+tail)

    def _assemble(self,field,gamma,*,second=False):
        field=np.asarray(field,float)
        if field.shape!=self.shape+(3,) or not np.all(np.isfinite(field)) or not np.isfinite(gamma):
            raise ValueError('finite three-component field with declared shape required')
        u=field.reshape((self.size,3))
        vectors=self.reference[None,:,:]+u[self.destination]-u[:,None,:]
        vectors[...,0]+=gamma*self.reference[:,2]
        actual_radius=np.hypot(vectors[...,1],vectors[...,2])
        flat=vectors.reshape((-1,3));count=len(flat)
        values=np.empty((count,17));gradients=np.empty((count,17,3))
        tensors=np.empty((count,17,3,3)) if second else None
        modes=0;envelope=0.
        for start in range(0,count,self.chunk_size):
            stop=min(count,start+self.chunk_size)
            result=self.kernel.evaluate(flat[start:stop],order=2 if second else 1)
            values[start:stop]=result['value'];gradients[start:stop]=result['gradient']
            if second:
                tensors[start:stop]=result['hessian']
            modes=max(modes,result['modes_used'])
            envelope=max(envelope,result['maximum_last_mode_envelope'])
        diff=values.reshape((self.size,self.neighbors,17))-self.reference_channels
        channels=diff.sum(axis=1)+self._scatter(diff*self.kernel.parity)
        rho=self.rows.rho_bulk+channels[:,1]
        if np.any(rho<=0) or not np.all(np.isfinite(rho)):
            raise ValueError('nonpositive/nonfinite actual site density; no repair applied')
        F=self.rows.embedding
        weights=np.empty_like(channels);weights[:,0]=.5
        weights[:,1]=F.first_derivative(rho)
        weights[:,2:]=2*self.rows.angular_weights*channels[:,2:]
        energy=.5*channels[:,0].sum()+np.sum(F.value(rho)-F.value(self.rows.rho_bulk))
        energy+=np.sum(self.rows.angular_weights*channels[:,2:]**2)
        bond_weights=weights[:,None,:]+weights[self.destination]*self.kernel.parity
        gradients=gradients.reshape((self.size,self.neighbors,17,3))
        bond_force=np.einsum('nrc,nrci->nri',bond_weights,gradients)
        force=-bond_force.sum(axis=1)+self._scatter(bond_force)
        affine=float(np.sum(bond_force[...,0]*self.reference[:,2]))
        out=dict(energy=float(energy),gradient=force.reshape(self.shape+(3,)),
                 affine_derivative=affine,minimum_site_density=float(rho.min()),
                 minimum_transverse_radius=float(actual_radius.min()),
                 channels=channels.reshape(self.shape+(17,)),density=rho.reshape(self.shape),
                 reciprocal_modes=modes,last_mode_envelope=envelope,transverse_ring=self.ring)
        if not second:
            return out,None
        # Contract the row Hessian after the actual site weights are known.
        # Keep only 3x3 per bond, not the full channel Hessian in repeated CG.
        local_hessian=np.einsum('nrc,nrcij->nrij',bond_weights,
                               tensors.reshape((self.size,self.neighbors,17,3,3)))
        Fsecond=F.second_derivative(rho)
        def apply(direction,gamma_direction=0.):
            v=np.asarray(direction,float)
            if v.shape!=self.shape+(3,) or not np.all(np.isfinite(v)) or not np.isfinite(gamma_direction):
                raise ValueError('finite vector Hessian direction required')
            v=v.reshape((self.size,3))
            dv=v[self.destination]-v[:,None,:]
            dv[...,0]+=gamma_direction*self.reference[:,2]
            tangent=np.einsum('nrci,nri->nrc',gradients,dv)
            dc=tangent.sum(axis=1)+self._scatter(tangent*self.kernel.parity)
            dw=np.zeros_like(dc);dw[:,1]=Fsecond*dc[:,1]
            dw[:,2:]=2*self.rows.angular_weights*dc[:,2:]
            bw=dw[:,None,:]+dw[self.destination]*self.kernel.parity
            dbond=np.einsum('nrc,nrci->nri',bw,gradients)+np.einsum('nrij,nrj->nri',local_hessian,dv)
            hv=-dbond.sum(axis=1)+self._scatter(dbond)
            hg=float(np.sum(dbond[...,0]*self.reference[:,2]))
            return hv.reshape(self.shape+(3,)),hg
        return out,apply

    def evaluate(self,field,gamma=0.,*,direction=None,gamma_direction=0.,fields=False):
        out,apply=self._assemble(field,gamma,second=direction is not None)
        if direction is not None:
            out['hessian_vector'],out['affine_hessian_vector']=apply(direction,gamma_direction)
        if not fields:
            out.pop('channels');out.pop('density')
        return out

    def linearize(self,field,gamma=0.):
        """Cache the actual analytic Hessian for repeated directional products."""
        return self._assemble(field,gamma,second=True)

    def pair_zero_mode_tail_bounds(self,field):
        """Conservative omitted LJ G=0 bounds; not a bound on all channels.

        Each square offset shell k contains 8k rows. Its radius >= alpha*k.
        Periodicity cancels the summed linear energy term. Taylor's remainder
        and the radial Hessian give an O(ring^-5) attractive energy/force bound.
        Exponential density/moments and nonzero modes still need refinement.
        """
        field=np.asarray(field,float)
        if field.shape!=self.shape+(3,) or not np.all(np.isfinite(field)):
            raise ValueError('finite vector field required')
        transverse=field[...,1:]-field[...,1:].mean(axis=(0,1))
        D=2*float(np.max(np.linalg.norm(transverse,axis=-1)))
        matrix=np.array([[self.rows.d,self.rows.d/3],[0.,self.rows.h]])
        alpha=float(np.linalg.svd(matrix,compute_uv=False).min())
        reduced=alpha-D/(self.ring+1)
        if reduced<=0:
            raise ValueError('tail bound cannot separate deformed rows from omitted neighborhood')
        p=self.rows.bulk.p;bound=0.
        for power,scale in ((6,p.sigma_lj**12),(3,p.sigma_lj**6)):
            c=math.sqrt(np.pi)*math.gamma(power-.5)/(self.rows.b*math.gamma(power))
            exponent=2*power+1
            bound+=4*p.epsilon*scale*(2*power)*(2*power-1)*c*8*reduced**(-exponent)*zeta(exponent-1,self.ring+1)
        return dict(pair_zero_mode_energy_bound_ev=self.size*D*D*bound/4,
                    pair_zero_mode_force_norm_bound=D*bound,
                    maximum_transverse_difference_bound=D,
                    all_channel_infinite_tail_certified=False)

    def layer_registry(self,field,gamma=0.):
        """Vector adjacent-layer shifts, not scalar n-defined plastic strain.

        Compare additional in-plane shifts with 0,tau,2tau modulo the triangular
        lattice. Small row-to-row dispersion is required before interpreting a
        layer as a rigid stacking fault. Scalar x well partitions can put an
        actual vector minimum exactly on x=b/2 and are not a vector mechanism.
        """
        field=np.asarray(field,float)
        if field.shape!=self.shape+(3,) or not np.all(np.isfinite(field)) or not np.isfinite(gamma):
            raise ValueError('finite vector field and affine shear required')
        bond=np.roll(field,-1,axis=1)-field
        bond[...,0]+=gamma*self.rows.h
        b,d=self.rows.b,self.rows.d;tau=np.array([b/2,d/3])
        def distance(shift,target):
            q=shift-target;n0=int(np.floor(q[1]/d+.5));values=[]
            for n in (n0-1,n0,n0+1):
                x=q[0]-n*b/2;m=int(np.floor(x/b+.5))
                values.append(np.hypot(x-m*b,q[1]-n*d))
            return float(min(values))
        output=[]
        for layer in range(self.shape[1]):
            phase=np.exp(2j*np.pi*bond[:,layer,0]/b)
            center=np.exp(1j*np.angle(phase.mean()))
            x=float(b*np.angle(center)/(2*np.pi));y=float(bond[:,layer,1].mean())
            distances=[distance(np.array([x,y]),k*tau) for k in (0,1,2)]
            margin=b*abs(bond[:,layer,0]/b-np.floor(bond[:,layer,0]/b)-.5)
            output.append(dict(layer=layer,additional_slip_x_over_L0=x,
                additional_slip_y_over_L0=y,opening_change_over_L0=float(bond[:,layer,2].mean()),
                x_phase_dispersion=float(np.max(abs(phase-center))),
                transverse_row_dispersion=float(np.max(np.ptp(bond[:,layer,1:],axis=0))),
                minimum_scalar_x_partition_margin=float(margin.min()),
                distance_to_perfect_extra_registry=distances[0],
                distance_to_extra_tau=distances[1],distance_to_extra_2tau=distances[2],
                closest_extra_registry_class=int(np.argmin(distances)),
                scalar_x_partition_is_not_vector_plastic_strain=True))
        return output

    def relax(self,initial,*,gamma=0.,shear_stress=None,force_tolerance=2e-7,
              max_iterations=500,polish_steps=4,callback=None):
        """Deterministic static relaxation; no topology/holding-force constraints.

        Optimizer stopping and actual force convergence are reported separately.
        Only three uniform translations are gauge fixed. No force clipping.
        """
        initial=np.asarray(initial,float).copy()
        if initial.shape!=self.shape+(3,) or not np.all(np.isfinite(initial)):
            raise ValueError('finite initial three-component field required')
        if (not np.isfinite(gamma) or not np.isfinite(force_tolerance) or force_tolerance<=0
                or max_iterations<1 or (shear_stress is not None and not np.isfinite(shear_stress))):
            raise ValueError('invalid static relaxation controls')
        initial-=initial.mean(axis=(0,1));stress_control=shear_stress is not None
        scale=self.rows.h*np.sqrt(self.size);n=3*self.size
        def unpack(x):
            return x[:n].reshape(self.shape+(3,)),x[-1]/scale if stress_control else gamma
        def pack_gradient(out):
            grad=out['gradient']-out['gradient'].mean(axis=(0,1))
            if not stress_control:
                return grad.ravel()
            return np.r_[grad.ravel(),(out['affine_derivative']-shear_stress*self.volume)/scale]
        def fun(x):
            field,strain=unpack(x);out=self.evaluate(field,strain)
            return out['energy']-(shear_stress*self.volume*strain if stress_control else 0.),pack_gradient(out)
        x0=np.r_[initial.ravel(),gamma*scale] if stress_control else initial.ravel()
        iterations=[0]
        def progress(x):
            iterations[0]+=1
            if callback is not None:
                field,strain=unpack(x);callback(iterations[0],field.copy(),float(strain))
        sol=minimize(fun,x0,jac=True,method='L-BFGS-B',callback=progress,
                     options=dict(gtol=force_tolerance/10,ftol=2e-15,maxiter=max_iterations,maxls=35,maxcor=24))
        x=sol.x.copy();polished=0
        for _ in range(polish_steps):
            field,strain=unpack(x);out,apply=self.linearize(field,strain)
            g=pack_gradient(out)
            if np.max(abs(g))<=force_tolerance:
                break
            def matvec(v):
                dv,dg=v[:n].reshape(self.shape+(3,)),v[-1]/scale if stress_control else 0.
                hv,hg=apply(dv,dg)
                hv=hv-hv.mean(axis=(0,1))+dv.mean(axis=(0,1))
                return np.r_[hv.ravel(),hg/scale] if stress_control else hv.ravel()
            operator=LinearOperator((len(x),len(x)),matvec=matvec,dtype=float)
            step,info=cg(operator,-g,rtol=3e-5,atol=force_tolerance/100,maxiter=350)
            if info!=0 or g@step>=0:
                break
            energy=out['energy']-(shear_stress*self.volume*strain if stress_control else 0.)
            accepted=False
            for damping in (1.,.5,.25,.125):
                trial=x+damping*step;e,trial_g=fun(trial)
                if e<=energy+3e-13*max(1.,abs(energy)) and np.linalg.norm(trial_g)<np.linalg.norm(g):
                    x=trial;polished+=1;accepted=True;break
            if not accepted:
                break
        field,strain=unpack(x);field=field.copy();field-=field.mean(axis=(0,1))
        out=self.evaluate(field,strain)
        residual=float(np.max(abs(out['gradient'])))
        affine_residual=abs(out['affine_derivative']-shear_stress*self.volume)/scale if stress_control else 0.
        out.update(displacement=field,affine_shear=float(strain),
                   internal_shear_stress=out['affine_derivative']/self.volume,
                   maximum_force_residual=residual,scaled_stress_residual=float(affine_residual),
                   force_converged=max(residual,affine_residual)<=force_tolerance,
                   optimizer_success=bool(sol.success),optimizer_message=str(sol.message),
                   iterations=sol.nit,evaluations=sol.nfev,newton_polish_steps=polished,
                   x_winding=self.scalar.winding(field[...,0]),
                   physical_time_available=False,transverse_cell_strains_fixed=True)
        return out

    def minimum_curvature(self,field,gamma=0.,*,stress_control=False,tolerance=1e-7):
        """Smallest nongauge Hessian mode; 3 translations excluded explicitly."""
        _,apply=self.linearize(field,gamma)
        scale=self.rows.h*np.sqrt(self.size);n=3*self.size+int(stress_control)
        def matvec(v):
            dv=v[:3*self.size].reshape(self.shape+(3,))
            hv,hg=apply(dv,v[-1]/scale if stress_control else 0.)
            hv=hv-hv.mean(axis=(0,1))+dv.mean(axis=(0,1))
            return np.r_[hv.ravel(),hg/scale] if stress_control else hv.ravel()
        operator=LinearOperator((n,n),matvec=matvec,dtype=float)
        values,vectors=eigsh(operator,k=4,which='SA',tol=tolerance,
                             v0=np.cos(np.arange(n)*.371),maxiter=1000)
        overlap=np.linalg.norm(vectors[:3*self.size].reshape((self.size,3,4)).sum(axis=0),axis=0)/np.sqrt(self.size)
        candidates=np.flatnonzero(overlap<.5)
        if not len(candidates):
            raise ArithmeticError('no material curvature mode isolated from translations')
        index=candidates[np.argmin(values[candidates])];vector=vectors[:,index]
        return dict(minimum_eigenvalue=float(values[index]),
                    eigen_residual=float(np.linalg.norm(matvec(vector)-values[index]*vector)),
                    translation_overlap=float(overlap[index]),translation_gauge_eigenvalue=1.,
                    affine_strain_included=stress_control)
