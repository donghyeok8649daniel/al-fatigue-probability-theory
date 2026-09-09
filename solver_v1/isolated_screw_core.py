"""Isolated straight screw, explicit Dirichlet far field, STATIC research.

The x atomic rows are infinite Poisson/Bessel sums. Only transverse *changes*
of site environments have a refined neighborhood. All sites whose energies
depend on a free row are included, including fixed boundary sites. This is
not the periodic opposite-pair problem, a finite loop, kinetics or a yield law.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, cg
from scipy.spatial import ConvexHull

from .nonlocal_interface_elasticity import screw_elastic_factor
from .vector_fcc_rows import VectorRowKernel


def row_environment(surface):
    """No scalar FFT approximation; supports independently ranged STF ranks."""
    bulk=surface.interface.bulk; b=float(bulk.geometry.b); h=float(bulk.a0)
    if not np.isclose(h/b,np.sqrt(2/3),rtol=2e-8,atol=0):
        raise ValueError('cubic background required; noncubic moments are not zero')
    return SimpleNamespace(surface=surface,bulk=bulk,b=b,h=h,d=np.sqrt(3)*b/2,
        rho_bulk=float(bulk.full_environment_density(h,0).value),embedding=bulk.embedding,
        angular_weights=np.repeat([surface.vector_amplitude_ev,
            surface.quadrupole_amplitude_ev,surface.amplitude_ev],[3,5,7]))


@dataclass(frozen=True)
class ScrewFarField:
    """Anti-plane anisotropic elasticity from the SAME energy's cubic tensor.

    Coordinates/translation b are in L0; tensor and applied traction must use
    the same units. For energy matching use tensor in eV/L0^3. Then mu*b^3/4pi
    is eV per repeat b*L0 of a straight line, NOT a finite activation energy.
    """
    tensor: np.ndarray
    b: float
    center: tuple[float,float]

    def __post_init__(self):
        screw_elastic_factor(self.tensor)
        if not np.isfinite(self.b) or self.b<=0 or not np.all(np.isfinite(self.center)):
            raise ValueError('positive Burgers translation and finite core center required')

    @property
    def matrix(self):
        return np.asarray(self.tensor)[0,1:,0,1:]

    @property
    def log_energy_coefficient(self):
        return screw_elastic_factor(self.tensor)*self.b**3/(4*np.pi)

    def displacement(self, positions, *, shear_traction=0., burgers_sign=1):
        q=np.asarray(positions,float)
        if q.ndim!=2 or q.shape[1]!=3 or not np.all(np.isfinite(q)):
            raise ValueError('finite row positions (n,3) required')
        if burgers_sign not in (-1,0,1) or not np.isfinite(shear_traction):
            raise ValueError('declared Burgers sign and finite traction required')
        A,B,D=self.matrix[0,0],self.matrix[0,1],self.matrix[1,1]
        p=(-B+1j*np.sqrt(A*D-B*B))/D
        yz=q[:,1:]-np.asarray(self.center)
        argument=yz[:,0]+p*yz[:,1]
        if burgers_sign and np.any(abs(argument)==0):
            raise ValueError('far-field branch center cannot coincide with an atomic row')
        # sigma_xy=0, sigma_xz=declared traction; not an independently fitted E.
        strain=np.linalg.solve(self.matrix,[0.,shear_traction])
        out=np.zeros_like(q)
        out[:,0]=burgers_sign*self.b/(2*np.pi)*np.angle(argument)+yz@strain
        return out


class IsolatedScrewCore:
    """Finite free disk, infinite-row environments and fixed elastic outside.

    Directed neighbor environments implement one-half pair counting and
    per-atom nonlinear embedding. The outer sites affected by free rows are
    energy sites too; omitting their F would give an incorrect core force.
    Frozen environment changes outside `ring` are NOT silently certified.
    """
    def __init__(self,surface,far_field,*,free_radius=3.,ring=3,tolerance=2e-12,
                 shear_traction=0.,burgers_sign=1,chunk_size=4096):
        self.rows=row_environment(surface); self.far_field=far_field
        if (not np.isfinite(free_radius) or free_radius<=0 or int(ring)!=ring or ring<1
                or int(chunk_size)!=chunk_size or chunk_size<1):
            raise ValueError('positive radius, integer neighborhood and chunk size required')
        if not np.isclose(far_field.b,self.rows.b,rtol=1e-12,atol=0):
            raise ValueError('far-field Burgers translation differs from the atomic row repeat')
        self.free_radius=float(free_radius); self.ring=int(ring); self.chunk_size=int(chunk_size)
        self.shear_traction=float(shear_traction); self.burgers_sign=burgers_sign
        self.kernel=VectorRowKernel(self.rows,tolerance=tolerance)
        offsets=[(j,l) for j,l in product(range(-ring,ring+1),repeat=2) if j or l]
        self.offsets=np.array(offsets,int); self.reference=self.positions(self.offsets)
        reference=self.kernel.evaluate(self.reference)
        self.reference_channels=reference['value']
        self.reference_bond_gradient=(.5*reference['gradient'][:,0,:]
            +self.rows.embedding.first_derivative(self.rows.rho_bulk)*reference['gradient'][:,1,:])
        transform=np.array([[self.rows.d,self.rows.d/3],[0.,self.rows.h]])
        extent=int(np.ceil((free_radius+np.linalg.norm(far_field.center))/
                           np.linalg.svd(transform,compute_uv=False).min()))+1
        free=set()
        for index in product(range(-extent,extent+1),repeat=2):
            if np.linalg.norm(self.positions(np.array([index]))[0,1:]-far_field.center)<free_radius:
                free.add(index)
        if not free:
            raise ValueError('free disk contains no rows')
        def expand(sites):
            return sites|{(j+dj,l+dl) for j,l in sites for dj,dl in offsets}
        energy=expand(free); all_sites=expand(energy)
        self.indices=np.array(sorted(all_sites),int); lookup={tuple(x):i for i,x in enumerate(self.indices)}
        self.free_ids=np.array([lookup[x] for x in sorted(free)],int)
        self.energy_ids=np.array([lookup[x] for x in sorted(energy)],int)
        self.destination=np.array([[lookup[(j+dj,l+dl)] for dj,dl in offsets]
                                   for j,l in sorted(energy)],int)
        self._flat_destination=self.destination.ravel()
        self.xyz=self.positions(self.indices)
        self.boundary=far_field.displacement(self.xyz,shear_traction=shear_traction,
                                             burgers_sign=burgers_sign)
        self.loop_ids=ConvexHull(self.xyz[:,1:]).vertices

    def positions(self,indices):
        j,l=np.asarray(indices).T; rows=self.rows
        return np.column_stack([rows.b*(j+l)/2,rows.d*(j+l/3),rows.h*l])

    @property
    def initial(self):
        return self.boundary[self.free_ids].copy()

    def full_field(self,free):
        free=np.asarray(free,float)
        if free.shape!=(len(self.free_ids),3) or not np.all(np.isfinite(free)):
            raise ValueError('finite vector displacement for every free row required')
        field=self.boundary.copy(); field[self.free_ids]=free
        return field

    def _scatter(self,values):
        flat=values.reshape((self.destination.size,-1))
        return np.stack([np.bincount(self._flat_destination,weights=flat[:,i],
                                    minlength=len(self.indices)) for i in range(flat.shape[1])],axis=1)

    def _assemble(self,free,*,second=False):
        field=self.full_field(free); ns,nr=self.destination.shape
        vectors=self.reference[None,:,:]+field[self.destination]-field[self.energy_ids,None,:]
        flat=vectors.reshape((-1,3)); values=np.empty((len(flat),17)); grads=np.empty((len(flat),17,3))
        tensors=np.empty((len(flat),17,3,3)) if second else None
        modes=0; envelope=0.
        for start in range(0,len(flat),self.chunk_size):
            stop=min(start+self.chunk_size,len(flat))
            r=self.kernel.evaluate(flat[start:stop],order=2 if second else 1)
            values[start:stop]=r['value']; grads[start:stop]=r['gradient']
            if second:
                tensors[start:stop]=r['hessian']
            modes=max(modes,r['modes_used']); envelope=max(envelope,r['maximum_last_mode_envelope'])
        channels=(values.reshape((ns,nr,17))-self.reference_channels).sum(axis=1)
        rho=self.rows.rho_bulk+channels[:,1]; F=self.rows.embedding
        if np.any(rho<=0) or not np.all(np.isfinite(rho)):
            raise ValueError('invalid actual site density; no clipping/repair')
        weights=np.empty_like(channels); weights[:,0]=.5; weights[:,1]=F.first_derivative(rho)
        weights[:,2:]=2*self.rows.angular_weights*channels[:,2:]
        site_energy=.5*channels[:,0]+F.value(rho)-F.value(self.rows.rho_bulk)
        site_energy+=np.sum(self.rows.angular_weights*channels[:,2:]**2,axis=1)
        # Exact Taylor-linear site-energy redistribution. Its summed free-row
        # derivative is zero in the perfect equilibrated lattice. Keep both
        # partitions: partial raw site sums include this boundary work and are
        # not directly the quadratic continuum strain energy in an annulus.
        linear_site_work=np.einsum('ri,nri->n',self.reference_bond_gradient,
                                  field[self.destination]-field[self.energy_ids,None,:])
        gradients=grads.reshape((ns,nr,17,3))
        bond=np.einsum('nc,nrci->nri',weights,gradients)
        gradient=self._scatter(bond); gradient[self.energy_ids]-=bond.sum(axis=1)
        out=dict(energy=float(site_energy.sum()),gradient=gradient[self.free_ids],
                 all_site_gradient=gradient,site_energy=site_energy,density=rho,
                 linear_reference_site_work=linear_site_work,
                 quadratic_remainder_site_energy=site_energy-linear_site_work,
                 minimum_density=float(rho.min()),minimum_transverse_radius=float(np.min(np.linalg.norm(vectors[...,1:],axis=-1))),
                 reciprocal_modes=modes,last_mode_envelope=envelope,
                 energy_sites=ns,free_sites=len(self.free_ids),transverse_ring=self.ring,
                 all_channel_infinite_tail_certified=False)
        if not second:
            return out,None
        local=np.einsum('nc,nrcij->nrij',weights,tensors.reshape((ns,nr,17,3,3)))
        fsecond=F.second_derivative(rho)
        def apply(direction):
            v=np.zeros_like(field); v[self.free_ids]=np.asarray(direction).reshape((-1,3))
            dv=v[self.destination]-v[self.energy_ids,None,:]
            dc=np.einsum('nrci,nri->nc',gradients,dv)
            dw=np.zeros_like(dc); dw[:,1]=fsecond*dc[:,1]
            dw[:,2:]=2*self.rows.angular_weights*dc[:,2:]
            dbond=np.einsum('nc,nrci->nri',dw,gradients)+np.einsum('nrij,nrj->nri',local,dv)
            result=self._scatter(dbond); result[self.energy_ids]-=dbond.sum(axis=1)
            return result[self.free_ids]
        return out,apply

    def evaluate(self,free):
        return self._assemble(free)[0]

    def linearize(self,free):
        return self._assemble(free,second=True)

    def boundary_winding(self,free):
        field=self.full_field(free)
        affine=self.far_field.displacement(self.xyz,shear_traction=self.shear_traction,burgers_sign=0)
        phase=2*np.pi*(field[self.loop_ids,0]-affine[self.loop_ids,0])/self.rows.b
        return float(np.sum(np.angle(np.exp(1j*(np.roll(phase,-1)-phase))))/(2*np.pi))

    def radial_energy(self,free,radii):
        out=self.evaluate(free)
        distance=np.linalg.norm(self.xyz[self.energy_ids,1:]-self.far_field.center,axis=1)
        return [dict(radius_over_L0=float(r),energy_eV_per_row_period=float(out['site_energy'][distance<r].sum()),
                     linear_reference_work_eV=float(out['linear_reference_site_work'][distance<r].sum()),
                     quadratic_remainder_energy_eV=float(out['quadratic_remainder_site_energy'][distance<r].sum()),
                     site_count=int(np.sum(distance<r)),inside_free_disk=bool(r<=self.free_radius)) for r in radii]

    def affine_antiplane_hessian(self):
        """Independent same-row-sum elastic check, energy per atomic volume.

        Two strain coordinates are u_x=gamma_y*y+gamma_z*z. This is not a
        new elastic fit and includes scalar F'' and separate-ranged STF terms.
        """
        r=self.kernel.evaluate(self.reference,order=2); R=self.reference[:,1:]
        first=np.einsum('nc,ni->ci',r['gradient'][:,:,0],R)
        second=np.einsum('nc,ni,nj->cij',r['hessian'][:,:,0,0],R,R)
        F=self.rows.embedding; rho=self.rows.rho_bulk
        result=(.5*second[0]+F.first_derivative(rho)*second[1]
                +F.second_derivative(rho)*np.outer(first[1],first[1])
                +np.einsum('c,ci,cj->ij',2*self.rows.angular_weights,first[2:],first[2:]))
        return result/(self.rows.b*self.rows.d*self.rows.h)

    def relax(self,initial=None,*,max_iterations=150,force_tolerance=2e-6,polish_steps=3,callback=None):
        if max_iterations<1 or not np.isfinite(force_tolerance) or force_tolerance<=0:
            raise ValueError('positive iterations and force tolerance required')
        shape=(len(self.free_ids),3); initial=self.initial if initial is None else initial
        self.full_field(initial); count=[0]
        def fun(x):
            out=self.evaluate(x.reshape(shape)); return out['energy'],out['gradient'].ravel()
        def progress(x):
            count[0]+=1
            if callback is not None:
                callback(count[0],x.reshape(shape).copy())
        result=minimize(fun,np.asarray(initial).ravel(),jac=True,method='L-BFGS-B',callback=progress,
                        options=dict(maxiter=max_iterations,gtol=force_tolerance/10,ftol=1e-14,maxls=35,maxcor=24))
        x=result.x.reshape(shape).copy(); polished=0
        for _ in range(polish_steps):
            out,apply=self.linearize(x); g=out['gradient'].ravel()
            if np.max(abs(g))<=force_tolerance:
                break
            operator=LinearOperator((len(g),len(g)),matvec=lambda v:apply(v.reshape(shape)).ravel(),dtype=float)
            step,info=cg(operator,-g,rtol=3e-5,atol=force_tolerance/100,maxiter=250)
            if info!=0 or g@step>=0:
                break
            accepted=False
            for factor in (1.,.5,.25,.125):
                trial=x+factor*step.reshape(shape); tested=self.evaluate(trial)
                if tested['energy']<=out['energy']+1e-13*max(1.,abs(out['energy'])) and np.linalg.norm(tested['gradient'])<np.linalg.norm(g):
                    x=trial; polished+=1; accepted=True; break
            if not accepted:
                break
        out=self.evaluate(x); residual=float(np.max(abs(out['gradient'])))
        return dict(field=x,energy=out['energy'],maximum_free_force=residual,
                    force_converged=bool(residual<=force_tolerance),optimizer_success=bool(result.success),
                    optimizer_message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),
                    polish_steps=polished,boundary_winding=self.boundary_winding(x),
                    minimum_density=out['minimum_density'],minimum_transverse_radius=out['minimum_transverse_radius'])
