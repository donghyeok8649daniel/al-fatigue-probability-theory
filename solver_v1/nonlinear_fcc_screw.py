"""Nonlinear infinite-row FCC energy, separate STATIC anti-plane research.

Each cross-section site represents an infinite atomic row parallel to e1.
All site environments are summed before nonlinear embedding/angular energies.
FFT evaluation is an exact convolution of the retained Poisson coefficients,
not a fitted PN law. Transverse and reciprocal omission require refinement.
Normal/transverse displacement and finite curved lines are NOT included here.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import permutations, product
import math

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, cg, eigsh

from .angular_environment_reference import traceless_second, traceless_third
from .discrete_fcc_screw import exponential_row_transform_derivatives
from .lattice_bessel import _mode_coefficient_and_deda


@lru_cache(maxsize=3)
def stf_basis(rank):
    """Orthonormal symmetric trace-free basis; preserves the full tensor norm."""
    size = 3**rank
    columns = []
    for vector in np.eye(size):
        raw = vector.reshape((3,)*rank)
        symmetric = sum(raw.transpose(p) for p in permutations(range(rank))) / math.factorial(rank)
        projected = (symmetric if rank == 1 else traceless_second(symmetric)
                     if rank == 2 else traceless_third(symmetric))
        columns.append(projected.ravel())
    values, vectors = np.linalg.eigh(np.array(columns).T)
    basis = vectors[:, values > .5]
    if basis.shape[1] != 2*rank+1:
        raise ArithmeticError("STF projection rank mismatch")
    return basis


def moment_mode(y, z, g, *, kappa, amplitude, b, rank):
    """Complex +m coefficient of sum exp(-k r) STF(r^rank), paired +/-m."""
    hd = exponential_row_transform_derivatives(np.hypot(y,z), g, kappa)
    raw = np.array([2*amplitude/b * 1j**indices.count(0) * hd[indices.count(0)]
                    * y**indices.count(1)*z**indices.count(2)
                    for indices in product(range(3), repeat=rank)]).T
    return raw @ stf_basis(rank)


class NonlinearScrewRows:
    """Infinite-row coefficients for arbitrary row phases at FIXED y,z.

    u_i and affine simple shear gamma specify x displacement u_i+gamma*z_i.
    Per-site channels: pair, density, STF rank1/2/3 in orthonormal coordinates.
    The whole nonlinear F(rho) is evaluated, never a frozen F'_bulk.
    G=0 row terms cancel in displacement differences in this anti-plane space.
    """
    def __init__(self, surface, *, tolerance=1e-11, max_ring=14,
                 fixed_ring=None, max_modes=64):
        self.surface = surface
        self.bulk = surface.interface.bulk
        self.b = float(self.bulk.geometry.b)
        self.h = float(self.bulk.a0)
        self.d = math.sqrt(3)*self.b/2
        if not np.isclose(self.h/self.b, math.sqrt(2/3), rtol=2e-8):
            raise ValueError("this implementation requires cubic Q_bulk=0")
        if tolerance <= 0 or max_ring < 5 or max_modes < 3:
            raise ValueError("invalid row-series controls")
        if fixed_ring is not None and not 1 <= fixed_ring < max_ring:
            raise ValueError("fixed ring needs an extra validation neighborhood")
        self.rho_bulk = float(self.bulk.full_environment_density(self.h,0.).value)
        self.embedding = self.bulk.embedding
        self.angular_weights = np.repeat([surface.vector_amplitude_ev,
            surface.quadrupole_amplitude_ev,surface.amplitude_ev], [3,5,7])
        j,l = np.meshgrid(np.arange(-max_ring,max_ring+1),
                          np.arange(-max_ring,max_ring+1),indexing="ij")
        keep = (j!=0)|(l!=0)
        j,l = j[keep],l[keep]
        y,z = self.d*(j+l/3),self.h*l
        radius = np.hypot(y,z)
        ring = np.maximum(abs(j),abs(l))
        coefficients = []; small = 0
        for m in range(1,max_modes+1):
            g = 2*np.pi*m/self.b
            pair = 4*self.bulk.p.epsilon*(self.bulk.p.sigma_lj**12*
                _mode_coefficient_and_deda(radius,6,m,self.b)[0]
                - self.bulk.p.sigma_lj**6*_mode_coefficient_and_deda(radius,3,m,self.b)[0])
            den = 2*self.bulk.density_params.C_rho/self.b * exponential_row_transform_derivatives(
                radius,g,self.bulk.density_params.kappa)[0]
            angular = [moment_mode(y,z,g,kappa=surface.angular.kappa,
                amplitude=surface.angular.amplitude,b=self.b,rank=r) for r in (1,2,3)]
            c = np.column_stack([pair,den,*angular]) * ((-1.)**(m*(j+l)))[:,None]
            coefficients.append(c)
            envelope = float(np.max(np.abs(c))*max(1.,g*g))
            small = small+1 if envelope < tolerance/100 else 0
            if small == 2:
                break
        else:
            raise ArithmeticError("nonlinear row reciprocal modes exhausted")
        coef = np.array(coefficients)  # mode,row,channel
        g = 2*np.pi*np.arange(1,len(coef)+1)/self.b
        # Phase-independent channel bounds, not a harmonic-only tail estimate.
        # Delta values <=2 sum|c|; derivatives use g and g^2 respectively.
        bounds = np.array([2*np.sum(abs(coef)*g[:,None,None]**order,axis=0)
                           for order in (0,1,2)])
        selected = fixed_ring
        if selected is None:
            for candidate in range(3,max_ring-1):
                if np.max(np.sum(bounds[:,ring>candidate],axis=1)) < tolerance:
                    selected = candidate
                    break
        if selected is None:
            raise ArithmeticError("nonlinear transverse series exhausted; increase max_ring")
        chosen = ring<=selected
        self.j,self.l = j[chosen],l[chosen]
        self.z = z[chosen]
        self.coefficients = coef[:,chosen]
        self.g = g
        self.diagnostics = dict(rings=selected,validation_rings=max_ring,
            reciprocal_modes=len(g),rows=int(sum(chosen)),last_mode_channel_envelope=envelope,
            omitted_validation_channel_bounds=np.sum(bounds[:,~chosen],axis=1).tolist(),
            tail_status="absolute channel bounds within extra rings; infinite remainder requires refinement")

    def row_changes(self, row, displacement, derivative=0):
        """Independent direct reciprocal evaluation of ONE row contribution."""
        if derivative not in (0,1,2):
            raise ValueError("derivative order must be 0,1,2")
        phase = np.exp(1j*self.g*displacement)
        factor = phase-1 if derivative==0 else (1j*self.g)**derivative*phase
        return np.real(factor @ self.coefficients[:,row])

    def periodic_cell(self, shape):
        return PeriodicScrewCell(self, shape)


class PeriodicScrewCell:
    """Periodic j,l cell with explicit affine shear and integer row relabeling.

    Energy is eV per line repeat b*L0 (NOT a finite activation event).
    Atomic sites, not a continuum numerical mesh, define this state space.
    """
    def __init__(self, rows, shape):
        self.rows=rows
        if len(shape)!=2 or any(int(n)!=n or n<4 for n in shape):
            raise ValueError("two periodic cross-section sizes >=4 required")
        self.shape=tuple(int(n) for n in shape)
        self.size=int(np.prod(self.shape))
        self._cached_gamma=None
        self.reference_channels = rows.coefficients.sum(axis=(0,1)).real

    def _operators(self,gamma):
        if self._cached_gamma == gamma:
            return self._cached_operators
        r=self.rows; operators=[]
        for g,c in zip(r.g,r.coefficients):
            current = c*np.exp(1j*g*gamma*r.z)[:,None]
            orders=[]
            for power in (0,1,2):
                stencil=np.zeros(self.shape+(17,),complex)
                np.add.at(stencil,((-r.j)%self.shape[0],(-r.l)%self.shape[1]),
                          current*(1j*g*r.z[:,None])**power)
                orders.append(np.fft.fft2(stencil,axes=(0,1)))
            operators.append(orders)
        self._cached_gamma=gamma; self._cached_operators=operators
        return operators

    @staticmethod
    def _apply(operator,values):
        return np.fft.ifft2(operator*np.fft.fft2(values,axes=(0,1)),axes=(0,1))

    def evaluate(self,u,gamma=0.,*,direction=None,gamma_direction=0.,fields=False):
        u=np.asarray(u,float)
        if u.shape!=self.shape or not np.all(np.isfinite(u)) or not np.isfinite(gamma):
            raise ValueError("finite displacement field with declared shape required")
        v=None if direction is None else np.asarray(direction,float)
        if not np.isfinite(gamma_direction):
            raise ValueError("finite affine Hessian direction required")
        if v is not None and (v.shape!=self.shape or not np.all(np.isfinite(v))):
            raise ValueError("finite Hessian direction with the same shape required")
        operators=self._operators(float(gamma))
        channels=np.zeros(self.shape+(17,)); tangent=np.zeros_like(channels)
        records=[]
        for g,(op,op1,op2) in zip(self.rows.g,operators):
            z=np.exp(1j*g*u)[...,None]; zbar=z.conj()
            az=self._apply(op,z)
            channels+=(zbar*az).real
            if v is not None:
                dz=1j*g*z*v[...,None]; db=dz.conj()
                da=self._apply(op,dz)+gamma_direction*self._apply(op1,z)
                tangent+=(db*az+zbar*da).real
                records.append((g,op,op1,op2,z,zbar,az,dz,db,da))
            else:
                records.append((g,op,op1,op2,z,zbar,az,None,None,None))
        channels-=self.reference_channels
        rho=self.rows.rho_bulk+channels[...,1]
        if np.min(rho)<=0:
            raise ValueError("nonpositive actual site density; no repair applied")
        F=self.rows.embedding
        weights=np.empty_like(channels); weights[...,0]=.5
        weights[...,1]=F.first_derivative(rho)
        weights[...,2:]=2*self.rows.angular_weights*channels[...,2:]
        embedding=F.value(rho)-F.value(self.rows.rho_bulk)
        energy=.5*channels[...,0].sum()+embedding.sum()+np.sum(
            self.rows.angular_weights*channels[...,2:]**2)
        grad=np.zeros(self.shape); grad_gamma=0.
        hv=np.zeros(self.shape); hv_gamma=0.
        if v is not None:
            dw=np.zeros_like(weights); dw[...,1]=F.second_derivative(rho)*tangent[...,1]
            dw[...,2:]=2*self.rows.angular_weights*tangent[...,2:]
        reverse=np.ix_((-np.arange(self.shape[0]))%self.shape[0],
                       (-np.arange(self.shape[1]))%self.shape[1])
        for g,op,op1,op2,z,zb,az,dz,db,da in records:
            adj=self._apply(op[reverse],weights*zb)
            grad+=np.real(1j*g*(z*adj-weights*zb*az)).sum(axis=-1)
            ag=self._apply(op1,z)
            grad_gamma+=float(np.sum(weights*(zb*ag).real))
            if v is not None:
                dadj=self._apply(op[reverse],dw*zb+weights*db)+gamma_direction*self._apply(op1[reverse],weights*zb)
                hv+=np.real(1j*g*(dz*adj+z*dadj-dw*zb*az-weights*db*az-weights*zb*da)).sum(axis=-1)
                dag=self._apply(op1,dz)+gamma_direction*self._apply(op2,z)
                hv_gamma+=float(np.sum(dw*(zb*ag).real+weights*(db*ag+zb*dag).real))
        out=dict(energy=float(energy),gradient=grad,affine_derivative=grad_gamma,
                 minimum_site_density=float(np.min(rho)))
        if fields:
            out.update(channels=channels,density=rho)
        if v is not None:
            out.update(hessian_vector=hv,affine_hessian_vector=hv_gamma)
        return out

    def direct_channels(self,u,gamma=0.):
        """Slow independently assembled real-index convolution, for validation."""
        u=np.asarray(u,float); output=np.zeros(self.shape+(17,))
        for row,(j,l,z) in enumerate(zip(self.rows.j,self.rows.l,self.rows.z)):
            delta=np.roll(u,(-int(j),-int(l)),axis=(0,1))-u+gamma*z
            for g,c in zip(self.rows.g,self.rows.coefficients[:,row]):
                output+=np.real(np.expm1(1j*g*delta)[...,None]*c)
        return output

    def winding(self,u):
        """Burgers winding per logical plaquette; integer row relabeling invariant.

        A b-valued step alone has zero winding and is NOT a resolved core.
        Aliasing at an exactly half-period bond is explicitly diagnosed.
        """
        phase=2*np.pi*np.asarray(u)/self.rows.b
        first=np.angle(np.exp(1j*(np.roll(phase,-1,axis=0)-phase)))
        second=np.angle(np.exp(1j*(np.roll(phase,-1,axis=1)-phase)))
        circulation=first+np.roll(second,-1,axis=0)-np.roll(first,-1,axis=1)-second
        return dict(charge=np.rint(circulation/(2*np.pi)).astype(int),
                    integer_residual=float(np.max(abs(circulation/(2*np.pi)-np.rint(circulation/(2*np.pi))))),
                    minimum_half_period_margin=float(np.pi-max(np.max(abs(first)),np.max(abs(second)))))

    def dipole_seed(self,separation_over_b):
        """Declared opposite screw pair, not an inferred Al defect population."""
        if not np.isfinite(separation_over_b) or separation_over_b<=0:
            raise ValueError("finite positive initial separation required")
        j,l=np.meshgrid(np.arange(self.shape[0])-(self.shape[0]-1)/2,
                        np.arange(self.shape[1])-(self.shape[1]-1)/2,indexing="ij")
        y=self.rows.d*(j+l/3); z=self.rows.h*l
        d=float(separation_over_b)*self.rows.b
        return self.rows.b/(2*np.pi)*(np.arctan2(z,y-d/2)-np.arctan2(z,y+d/2))

    def registry_shear(self,u,gamma=0.):
        """Kinematic winding/intrabond decomposition of PERIODIC cell shear.

        This is not the production chi axial bridge. Integer row relabeling
        changes individual bonds but not their summed registry shear. Existing
        defect content is an initial condition, not newly produced plasticity.
        """
        bond=np.roll(u,-1,axis=1)-u+gamma*self.rows.h
        index=np.floor(bond/self.rows.b+.5)
        intra=bond-self.rows.b*index
        registry=float(self.rows.b*np.mean(index)/self.rows.h)
        intrabond=float(np.mean(intra)/self.rows.h)
        return dict(affine_shear=float(gamma),registry_shear=registry,intrabond_shear=intrabond,
                    decomposition_residual=float(gamma-registry-intrabond))

    def relax(self,initial,*,gamma=0.,shear_stress=None,force_tolerance=2e-8,max_iterations=1200):
        """Local STATIC relaxation, no optimizer clock/velocity.

        shear_stress=None fixes affine strain gamma. Otherwise stress is in
        eV/L0^3 and gamma is another unknown: G=E-tau*V_cell*gamma. Gamma=0
        is NOT zero stress for a cell containing a dislocation dipole.
        No winding constraints/rates are inserted. Check stability separately.
        """
        start=np.asarray(initial,float).copy()
        if start.shape!=self.shape or not np.all(np.isfinite(start)):
            raise ValueError("finite initial displacement with the declared shape required")
        if (not np.isfinite(gamma) or not np.isfinite(force_tolerance) or force_tolerance<=0
                or max_iterations<1 or (shear_stress is not None and not np.isfinite(shear_stress))):
            raise ValueError("invalid static relaxation controls")
        start-=start.mean()
        stress_control=shear_stress is not None
        scale=self.rows.h*np.sqrt(self.size)
        volume=self.size*self.rows.b*self.rows.d*self.rows.h
        def unpack(x):
            return x[:self.size].reshape(self.shape),x[-1]/scale if stress_control else gamma
        def packed(x,direction=None):
            field,strain=unpack(x)
            v,dg=(None,0.) if direction is None else (direction[:self.size].reshape(self.shape),
                                direction[-1]/scale if stress_control else 0.)
            result=self.evaluate(field,strain,direction=v,gamma_direction=dg)
            g=result['gradient']; g-=g.mean()
            energy=result['energy']
            gradient=g.ravel()
            if stress_control:
                energy-=shear_stress*volume*strain
                gradient=np.r_[gradient,(result['affine_derivative']-shear_stress*volume)/scale]
            if direction is not None:
                h=result['hessian_vector']; h-=h.mean()
                # Only the translation gauge is shifted, not material modes.
                h+=v.mean()
                hv=h.ravel()
                if stress_control:
                    hv=np.r_[hv,result['affine_hessian_vector']/scale]
                return hv
            return energy,gradient
        def fun(x):
            return packed(x)
        x0=np.r_[start.ravel(),gamma*scale] if stress_control else start.ravel()
        sol=minimize(fun,x0,jac=True,method="L-BFGS-B",
            options=dict(gtol=force_tolerance/10,ftol=1e-15,maxiter=max_iterations,maxls=40,maxcor=20))
        x=sol.x.copy(); polish=0
        # Objective stagnation is not force convergence. Analytic Newton-CG,
        # restricted to descent steps, removes roundoff-limited L-BFGS residual.
        for polish in range(5):
            energy,gradient=packed(x)
            if max(abs(gradient))<=force_tolerance:
                break
            op=LinearOperator((len(x),len(x)),matvec=lambda v:packed(x,v),dtype=float)
            step,info=cg(op,-gradient,rtol=1e-5,atol=force_tolerance/100,maxiter=400)
            if info!=0 or gradient@step>=0:
                break
            accepted=False
            for damping in (1.,.5,.25,.125):
                trial=x+damping*step
                e,g=packed(trial)
                if e<=energy+2e-13*max(1.,abs(energy)) and np.linalg.norm(g)<np.linalg.norm(gradient):
                    x=trial; accepted=True; break
            if not accepted:
                break
        state,strain=unpack(x);state=state.copy();state-=state.mean()
        result=self.evaluate(state,strain)
        residual=float(np.max(abs(result['gradient'])))
        mechanical=result['affine_derivative']/volume
        scaled_stress_residual=abs(result['affine_derivative']-shear_stress*volume)/scale if stress_control else 0.
        return dict(displacement=state,energy=result['energy'],affine_derivative=result['affine_derivative'],
                    affine_shear=float(strain),internal_shear_stress=float(mechanical),
                    maximum_force_residual=residual,scaled_stress_residual=float(scaled_stress_residual),
                    force_converged=max(residual,scaled_stress_residual)<=force_tolerance,
                    optimizer_success=bool(sol.success),optimizer_message=str(sol.message),iterations=sol.nit,
                    newton_polish_steps=polish,winding=self.winding(state),physical_time_available=False)

    def minimum_curvature(self,u,gamma=0.,*,stress_control=False,tolerance=1e-7):
        """Smallest static Hessian eigenvalue with the translation gauge removed.

        At imposed stress the affine strain is included with the same scaling
        as relax(). This is local stability, not a global or nucleation proof.
        """
        scale=self.rows.h*np.sqrt(self.size)
        n=self.size+int(stress_control)
        def apply(v):
            field=v[:self.size].reshape(self.shape)
            out=self.evaluate(u,gamma,direction=field,
                              gamma_direction=v[-1]/scale if stress_control else 0.)
            h=out['hessian_vector']; h=h-h.mean()+field.mean()
            return np.r_[h.ravel(),out['affine_hessian_vector']/scale] if stress_control else h.ravel()
        op=LinearOperator((n,n),matvec=apply,dtype=float)
        # Deterministic, nonuniform initial vector (not stochastic sampling).
        value,vector=eigsh(op,k=1,which="SA",tol=tolerance,
                          v0=np.cos(np.arange(n)*.371),maxiter=800)
        overlap=abs(np.sum(vector[:self.size,0]))/np.sqrt(self.size)
        if overlap>.5:
            # If every material mode is stiffer than the gauge shift, the
            # lowest returned value would be the artificial translation mode.
            # Keep the true nongauge eigenpair, never report that shift as stiffness.
            values,vectors=eigsh(op,k=2,which="SA",tol=tolerance,
                                v0=np.cos(np.arange(n)*.371),maxiter=800)
            overlaps=abs(np.sum(vectors[:self.size],axis=0))/np.sqrt(self.size)
            candidates=np.flatnonzero(overlaps<.5)
            if not len(candidates):
                raise ArithmeticError("failed to isolate a nongauge curvature mode")
            index=candidates[np.argmin(values[candidates])]
            value,vector=values[index:index+1],vectors[:,index:index+1]
            overlap=overlaps[index]
        residual=np.linalg.norm(apply(vector[:,0])-value[0]*vector[:,0])
        return dict(minimum_eigenvalue=float(value[0]),eigen_residual=float(residual),
                    translation_gauge_eigenvalue=1.,translation_overlap=float(overlap),
                    affine_strain_included=stress_control)
