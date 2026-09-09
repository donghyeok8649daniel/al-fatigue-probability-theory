"""Static finite-strain validation at fixed analytic material parameters.

Isotropic FCC dilation changes atom positions, NOT L0, density normalization,
decays, energy coefficients or the radial-quadrupole gauge. A 300 K MD box
may be tested with a static potential Hessian but is not an anharmonic free
energy or a kinetic calibration. Production energy/PDE defaults are untouched.
"""
from dataclasses import replace
from types import SimpleNamespace

import numpy as np

from .angular_environment_reference import traceless_second,traceless_third
from .fcc111_full_energy import FullFCC111StackEnergy
from .range_resolved_material import build_range_surface
from .static_bulk_stability import neighbors_in_plane_frame
from .symmetry_resolved_material import quadrupole_channel_gauge
from .tail_constrained_material import TailBulkCoefficientBasis,FCCTailEnvelope,normalized_amplitude


class IsotropicBulkBasis(TailBulkCoefficientBasis):
    """Unweighted H(q) columns [eV/L0²], for the existing quartic family.

    At current x=sum rho/rho_ref, define R=sum(1-cos) Hess x_bond, g=sum sin grad x_bond.
    The scalar columns become
      H_A=-R/sqrt(x)+gg^T/(4*x^(3/2)), H_B=2R,
      H_C=4(x-1)R+2gg^T.
    Angular gauges remain those of the original reference. Q1,Q2,Q3 vanish
    in ANY cubic dilation, so the quadratic-moment formula is still exact.
    The quartic invariant has no harmonic term; the mixed (x-1)||Q3||² term
    DOES have one away from x=1: H_cross=(x-1) H_D3.
    """
    def __init__(self,decays,*,stretch=1.,radius=12.,include_cross=False):
        decays=tuple(map(float,decays))
        if len(decays)!=3 or np.any(~np.isfinite(decays)) or min(decays)<=0:
            raise ValueError('three fixed positive reference decays required')
        if not np.isfinite(stretch) or stretch<=0:raise ValueError('positive cubic dilation required')
        self.decays=decays;self.stretch=float(stretch);self.radius=float(radius)
        self.include_cross=bool(include_cross)
        reference=build_range_surface(*decays,np.ones(8)).face.bulk
        self.reference_density=reference.embedding.rho_ref
        # Clone ONLY geometry; no mutation of the material or reference gauge.
        self.current_bulk=FullFCC111StackEnergy(replace(reference.p,b=reference.p.b*stretch),
            density_params=reference.density_params,embedding=reference.embedding,
            stack_config=reference.stack_config,find_equilibrium=False)
        self.geometry=self.current_bulk.geometry
        current=SimpleNamespace(geometry=self.geometry,a0=self.geometry.h111)
        self.R=neighbors_in_plane_frame(current,radius)
        self.r=np.linalg.norm(self.R,axis=1);self.e=self.R/self.r[:,None]
        self.rr=self.e[:,:,None]*self.e[:,None,:]
        self.tail=FCCTailEnvelope(radius,self.geometry.b)
        self.pair=[self.r[:,None,None]**-14*(168*self.rr-12*np.eye(3)),
                   self.r[:,None,None]**-8*(6*np.eye(3)-48*self.rr)]
        k=decays[0];self.rho_weight=normalized_amplitude(k)*np.exp(-k*self.r)
        self.rho_gradient=-k*self.rho_weight[:,None]*self.e
        self.rho_hessian=self.rho_weight[:,None,None]*(k*k*self.rr
            -k/self.r[:,None,None]*(np.eye(3)-self.rr))
        density=self.current_bulk.full_environment_density(self.geometry.h111,0.)
        self.x=float(density.value/self.reference_density)
        self.density_series_tail=float(density.estimated_tail_absolute/self.reference_density)
        if stretch==1.:
            # The numerator IS the defining reference density in this exact
            # limit; its ratio is identically one, not two independent errors.
            self.x=1.;self.density_series_tail=0.
        if self.x<=0:raise ArithmeticError('positive full environment required')
        self.moments=[(3,decays[1]),(1,decays[1]),(2,decays[2])]
        self.gradients=[self._moment_gradient(rank,k) for rank,k in self.moments]
        self.gauge=quadrupole_channel_gauge(decays[1],decays[2])
        self.extra_gradient=(self._moment_gradient(2,decays[1])
            -self.gauge['eta']*self.gradients[-1])/self.gauge['normalization']

    def evaluate(self,q_cubic):
        columns,errors=super().evaluate(q_cubic);x=self.x
        original=columns.copy();old_errors=errors.copy()
        columns[2]=original[2]/np.sqrt(x)+.125*(x**-1.5-x**-.5)*original[4]
        columns[4]=original[4]+2*(x-1)*original[3]
        errors[2]=old_errors[2]/np.sqrt(x)+.125*abs(x**-1.5-x**-.5)*old_errors[4]
        errors[4]=old_errors[4]+2*abs(x-1)*old_errors[3]
        dx=self.density_series_tail;lower=x-dx
        if lower<=0:raise ArithmeticError('density series does not resolve a positive environment')
        Rnorm=(np.linalg.norm(original[3],2)+old_errors[3])/2
        ggnorm=(np.linalg.norm(original[4],2)+old_errors[4])/2
        errors[2]+=dx*(.5*lower**-1.5*Rnorm+.375*lower**-2.5*ggnorm)
        errors[4]+=4*dx*Rnorm
        q=self.wavevector(q_cubic);qn=np.linalg.norm(q)
        variation=np.einsum('n,nij->ij',np.sin(self.R@q),self.extra_gradient)
        def error(k):
            def integral(power):
                return normalized_amplitude(k)*(2*self.tail.exponential(k,1+power)
                                                +k*self.tail.exponential(k,2+power))
            return min(integral(0),qn*integral(1))
        dl=(error(self.decays[1])+abs(self.gauge['eta'])*error(self.decays[2]))/self.gauge['normalization']
        eextra=2*(2*np.linalg.norm(variation,2)*dl+dl*dl)
        columns=np.concatenate([columns,(2*variation@variation.T)[None],np.zeros((1,3,3))])
        errors=np.r_[errors,eextra,0.]
        if self.include_cross:
            columns=np.concatenate([columns,((x-1)*columns[5])[None]])
            errors=np.r_[errors,abs(x-1)*errors[5]+dx*(np.linalg.norm(columns[5],2)+errors[5])]
        # Neighbor envelopes are analytic. The separate reciprocal density
        # truncation estimate is propagated, exposed and independently checked;
        # it is not promoted to a proof of the entire Brillouin zone.
        return columns,errors

    def direct_sinusoidal_energy(self,coefficients,*,planes,mode,polarization_plane,
                                  amplitude,validation_radius,saturation=0.):
        """Independent finite direct neighbor validation of the same site law.

        NOT a replacement cutoff potential. Radius/amplitude refinement and
        the analytical infinite-neighbor tail comparison remain mandatory.
        Each displaced site sums ALL its included neighbors before F/invariants.
        """
        c=np.asarray(coefficients,float);v=np.asarray(polarization_plane,float)
        if (c.shape!=(11 if self.include_cross else 10,) or np.any(~np.isfinite(c))
                or int(planes)!=planes or planes<3 or int(mode)!=mode or not 0<mode<planes
                or v.shape!=(3,) or not np.isclose(v@v,1.,rtol=1e-12,atol=1e-12)
                or np.any(~np.isfinite(v)) or not np.isfinite(amplitude)
                or not np.isfinite(saturation) or saturation<0):
            raise ValueError('fixed finite coefficients, unit polarization and periodic mode required')
        bulk=SimpleNamespace(geometry=self.geometry,a0=self.geometry.h111)
        R=neighbors_in_plane_frame(bulk,validation_radius)
        layer=np.rint(R[:,2]/self.geometry.h111).astype(int)
        phase=2*np.pi*mode/planes;values=[]
        for site in range(planes):
            shift=amplitude*(np.cos(phase*(site+layer))-np.cos(phase*site))
            positions=R+shift[:,None]*v;r=np.linalg.norm(positions,axis=1)
            weights=[normalized_amplitude(k)*np.exp(-k*r) for k in self.decays]
            x=weights[0].sum()
            Q1=weights[1]@positions
            Q3=traceless_third(np.einsum('n,ni,nj,nk->ijk',weights[1],positions,positions,positions))
            oddQ2=traceless_second(np.einsum('n,ni,nj->ij',weights[1],positions,positions))
            Q2=traceless_second(np.einsum('n,ni,nj->ij',weights[2],positions,positions))
            QE=(oddQ2-self.gauge['eta']*Q2)/self.gauge['normalization']
            I=float(np.sum(Q3*Q3))
            energy=(.5*np.sum(c[0]*r**-12-c[1]*r**-6)
                -c[2]*np.sqrt(x)+c[3]*(x-1)+c[4]*(x-1)**2
                +c[5]*I+c[6]*(Q1@Q1)+c[7]*np.sum(Q2*Q2)+c[8]*np.sum(QE*QE)
                +c[9]*I*I/(1+saturation*I))
            if self.include_cross:energy+=c[10]*(x-1)*I
            values.append(energy)
        return float(np.mean(values))
