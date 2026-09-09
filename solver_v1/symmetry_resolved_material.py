"""One-amplitude, nested analytic Eg radial-quadrupole research extension.

At cubic FCC, the strain derivative of a radial STF rank2 environment splits
into Eg (tetragonal) and T2g (cubic off-diagonal shear) irreducible channels.
Use the ALREADY PRESENT odd/even decay scales, not a new fitted length. If
Q_k is the full per-site quadrupole and B_k its signed T2g response, define

 Q_E = [Q_odd - (B_odd/B_even) Q_even] / N,
 N^2 = 2 ||d_gamma(Q_odd-eta Q_even)||^2.

gamma is the existing engineering simple shear of the [111] plane frame.
This normalization fixes only an amplitude gauge. The combination has zero
T2g/hydro response but an independent Eg response. Add D_E sum_i ||Q_E,i||^2,
D_E>=0, AFTER each site's full environment sum. Only ONE new energy amplitude
is fitted; all directions/normalizations derive from FCC symmetry, not yield.
Every plane coefficient/derivative is the existing analytic exponential
Poisson sum. The common-range/coalescent gauge is explicitly rejected rather
than divided by zero. The unchanged family is recovered exactly at D_E=0.
"""
from functools import lru_cache
from types import SimpleNamespace

import numpy as np

from .angular_environment_reference import AngularInterfaceInvariant
from .range_resolved_material import build_range_surface, RangeObservationCache
from .vector_interface_reference import FullRegistryInterface, VectorMomentPlane, VectorInterfaceEvaluation, _embedding_jet
from .tail_constrained_material import TailBulkCoefficientBasis, normalized_amplitude


@lru_cache(maxsize=192)
def signed_quadrupole_response(decay):
    model=build_range_surface(3.,4.,float(decay),np.ones(8),tolerance=2e-13)
    invariant=model.surface.quadrupole
    derivative=np.zeros((2,3,3)); small=0
    for layer in range(1,invariant.max_layers+1):
        values,_,_=invariant._plane(layer*invariant.h, invariant.interface._baseline_delta(layer))
        term=2*layer*invariant.h*values[1:3]
        derivative+=term
        small=small+1 if np.max(abs(term))<2e-13 else 0
        if small>=4: break
    else: raise ArithmeticError('signed cubic response layer sum did not converge')
    normal_shape=np.diag([-1/3,-1/3,2/3])
    signed=float(np.sum(derivative[0]*normal_shape)/np.sum(normal_shape**2))
    if np.max(abs(derivative[0]-signed*normal_shape))>1e-10*max(1.,abs(signed)):
        raise ArithmeticError('cubic normal response is not the T2g symmetry channel')
    return signed,derivative


def quadrupole_channel_gauge(odd,even):
    Bs,Ds=signed_quadrupole_response(float(odd))
    Bl,Dl=signed_quadrupole_response(float(even))
    if abs(Bl)<1e-12: raise ValueError('unresolved T2g denominator')
    eta=Bs/Bl; combined=Ds-eta*Dl
    norm=float(np.sqrt(2*np.sum(combined[1]**2)))
    if norm<1e-9: raise ValueError('coalescent radial quadrupole gauge is not identifiable')
    return dict(eta=float(eta),normalization=norm,
        normal_response_residual=float(np.max(abs(combined[0]))/norm),
        bulk_column=np.r_[0.,0.,0.,2*np.sum(combined[0]**2)/norm**2,1.])


class CombinedQuadrupolePlane:
    def __init__(self, face, odd, even, *, tolerance=2e-11):
        self.gauge=quadrupole_channel_gauge(odd,even)
        self.eta=self.gauge['eta'];self.normalization=self.gauge['normalization']
        internal=tolerance*self.normalization/(2*(1+abs(self.eta)))
        self.parts=[VectorMomentPlane(AngularInterfaceInvariant(face,rank=2,
            angular_decay=k,tolerance=internal)) for k in (odd,even)]
        self.invariant=SimpleNamespace(rank=2,tolerance=tolerance,max_layers=128)

    def evaluate(self,d,delta):
        first=self.parts[0].evaluate(d,delta);second=self.parts[1].evaluate(d,delta)
        return ((first[0]-self.eta*second[0])/self.normalization,
                max(first[1],second[1]),
                (first[2]+abs(self.eta)*second[2])/self.normalization)


class SymmetryResolvedInterface:
    """Static research surface only, with explicit non-colliding components."""
    def __init__(self,decays,coefficients,*,tolerance=2e-11):
        c=np.asarray(coefficients,float)
        if c.shape!=(9,) or np.any(~np.isfinite(c)) or c[8]<0:
            raise ValueError('nine finite coefficients and D_E>=0 required')
        self.decays=tuple(decays);self.coefficients=c.copy()
        self.base=build_range_surface(*decays,c[:8],tolerance=tolerance)
        self.h=self.base.h;self.face=self.base.face;self.geometry=self.base.geometry
        self.moment=CombinedQuadrupolePlane(self.face,decays[1],decays[2],tolerance=tolerance)
        # A separate baseline cache: the old rank2 and new combined rank2 must
        # NEVER share the legacy ('angular',rank,layer) keys.
        self.aux=FullRegistryInterface(self.base.surface)

    def extra_jet(self,state):
        return self.aux._angular(np.asarray(state),self.moment)

    def evaluate(self,state):
        value=self.base.evaluate(state)
        if self.coefficients[8]==0: return value
        jet,*diag=self.extra_jet(state)
        components=dict(value.components)
        components['angular_Eg_radial']=self.coefficients[8]*jet
        diagnostics=dict(value.diagnostics);diagnostics['angular_Eg_radial']=diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()),components,diagnostics)


class SymmetryObservationCache:
    def __init__(self,observations):
        self.observations=tuple(observations);self.legacy=RangeObservationCache(observations)

    @lru_cache(maxsize=128)
    def extra(self,odd,even):
        probe=SymmetryResolvedInterface((3.,odd,even),np.ones(9))
        bulk=probe.moment.gauge['bulk_column'];jets={}
        for obs in self.observations:
            if obs.bulk_index is None and obs.state not in jets:
                jets[obs.state]=probe.extra_jet(obs.state)[0]
        return np.array([bulk[o.bulk_index] if o.bulk_index is not None else
                         np.asarray(o.jet_weights)@jets[o.state] for o in self.observations])

    def matrix(self,decays):
        return np.column_stack([self.legacy.matrix(decays),self.extra(float(decays[1]),float(decays[2]))])


class SymmetryTailBulkBasis(TailBulkCoefficientBasis):
    def __init__(self,decays,*,radius=12.):
        super().__init__(decays,radius=radius)
        self.gauge=quadrupole_channel_gauge(decays[1],decays[2])
        self.extra_gradient=(self._moment_gradient(2,decays[1])
            -self.gauge['eta']*self.gradients[-1])/self.gauge['normalization']

    def evaluate(self,q_cubic):
        columns,errors=super().evaluate(q_cubic)
        q=self.wavevector(q_cubic);qn=np.linalg.norm(q)
        variation=np.einsum('n,nij->ij',np.sin(self.R@q),self.extra_gradient)
        def error(k):
            def integral(power):
                return normalized_amplitude(k)*(2*self.tail.exponential(k,1+power)
                                                 +k*self.tail.exponential(k,2+power))
            return min(integral(0),qn*integral(1))
        dl=(error(self.decays[1])+abs(self.gauge['eta'])*error(self.decays[2]))/self.gauge['normalization']
        extra_error=2*(2*np.linalg.norm(variation,2)*dl+dl*dl)
        return (np.concatenate([columns,(2*variation@variation.T)[None]],axis=0),
                np.r_[errors,extra_error])


class CubicDensityBasis:
    """The next one-coefficient embedding ablation: H(x)=(x-1)^3.

    H(1)=H'(1)=H''(1)=0, H(0)=-1. Its -1 isolated-atom contribution MUST enter
    cohesion. No change to the density kernel or harmonic force constants at
    the reference density. Signed coefficient requires total-energy validation;
    it is not asserted to be a globally convex embedding function.
    """
    def __init__(self,rho_ref): self.rho_ref=float(rho_ref)
    def value(self,rho): return (np.asarray(rho)/self.rho_ref-1)**3
    def first_derivative(self,rho): return 3*(np.asarray(rho)/self.rho_ref-1)**2/self.rho_ref
    def second_derivative(self,rho): return 6*(np.asarray(rho)/self.rho_ref-1)/self.rho_ref**2


class CubicSymmetryInterface(SymmetryResolvedInterface):
    def __init__(self,decays,coefficients,*,tolerance=2e-11):
        c=np.asarray(coefficients,float)
        if c.shape!=(10,) or np.any(~np.isfinite(c)):
            raise ValueError('ten finite research coefficients required')
        super().__init__(decays,c[:9],tolerance=tolerance)
        self.coefficients=c.copy()
        self.cubic_amplitude=float(c[9])
        self.cubic_basis=CubicDensityBasis(self.face.bulk.embedding.rho_ref)

    def evaluate(self,state):
        value=super().evaluate(state)
        if self.cubic_amplitude==0: return value
        changes,*diag=self.base._density(np.asarray(state,float))
        jet=_embedding_jet(changes,self.face.rho_bulk,self.cubic_basis)
        components=dict(value.components);components['embedding_cubic_density']=self.cubic_amplitude*jet
        diagnostics=dict(value.diagnostics);diagnostics['embedding_cubic_density']=diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()),components,diagnostics)


class CubicSymmetryObservationCache(SymmetryObservationCache):
    @lru_cache(maxsize=128)
    def cubic(self,scalar):
        probe=build_range_surface(scalar,4.,5.,np.ones(8));face=probe.face
        emb=CubicDensityBasis(face.bulk.embedding.rho_ref)
        bulk=np.array([0.,-1.,0.,0.,0.]);jets={}
        for o in self.observations:
            if o.bulk_index is None and o.state not in jets:
                jets[o.state]=_embedding_jet(probe._density(np.asarray(o.state))[0],face.rho_bulk,emb)
        return np.array([bulk[o.bulk_index] if o.bulk_index is not None else
                         np.asarray(o.jet_weights)@jets[o.state] for o in self.observations])

    def matrix(self,decays):
        return np.column_stack([super().matrix(decays),self.cubic(float(decays[0]))])


class CubicSymmetryTailBulkBasis(SymmetryTailBulkBasis):
    def evaluate(self,q_cubic):
        columns,tail=super().evaluate(q_cubic)
        return np.concatenate([columns,np.zeros((1,3,3))]),np.r_[tail,0.]
