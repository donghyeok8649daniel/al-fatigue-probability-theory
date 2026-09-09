"""Full (opening, registry-x, registry-y) STATIC interface audit.

This exposes the omitted registry direction of the existing analytic candidate;
it does not change its parameters, production scalar PDE or kinetic status.
Coordinates are reduced by the caller's fixed L0; energy is eV/interface cell.
The separately named Mishin adapter is TARGET/VALIDATION ONLY.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from .angular_environment_reference import (
    _vector_fourier, _second_fourier, _tensor_fourier,
)
from .fcc111_active_interface import _power_mean_difference
from .fcc111_lattice_sum import (
    plane_power_sum_reciprocal, plane_exponential_sum_reciprocal,
    triangular_reciprocal_shells,
)
from .nonlinear_fcc_screw import stf_basis


# Internal analytic jet: value, a,x,y, aa,ax,ay,xx,xy,yy.
HESSIAN_INDICES = np.array([[4, 5, 6], [5, 7, 8], [6, 8, 9]])


@dataclass(frozen=True)
class VectorInterfaceEvaluation:
    energy: float
    gradient: np.ndarray
    hessian: np.ndarray
    components: dict
    diagnostics: dict

    @classmethod
    def from_jet(cls, jet, components=None, diagnostics=None):
        if np.asarray(jet).shape != (10,) or not np.all(np.isfinite(jet)):
            raise FloatingPointError(f'nonfinite/invalid interface jet; components={components}; diagnostics={diagnostics}')
        return cls(float(jet[0]), np.array(jet[1:4]),
                   np.array(jet[HESSIAN_INDICES]), components or {}, diagnostics or {})


def _kernel_jet(value):
    return np.r_[value.value, value.d_d, value.grad_delta, value.d2_dd,
                 value.mixed_d_delta, value.hess_delta[0], value.hess_delta[1, 1]]


def _square_jet(values):
    """Jet of sum_depth ||moment||^2 AFTER each site's environment sum."""
    t = values[:, 0]
    g = values[:, 1:4]
    h = values[:, HESSIAN_INDICES]
    grad = 2*np.einsum('rc,ric->i', t, g)
    hess = 2*(np.einsum('ric,rjc->ij', g, g) + np.einsum('rc,rijc->ij', t, h))
    return np.r_[np.sum(t*t), grad, hess[0], hess[1, 1:], hess[2, 2]]


def _embedding_jet(changes, bulk_density, embedding):
    density = bulk_density + changes[:, 0]
    if np.any(density <= 0) or np.any(~np.isfinite(density)):
        raise FloatingPointError('nonpositive/nonfinite per-atom environment density')
    first = np.asarray(embedding.first_derivative(density))
    second = np.asarray(embedding.second_derivative(density))
    grad = np.einsum('r,ri->i', first, changes[:, 1:4])
    hess = (np.einsum('r,ri,rj->ij', second, changes[:, 1:4], changes[:, 1:4])
            + np.einsum('r,rij->ij', first, changes[:, HESSIAN_INDICES]))
    value = np.sum(embedding.value(density) - embedding.value(bulk_density))
    return 2*np.r_[value, grad, hess[0], hess[1, 1:], hess[2, 2]]


class VectorMomentPlane:
    """Analytic 2D Fourier derivatives, all three state directions together.

    Cached coefficients depend on actual d, not on a frozen transverse radius.
    Full reciprocal shells and their absolute derivative envelopes are retained.
    The last-shell diagnostic is empirical, not a rigorous infinite-tail bound.
    """
    def __init__(self, invariant):
        self.invariant = invariant

    @lru_cache(maxsize=128)
    def _coefficients(self, d):
        obj = self.invariant
        transform = {1: _vector_fourier, 2: _second_fourier, 3: _tensor_fourier}[obj.rank]
        shells = [np.zeros((1, 2)), *(s.vectors for s in
                  triangular_reciprocal_shells(obj.geometry, 32))]
        retained = []; vectors_all = []; small = 0
        for count, vectors in enumerate(shells):
            c0, c1, c2 = (transform(d, vectors, obj.kappa, n)@stf_basis(obj.rank)
                         for n in (0, 1, 2))
            gx, gy = vectors.T
            cs = np.stack([c0, c1, 1j*gx[:, None]*c0, 1j*gy[:, None]*c0,
                           c2, 1j*gx[:, None]*c1, 1j*gy[:, None]*c1,
                           -gx[:, None]**2*c0, -(gx*gy)[:, None]*c0,
                           -gy[:, None]**2*c0], axis=1)
            cs *= obj.amplitude/obj.geometry.atomic_cell_area
            last = float(np.max(np.sum(np.abs(cs), axis=0)))
            retained.append(cs); vectors_all.append(vectors)
            small = small+1 if count > 0 and last < obj.tolerance*.1 else 0
            if small >= 3:
                return np.concatenate(vectors_all), np.concatenate(retained), count, last
        raise RuntimeError('vector moment reciprocal shell envelope did not converge')

    def evaluate(self, d, delta):
        vectors, coefficients, shells, last = self._coefficients(float(d))
        values = np.real(np.einsum('g,gic->ic', np.exp(1j*(vectors@delta)), coefficients))
        return values, shells, last

    def clear_cache(self):
        self._coefficients.cache_clear()


class FullRegistryInterface:
    """Same rigid half-crystal candidate with explicit vector registry u.

    q=(a,u_x,u_y). EAM remains PER ATOM, with the sum over both symmetric
    half-crystals. The power-kernel G=0 layer sum remains exact Hurwitz zeta.
    No numerical differentiation or new phenomenological energy is used.
    """
    def __init__(self, surface):
        self.surface = surface
        self.face = surface.interface
        self.geometry = self.face.geometry
        self.h = surface.h
        self.moments = [(surface.amplitude_ev, VectorMomentPlane(surface.angular))]
        if surface.vector is not None:
            self.moments.append((surface.vector_amplitude_ev, VectorMomentPlane(surface.vector)))
        if surface.quadrupole is not None:
            self.moments.append((surface.quadrupole_amplitude_ev, VectorMomentPlane(surface.quadrupole)))
        self._baselines = {}

    def _pair(self, q, exponent):
        a, x, y = q; face = self.face
        mean = _power_mean_difference(exponent=exponent, h=self.h,
                                     delta_a=a-self.h, area=self.geometry.atomic_cell_area)
        total = np.zeros(10); total[[0, 1, 4]] = mean
        small = 0; maximum_shells = 0; last = math.inf
        for k in range(1, face.max_layers+1):
            delta = self.geometry.abc_shift(k)
            def plane(d, shift):
                v = plane_power_sum_reciprocal(d, shift, p=exponent,
                    geometry=self.geometry, config=face.bulk.stack_config.reciprocal)
                jet = _kernel_jet(v)
                alpha = 2-2*exponent
                zero = math.pi/(self.geometry.atomic_cell_area*(exponent-1))*d**alpha
                jet[[0, 1, 4]] -= [zero, alpha*zero/d, alpha*(alpha-1)*zero/d**2]
                return jet, v.shells_used
            key = ('pair', exponent, k)
            if key not in self._baselines:
                self._baselines[key] = plane(k*self.h, delta)[0][0]
            jet, shells = plane(a+(k-1)*self.h, delta+np.array([x, y]))
            jet[0] -= self._baselines[key]
            total += k*jet
            maximum_shells = max(maximum_shells, shells)
            last = float(np.max(np.abs(k*jet)))
            small = small+1 if last < face.tolerance else 0
            if small >= face.consecutive_small_layers:
                return total, k, maximum_shells, last
        raise RuntimeError('vector interface pair layers did not converge')

    def _density(self, q):
        a, x, y = q; face = self.face; den = face.bulk.density_params
        changes = []; small = 0; maximum_shells = 0
        for k in range(1, face.max_layers+1):
            delta = self.geometry.abc_shift(k)
            def plane(d, shift):
                v = plane_exponential_sum_reciprocal(d, shift, kappa=den.kappa,
                    amplitude=den.C_rho, geometry=self.geometry,
                    config=face.bulk.stack_config.reciprocal)
                return _kernel_jet(v), v.shells_used
            key = ('density', k)
            if key not in self._baselines:
                self._baselines[key] = plane(k*self.h, delta)[0][0]
            jet, shells = plane(a+(k-1)*self.h, delta+np.array([x, y]))
            jet[0] -= self._baselines[key]; changes.append(jet)
            maximum_shells = max(maximum_shells, shells)
            last = float(np.max(np.abs(jet)))
            small = small+1 if last < face.tolerance*max(1., face.rho_bulk) else 0
            if small >= face.consecutive_small_layers:
                return np.cumsum(changes[::-1], axis=0)[::-1], k, maximum_shells, last
        raise RuntimeError('vector interface density layers did not converge')

    def _angular(self, q, moment):
        a, x, y = q; obj = moment.invariant
        changes = []; small = 0; maximum_shells = 0
        for k in range(1, obj.max_layers+1):
            delta = self.geometry.abc_shift(k)
            key = ('angular', obj.rank, k)
            if key not in self._baselines:
                self._baselines[key] = moment.evaluate(k*self.h, delta)[0][0].copy()
            jet, shells, envelope = moment.evaluate(a+(k-1)*self.h, delta+np.array([x, y]))
            jet[0] -= self._baselines[key]; changes.append(jet)
            maximum_shells = max(maximum_shells, shells)
            last = float(np.max(np.abs(jet)))
            small = small+1 if last < obj.tolerance else 0
            if small >= 4:
                depths = np.cumsum(changes[::-1], axis=0)[::-1]
                return 2*_square_jet(depths), k, maximum_shells, last
        raise RuntimeError('vector angular interface neighborhoods did not converge')

    def evaluate(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (3,) or np.any(~np.isfinite(q)) or q[0] <= 0:
            raise ValueError('finite q=(positive opening, registry-x, registry-y) required')
        components = {}; diagnostics = {}
        p = self.face.bulk.p
        for exponent, coefficient in ((6, 4*p.epsilon*p.sigma_lj**12),
                                      (3, -4*p.epsilon*p.sigma_lj**6)):
            jet, *diag = self._pair(q, exponent)
            components['pair_'+str(exponent)] = coefficient*jet
            diagnostics['pair_'+str(exponent)] = diag
        changes, *diag = self._density(q)
        components['embedding'] = _embedding_jet(changes, self.face.rho_bulk, self.face.bulk.embedding)
        diagnostics['density'] = diag
        for amplitude, moment in self.moments:
            jet, *diag = self._angular(q, moment)
            name = 'angular_'+str(moment.invariant.rank)
            components[name] = amplitude*jet; diagnostics[name] = diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()), components, diagnostics)


class MishinVectorInterfaceReference:
    """Independent published-EAM SOURCE on identical geometry, never PDE input.

    Uses the declared interpolation and cutoff of that reference potential.
    Physical source distances are angstrom, converted to the caller's fixed L0.
    """
    def __init__(self, reference, length_scale_angstrom):
        if reference.interpolation != 'cubic' or length_scale_angstrom <= 0:
            raise ValueError('positive fixed L0 and cubic SOURCE interpolation required')
        self.reference = reference
        self.length = float(length_scale_angstrom)
        self.h = reference.h/self.length

    def _plane(self, d, shift):
        obj = self.reference
        xy = obj.R + shift
        v = np.column_stack([np.full(len(xy), d), xy])
        radii = np.linalg.norm(v, axis=1)
        select = (radii > 0) & (radii < obj.r[-1])
        r = radii[select]; n = v[select]/r[:, None]
        z = obj._rphi(r); zp = obj._rphi(r, 1); zpp = obj._rphi(r, 2)
        def pack(value, first, second):
            g = np.einsum('r,ri->i', first, n)
            h = (np.einsum('r,ri,rj->ij', second-first/r, n, n)
                 + np.eye(3)*np.sum(first/r))
            result = np.r_[np.sum(value), g, h[0], h[1, 1:], h[2, 2]]
            if not np.all(np.isfinite(result)):
                raise FloatingPointError(f'SOURCE plane jet nonfinite: d={d}, shift={shift}, r={r}, '
                    f'value={value}, first={first}, second={second}, n={n}, result={result}')
            return result
        return (pack(z/r, zp/r-z/r**2, zpp/r-2*zp/r**2+2*z/r**3),
                pack(obj._rho(r), obj._rho(r, 1), obj._rho(r, 2)))

    def evaluate(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (3,) or np.any(~np.isfinite(q)) or q[0] <= 0:
            raise ValueError('finite three-coordinate interface state required')
        obj = self.reference; a, x, y = q*self.length
        count = int(np.ceil((obj.cutoff+abs(a-obj.h))/obj.h))+1
        pair = np.zeros(10); changes = []
        for k in range(1, count+1):
            delta = obj.geometry.abc_shift(k)
            bp, brho = obj.plane(k*obj.h, delta)
            pj, rj = self._plane(a+(k-1)*obj.h, delta+np.array([x, y]))
            pj[0] -= bp; rj[0] -= brho
            pair += k*pj; changes.append(rj)
        depths = np.cumsum(changes[::-1], axis=0)[::-1]
        class SourceEmbedding:
            value = staticmethod(obj.F)
            first_derivative = staticmethod(lambda r: obj._F(r, 1))
            second_derivative = staticmethod(lambda r: obj._F(r, 2))
        emb = _embedding_jet(depths, obj._rho_bulk, SourceEmbedding)
        scale = np.r_[1., np.full(3, self.length), np.full(6, self.length**2)]
        return VectorInterfaceEvaluation.from_jet((pair+emb)*scale,
            {'source_pair': pair*scale, 'source_embedding': emb*scale},
            {'source_sha256': obj.sha256, 'source_only': True, 'layers': count})
