"""Fourier--Bessel representation of unchanged pure-Si SW (111) planes.

Research reference only: no MACE conversion, material fit, finite-T PMF,
first-initiation basin, mobility, physical clock, or production registration.
J0/J1/J2 Hankel integrals retain the original compact SW radial functions.
All moments are summed per site BEFORE the angular quadratic is evaluated.
Units: Angstrom, eV; interface energy per primitive surface cell.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product

import numpy as np
from scipy.special import jv, roots_legendre

from .fcc111_geometry import fcc111_geometry
from .fcc111_lattice_sum import triangular_reciprocal_shells
from .silicon_environment_research import Jet, _parameters, _site_energy


@dataclass(frozen=True)
class PlaneMoments:
    """pair, rho, sum(w^2), vector[3], STF tensor[9]; derivatives in (d,x,y)."""
    value: np.ndarray
    gradient: np.ndarray
    hessian: np.ndarray

    def mapped(self, sign):
        return PlaneMoments(self.value, self.gradient*sign, self.hessian*sign**2)

    def jets(self):
        return [Jet(v, g, h) for v, g, h in zip(self.value, self.gradient, self.hessian)]


def _zero():
    return PlaneMoments(np.zeros(15), np.zeros((15, 3)), np.zeros((15, 3, 3)))


def _combine(moments):
    result = _zero()
    for m in moments:
        result = PlaneMoments(result.value+m.value, result.gradient+m.gradient,
                              result.hessian+m.hessian)
    return result


def site_energy(moments, parameters):
    """Original SW centered-site identity, including neighbor self subtraction."""
    p = _parameters(parameters)
    pair, rho, diagonal, *other = moments.jets()
    vector, tensor = other[:3], other[3:]
    c = p['costheta0']
    angular = (.5*sum(t*t for t in tensor) - c*sum(v*v for v in vector)
               + (1/6+c*c/2)*rho*rho - .5*(1-c)**2*diagonal)
    return .5*pair+p['lambda']*p['epsilon']*angular


def _point_moments(vector, p):
    r = np.linalg.norm(vector)
    cutoff = p['a']*p['sigma']
    if r <= 1e-14:
        raise ValueError('coincident neighbor')
    if r >= cutoff:
        return np.zeros(15)
    z = p['sigma']/(r-cutoff)
    w = np.exp(p['gamma']*z)
    pair = p['A']*p['epsilon']*(p['B']*(p['sigma']/r)**p['p']
                               -(p['sigma']/r)**p['q'])*np.exp(z)
    n = vector/r
    return np.r_[pair, w, w*w, w*n, (w*(np.outer(n, n)-np.eye(3)/3)).ravel()]


class SWPlaneBessel:
    """Complete reciprocal metric disk and independently variable quadrature.

    shell_index and nodes are resolution controls, NOT certified tail bounds.
    d=0 is excluded: the pair kernel has a nonintegrable self singularity.
    Same-plane neighbors are handled separately by an exact finite SW sum.
    """
    def __init__(self, parameters, lattice, *, shell_index=96, nodes=512):
        self.p = _parameters(parameters)
        self.geometry = fcc111_geometry(lattice)
        if int(shell_index) != shell_index or shell_index < 2:
            raise ValueError('integer shell_index >= 2 required')
        if int(nodes) != nodes or nodes < 32:
            raise ValueError('integer quadrature nodes >= 32 required')
        self.shell_index, self.nodes = int(shell_index), int(nodes)
        self.shells = triangular_reciprocal_shells(self.geometry, self.shell_index)
        self.waves = np.r_[0., [s.magnitude for s in self.shells]]
        self.abscissa, self.weights = roots_legendre(self.nodes)

    @lru_cache(maxsize=64)
    def _coefficients(self, distance):
        d = Jet.variable(float(distance), 1, 0)
        cutoff = self.p['a']*self.p['sigma']
        maximum = np.sqrt(cutoff**2-distance**2)
        rho = maximum*(self.abscissa+1)/2
        weight = np.pi*maximum*self.weights*rho
        channels = np.empty((8, 3, self.nodes))
        p = self.p
        for k, rxy in enumerate(rho):
            r = (d*d+rxy*rxy)**.5
            z = p['sigma']/(r-cutoff)
            w = (p['gamma']*z).exp()
            pair = p['A']*p['epsilon']*(p['B']*(p['sigma']/r)**p['p']
                                      -(p['sigma']/r)**p['q'])*z.exp()
            t, u = rxy/r, d/r
            values = [pair, w, w*w, w*u, w*t, w*t*t, w*t*u, w*(u*u-1/3)]
            for i, value in enumerate(values):
                channels[i, :, k] = [value.value, value.gradient[0], value.hessian[0, 0]]
        argument = self.waves[:, None]*rho
        bessel = [jv(order, argument) for order in range(3)]
        # I_phi0,I_w0,I_w2_0,I_wu0,I_wt1,I_wt2_0,I_wt2_2,I_wtu1,I_zz0.
        recipe = [(0,0),(1,0),(2,0),(3,0),(4,1),(5,0),(5,2),(6,1),(7,0)]
        output = np.empty((len(self.waves), 3, 9))
        for i, (channel, order) in enumerate(recipe):
            output[:, :, i] = bessel[order]@(channels[channel]*weight).T
        return output

    def evaluate(self, distance, delta):
        d, delta = float(distance), np.asarray(delta, float)
        if not np.isfinite(d) or d == 0 or delta.shape != (2,) or not np.all(np.isfinite(delta)):
            raise ValueError('nonzero finite d and finite two-dimensional delta required')
        if abs(d) >= self.p['a']*self.p['sigma']:
            return _zero()
        return self._evaluate_cached(d,float(delta[0]),float(delta[1]))

    @lru_cache(maxsize=128)
    def _evaluate_cached(self,d,dx,dy):
        delta = np.array([dx,dy])
        total = _zero()
        coefficients = self._coefficients(d)
        groups = [(np.zeros((1,2)), coefficients[0])]
        groups.extend((s.vectors, c) for s, c in zip(self.shells, coefficients[1:]))
        for vectors, c in groups:
            magnitude = np.linalg.norm(vectors, axis=1)
            direction = np.divide(vectors, magnitude[:, None], out=np.zeros_like(vectors),
                                  where=magnitude[:, None] != 0)
            # First dimension indexes d derivative order (0,1,2).
            f = np.zeros((3, len(vectors), 15), complex)
            f[:, :, :3] = c[:, None, :3]
            f[:, :, 3:5] = -1j*c[:, None, 4, None]*direction[None]
            f[:, :, 5] = c[:, None, 3]
            tensor = np.zeros((3, len(vectors), 3, 3), complex)
            for a in range(2):
                for b in range(2):
                    tensor[:, :, a, b] = (-c[:, None, 6]*direction[None,:,a]*direction[None,:,b]
                        + ((c[:, None, 5]+c[:, None, 6])/2-c[:, None, 1]/3 if a == b else 0))
                tensor[:, :, a, 2] = -1j*c[:, None, 7]*direction[None,:,a]
                tensor[:, :, 2, a] = tensor[:, :, a, 2]
            tensor[:, :, 2, 2] = c[:, None, 8]
            f[:, :, 6:] = tensor.reshape(3, len(vectors), 9)
            phase = np.exp(1j*(vectors@delta))/self.geometry.atomic_cell_area
            value = np.einsum('vf,v->f', f[0], phase)
            gradient = np.zeros((15,3), complex)
            hessian = np.zeros((15,3,3), complex)
            gradient[:,0] = np.einsum('vf,v->f', f[1], phase)
            hessian[:,0,0] = np.einsum('vf,v->f', f[2], phase)
            for a in range(2):
                gradient[:,a+1] = np.einsum('vf,v->f', f[0], phase*1j*vectors[:,a])
                hessian[:,0,a+1] = np.einsum('vf,v->f', f[1], phase*1j*vectors[:,a])
                hessian[:,a+1,0] = hessian[:,0,a+1]
                for b in range(2):
                    hessian[:,a+1,b+1] = np.einsum('vf,v->f', f[0], -phase*vectors[:,a]*vectors[:,b])
            scale = max(1., np.max(abs(value)), np.max(abs(gradient)), np.max(abs(hessian)))
            if max(np.max(abs(value.imag)), np.max(abs(gradient.imag)),
                   np.max(abs(hessian.imag))) > 1e-11*scale:
                raise FloatingPointError('reciprocal symmetry left a nonreal moment')
            total = _combine([total, PlaneMoments(value.real, gradient.real, hessian.real)])
        for array in (total.value,total.gradient,total.hessian):
            array.setflags(write=False)
        return total

    def direct_vectors(self, distance, delta, *, exclude_origin=False):
        delta = np.asarray(delta, float)
        cutoff = self.p['a']*self.p['sigma']
        if abs(distance) >= cutoff:
            return np.empty((0,3))
        direct = np.column_stack((self.geometry.a1, self.geometry.a2))
        span = int(np.ceil((cutoff+np.linalg.norm(delta))/np.linalg.svd(direct, compute_uv=False)[-1]))+1
        vectors = []
        for m, n in product(range(-span, span+1), repeat=2):
            v = np.r_[self.geometry.lattice_vector(m,n)+delta, distance]
            if exclude_origin and m == 0 and n == 0:
                continue
            if np.linalg.norm(v) < cutoff:
                vectors.append(v)
        return np.asarray(vectors).reshape(-1,3)

    def intraplane(self):
        vectors = self.direct_vectors(0., [0.,0.], exclude_origin=True)
        m = _zero()
        for v in vectors:
            m = PlaneMoments(m.value+_point_moments(v,self.p), m.gradient, m.hessian)
        return m


class SW111BesselInterface:
    """Unrelaxed pure-Si diamond (111) shuffle/glide coherent interface.

    Diamond planes: z=(l+3b/4)h111, delta=l*tau, b=0,1.
    Retained coordinates: (additional opening, full lateral x, full lateral y).
    The per-site subtraction covers every changed center on both sides.
    This static rigid interface is NOT a first crack initiation calculation.
    """
    def __init__(self, parameters, lattice, *, cut_kind='shuffle', shell_index=96, nodes=512):
        if cut_kind not in ('shuffle', 'glide'):
            raise ValueError('cut_kind must be shuffle or glide')
        self.plane = SWPlaneBessel(parameters, lattice, shell_index=shell_index, nodes=nodes)
        self.p, self.geometry = self.plane.p, self.plane.geometry
        self.cut_kind = cut_kind
        self.cut = (.375 if cut_kind == 'shuffle' else .875)*self.geometry.h111
        self.gap = (.75 if cut_kind == 'shuffle' else .25)*self.geometry.h111

    def _planes(self, opening):
        radius = 2*self.p['a']*self.p['sigma']+abs(opening)+2*self.geometry.h111
        span = int(np.ceil(radius/self.geometry.h111))+2
        return [(l,b,(l+.75*b)*self.geometry.h111,l*self.geometry.tau)
                for l in range(-span,span+1) for b in range(2)]

    def _center_moments(self, center, planes, q, active):
        _, _, zc, dc = center
        side = zc > self.cut
        parts = []
        for plane in planes:
            _, _, z, delta = plane
            sign = int(z > self.cut)-int(side) if active else 0
            d = z-zc+sign*q[0]
            shift = delta-dc+sign*q[1:]
            if plane[:2] == center[:2]:
                parts.append(self.plane.intraplane())
            else:
                parts.append(self.plane.evaluate(d,shift).mapped(sign))
        return _combine(parts)

    def evaluate(self, state, *, method='bessel'):
        q = np.asarray(state,float)
        if q.shape != (3,) or not np.all(np.isfinite(q)) or q[0]+self.gap <= 0:
            raise ValueError('finite three-coordinate state must keep planes ordered')
        if method not in ('bessel','direct'):
            raise ValueError('method must be bessel or direct')
        planes, cutoff = self._planes(q[0]), self.p['a']*self.p['sigma']
        result = Jet.constant(0.,3)
        for center in planes:
            opposite = [p for p in planes if (p[2] > self.cut) != (center[2] > self.cut)]
            sign = 1 if center[2] < self.cut else -1
            if all(abs(p[2]-center[2]) >= cutoff and
                   abs(p[2]-center[2]+sign*q[0]) >= cutoff for p in opposite):
                continue
            if method == 'bessel':
                now = site_energy(self._center_moments(center,planes,q,True),self.p)
                bulk = site_energy(self._center_moments(center,planes,np.zeros(3),False),self.p)
            else:
                now = self._direct_site(center,planes,q,True)
                bulk = self._direct_site(center,planes,np.zeros(3),False)
            result += now-bulk
        return result

    def _direct_site(self, center, planes, q, active):
        coordinates = [Jet.variable(v,3,i) for i,v in enumerate(q)]
        vectors = []
        for plane in planes:
            sign = int(plane[2] > self.cut)-int(center[2] > self.cut) if active else 0
            dz, delta = plane[2]-center[2], plane[3]-center[3]
            current_d, current_delta = dz+sign*q[0], delta+sign*q[1:]
            for v in self.plane.direct_vectors(current_d,current_delta,
                                                exclude_origin=plane[:2] == center[:2]):
                base = v-sign*q[[1,2,0]]
                vectors.append([base[0]+sign*coordinates[1],
                                base[1]+sign*coordinates[2],
                                base[2]+sign*coordinates[0]])
        return _site_energy(vectors,self.p,'direct')
