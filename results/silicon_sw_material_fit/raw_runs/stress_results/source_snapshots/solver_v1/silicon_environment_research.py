"""Static SW representation/coordinate audit; never a production material.

The same per-site energy is evaluated by explicit triples or scalar/vector/STF
moments. Second-order forward jets differentiate the expressions, not numerical
differences. The finite cutoff is the SW reference's own definition; no existing
infinite LJ/Bessel evaluator, mobility, or probability generator is replaced.
All lengths are Angstrom and energies eV. No fitted or default Si parameters.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass
class Jet:
    value: float
    gradient: np.ndarray
    hessian: np.ndarray

    @classmethod
    def constant(cls, value, dimension):
        return cls(float(value), np.zeros(dimension), np.zeros((dimension, dimension)))

    @classmethod
    def variable(cls, value, dimension, index):
        result = cls.constant(value, dimension)
        result.gradient[index] = 1.
        return result

    def _coerce(self, other):
        return other if isinstance(other, Jet) else Jet.constant(other, len(self.gradient))

    def __add__(self, other):
        other = self._coerce(other)
        return Jet(self.value + other.value, self.gradient + other.gradient,
                   self.hessian + other.hessian)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, -self.gradient, -self.hessian)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) + (-self)

    def __mul__(self, other):
        other = self._coerce(other)
        return Jet(self.value * other.value,
                   self.gradient * other.value + self.value * other.gradient,
                   self.hessian * other.value + self.value * other.hessian
                   + np.outer(self.gradient, other.gradient)
                   + np.outer(other.gradient, self.gradient))

    __rmul__ = __mul__

    def compose(self, value, first, second):
        return Jet(value, first * self.gradient,
                   first * self.hessian + second * np.outer(self.gradient, self.gradient))

    def __pow__(self, power):
        if power == 0:
            return Jet.constant(1., len(self.gradient))
        if power == 1:
            return self
        return self.compose(self.value**power, power * self.value**(power-1),
                            power * (power-1) * self.value**(power-2))

    def __truediv__(self, other):
        return self * self._coerce(other)**-1

    def __rtruediv__(self, other):
        return self._coerce(other) * self**-1

    def exp(self):
        value = np.exp(self.value)
        return self.compose(value, value, value)


def _value(x):
    return x.value if isinstance(x, Jet) else float(x)


def _exp(x):
    return x.exp() if isinstance(x, Jet) else np.exp(x)


def _dot(first, second):
    return sum(a*b for a, b in zip(first, second))


def _parameters(parameters):
    keys = 'epsilon sigma a lambda gamma costheta0 A B p q'.split()
    p = {key: float(parameters[key]) for key in keys}
    if (not np.all(np.isfinite(list(p.values())))
            or any(p[key] <= 0 for key in 'epsilon sigma a gamma A B p'.split())
            or p['lambda'] < 0 or p['q'] < 0 or abs(p['costheta0']) > 1):
        raise ValueError('finite admissible SW reference parameters required')
    return p


def _site_energy(vectors, p, representation):
    """Half pair sum plus centered unordered triples; all per-site terms kept."""
    if representation not in ('direct', 'moments'):
        raise ValueError('representation must be direct or moments')
    cutoff = p['a'] * p['sigma']
    pair, neighbors = 0., []
    for vector in vectors:
        square = _dot(vector, vector)
        if _value(square) <= 1e-20:
            raise ValueError('coincident atoms in reference environment')
        radius = square**.5
        if _value(radius) >= cutoff:
            continue
        z = p['sigma'] / (radius-cutoff)
        radial = p['B'] * (p['sigma']/radius)**p['p'] - (p['sigma']/radius)**p['q']
        pair += .5 * p['A'] * p['epsilon'] * radial * _exp(z)
        neighbors.append((_exp(p['gamma']*z), [x/radius for x in vector]))
    c = p['costheta0']
    if representation == 'direct':
        angular = sum(w*v*(_dot(n, m)-c)**2
                      for j, (w, n) in enumerate(neighbors)
                      for v, m in neighbors[j+1:])
    else:
        rho = sum(w for w, _ in neighbors)
        diagonal = sum(w*w for w, _ in neighbors)
        vector = [sum(w*n[a] for w, n in neighbors) for a in range(3)]
        tensor = [[sum(w*(n[a]*n[b]-(1/3 if a == b else 0))
                       for w, n in neighbors) for b in range(3)] for a in range(3)]
        angular = (.5*sum(x*x for row in tensor for x in row)
                   - c*_dot(vector, vector) + (1/6+c*c/2)*rho*rho
                   - .5*(1-c)**2*diagonal)
    return pair + p['lambda'] * p['epsilon'] * angular


def _state(values, dimension, derivatives):
    q = np.asarray(values, float)
    if q.shape != (dimension,) or not np.all(np.isfinite(q)):
        raise ValueError(f'{dimension} finite coordinates required')
    return ([Jet.variable(x, dimension, i) for i, x in enumerate(q)]
            if derivatives else list(q))


def _finish(result, dimension, derivatives):
    if not derivatives:
        if not np.isfinite(result):
            raise FloatingPointError('nonfinite reference energy')
        return float(result)
    if not isinstance(result, Jet):
        result = Jet.constant(result, dimension)
    if not (np.isfinite(result.value) and np.all(np.isfinite(result.gradient))
            and np.all(np.isfinite(result.hessian))):
        raise FloatingPointError('nonfinite reference jet')
    return result


def environment_jet(vectors, parameters, representation='moments'):
    """Cartesian derivatives of one centered site's energy, no bulk subtraction."""
    vectors = np.asarray(vectors, float)
    if vectors.ndim != 2 or vectors.shape[1] != 3:
        raise ValueError('neighbor vectors must have shape (n,3)')
    flat = _state(vectors.ravel(), vectors.size, True)
    result = _site_energy([flat[i:i+3] for i in range(0, len(flat), 3)],
                          _parameters(parameters), representation)
    return _finish(result, vectors.size, True)


class DiamondCell:
    """Periodic two-atom primitive cell including relative basis displacement.

    State: exx, eyy, ezz, gamma_yz, gamma_xz, gamma_xy, ux, uy, uz.
    F=I+engineering symmetric strain. u is an additional LAB-frame displacement
    of basis atom 1 relative to the affinely transformed basis. Energy per cell.
    """
    def __init__(self, parameters, lattice, image_shell=2):
        self.p = _parameters(parameters)
        self.lattice = float(lattice)
        if not np.isfinite(lattice) or lattice <= 0:
            raise ValueError('positive finite lattice parameter required')
        if int(image_shell) != image_shell or image_shell < 1:
            raise ValueError('positive integer image shell required')
        self.shell = int(image_shell)
        self.cell = lattice/2*np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0.]])
        self.volume = abs(float(np.linalg.det(self.cell)))
        basis = np.array([[0., 0., 0.], [lattice/4]*3])
        self.environments = []
        for i in range(2):
            sites = []
            for t in product(range(-self.shell, self.shell+1), repeat=3):
                for j in range(2):
                    if i == j and t == (0, 0, 0):
                        continue
                    sites.append((np.asarray(t)@self.cell+basis[j]-basis[i], j-i))
            self.environments.append(sites)

    def evaluate(self, state, representation='moments', derivatives=True):
        q = _state(state, 9, derivatives)
        numeric = np.asarray(state, float)
        f_numeric = np.eye(3)+np.array([[numeric[0], numeric[5]/2, numeric[4]/2],
            [numeric[5]/2, numeric[1], numeric[3]/2],
            [numeric[4]/2, numeric[3]/2, numeric[2]]])
        if np.linalg.det(f_numeric) <= 0:
            raise ValueError('orientation-preserving cell required')
        # Any omitted integer translate has |t| >= shell+1. Account for the
        # basis offset and internal shift in a conservative distance bound.
        basis_delta = self.lattice/4*np.ones(3)
        omitted_lower = (np.linalg.svd(self.cell@f_numeric.T, compute_uv=False)[-1]
                         * (self.shell+1) - np.linalg.norm(basis_delta@f_numeric.T)
                         - np.linalg.norm(numeric[6:]))
        if omitted_lower <= self.p['a']*self.p['sigma']:
            raise ValueError('image shell cannot certify the finite SW cutoff; increase shell')
        f = [[1+q[0], q[5]/2, q[4]/2], [q[5]/2, 1+q[1], q[3]/2],
             [q[4]/2, q[3]/2, 1+q[2]]]
        result = 0.
        cutoff = self.p['a']*self.p['sigma']
        for sites in self.environments:
            vectors = []
            for base, sign in sites:
                current = base@f_numeric.T+sign*numeric[6:]
                if np.linalg.norm(current) >= cutoff:
                    continue
                vectors.append([sum(f[a][b]*base[b] for b in range(3))+sign*q[6+a]
                                for a in range(3)])
            result += _site_energy(vectors, self.p, representation)
        return _finish(result, 9, derivatives)


class RigidInterface:
    """Three-coordinate (opening,u1,u2) audit of an existing finite SW slab.

    Takes the prior independent slab geometry, but constructs neighbor vectors
    and evaluates energy/derivatives here. Both half-crystals' per-site changes
    are summed; a zero-bulk angular background is never presumed.
    """
    def __init__(self, slab, image_shell=3):
        self.p = _parameters(slab.p)
        self.positions = np.asarray(slab.positions).copy()
        self.cell = np.asarray(slab.cell).copy()
        self.upper = np.asarray(slab.upper, int).copy()
        self.frame = np.array([slab.normal, slab.slip, np.cross(slab.normal, slab.slip)])
        self.cells, self.area, self.period = slab.cells, slab.area, slab.period
        if int(image_shell) != image_shell or image_shell < 1:
            raise ValueError('positive integer image shell required')
        self.shell = int(image_shell)
        translations = np.array(list(product(range(-self.shell, self.shell+1), repeat=2)))
        images = translations@self.cell[:2]
        replicas = (images[:, None, :]+self.positions).reshape(-1, 3)
        upper = np.tile(self.upper, len(images))
        self.environments = []
        for i, center in enumerate(self.positions):
            vectors = replicas-center
            keep = np.linalg.norm(vectors, axis=1) > 1e-9
            self.environments.append((vectors[keep], upper[keep]-self.upper[i]))

    def evaluate(self, state, representation='moments', derivatives=True):
        q = _state(state, 3, derivatives)
        displacement = np.asarray(state, float)@self.frame
        # Bound omitted translations independently of the retained neighbor list.
        planar = self.positions@self.frame[1:].T
        diameter = np.max(np.linalg.norm(planar[:, None, :]-planar, axis=2))
        omitted_lower = (np.linalg.svd(self.cell[:2], compute_uv=False)[-1]*(self.shell+1)
                         - diameter - np.linalg.norm(displacement@self.frame[1:].T))
        if omitted_lower <= self.p['a']*self.p['sigma']:
            raise ValueError('image shell cannot certify the finite SW cutoff; increase shell')
        shift = [sum(q[b]*self.frame[b, a] for b in range(3)) for a in range(3)]
        cutoff = self.p['a']*self.p['sigma']
        result = 0.
        for bases, signs in self.environments:
            active = ((np.linalg.norm(bases, axis=1) < cutoff)
                      | (np.linalg.norm(bases+signs[:, None]*displacement, axis=1) < cutoff))
            bases, signs = bases[active], signs[active]
            if not np.any(signs):
                continue
            vectors = [[base[a]+sign*shift[a] for a in range(3)]
                       for base, sign in zip(bases, signs)]
            result += (_site_energy(vectors, self.p, representation)
                       - _site_energy(bases, self.p, representation))
        return _finish(result/self.cells, 3, derivatives)


def relaxed_hessian(jet, retained):
    """Static Schur complement, only at an internally stable stationary state.

    Caller must verify internal stationarity. This is neither finite-T PMF nor
    mobility, and does not certify dynamical elimination of internal variables.
    """
    retained = np.asarray(retained, int)
    if (retained.ndim != 1 or not len(retained) or len(set(retained)) != len(retained)
            or np.any(retained < 0) or np.any(retained >= len(jet.gradient))):
        raise ValueError('distinct valid retained coordinate indices required')
    eliminated = np.array([i for i in range(len(jet.gradient)) if i not in retained])
    if not len(eliminated):
        return jet.hessian[np.ix_(retained, retained)].copy(), np.empty((0, len(retained)))
    hu = jet.hessian[np.ix_(eliminated, eliminated)]
    if np.linalg.eigvalsh(hu)[0] <= 0:
        raise ValueError('eliminated directions must have positive curvature')
    response = -np.linalg.solve(hu, jet.hessian[np.ix_(eliminated, retained)])
    effective = (jet.hessian[np.ix_(retained, retained)]
                 + jet.hessian[np.ix_(retained, eliminated)]@response)
    return effective, response
