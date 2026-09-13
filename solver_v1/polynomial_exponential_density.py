"""Positive quadratic-envelope exponential: isolated RESEARCH primitive.

f(r)=C exp(-k r)[(r-c)^2+w^2], C,k,w>0. This does not replace any
production density. Unlike a positive sum of pure exponentials it can have
negative logarithmic curvature. Source EAM density is auxiliary/gauge-dependent,
so matching this primitive alone is NOT material calibration.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=12)
def transform_terms(k_order, d_order):
    """Exact derivatives of k exp(-d Q)(Q^-3+d Q^-2), Q²=k²+G².

Term keys are powers of k, d, and inverse Q. No finite differences.
The common factor 2*pi*exp(-d Q) is applied after this polynomial recurrence.
"""
    terms = {(1, 0, 3): 1., (1, 1, 2): 1.}
    for _ in range(k_order):
        new = {}
        for (h, j, p), value in terms.items():
            for key, coefficient in (((h-1, j, p), h), ((h+1, j+1, p+1), -1),
                                     ((h+1, j, p+2), -p)):
                if coefficient:
                    new[key] = new.get(key, 0.)+coefficient*value
        terms = new
    for _ in range(d_order):
        new = {}
        for (h, j, p), value in terms.items():
            for key, coefficient in (((h, j-1, p), j), ((h, j, p-1), -1)):
                if coefficient:
                    new[key] = new.get(key, 0.)+coefficient*value
        terms = new
    return tuple((h, j, p, v) for (h, j, p), v in terms.items() if v)


@dataclass(frozen=True)
class QuadraticEnvelopeDensity:
    """All lengths use one explicit caller-selected coordinate unit.

If r is dimensional, C has density/length², k has 1/length, c,w have length.
For reduced r, these become dimensionless shape/gauge parameters.
"""
    C: float
    k: float
    center: float
    width: float

    def __post_init__(self):
        values = np.array([self.C, self.k, self.center, self.width])
        if np.any(~np.isfinite(values)) or min(self.C, self.k, self.width) <= 0:
            raise ValueError('finite parameters with positive amplitude/decay/width required')

    def radial(self, r, order=0):
        r = np.asarray(r, float)
        if np.any(~np.isfinite(r)) or np.any(r < 0) or order not in (0, 1, 2):
            raise ValueError('nonnegative finite radius and derivative order 0,1,2 required')
        y = r-self.center
        polynomial = y*y+self.width*self.width
        value = (polynomial if order == 0 else 2*y-self.k*polynomial if order == 1
                 else 2-4*self.k*y+self.k*self.k*polynomial)
        return self.C*np.exp(-self.k*r)*value

    def plane_transform(self, d, G, d_order=0):
        """2D Fourier transform and exact first/second normal derivatives.

T_quad=C(∂k²+2c∂k+c²+w²) T_exp; normal derivatives commute.
This is an exact analytic identity, evaluated in floating point.
"""
        G = np.asarray(G, float)
        if not np.isfinite(d) or d < 0 or np.any(~np.isfinite(G)) or np.any(G < 0) or d_order not in (0, 1, 2):
            raise ValueError('finite nonnegative separation/wavenumber and order 0,1,2 required')
        Q = np.sqrt(self.k*self.k+G*G)
        terms = []
        for order, coefficient in enumerate((self.center**2+self.width**2, 2*self.center, 1.)):
            terms.extend(coefficient*v*self.k**h*d**j*Q**(-p)
                         for h, j, p, v in transform_terms(order, d_order))
        result = self.C*2*np.pi*np.exp(-d*Q)*sum(terms)
        if np.any(~np.isfinite(result)):
            raise ArithmeticError('nonfinite Fourier transform')
        return result


def positive_exponential_log_curvature(r, amplitudes, decays):
    """(log sum A exp(-k r))'' = Var_weighted(k) >=0, a shape identity."""
    a, k = np.asarray(amplitudes, float), np.asarray(decays, float)
    if (a.ndim != 1 or not len(a) or k.shape != a.shape or np.any(a <= 0)
            or np.any(k <= 0) or np.any(~np.isfinite(a+k)) or not np.isfinite(r)):
        raise ValueError('positive matching exponential coefficients/decays required')
    log_weight = np.log(a)-k*r
    weight = np.exp(log_weight-log_weight.max()); weight /= weight.sum()
    mean = weight@k
    return float(weight@((k-mean)**2))


def plane_sum(kernel, d, delta, geometry, *, tolerance=1e-12, max_shell_index=32):
    """Infinite 2D Poisson identity, semi-analytic reciprocal truncation.

The last-shell envelope is empirical, NOT a rigorous full-stack tail bound.
No full-crystal embedding, material fit or production registration is implied.
"""
    from .fcc111_lattice_sum import triangular_reciprocal_shells, PlaneKernelResult
    delta = np.asarray(delta, float)
    if (not np.isfinite(d) or d <= 0 or delta.shape != (2,) or np.any(~np.isfinite(delta))
            or not np.isfinite(tolerance) or tolerance <= 0 or max_shell_index < 2):
        raise ValueError('positive separation/tolerance and finite registry required')
    jet = np.zeros(10)
    jet[[0, 1, 4]] = [kernel.plane_transform(d, 0., n) for n in (0, 1, 2)]
    jet /= geometry.atomic_cell_area
    small, count, vectors_used = 0, 0, 1
    for shell in triangular_reciprocal_shells(geometry, max_shell_index):
        vectors = shell.vectors
        gx, gy = vectors.T
        value, normal, second = [kernel.plane_transform(d, shell.magnitude, n)/geometry.atomic_cell_area
                                  for n in (0, 1, 2)]
        cs = np.array([np.full(len(gx), value), np.full(len(gx), normal),
                       1j*gx*value, 1j*gy*value, np.full(len(gx), second),
                       1j*gx*normal, 1j*gy*normal, -gx*gx*value, -gx*gy*value, -gy*gy*value])
        jet += np.real(cs@np.exp(1j*(vectors@delta)))
        last = float(np.max(np.sum(abs(cs), axis=1)))
        count += 1; vectors_used += len(vectors)
        small = small+1 if last < .1*tolerance else 0
        if small >= 3:
            return PlaneKernelResult(jet[0], jet[1], jet[2:4], jet[4], jet[5:7],
                np.array([[jet[7], jet[8]], [jet[8], jet[9]]]), count, vectors_used, last)
    raise ArithmeticError('quadratic exponential reciprocal sum did not converge')
