"""One-shape nonlinear even-environment ablation, NOT a production model.

Keep every LJ/Poisson/Bessel neighbor sum and the per-site counting. Replace
only the two quadratic rank2 site invariants by
    f(I; alpha, N2) = I/(1 + alpha*I/N2), I=Q:Q.
N2 is the fixed reference engineering-shear Hessian of Q:Q; the combined Eg
channel already has N2=1. This fixes an amplitude gauge, not a new length.
alpha>=0 is ONE dimensionless material-shape hypothesis, not a clipping rule.
f=I+O(I^2) preserves ALL harmonic response at Q=0, including finite-q and
cubic dilations. It does not preserve nonlinear response or imply stability.
"""
from functools import lru_cache

import numpy as np

from .even_moment_calibration import quadrupole_bulk_curvatures
from .quartic_angular_material import (
    QuarticSymmetryInterface, QuarticSymmetryObservationCache,
    QuarticSymmetryTailBulkBasis,
)
from .vector_interface_reference import HESSIAN_INDICES, VectorInterfaceEvaluation


def saturated_site_norm_jet(values, *, alpha, reference_curvature):
    """Analytic sum of per-atom f(I); never f(sum I). No force clipping."""
    z = np.asarray(values, float)
    if (z.ndim != 3 or z.shape[1] != 10 or z.shape[2] == 0
            or np.any(~np.isfinite(z)) or not np.isfinite(alpha) or alpha < 0
            or not np.isfinite(reference_curvature) or reference_curvature <= 0):
        raise ValueError('finite site jets, alpha>=0 and positive fixed gauge required')
    Q, dQ, ddQ = z[:, 0], z[:, 1:4], z[:, HESSIAN_INDICES]
    I = np.einsum('nc,nc->n', Q, Q)
    dI = 2*np.einsum('nc,nic->ni', Q, dQ)
    ddI = 2*(np.einsum('nic,njc->nij', dQ, dQ)
              + np.einsum('nc,nijc->nij', Q, ddQ))
    k = alpha/reference_curvature
    denominator = 1 + k*I
    first, second = denominator**-2, -2*k*denominator**-3
    grad = np.einsum('n,ni->i', first, dI)
    H = (np.einsum('n,ni,nj->ij', second, dI, dI)
         + np.einsum('n,nij->ij', first, ddI))
    return np.r_[np.sum(I/denominator), grad, H[0], H[1, 1:], H[2, 2]]


class CubicInterfaceMomentSites:
    """Full per-site moment change about a ZERO-moment cubic reference.

    This helper is not valid for a noncubic bulk with nonzero Q_bulk. All
    layers are summed before taking a norm. The symmetric upper/lower site
    factor is applied by the caller, not within the tensor sum.
    """
    def __init__(self, base, moment):
        self.base, self.moment = base, moment
        self.baselines = {}

    @lru_cache(maxsize=192)
    def site_jets(self, state):
        a, x, y = state
        obj = self.moment.invariant
        terms, small, maximum_shell = [], 0, 0
        for layer in range(1, obj.max_layers+1):
            delta = self.base.geometry.abc_shift(layer)
            if layer not in self.baselines:
                self.baselines[layer] = self.moment.evaluate(
                    layer*self.base.h, delta)[0][0].copy()
            jet, shells, envelope = self.moment.evaluate(
                a+(layer-1)*self.base.h, delta+np.array([x, y]))
            jet[0] -= self.baselines[layer]
            terms.append(jet)
            maximum_shell = max(maximum_shell, shells)
            last = float(np.max(abs(jet)))
            small = small+1 if last < obj.tolerance else 0
            if small >= 4:
                sites = np.cumsum(np.array(terms)[::-1], axis=0)[::-1]
                return sites, dict(layers=layer, shells=maximum_shell,
                                   last_layer_jet=last, last_shell_envelope=envelope,
                                   empirical_truncation=True)
        raise ArithmeticError('even-moment interface neighborhood did not converge')


class SaturatedQuadrupoleInterface(QuarticSymmetryInterface):
    """Nested analytic static reference; no production registration."""
    def __init__(self, shape, coefficients, *, tolerance=2e-11):
        shape = np.asarray(shape, float)
        if (shape.shape != (4,) or np.any(~np.isfinite(shape))
                or np.min(shape[:3]) <= 0 or shape[3] < 0):
            raise ValueError('three positive decays and alpha>=0 required')
        super().__init__(shape[:3], coefficients, tolerance=tolerance)
        self.alpha = float(shape[3])
        self.even_reference_curvature = float(
            quadrupole_bulk_curvatures(self.base.surface.quadrupole)[2])
        if self.even_reference_curvature <= 0:
            raise ArithmeticError('unresolved even-moment reference gauge')
        even = next(m for _, m in self.base.moments if m.invariant.rank == 2)
        self.even_sites = CubicInterfaceMomentSites(self.base, even)
        self.eg_sites = CubicInterfaceMomentSites(self.base, self.moment)

    def even_jets(self, state):
        state = tuple(map(float, state))
        even, ed = self.even_sites.site_jets(state)
        eg, gd = self.eg_sites.site_jets(state)
        return (2*saturated_site_norm_jet(even, alpha=self.alpha,
                                         reference_curvature=self.even_reference_curvature),
                2*saturated_site_norm_jet(eg, alpha=self.alpha, reference_curvature=1.),
                dict(even=ed, Eg=gd, alpha=self.alpha,
                     even_reference_curvature=self.even_reference_curvature))

    def evaluate(self, state):
        result = super().evaluate(state)
        if self.alpha == 0:
            return result  # exact unchanged-model compatibility
        even, eg, diag = self.even_jets(state)
        components, diagnostics = dict(result.components), dict(result.diagnostics)
        components['angular_2'] = self.coefficients[7]*even
        components['angular_Eg_radial'] = self.coefficients[8]*eg
        diagnostics['quadrupole_saturation'] = diag
        return VectorInterfaceEvaluation.from_jet(sum(components.values()), components, diagnostics)

    def direct_even_energy(self, state, *, radius, layers):
        """Independent finite real-space validation, NEVER the energy generator."""
        from .angular_environment_reference import traceless_second
        if (int(radius) != radius or radius < 1 or int(layers) != layers or layers < 1):
            raise ValueError('positive integer validation radius and layer count required')
        a, x, y = map(float, state)
        m, n = np.meshgrid(np.arange(-radius, radius+1), np.arange(-radius, radius+1), indexing='ij')
        R = m.ravel()[:, None]*self.geometry.a1 + n.ravel()[:, None]*self.geometry.a2
        invariants = [p.invariant for p in self.moment.parts]
        changes = [[], []]
        for k in range(1, layers+1):
            delta = self.geometry.abc_shift(k)
            for j, obj in enumerate(invariants):
                values = []
                for d, shift in ((a+(k-1)*self.h, delta+[x, y]), (k*self.h, delta)):
                    pos = np.column_stack([R+shift, np.full(len(R), d)])
                    w = obj.amplitude*np.exp(-obj.kappa*np.linalg.norm(pos, axis=1))
                    values.append(traceless_second(np.einsum('n,ni,nj->ij', w, pos, pos)))
                changes[j].append(values[0]-values[1])
        odd, even = [np.cumsum(np.array(c)[::-1], axis=0)[::-1] for c in changes]
        eg = (odd-self.moment.eta*even)/self.moment.normalization
        I2, IE = [np.sum(Q*Q, axis=(1, 2)) for Q in (even, eg)]
        return np.array([2*np.sum(I2/(1+self.alpha*I2/self.even_reference_curvature)),
                         2*np.sum(IE/(1+self.alpha*IE))])


class SaturatedQuadrupoleCache(QuarticSymmetryObservationCache):
    @lru_cache(maxsize=128)
    def nonlinear_even(self, odd, even, alpha):
        model = SaturatedQuadrupoleInterface((3., odd, even, alpha), np.ones(10))
        jets = {}
        for o in self.observations:
            if o.bulk_index is None and o.state not in jets:
                jets[o.state] = model.even_jets(o.state)[:2]
        # Cubic harmonic columns are copied from the legacy matrix unchanged.
        return np.array([[0., 0.] if o.bulk_index is not None else
                         [np.asarray(o.jet_weights)@v for v in jets[o.state]]
                         for o in self.observations])

    def matrix(self, shape):
        shape = np.asarray(shape, float)
        if (shape.shape != (4,) or np.any(~np.isfinite(shape)) or shape[3] < 0):
            raise ValueError('four finite shape parameters and alpha>=0 required')
        matrix = super().matrix(shape[:3]).copy()
        if shape[3] == 0:
            return matrix
        changed = self.nonlinear_even(*map(float, shape[1:]))
        rows = np.array([o.bulk_index is None for o in self.observations])
        matrix[rows, 7:9] = changed[rows]
        return matrix


# Q2=QE=0 in every cubic dilation: f'(0)=1 preserves the existing complete
# harmonic finite-q matrix, its coefficient linearity, and its tail bound.
SaturatedQuadrupoleBulk = QuarticSymmetryTailBulkBasis
