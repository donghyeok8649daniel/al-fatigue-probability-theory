"""One-parameter density-normalized vector environment, static research only.

E1_i = D1 |Q1_i|^2 / (1-z+z*x_i), x_i=rho_i/rho_ref, 0<=z<=1.
All infinite per-atom sums precede this function. Old models are not modified.
See COORDINATION_SCREENING_V19.md for units, admissibility and limitations.
"""
from functools import lru_cache

import numpy as np

from .quadrupole_saturation import CubicInterfaceMomentSites
from .rank_one_range_material import RankOneRangeInterface, RankOneRangeCache, RankOneRangeBulk
from .vector_interface_reference import HESSIAN_INDICES, VectorInterfaceEvaluation


def screening_factor(x, screening, *, law='rational'):
    x = np.asarray(x, float)
    z = float(screening)
    if (law not in ('rational', 'power') or not np.isfinite(z)
            or not (-1 if law == 'power' else 0) <= z <= 1 or np.any(~np.isfinite(x))
            or np.any(x <= 0)):
        raise ValueError('positive site densities and valid explicit screening law/shape required')
    if law == 'power':
        g = np.exp(z*np.log(x))
        return g, z*g/x, z*(z-1)*g/x**2
    d = (1-z)+z*x
    return 1/d, -z/d**2, 2*z*z/d**3


def screened_site_jet(density_changes, moment_jets, screening, *, law='rational'):
    """Sum g(1+delta_x_i)||Q_i||² with exact a/x/y first/second jets."""
    density = np.asarray(density_changes, float)
    moments = np.asarray(moment_jets, float)
    if (density.ndim != 2 or density.shape[1] != 10 or moments.ndim != 3
            or moments.shape[1:] != (10, 3) or not len(density) or not len(moments)
            or np.any(~np.isfinite(density)) or np.any(~np.isfinite(moments))):
        raise ValueError('finite nonempty scalar and rank1 per-site jets required')
    count = max(len(density), len(moments))
    x = np.zeros((count, 10)); x[:len(density)] = density; x[:, 0] += 1
    t = np.zeros((count, 10, 3)); t[:len(moments)] = moments
    Q, Qp, Qpq = t[:, 0], t[:, 1:4], t[:, HESSIAN_INDICES]
    I = np.einsum('nc,nc->n', Q, Q)
    Ip = 2*np.einsum('nc,nic->ni', Q, Qp)
    Ipq = 2*(np.einsum('nic,njc->nij', Qp, Qp)+np.einsum('nc,nijc->nij', Q, Qpq))
    xp, xpq = x[:, 1:4], x[:, HESSIAN_INDICES]
    g, gx, gxx = screening_factor(x[:, 0], screening, law=law)
    gradient = np.einsum('n,ni->i', g, Ip)+np.einsum('n,n,ni->i', gx, I, xp)
    H = (np.einsum('n,nij->ij', g, Ipq)
         + np.einsum('n,ni,nj->ij', gx, xp, Ip)
         + np.einsum('n,ni,nj->ij', gx, Ip, xp)
         + np.einsum('n,n,nij->ij', gx, I, xpq)
         + np.einsum('n,n,ni,nj->ij', gxx, I, xp, xp))
    return np.r_[g@I, gradient, H[0], H[1, 1:], H[2, 2]]


class CoordinationScreenedInterface(RankOneRangeInterface):
    def __init__(self, shape, coefficients, *, tolerance=2e-11, law='rational'):
        shape = np.asarray(shape, float)
        if shape.shape != (6,):
            raise ValueError('five inherited shapes plus explicit screening required')
        screening_factor(1., shape[5], law=law)
        # At the inverse-density endpoint a vanishing density must not create
        # a growing isolated-neighbor energy. Refuse, do not clip parameters.
        if law == 'rational' and shape[5] == 1 and 2*shape[4] <= shape[0]:
            raise ValueError('inverse-density isolated-neighbor limit requires 2*k1>k_scalar')
        if law == 'power' and 2*shape[4]+min(shape[5], 0)*shape[0] <= 0:
            raise ValueError('power-law isolated-neighbor energy must decay')
        super().__init__(shape[:5], coefficients, tolerance=tolerance)
        self.law = law
        self.screening = float(shape[5])
        self.screened_shape = shape.copy()
        moment = next(m for _, m in self.base.moments if m.invariant.rank == 1)
        self.vector_sites = CubicInterfaceMomentSites(self.base, moment)

    @lru_cache(maxsize=192)
    def site_inputs(self, state):
        density, *diag = self.base._density(np.asarray(state, float))
        sites, moment_diag = self.vector_sites.site_jets(state)
        return (density/self.face.rho_bulk, sites,
                dict(density=diag, vector=moment_diag, empirical_truncation=True))

    def screened_jet(self, state, screening=None):
        density, sites, diag = self.site_inputs(tuple(map(float, state)))
        z = self.screening if screening is None else screening
        return 2*screened_site_jet(density, sites, z, law=self.law), diag

    def evaluate(self, state):
        result = super().evaluate(state)
        if self.screening == 0:
            return result
        jet, diag = self.screened_jet(state)
        components, diagnostics = dict(result.components), dict(result.diagnostics)
        components['angular_1'] = self.coefficients[6]*jet
        diagnostics['coordination_screening'] = dict(screening=self.screening, law=self.law, **diag)
        return VectorInterfaceEvaluation.from_jet(sum(components.values()), components, diagnostics)

    def direct_screened_energy(self, state, *, radius, layers):
        """Independent finite real-space validation only; per-site normalization."""
        if int(radius) != radius or radius < 1 or int(layers) != layers or layers < 1:
            raise ValueError('positive integer direct validation extents required')
        a, x, y = map(float, state)
        m, n = np.meshgrid(np.arange(-radius, radius+1), np.arange(-radius, radius+1), indexing='ij')
        R = m.ravel()[:, None]*self.geometry.a1+n.ravel()[:, None]*self.geometry.a2
        den = self.face.bulk.density_params
        vector = self.base.surface.vector
        scalar_changes, vector_changes = [], []
        for layer in range(1, layers+1):
            delta = self.geometry.abc_shift(layer)
            values = []
            for d, shift in ((a+(layer-1)*self.h, delta+[x, y]), (layer*self.h, delta)):
                pos = np.column_stack([R+shift, np.full(len(R), d)])
                r = np.linalg.norm(pos, axis=1)
                values.append((np.sum(den.C_rho*np.exp(-den.kappa*r))/self.face.rho_bulk,
                               vector.amplitude*np.exp(-vector.kappa*r)@pos))
            scalar_changes.append(values[0][0]-values[1][0])
            vector_changes.append(values[0][1]-values[1][1])
        density = 1+np.cumsum(scalar_changes[::-1])[::-1]
        Q = np.cumsum(np.asarray(vector_changes)[::-1], axis=0)[::-1]
        return float(2*np.sum(screening_factor(density, self.screening, law=self.law)[0]*np.sum(Q*Q, axis=1)))


class CoordinationScreenedCache(RankOneRangeCache):
    def __init__(self, observations, *, law='rational'):
        super().__init__(observations)
        screening_factor(1., 0., law=law)
        self.law = law

    @lru_cache(maxsize=24)
    def site_model(self, scalar, rank1):
        return CoordinationScreenedInterface((scalar, 6., 8., 0., rank1, 0.), np.ones(10), law=self.law)

    def matrix(self, shape):
        shape = np.asarray(shape, float)
        if shape.shape != (6,):
            raise ValueError('six explicit material shape parameters required')
        screening_factor(1., shape[5], law=self.law)
        if ((self.law == 'rational' and shape[5] == 1 and 2*shape[4] <= shape[0])
                or (self.law == 'power' and 2*shape[4]+min(shape[5], 0)*shape[0] <= 0)):
            raise ValueError('invalid isolated-neighbor endpoint limit')
        matrix = super().matrix(shape[:5]).copy()
        if shape[5] == 0:
            return matrix
        model = self.site_model(float(shape[0]), float(shape[4]))
        for i, o in enumerate(self.observations):
            if o.bulk_index is None:
                matrix[i, 6] = np.asarray(o.jet_weights)@model.screened_jet(o.state, shape[5])[0]
        return matrix


class CoordinationScreenedBulk(RankOneRangeBulk):
    def __init__(self, shape, *, radius=12., stretch=1., law='rational'):
        shape = np.asarray(shape, float)
        if shape.shape != (6,):
            raise ValueError('full screened material shape required')
        screening_factor(1., shape[5], law=law)
        if ((law == 'rational' and shape[5] == 1 and 2*shape[4] <= shape[0])
                or (law == 'power' and 2*shape[4]+min(shape[5], 0)*shape[0] <= 0)):
            raise ValueError('invalid isolated-neighbor endpoint limit')
        super().__init__(shape[:5], radius=radius, stretch=stretch)
        self.law = law
        self.screening = float(shape[5])

    def evaluate(self, q_cubic):
        columns, errors = super().evaluate(q_cubic)
        factor = float(screening_factor(self.x, self.screening, law=self.law)[0])
        lower = self.x-self.density_series_tail
        derivative_bound = abs(float(screening_factor(lower, self.screening, law=self.law)[1]))
        errors[6] = (factor*errors[6]+derivative_bound*self.density_series_tail*
                     (np.linalg.norm(columns[6], 2)+errors[6]))
        columns[6] *= factor
        return columns, errors

    def direct_sinusoidal_energy(self, coefficients, *, planes, mode, polarization_plane,
                                amplitude, validation_radius, saturation=0., quadrupole_saturation=0.):
        """Same site law in independent displaced-neighbor harmonic validation."""
        from types import SimpleNamespace
        from .static_bulk_stability import neighbors_in_plane_frame
        from .tail_constrained_material import normalized_amplitude
        original = super().direct_sinusoidal_energy(coefficients, planes=planes, mode=mode,
            polarization_plane=polarization_plane, amplitude=amplitude,
            validation_radius=validation_radius, saturation=saturation,
            quadrupole_saturation=quadrupole_saturation)
        if self.screening == 0:
            return original
        R = neighbors_in_plane_frame(SimpleNamespace(geometry=self.geometry, a0=self.geometry.h111),
                                     validation_radius)
        layer = np.rint(R[:, 2]/self.geometry.h111).astype(int)
        phase = 2*np.pi*mode/planes; correction = []
        for site in range(planes):
            shift = amplitude*(np.cos(phase*(site+layer))-np.cos(phase*site))
            positions = R+shift[:, None]*np.asarray(polarization_plane)
            r = np.linalg.norm(positions, axis=1)
            x = np.sum(normalized_amplitude(self.decays[0])*np.exp(-self.decays[0]*r))
            Q = normalized_amplitude(self.rank1_decay)*np.exp(-self.rank1_decay*r)@positions
            g = float(screening_factor(x, self.screening, law=self.law)[0])
            correction.append(coefficients[6]*(g-1)*(Q@Q))
        return float(original+np.mean(correction))
