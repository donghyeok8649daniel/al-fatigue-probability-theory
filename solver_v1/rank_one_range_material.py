"""One independent rank1 radial range, separate analytic static research.

The normalized exponential/LJ/Poisson sums and per-atom energy are unchanged.
k1=k_odd recovers v16. Homogeneous bulk columns are unchanged by inversion;
nonuniform Bloch/interface response IS changed and must be revalidated.
"""
from functools import lru_cache

import numpy as np

from .angular_environment_reference import AngularInterfaceInvariant
from .isotropic_bulk_validation import IsotropicBulkBasis
from .quadrupole_saturation import SaturatedQuadrupoleInterface, SaturatedQuadrupoleCache
from .range_resolved_material import build_range_surface
from .vector_interface_reference import VectorMomentPlane


class RankOneRangeInterface(SaturatedQuadrupoleInterface):
    def __init__(self, shape, coefficients, *, tolerance=2e-11):
        shape = np.asarray(shape, float)
        if (shape.shape != (5,) or np.any(~np.isfinite(shape)) or shape[4] <= 0):
            raise ValueError('v16 shape plus finite positive rank1 decay required')
        super().__init__(shape[:4], coefficients, tolerance=tolerance)
        self.full_shape = shape.copy()
        self.rank1_decay = float(shape[4])
        if self.rank1_decay != shape[1]:
            moment = AngularInterfaceInvariant(self.face, rank=1,
                angular_decay=self.rank1_decay, tolerance=tolerance)
            self.base.surface.vector = moment
            self.base.moments = [(amplitude, VectorMomentPlane(moment) if old.invariant.rank == 1 else old)
                                 for amplitude, old in self.base.moments]
            # This is a fresh private model; do not reuse a tied-range baseline.
            if any(key[0] == 'angular' and key[1] == 1 for key in self.base._baselines):
                raise AssertionError('new rank1 model unexpectedly has stale site baselines')

    def direct_rank1_energy(self, state, *, radius, layers):
        """Independent real-space per-site energy, not the canonical evaluator."""
        if int(radius) != radius or radius < 1 or int(layers) != layers or layers < 1:
            raise ValueError('positive integer independent radius and layer count required')
        a, x, y = map(float, state)
        m, n = np.meshgrid(np.arange(-radius, radius+1), np.arange(-radius, radius+1), indexing='ij')
        R = m.ravel()[:, None]*self.geometry.a1+n.ravel()[:, None]*self.geometry.a2
        moment = self.base.surface.vector
        changes = []
        for k in range(1, layers+1):
            delta = self.geometry.abc_shift(k)
            values = []
            for d, shift in ((a+(k-1)*self.h, delta+[x, y]), (k*self.h, delta)):
                pos = np.column_stack([R+shift, np.full(len(R), d)])
                weight = moment.amplitude*np.exp(-moment.kappa*np.linalg.norm(pos, axis=1))
                values.append(weight@pos)
            changes.append(values[0]-values[1])
        sites = np.cumsum(np.array(changes)[::-1], axis=0)[::-1]
        return float(2*np.sum(sites*sites))


class RankOneRangeCache(SaturatedQuadrupoleCache):
    @lru_cache(maxsize=128)
    def vector_column(self, decay):
        model = build_range_surface(3., decay, 5., np.ones(8))
        moment = next(m for _, m in model.moments if m.invariant.rank == 1)
        return np.array([0. if o.bulk_index is not None else
            np.asarray(o.jet_weights)@model._angular(o.state, moment)[0]
            for o in self.observations])

    def matrix(self, shape):
        shape = np.asarray(shape, float)
        if (shape.shape != (5,) or np.any(~np.isfinite(shape)) or shape[4] <= 0):
            raise ValueError('v16 shape plus positive independent rank1 decay required')
        matrix = super().matrix(shape[:4]).copy()
        if shape[4] != shape[1]:
            matrix[:, 6] = self.vector_column(float(shape[4]))
        return matrix


class RankOneRangeBulk(IsotropicBulkBasis):
    """Same analytical coefficient/tail operator, now with the correct k1."""
    def __init__(self, shape, *, radius=12., stretch=1.):
        shape = np.asarray(shape, float)
        if (shape.shape != (5,) or np.any(~np.isfinite(shape))
                or np.any(shape[[0, 1, 2, 4]] <= 0) or shape[3] < 0):
            raise ValueError('full five-component material shape required; do not discard k1')
        super().__init__(shape[:3], radius=radius, stretch=stretch, rank1_decay=shape[4])


def shape_for_bulk(shape, definition):
    """Do not silently discard an independent radial shape in old helpers."""
    return shape if (definition.get('density_mixture_extension')
                     or definition.get('rank_one_range_extension')) else shape[:3]
