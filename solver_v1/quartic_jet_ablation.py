"""Explicitly separate rational-quartic static ablation, NOT an adopted model.

The rational function was already derived/tested in v14. Here it is tested
explicitly with the later full environmental family, after the v23 fixed-
shape incompatibility audit. It is not silently combined into old defaults.
The existing v22 core bridge refuses this nonzero ablation until supported.
"""
from functools import lru_cache

import numpy as np

from .coordination_screening import CoordinationScreenedInterface
from .quartic_angular_material import QuarticMoment, per_site_quartic_jet


class RationalQuarticStaticInterface(CoordinationScreenedInterface):
    def __init__(self, shape, coefficients, *, quartic_saturation, tolerance=2e-11):
        if not np.isfinite(quartic_saturation) or quartic_saturation < 0:
            raise ValueError('explicit finite nonnegative rational-quartic shape required')
        super().__init__(shape, coefficients, tolerance=tolerance, law='power')
        self.quartic = QuarticMoment(self.base, saturation=float(quartic_saturation))
        self.static_ablation_only = True


class RationalQuarticJetColumn:
    def __init__(self, model, observations):
        self.moment = QuarticMoment(model.base)
        self.observations = list(observations)

    @lru_cache(maxsize=256)
    def depths(self, state):
        return self.moment.site_jets(state)[0]

    def reference_invariant(self, state):
        q = self.depths(tuple(state))[:, 0]
        return float(np.max(np.sum(q*q, axis=tuple(range(1, q.ndim)))))

    def column(self, saturation):
        if not np.isfinite(saturation) or saturation < 0:
            raise ValueError('nonnegative finite rational shape required')
        jets = {}
        out = np.zeros(len(self.observations))
        for i, o in enumerate(self.observations):
            if o.state is None:
                continue  # Q3=0 at every centrosymmetric affine bulk state
            if o.state not in jets:
                jets[o.state] = 2*per_site_quartic_jet(self.depths(o.state), saturation=saturation)
            out[i] = np.asarray(o.jet_weights)@jets[o.state]
        return out
