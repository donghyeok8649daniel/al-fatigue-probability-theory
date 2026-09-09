"""One-parameter angular-range extension, separate from production energies.

Ranks1/3 share k_odd; rank2 may use k_even. Scalar density, LJ functional form,
per-atom embedding and STF normalization are unchanged. Equal ranges recover
the existing eight-coefficient model exactly. This is an analytic research
family, not a calibrated MEAM potential or a strength law.
"""
from functools import lru_cache

import numpy as np

from .angular_environment_reference import AngularInterfaceInvariant
from .even_moment_calibration import even_basis, quadrupole_bulk_curvatures
from .vector_interface_reference import VectorMomentPlane
from .vector_material_calibration import (
    VectorCoefficientBasis, build_coefficient_surface, exact_constraint_tangent,
    coefficient_identifiability,
)


def build_range_surface(scalar_decay, odd_decay, quadrupole_decay, coefficients,
                        *, tolerance=2e-11):
    """Construct a fresh full-vector static surface; no caller mutation."""
    decay = np.asarray([scalar_decay, odd_decay, quadrupole_decay], float)
    if np.any(~np.isfinite(decay)) or np.min(decay) <= 0:
        raise ValueError('three positive finite microscopic decays required')
    model = build_coefficient_surface(scalar_decay, odd_decay, coefficients,
                                      tolerance=tolerance)
    even = AngularInterfaceInvariant(model.face, tolerance=tolerance,
                                    angular_decay=quadrupole_decay, rank=2)
    model.surface.quadrupole = even
    model.moments = [(amplitude, moment) for amplitude, moment in model.moments
                     if moment.invariant.rank != 2]
    model.moments.append((float(coefficients[7]), VectorMomentPlane(even)))
    return model


class RangeResolvedBasis(VectorCoefficientBasis):
    """Exact coefficient jet, with one additional radial shape coordinate."""
    def __init__(self, scalar_decay, odd_decay, quadrupole_decay,
                 *, tolerance=2e-11):
        self.scalar_decay = float(scalar_decay)
        self.angular_decay = float(odd_decay)  # legacy bulk helper convention
        self.quadrupole_decay = float(quadrupole_decay)
        self.model = build_range_surface(scalar_decay, odd_decay, quadrupole_decay,
                                        np.ones(8), tolerance=tolerance)
        self.h = self.model.h


def range_observation_matrix(basis, observations):
    """Same observations and gauge, with only the rank2 radial column changed."""
    bulk = even_basis(basis.scalar_decay, basis.angular_decay, convex=True)[:5].copy()
    bulk[:, 7] = np.r_[0., 0., quadrupole_bulk_curvatures(basis.model.surface.quadrupole)]
    return np.asarray([
        bulk[o.bulk_index] if o.bulk_index is not None else
        np.asarray(o.jet_weights)@basis.jet(o.state) for o in observations])


class RangeObservationCache:
    """Cache exact channel columns, not a tabulated state-energy surrogate.

    Each observed energy/jet is still an infinite-series evaluation. Only rank2
    depends on k_even; separating this exact dependency accelerates profiling.
    Observations are bound to this object and never taken from a global cache.
    """
    def __init__(self, observations):
        self.observations = tuple(observations)

    @lru_cache(maxsize=192)
    def legacy(self, scalar_decay, odd_decay):
        from .vector_material_calibration import observation_matrix
        return observation_matrix(VectorCoefficientBasis(scalar_decay, odd_decay),
                                  self.observations)

    @lru_cache(maxsize=192)
    def even(self, quadrupole_decay):
        # Scalar/odd values here cannot affect this separately normalized rank2
        # column. The equal-range reference is a convenient exact evaluator.
        basis = RangeResolvedBasis(2., 4., quadrupole_decay)
        bulk = np.r_[0., 0., quadrupole_bulk_curvatures(basis.model.surface.quadrupole)]
        moment = basis.model.moments[-1][1]
        return np.asarray([
            bulk[o.bulk_index] if o.bulk_index is not None else
            np.asarray(o.jet_weights)@basis.model._angular(o.state, moment)[0]
            for o in self.observations])

    def matrix(self, decays):
        scalar, odd, even = map(float, decays)
        matrix = self.legacy(scalar, odd).copy()
        matrix[:, 7] = self.even(even)
        return matrix


class RangeBulkValidation:
    """Independent direct finite-q check; do not double-count scalar EAM.

    The two evaluators share the same pair/density energy. Retain that energy
    ONCE and add only the separate-ranged quadrupole angular contribution.
    Radius refinement remains required; this is not a cutoff production model.
    """
    def __init__(self,model,*,cutoff=8.):
        from .static_bulk_stability import StaticBulkHessian
        surface=model.surface; bulk=model.face.bulk
        self.odd=StaticBulkHessian(bulk,cutoff=cutoff,D3=surface.amplitude_ev,
            D1=surface.vector_amplitude_ev,angular_decay=surface.angular.kappa)
        self.even=StaticBulkHessian(bulk,cutoff=cutoff,D2=surface.quadrupole_amplitude_ev,
            angular_decay=surface.quadrupole.kappa)

    def evaluate(self,q):
        odd=self.odd.evaluate(q); even=self.even.evaluate(q)
        matrix=odd['matrix']+even['angular']
        return dict(matrix=matrix,eigenvalues=np.linalg.eigvalsh(matrix))

    def crystallographic_wavevector(self,fractional):
        return self.odd.crystallographic_wavevector(fractional)


def radial_identifiability(decays, coefficients, observations):
    """Nine local directions after two exact constraints; not confidence.

    Compare the tied-log-angular-range direction (8 columns) with the added
    contrast. It is the nested common-range tangent ONLY on k_even=k_odd.
    Finite log steps are refined independently; no yield target.
    """
    from .yield_elastic_metric import cubic_metric_problem
    cache = RangeObservationCache(observations)
    base, obs = cubic_metric_problem(cache.matrix(decays), observations)
    info = coefficient_identifiability(base, obs, coefficients)
    selected = [i for i, o in enumerate(obs) if o.role == 'fit']
    scales = np.array([o.scale for o in obs])
    target = np.array([o.target for o in obs])
    z = np.asarray(coefficients)[2:]
    refinements = []
    for step in (4e-4, 2e-4):
        columns = []
        for axis in range(3):
            predictions = []
            for sign in (1, -1):
                shifted = np.array(decays, float)
                shifted[axis] *= np.exp(sign*step)
                matrix, _ = cubic_metric_problem(cache.matrix(shifted), observations)
                offset, tangent = exact_constraint_tangent(matrix, target)
                predictions.append(matrix@(offset+tangent@z))
            columns.append((predictions[0]-predictions[1])[selected]/(2*step*scales[selected]))
        refinements.append(np.asarray(columns).T)
    full = np.column_stack([info['jacobian'], refinements[-1]])
    norms = np.linalg.norm(full, axis=0)
    _, singular, vt = np.linalg.svd(full, full_matrices=False)
    floor = np.finfo(float).eps*max(full.shape)*singular[0]
    # At a common-range point, dkodd=dkeven is the original radial direction.
    nested = np.column_stack([full[:, :7], full[:, 7]+full[:, 8]])
    nested_singular = np.linalg.svd(nested, compute_uv=False)
    info.update(full_jacobian=full, full_singular_values=singular,
                full_rank=int(np.sum(singular > floor)),
                full_condition_number=float(singular[0]/singular[-1]) if singular[-1] > floor else None,
                full_right_vectors=vt, column_cosines=full.T@full/(norms[:, None]*norms[None, :]),
                radial_difference_change=float(np.max(abs(refinements[0]-refinements[1]))),
                tied_log_angular_range_singular_values=nested_singular,
                on_common_range_submanifold=bool(np.isclose(decays[1],decays[2],rtol=1e-13,atol=0)),
                full_radial_identifiability_checked=True,
                confidence_intervals_available=False,
                notes='Local scaled sensitivities; active coefficient bounds and material discrepancy preclude confidence claims.')
    return info
