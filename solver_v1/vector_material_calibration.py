"""Matched vector-interface/material calibration, separate from specimen strength.

No new energy family: LJ + scalar embedding + per-site STF ranks 1/2/3.
At fixed radial decays every energy/force/Hessian entry is linear in eight
coefficients. This is an exact coefficient decomposition, NOT a surrogate or
tabulation of the candidate energy. SOURCE EAM supplies observations only.
Experimental yield, source lengths, mobilities and fatigue are never fit here.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq

from .angular_environment_reference import AngularInterfaceResearchSurface
from .even_moment_calibration import even_basis
from .fcc111_active_interface import FCC111ActiveInterface
from .full_fcc_calibration_audit import IDEAL_H, NAMES, independent_bulk_targets
from .interface_static_scenarios import InterfaceUnits
from .joint_fcc_interface_calibration import JointFit, MinimalConvexEmbedding, build_joint_model
from .run_constrained_odd_calibration import constrained_matrix_fit
from .vector_interface_reference import FullRegistryInterface, _embedding_jet, VectorInterfaceEvaluation
from .vector_registry_audit import stationary_state


COEFFICIENTS = ('u', 'v', 'A', 'B', 'C', 'D3', 'D1', 'D2')
LENGTH_M = 4.05/math.sqrt(2)*1e-10
UNITS = InterfaceUnits(LENGTH_M, math.sqrt(3)/2*LENGTH_M**2)
SOURCE_CITATION = 'Mishin et al., Phys. Rev. B 59 (1999) 3393, doi:10.1103/PhysRevB.59.3393; NIST Al99'


def build_coefficient_surface(scalar_decay, angular_decay, coefficients, *, tolerance=2e-11):
    c = np.asarray(coefficients, float)
    if c.shape != (8,) or np.any(~np.isfinite(c)) or np.any(c[[0, 1]] <= 0):
        raise ValueError('eight finite coefficients and positive LJ pair required')
    fit = JointFit('vector_material_research', float(scalar_decay), c[:5], 0., None, None)
    bulk = build_joint_model(fit)
    # Fixed crystallographic reference: only fits satisfying the exact zero-force
    # equality can use it as equilibrium. Construction is NOT certification.
    bulk.a0 = IDEAL_H
    face = FCC111ActiveInterface(bulk, tolerance=tolerance/10)
    surface = AngularInterfaceResearchSurface(face, c[5], angular_decay=float(angular_decay),
        vector_amplitude_ev=c[6], quadrupole_amplitude_ev=c[7], tolerance=tolerance)
    return FullRegistryInterface(surface)


class VectorCoefficientBasis:
    """Analytic ten-component jet times [u,v,A,B,C,D3,D1,D2]."""
    def __init__(self, scalar_decay, angular_decay, *, tolerance=2e-11):
        self.scalar_decay, self.angular_decay = float(scalar_decay), float(angular_decay)
        if (not np.all(np.isfinite([self.scalar_decay, self.angular_decay, tolerance]))
                or min(self.scalar_decay, self.angular_decay, tolerance) <= 0):
            raise ValueError('positive decays/tolerance required')
        self.model = build_coefficient_surface(scalar_decay, angular_decay,
            [1., 1., 1., 1., 1., 1., 1., 1.], tolerance=tolerance)
        self.h = self.model.h

    @lru_cache(maxsize=192)
    def jet(self, state):
        q = np.asarray(state, float)
        if q.shape != (3,) or np.any(~np.isfinite(q)) or q[0] <= 0:
            raise ValueError('finite (positive a,x,y) required')
        m = self.model
        columns = [m._pair(q, 6)[0], -m._pair(q, 3)[0]]
        changes = m._density(q)[0]
        reference = m.face.bulk.embedding.rho_ref
        for c in np.eye(3):
            columns.append(_embedding_jet(changes, m.face.rho_bulk,
                MinimalConvexEmbedding(*c, reference)))
        columns.extend(m._angular(q, moment)[0] for _, moment in m.moments)
        result = np.column_stack(columns)
        if result.shape != (10, 8) or np.any(~np.isfinite(result)):
            raise FloatingPointError('invalid coefficient jet')
        result.setflags(write=False)
        return result

    def evaluate(self, state, coefficients):
        c = np.asarray(coefficients, float)
        if c.shape != (8,) or np.any(~np.isfinite(c)):
            raise ValueError('eight finite coefficients required')
        return VectorInterfaceEvaluation.from_jet(self.jet(tuple(state))@c)


@dataclass(frozen=True)
class MaterialObservation:
    name: str
    target: float
    scale: float
    units: str
    role: str
    state: tuple | None = None
    jet_weights: tuple | None = None
    bulk_index: int | None = None

    def __post_init__(self):
        if not np.isfinite(self.target) or not np.isfinite(self.scale) or self.scale <= 0:
            raise ValueError('finite target and positive normalization scale required')
        if self.role not in ('exact', 'fit', 'heldout'):
            raise ValueError('explicit observation role required')


def matched_observations(source):
    """Declared 0-K source conditions; no experimental yield in this dataset.

    5% bulk and 10% interface energy/curvature model-discrepancy scales are
    NOT measurement uncertainties. Zero-force tolerance corresponds to .25 GPa.
    Reverse barrier = saddle minus fault: do NOT count it as independent data.
    """
    bulk, scales = independent_bulk_targets()
    units = ('eV/strain', 'eV/atom', 'eV/strain^2', 'eV/strain^2', 'eV/strain^2')
    rows = [MaterialObservation(n, float(y), float(z), unit,
        'exact' if i < 2 else 'fit', bulk_index=i)
        for i, (n, y, z, unit) in enumerate(zip(NAMES, bulk, scales, units))]
    tau = np.array([.5, math.sqrt(3)/6])
    fault = stationary_state(source, np.r_[1.025*IDEAL_H, tau], expected_index=0)
    saddle = stationary_state(source, np.r_[1.04*IDEAL_H, .74*tau], expected_index=1)
    if not fault['valid'] or not saddle['valid']:
        raise ValueError('source full-vector stationary state is not validated')
    states = {'perfect': np.array([IDEAL_H, 0., 0.]), 'fault': fault['q'], 'saddle': saddle['q']}
    def append(name, state, component, role='fit', weights=None):
        v = source.evaluate(state)
        jet = np.r_[v.energy, v.gradient, v.hessian[0], v.hessian[1, 1:], v.hessian[2, 2]]
        w = np.eye(10)[component] if weights is None else np.asarray(weights, float)
        target = float(w@jet)
        scale = float(UNITS.traction_mpa_to_force(250.)) if component in (1, 2, 3) else .10*abs(target)
        unit = 'eV/cell' if component == 0 else ('eV/L0' if component < 4 else 'eV/L0^2')
        rows.append(MaterialObservation(name, target, scale, unit, role, tuple(state), tuple(w)))
    append('perfect_Haa', states['perfect'], 4)
    append('perfect_Hxx', states['perfect'], 7)
    for name in ('saddle', 'fault'):
        append(name+'_energy_at_source_state', states[name], 0)
        append(name+'_normal_force_at_source_state', states[name], 1)
    direction = tau/np.linalg.norm(tau)
    w = np.zeros(10); w[2:4] = direction
    append('saddle_path_force_at_source_state', states['saddle'], 2, weights=w)
    for ratio in (1.1, 1.5, 2., 40.):
        # At 40h source cross interactions are exactly beyond SOURCE cutoff.
        # Analytic candidate remains an infinite sum at the same finite a.
        append(f'opening_{ratio:g}h_energy', [ratio*IDEAL_H, 0., 0.], 0)
    for label, state in [('direct_quarter', [IDEAL_H, .25, 0.]),
                         ('direct_half', [IDEAL_H, .5, 0.]),
                         ('oblique', [1.08*IDEAL_H, .21, -.12]),
                         ('opening_1.25h', [1.25*IDEAL_H, 0., 0.]),
                         ('opening_3h', [3*IDEAL_H, 0., 0.])]:
        append(label+'_energy', state, 0, 'heldout')
    for name in ('saddle', 'fault'):
        for comp in (4, 7):
            append(name+('_Haa' if comp == 4 else '_Hxx'), states[name], comp, 'heldout')
    return rows, states


def observation_matrix(basis, observations):
    bulk = even_basis(basis.scalar_decay, basis.angular_decay, convex=True)[:5]
    rows = []
    for obs in observations:
        rows.append(bulk[obs.bulk_index] if obs.bulk_index is not None else
                    np.asarray(obs.jet_weights)@basis.jet(obs.state))
    return np.asarray(rows)


def fit_coefficients(matrix, observations, *, opening_inequalities=None):
    selected = [i for i, o in enumerate(observations) if o.role != 'heldout']
    if selected[:2] != [0, 1] or any(observations[i].role != 'exact' for i in (0, 1)):
        raise ValueError('first two rows must be exact bulk geometry/cohesion constraints')
    target = np.array([o.target for o in observations])
    scales = np.array([o.scale for o in observations])
    result = constrained_matrix_fit(matrix[selected], target[selected], scales[selected],
        coefficient_order=COEFFICIENTS, lower_rest=[0., -np.inf, 0., 0., -np.inf, -np.inf],
        nonnegative_observation_rows=opening_inequalities)
    result['selected_rows'] = selected
    return result


def exact_constraint_tangent(matrix, target):
    """c=offset+T z removes force/cohesion equalities before identifiability."""
    block = np.asarray(matrix)[:2, :2]
    offset = np.r_[np.linalg.solve(block, np.asarray(target)[:2]), np.zeros(6)]
    tangent = np.vstack((-np.linalg.solve(block, np.asarray(matrix)[:2, 2:]), np.eye(6)))
    return offset, tangent


def coefficient_identifiability(matrix, observations, coefficients):
    selected = [i for i, o in enumerate(observations) if o.role == 'fit']
    _, tangent = exact_constraint_tangent(matrix, [o.target for o in observations])
    scales = np.array([observations[i].scale for i in selected])
    # Parameter units: 1 eV for small amplitudes, actual amplitude otherwise.
    # These are coordinate normalizations, not a probability prior.
    normalization = np.maximum(1., abs(np.asarray(coefficients)[2:]))
    jac = (matrix[selected]@tangent)*normalization/scales[:, None]
    u, singular, vt = np.linalg.svd(jac, full_matrices=False)
    floor = np.finfo(float).eps*max(jac.shape)*singular[0]
    rank = int(np.sum(singular > floor))
    return dict(jacobian=jac, singular_values=singular, rank=rank,
        condition_number=float(singular[0]/singular[-1]) if singular[-1] > floor else None,
        normalization=normalization, right_vectors=vt,
        full_radial_identifiability_checked=False, statistical_confidence_claimed=False)


def opening_traction_extrema(model, count, *, upper_ratio=12.):
    """All sign-bracketed W_aa roots at zero registry on a declared interval.

    Independent count/tolerance refinement is required. This is not a proof
    excluding double roots, all registry directions, or features beyond the
    interval. In particular, coarse nonnegative W_a samples are insufficient.
    """
    if count < 5 or not np.isfinite(upper_ratio) or upper_ratio <= 2:
        raise ValueError('resolved bracket count and opening domain required')
    grid=np.unique(np.r_[np.linspace(1.00001,2.,count),
                         np.geomspace(2.,upper_ratio,count//2)])*model.h
    curvature=np.array([model.evaluate((a,0.,0.)).hessian[0,0] for a in grid])
    if np.any(~np.isfinite(curvature)):
        raise FloatingPointError('nonfinite opening curvature')
    points=[]
    for left,right,cl,cr in zip(grid[:-1],grid[1:],curvature[:-1],curvature[1:]):
        if cl==0:
            a=left
        elif cl*cr<0:
            a=brentq(lambda a:model.evaluate((a,0.,0.)).hessian[0,0],left,right,xtol=2e-12)
        else:
            continue
        value=model.evaluate((a,0.,0.))
        points.append(dict(a_over_h=float(a/model.h),force_eV_L0=float(value.gradient[0]),
                           Haa_eV_L0sq=float(value.hessian[0,0])))
    if curvature[-1]==0:
        value=model.evaluate((grid[-1],0.,0.))
        points.append(dict(a_over_h=float(upper_ratio),force_eV_L0=float(value.gradient[0]),
                           Haa_eV_L0sq=float(value.hessian[0,0])))
    return points
