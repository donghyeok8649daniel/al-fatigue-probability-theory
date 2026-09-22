"""Internal diamond-basis relaxation for a continuum-loaded Si atomic boundary.

The 0 K optical correction is computed from the SAME atomistic engine as the
bulk elastic constants. It is a first-order boundary expansion, not a new
potential, a finite-temperature PMF, or flexible crack-tip equilibration.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .silicon_specimen_research import AtomicCrackBoundary, SI_111_FRAME


def internal_strain_response(engine, lattice_A, *, step=1e-4):
    """du_B-A/d(engineering cubic strain) in A, from H_uu du + B de = 0."""
    if not np.isfinite(lattice_A) or lattice_A <= 0 or not np.isfinite(step) or step <= 0:
        raise ValueError('positive finite lattice and step required')
    cell = lattice_A*.5*np.array([[0., 1., 1.], [1., 0., 1.], [1., 1., 0.]])
    positions = lattice_A*np.array([[0., 0., 0.], [.25, .25, .25]])
    h = []
    for direction in np.eye(3):
        plus, minus = positions.copy(), positions.copy()
        plus[1] += step*direction
        minus[1] -= step*direction
        h.append((engine.evaluate(plus, cell).gradient[1]-engine.evaluate(minus, cell).gradient[1])/(2*step))
    h = np.stack(h, axis=1)
    if np.linalg.eigvalsh((h+h.T)/2).min() <= 0 or np.max(abs(h-h.T)) > 1e-7:
        raise ValueError('optical bulk Hessian is not stable/symmetric')
    modes = []
    for i, j in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]:
        mode = np.zeros((3, 3))
        mode[i, j] += .5
        mode[j, i] += .5
        modes.append(mode)
    coupling = []
    for mode in modes:
        fp, fm = np.eye(3)+step*mode, np.eye(3)-step*mode
        coupling.append((engine.evaluate(positions@fp.T, cell@fp.T).gradient[1]
                         -engine.evaluate(positions@fm.T, cell@fm.T).gradient[1])/(2*step))
    coupling = np.stack(coupling, axis=1)
    response = np.linalg.solve(h, -coupling)
    return response, dict(step=step, hessian_eV_A2=h.tolist(), coupling_eV_A=coupling.tolist(),
                           response_A_per_strain=response.tolist(),
                           optical_eigenvalues_eV_A2=np.linalg.eigvalsh((h+h.T)/2).tolist(),
                           engineering_order=['xx', 'yy', 'zz', 'yz', 'xz', 'xy'])


@dataclass
class InternallyRelaxedBoundary(AtomicCrackBoundary):
    internal_response: np.ndarray
    sublattice_sign: np.ndarray

    def displaced(self, elasticity, K_MPa_sqrt_m):
        positions = self.reference.copy()
        u, gradient, _ = elasticity.field(positions[:, :2], K_MPa_sqrt_m)
        strain_local = np.zeros((len(positions), 3, 3))
        strain_local[:, :2, :2] = .5*(gradient+gradient.swapaxes(-1, -2))
        strain_cubic = np.einsum('ai,nij,jb->nab', SI_111_FRAME.T, strain_local, SI_111_FRAME)
        engineering = np.stack([strain_cubic[:, 0, 0], strain_cubic[:, 1, 1], strain_cubic[:, 2, 2],
                                2*strain_cubic[:, 1, 2], 2*strain_cubic[:, 0, 2], 2*strain_cubic[:, 0, 1]], axis=1)
        relative = (engineering@self.internal_response.T)@SI_111_FRAME.T
        positions[:, :2] += u
        positions += .5*self.sublattice_sign[:, None]*relative
        return positions


def with_internal_relaxation(boundary, lattice_A, response):
    """Identify A/B planes using the known centered (111) shuffle cut.

    The big shuffle gap is A at -3h/8 to B at +3h/8, modulo h=a/sqrt(3).
    Average internal displacement is zero: A gets -u/2, B gets +u/2.
    """
    response = np.asarray(response, float)
    if response.shape != (3, 6) or not np.all(np.isfinite(response)):
        raise ValueError('finite (3,6) optical response required')
    phase = np.mod(boundary.reference[:, 1]/(lattice_A/np.sqrt(3)), 1.)
    a = np.isclose(phase, .625, rtol=0, atol=1e-7)
    b = np.isclose(phase, .375, rtol=0, atol=1e-7)
    if not np.all(a | b) or not np.any(a) or not np.any(b):
        raise ValueError('unrecognized diamond sublattice convention')
    return InternallyRelaxedBoundary(**vars(boundary), internal_response=response,
                                     sublattice_sign=np.where(b, 1., -1.))
