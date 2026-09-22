"""Fixed-volume chemical concentration and constrained cleavage diagnostics.

Research utilities only. Neutral atom counts do not specify carrier density.
Rigid separation energies are neither finite-T PMFs nor specimen strengths.
All cells and positions use ASE row-vector convention, in Angstrom/eV.
"""
from __future__ import annotations

import numpy as np


def chemical_concentration_cm3(dopant_count, *, reference_material_volume_A3):
    """Use undeformed MATERIAL volume; opening a vacuum gap cannot dilute B/P."""
    if isinstance(dopant_count, (bool, np.bool_)) or not isinstance(dopant_count, (int, np.integer)) or dopant_count < 0:
        raise ValueError('nonnegative integer chemical atom count required')
    volume = float(reference_material_volume_A3)
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError('positive finite reference material volume required')
    return float(dopant_count)*1e24/volume


def constrained_energy_derivative(positions, forces, cell, stress, *, position_tangent, cell_tangent):
    """dE/dq from Cartesian force and tensile-positive Cauchy stress.

    For row cell C, affine velocities are r C^-1 C'. With arbitrary prescribed
    r'(q), E' = V sigma : (C^-1 C') - sum_i F_i . (r'_i-r_i C^-1 C').
    This includes the non-affine correction omitted by simply using V sigma.
    Cell/stress convention is checked independently in the MACE execution.
    """
    r, f, c, s, dr, dc = [np.asarray(x, float) for x in
                          (positions, forces, cell, stress, position_tangent, cell_tangent)]
    if (r.ndim != 2 or r.shape[1] != 3 or f.shape != r.shape or dr.shape != r.shape
            or c.shape != (3, 3) or s.shape != (3, 3) or dc.shape != (3, 3)
            or not all(np.all(np.isfinite(x)) for x in (r, f, c, s, dr, dc))):
        raise ValueError('finite compatible atom, cell and tensor arrays required')
    volume = float(np.linalg.det(c))
    if volume <= 0 or not np.allclose(s, s.T, atol=1e-10, rtol=1e-10):
        raise ValueError('right-handed nonsingular cell and symmetric Cauchy stress required')
    tangent = np.linalg.solve(c, dc)
    return float(volume*np.sum(s*tangent)-np.sum(f*(dr-r@tangent)))


def deterministic_substitution_order(positions, cell, *, arrangement):
    """Nested selected configurations, NOT a random-alloy ensemble.

    The cleavage is the periodic z boundary. 'interface_cluster' grows around
    one interface site. 'interface_spread' uses farthest sites in the two
    adjacent planes. 'bulk_spread' starts deepest inside and spreads in 3D.
    Ties use stable original atom indices; rotations need not preserve this
    placement rule, so physical symmetry tests transform the selected atoms.
    """
    r, c = np.asarray(positions, float), np.asarray(cell, float)
    if (r.ndim != 2 or r.shape[1] != 3 or not len(r) or c.shape != (3, 3)
            or not np.all(np.isfinite(r)) or not np.all(np.isfinite(c))
            or not np.allclose(c, np.diag(np.diag(c)), atol=1e-12, rtol=0)
            or np.any(np.diag(c) <= 0)):
        raise ValueError('finite positions and positive orthorhombic cell required')
    length = np.diag(c)
    p = r % length
    depth = np.minimum(p[:, 2], length[2]-p[:, 2])
    delta = p[:, None, :]-p[None, :, :]
    delta -= np.rint(delta/length)*length
    distance2 = np.sum(delta**2, axis=-1)
    xy = p[:, :2]-length[:2]/2
    xy -= np.rint(xy/length[:2])*length[:2]
    xy2 = np.sum(xy**2, axis=1)
    surface = np.flatnonzero(depth <= depth.min()+1e-6)
    if arrangement == 'interface_cluster':
        seed = int(surface[np.argmin(xy2[surface])])
        return np.lexsort((np.arange(len(p)), distance2[seed]))
    if arrangement == 'interface_spread':
        pool = surface
        seed = int(surface[np.argmin(xy2[surface])])
    elif arrangement == 'bulk_spread':
        pool = np.arange(len(p))
        deep = np.flatnonzero(depth >= depth.max()-1e-6)
        seed = int(deep[np.argmin(xy2[deep])])
    else:
        raise ValueError('unknown deterministic arrangement')
    order = [seed]
    remaining = [int(i) for i in pool if i != seed]
    while remaining:
        scores = distance2[np.ix_(remaining, order)].min(axis=1)
        chosen = int(np.argmax(scores))
        order.append(remaining.pop(chosen))
    return np.asarray(order, int)


def rigid_opening_geometry(positions, cell, opening_A):
    """One gap at the periodic z boundary; atom positions stay fixed.

    Each period creates two exposed surfaces at ONE interface, so Wsep is
    Delta E / area, not Delta E / (2 area). Input atoms must remain in cell.
    """
    r, c = np.array(positions, float, copy=True), np.array(cell, float, copy=True)
    q = float(opening_A)
    if (r.ndim != 2 or r.shape[1] != 3 or c.shape != (3, 3)
            or not np.all(np.isfinite(r)) or not np.all(np.isfinite(c))
            or not np.isfinite(q) or q < 0
            or not np.allclose(c, np.diag(np.diag(c)), atol=1e-12, rtol=0)
            or np.any(np.diag(c) <= 0)
            or np.any(r < -1e-9) or np.any(r >= np.diag(c)+1e-9)):
        raise ValueError('wrapped positions, positive orthorhombic cell and nonnegative opening required')
    c[2, 2] += q
    dc = np.zeros((3, 3)); dc[2, 2] = 1.
    return r, c, np.zeros_like(r), dc
