"""Intact finite Si tensile specimens; research geometry and diagnostics only.

No bonds/atoms are removed to seed a crack. The free external surfaces are bare,
unpassivated and initially ideal. A graph disconnection is a geometric diagnostic,
not a calibrated definition of the first crack, nor an initiation probability.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def intact_prism(*, repeats=(3, 3, 10), lattice_A=5.470992313871248, grip_width_A=6.):
    """Diamond: x=[1,-1,0], y=[0,0,-1], tensile z=[1,1,0].

    Nominal cross section is the original periodic-cell rectangle before cutting
    free surfaces. This convention does not assign a unique atomic surface area.
    Both grips translate rigidly; their internal bonds are not strained by loading.
    """
    from ase.lattice.cubic import Diamond
    repeats = np.asarray(repeats)
    if (repeats.shape != (3,) or np.any(repeats < 1)
            or np.any(repeats != np.rint(repeats))
            or not np.isfinite(lattice_A) or lattice_A <= 0
            or not np.isfinite(grip_width_A) or grip_width_A <= 0):
        raise ValueError('positive integer repeats and positive finite lengths required')
    atoms = Diamond('Si', directions=[[1, -1, 0], [0, 0, -1], [1, 1, 0]],
                    size=tuple(repeats.astype(int)), latticeconstant=lattice_A, pbc=True)
    lengths = np.diag(atoms.cell).copy()
    # Nonzero unit-cell phase avoids ambiguous atoms exactly on a cut boundary.
    phase = lattice_A*np.array([1/np.sqrt(2), 1., 1/np.sqrt(2)])/8
    atoms.positions[:] = (atoms.positions-phase) % lengths
    atoms.positions[:] -= lengths/2
    atoms.pbc = False
    z = atoms.positions[:, 2]
    low, high = z.min(), z.max()
    lower, upper = z < low+grip_width_A, z > high-grip_width_A
    free = ~(lower | upper)
    if np.any(lower & upper) or np.count_nonzero(free) < 4:
        raise ValueError('grips overlap or leave no resolved free gauge')
    gauge = float(z[upper].mean()-z[lower].mean())
    return atoms, dict(lower=lower, upper=upper, free=free, gauge_length_A=gauge,
        area_A2=float(lengths[0]*lengths[1]), nominal_lengths_A=lengths,
        directions_cubic=np.array([[1, -1, 0], [0, 0, -1], [1, 1, 0]]))


def prescribed_grips(reference, *, lower, upper, extension_A):
    r = np.array(reference, float, copy=True)
    lower, upper = np.asarray(lower, bool), np.asarray(upper, bool)
    if (r.ndim != 2 or r.shape[1] != 3 or lower.shape != (len(r),)
            or upper.shape != lower.shape or np.any(lower & upper)
            or not np.any(lower) or not np.any(upper)
            or not np.all(np.isfinite(r)) or not np.isfinite(extension_A)):
        raise ValueError('finite coordinates, extension and disjoint nonempty grips required')
    r[lower, 2] -= extension_A/2
    r[upper, 2] += extension_A/2
    return r


def grip_observables(positions, forces, *, lower, upper, area_A2):
    """Exact partial dU/d(extension) for rigid, opposite grip translations.

    The total derivative on a relaxed branch equals this only when free forces
    vanish. This is not the vacuum-volume-normalized ASE stress.
    """
    r, f = np.asarray(positions, float), np.asarray(forces, float)
    lower, upper = np.asarray(lower, bool), np.asarray(upper, bool)
    if (r.ndim != 2 or r.shape[1] != 3 or f.shape != r.shape
            or lower.shape != (len(r),) or upper.shape != lower.shape
            or np.any(lower & upper) or not np.any(lower) or not np.any(upper)
            or not np.all(np.isfinite(r)) or not np.all(np.isfinite(f))
            or not np.isfinite(area_A2) or area_A2 <= 0):
        raise ValueError('finite geometry/force, positive area and disjoint grips required')
    fixed, free = lower | upper, ~(lower | upper)
    reaction = -f[fixed]
    conjugate = float((f[lower, 2].sum()-f[upper, 2].sum())/2)
    return dict(conjugate_force_eV_A=conjugate,
        nominal_stress_GPa=conjugate/area_A2*160.2176634,
        lower_reaction_eV_A=(-f[lower].sum(axis=0)).tolist(),
        upper_reaction_eV_A=(-f[upper].sum(axis=0)).tolist(),
        reaction_force_residual_eV_A=float(np.linalg.norm(reaction.sum(axis=0))),
        reaction_torque_eV=float(np.linalg.norm(np.cross(r[fixed], reaction).sum(axis=0))),
        total_internal_force_eV_A=float(np.linalg.norm(f.sum(axis=0))),
        total_internal_torque_eV=float(np.linalg.norm(np.cross(r, f).sum(axis=0))),
        free_force_max_eV_A=float(np.linalg.norm(f[free], axis=1).max()) if np.any(free) else 0.)


def connectivity_diagnostics(positions, *, lower, upper, cutoffs_A=(2.8, 3.1, 3.4, 3.7)):
    """Cutoff sensitivity, never a standalone first-crack detector."""
    r = np.asarray(positions, float)
    lower, upper = np.asarray(lower, bool), np.asarray(upper, bool)
    cutoffs = np.asarray(cutoffs_A, float)
    if (r.ndim != 2 or r.shape[1] != 3 or lower.shape != (len(r),)
            or upper.shape != lower.shape or not np.any(lower) or not np.any(upper)
            or not np.all(np.isfinite(r)) or cutoffs.ndim != 1
            or not len(cutoffs) or not np.all(np.isfinite(cutoffs)) or np.any(cutoffs <= 0)):
        raise ValueError('finite geometry, nonempty grips and positive cutoffs required')
    tree, rows = cKDTree(r), []
    for cutoff in cutoffs:
        pairs = tree.query_pairs(float(cutoff), output_type='ndarray')
        graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(r), len(r)))
        count, labels = connected_components(graph, directed=False)
        bridge = bool(np.intersect1d(np.unique(labels[lower]), np.unique(labels[upper])).size)
        rows.append(dict(cutoff_A=float(cutoff), components=int(count), grip_connected=bridge,
            largest_component_atoms=int(np.bincount(labels).max()), pair_count=len(pairs)))
    return rows
