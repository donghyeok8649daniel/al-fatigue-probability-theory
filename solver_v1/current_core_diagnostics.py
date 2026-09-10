"""Unscaled geometric observations of a straight-row static core.

The compact x-phase winding locates integer projected circulation on lattice
triangles; it is NOT a complete vector/partial-dislocation classifier. An
adjacent-plane profile compares logical neighboring rows with the actual ABC
offset, not an invented continuum interpolation or a fitted core radius.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import Delaunay


def reference_configuration_seed(core, metadata, state, *, length_scale_m):
    """Copy coordinates ONLY as a declared zero-load multistart initial guess.

    Different material coefficients are allowed here, unlike continuation.
    Geometry, units, logical free rows and topology must be identical. The
    target's own exterior, energy, elastic tensor and final forces remain.
    """
    if (core.shear_traction != 0 or metadata['shear_traction_MPa'] != 0
            or metadata['burgers_sign'] != core.burgers_sign
            or not np.isclose(metadata['length_scale_m'],length_scale_m,rtol=2e-13,atol=0)
            or not np.allclose([metadata['b_reduced'],metadata['h_reduced']],
                               [core.rows.b,core.rows.h],rtol=0,atol=1e-13)
            or not np.allclose(metadata['center_over_L0'],core.far_field.center,rtol=0,atol=1e-13)):
        raise ValueError('configuration seed requires identical zero-load geometry, units and topology')
    indices=[tuple(i) for i in core.indices[core.free_ids]]
    if set(state)!=set(indices):
        raise ValueError('configuration seed requires the same logical free domain')
    result=np.asarray([state[i] for i in indices],float)
    if result.shape!=core.initial.shape or np.any(~np.isfinite(result)):
        raise ValueError('finite three-component source coordinates required')
    return result.copy()


def translated_core_seed(core, field, j_steps):
    """Initialize crystallographic glide; the exterior is NOT translated.

    j -> j+1 is a2=(b/2,sqrt(3)b/2,0), a valid FCC in-plane lattice vector.
    The infinite line is parallel to x, so the core's glide displacement is
    sqrt(3)b/2 along y. Burgers direction x and glide direction y differ.
    Move only the defect field; keep the homogeneous affine shear unchanged.
    This is an initial guess, never a forced trajectory or a pinning law.
    """
    if int(j_steps) != j_steps:
        raise ValueError('integer crystallographic row translation required')
    j_steps = int(j_steps)
    if abs(j_steps)*core.rows.d >= core.free_radius/2:
        raise ValueError('translation must remain well inside the current free disk')
    full = core.full_field(field)
    affine = core.far_field.displacement(core.xyz,
        shear_traction=core.shear_traction, burgers_sign=0)
    lookup = {tuple(index): i for i, index in enumerate(core.indices)}
    result = np.empty_like(field)
    for number, (index, target) in enumerate(zip(core.indices[core.free_ids], core.free_ids)):
        source_index = (int(index[0])-j_steps, int(index[1]))
        if source_index not in lookup:
            raise ValueError('source row not represented; enlarge the neighborhood')
        source = lookup[source_index]
        result[number] = full[source]-affine[source]+affine[target]
    return result


def twofold_core_partner(core, field):
    """Exact zero-load FCC twofold partner about the straight [1,-1,0] line.

    The centered disk maps (j,l)->(-j,1-l). Subtracting the same Volterra
    reference handles the half-row x-phase shift without choosing arbitrary
    unwrapping. The residual transforms as (wx,wy,wz)->(wx,-wy,-wz).
    This is a lattice symmetry, NOT a prescribed physical core motion.
    """
    expected = (core.rows.d/6,core.rows.h/2)
    if core.shear_traction != 0 or not np.allclose(core.far_field.center,expected,atol=1e-13,rtol=0):
        raise ValueError('twofold partner requires the declared symmetric zero-load boundary')
    full = core.full_field(field)
    residual = full-core.boundary
    lookup = {tuple(index):i for i,index in enumerate(core.indices)}
    result = np.empty_like(field)
    free = set(core.free_ids)
    for k,(j,layer) in enumerate(core.indices[core.free_ids]):
        target = core.free_ids[k]
        partner = lookup.get((-j,1-layer))
        if partner not in free:
            raise ValueError('free domain is not closed under the stated FCC symmetry')
        result[k] = core.boundary[target]+residual[partner]*np.array([1.,-1.,-1.])
    return result


def phase_winding_cells(core, field):
    full = core.full_field(field)
    affine = core.far_field.displacement(core.xyz,
        shear_traction=core.shear_traction, burgers_sign=0)
    phase = 2*np.pi*(full[:, 0]-affine[:, 0])/core.rows.b
    triangles = Delaunay(core.xyz[:, 1:]).simplices.copy()
    points = core.xyz[triangles, 1:]
    d1, d2 = points[:, 1]-points[:, 0], points[:, 2]-points[:, 0]
    reverse = d1[:, 0]*d2[:, 1]-d1[:, 1]*d2[:, 0] < 0
    triangles[reverse] = triangles[reverse][:, [0, 2, 1]]
    p = phase[triangles]
    winding = np.angle(np.exp(1j*(np.roll(p, -1, axis=1)-p))).sum(axis=1)/(2*np.pi)
    residual = float(np.max(abs(winding-np.rint(winding))))
    if residual > 1e-10:
        raise ArithmeticError('phase circuit failed integer closure')
    free = set(core.free_ids)
    rows = []
    for index in np.flatnonzero(abs(winding) > .5):
        ids = triangles[index]
        reference = core.xyz[ids, 1:].mean(axis=0)
        displaced = (core.xyz[ids, 1:]+full[ids, 1:]).mean(axis=0)
        rows.append(dict(winding=float(winding[index]), y_reference=float(reference[0]),
            z_reference=float(reference[1]), y_displaced=float(displaced[0]),
            z_displaced=float(displaced[1]), all_rows_free=all(i in free for i in ids)))
    return dict(cells=rows, total_winding=float(winding.sum()),
        maximum_integer_closure_residual=residual,
        projected_integer_circulation_only=True, partial_core_width_inferred=False)


def adjacent_registry_profile(core, field, *, lower_layer=0):
    full = core.full_field(field)
    affine = core.far_field.displacement(core.xyz,
        shear_traction=core.shear_traction, burgers_sign=0)
    lookup = {tuple(index): i for i, index in enumerate(core.indices)}
    free = set(core.free_ids)
    rows = []
    for j, layer in core.indices:
        if layer != lower_layer or (j, layer+1) not in lookup:
            continue
        lo, hi = lookup[j, layer], lookup[j, layer+1]
        raw = full[hi]-full[lo]
        defect = raw-(affine[hi]-affine[lo])
        rows.append(dict(j=int(j), lower_layer=int(layer),
            y_over_L0=float((core.xyz[lo, 1]+core.xyz[hi, 1])/2),
            raw_delta_x=float(raw[0]), raw_delta_y=float(raw[1]), raw_delta_z=float(raw[2]),
            affine_subtracted_delta_x=float(defect[0]),
            affine_subtracted_delta_y=float(defect[1]), affine_subtracted_delta_z=float(defect[2]),
            both_rows_free=lo in free and hi in free))
    return rows


def inner_field_difference(core, first, second, *, radius):
    """Same logical grid and boundary, phase-equivalent x comparison only."""
    if not np.isfinite(radius) or radius <= 0 or radius > core.free_radius:
        raise ValueError('comparison radius must be within the free disk')
    a, b = core.full_field(first), core.full_field(second)
    ids = core.free_ids
    mask = np.linalg.norm(core.xyz[ids, 1:]-core.far_field.center, axis=1) < radius
    delta = b[ids[mask]]-a[ids[mask]]
    delta[:, 0] -= core.rows.b*np.floor(delta[:, 0]/core.rows.b+.5)
    if not len(delta):
        raise ValueError('no rows inside the comparison radius')
    return dict(rows=len(delta), maximum_vector_change_over_L0=float(np.max(np.linalg.norm(delta, axis=1))),
        rms_vector_change_over_L0=float(np.sqrt(np.mean(np.sum(delta*delta, axis=1)))))
