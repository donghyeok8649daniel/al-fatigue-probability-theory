"""Quasistatic, section-averaged axial equilibrium coupled to the local PDE.

The prescribed scalar stress is a UNIFORM AXIAL END TRACTION, so N(t) =
traction(t) * loaded_end_area. Sections carry sigma(x,t) = N(t)/A(x).
This is a straight, centered bar reduction, not a 3D stress-concentration
solver. Triangle surface area and correlation area never scale local stress
or local probability. See AXIAL_SPECIMEN.md for the assumptions and tests.
"""
from dataclasses import dataclass, replace

import numpy as np

from .surface_setup import mesh_id
from .solver_adapter import run_ui_analysis


@dataclass(frozen=True)
class AxialSections:
    mesh_sha256: str
    origin_mm: np.ndarray
    axis: np.ndarray
    edges_mm: np.ndarray
    positions_mm: np.ndarray
    areas_mm2: np.ndarray
    loaded_area_mm2: float
    loaded_faces: np.ndarray
    opposite_faces: np.ndarray
    face_sections: np.ndarray

    @property
    def stress_factors(self):
        return self.loaded_area_mm2 / self.areas_mm2


def _section(vertices, faces, axis, origin, station, tolerance):
    """Intersect canonical mesh edges, then integrate one closed polygon.

    Sharing edge IDs avoids tolerance-based welding of intersection points.
    Stations must avoid vertices (handled by prepare_axial_sections).
    """
    relative = vertices - origin
    depth = relative @ axis - station
    edges, inverse = np.unique(np.sort(np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1),
        axis=0, return_inverse=True)
    crossing = (depth[edges[:, 0]] < 0) != (depth[edges[:, 1]] < 0)
    face_edges = inverse.reshape(3, -1).T
    cut_faces = face_edges[crossing[face_edges].sum(axis=1) == 2]
    if not len(cut_faces):
        raise ValueError('axial.section')
    segments = cut_faces[crossing[cut_faces]].reshape(-1, 2)
    ids = np.flatnonzero(crossing)
    a, b = relative[edges[ids, 0]], relative[edges[ids, 1]]
    fraction = -depth[edges[ids, 0]] / (depth[edges[ids, 1]] - depth[edges[ids, 0]])
    points = a + fraction[:, None] * (b - a)
    lookup = {int(edge): i for i, edge in enumerate(ids)}
    adjacent = {int(edge): [] for edge in ids}
    for a, b in segments:
        adjacent[int(a)].append(int(b)); adjacent[int(b)].append(int(a))
    if any(len(neighbors) != 2 for neighbors in adjacent.values()):
        raise ValueError('axial.section')
    first = int(ids[0]); order = [first]; previous = None; current = first
    while True:
        neighbors = adjacent[current]
        following = neighbors[0] if neighbors[0] != previous else neighbors[1]
        if following == first:
            break
        if len(order) >= len(ids):
            raise ValueError('axial.section')
        order.append(following); previous, current = current, following
    if len(order) != len(ids):  # holes / separate components are outside this reduction
        raise ValueError('axial.section')
    polygon = points[[lookup[i] for i in order]]
    e1 = np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    xy = np.column_stack((polygon @ e1, polygon @ e2))
    following = np.roll(xy, -1, axis=0)
    cross = xy[:, 0]*following[:, 1] - following[:, 0]*xy[:, 1]
    signed_area = cross.sum()/2
    if not np.isfinite(signed_area) or abs(signed_area) <= tolerance**2:
        raise ValueError('axial.section')
    centroid = ((xy + following)*cross[:, None]).sum(axis=0)/(6*signed_area)
    return abs(float(signed_area)), centroid


def prepare_axial_sections(mesh, loaded_faces, count=16):
    """Require a complete planar end, an opposite end, and centered sections.

    The opposite end supplies the equal reaction. There are no lateral loads,
    distributed axial loads, bending, torsion, inertia or evolving geometry.
    The sampling grid is independent of triangle sizes and can be refined.
    """
    if mesh is None:
        raise ValueError('axial.mesh')
    if isinstance(count, bool) or int(count) != count or not 2 <= count <= 64:
        raise ValueError('axial.count_error')
    try:
        _ = mesh.enclosed_volume_mm3
    except ValueError as exc:
        raise ValueError('axial.mesh') from exc
    ids = np.asarray(loaded_faces)
    if (ids.ndim != 1 or ids.dtype.kind not in 'iu' or not ids.size or ids.min() < 0
            or ids.max() >= len(mesh.faces) or len(np.unique(ids)) != len(ids)):
        raise ValueError('axial.end')
    triangles = mesh.vertices[mesh.faces]
    areas = mesh.areas
    cap = triangles[ids]
    origin = np.average(cap.mean(axis=1), weights=areas[ids], axis=0)
    area_vectors = np.cross(cap[:, 1]-cap[:, 0], cap[:, 2]-cap[:, 0])/2
    axis = area_vectors.sum(axis=0)
    magnitude = np.linalg.norm(axis)
    if magnitude <= 1e-12*areas[ids].sum():
        raise ValueError('axial.end')
    axis /= magnitude
    depth = (mesh.vertices-origin) @ axis
    if depth.mean() < 0:
        axis = -axis; depth = -depth
    scale = float(np.ptp(mesh.vertices, axis=0).max())
    tolerance = max(scale*1e-9, 1e-11)
    if np.max(np.abs(depth[mesh.faces[ids]])) > tolerance or depth.min() < -tolerance:
        raise ValueError('axial.end')
    length = float(depth.max())
    cap_mask = np.max(np.abs(depth[mesh.faces]), axis=1) <= tolerance
    opposite = np.max(np.abs(depth[mesh.faces]-length), axis=1) <= tolerance
    if set(np.flatnonzero(cap_mask)) != set(map(int, ids)) or not np.any(opposite):
        raise ValueError('axial.end')
    opposite_center = np.average(triangles[opposite].mean(axis=1),
                                 weights=areas[opposite], axis=0)-origin
    if np.linalg.norm(opposite_center - (opposite_center @ axis)*axis) > 10*tolerance:
        raise ValueError('axial.centered')
    edges = np.linspace(0, length, int(count)+1)
    stations = (edges[:-1]+edges[1:])/2
    section_areas = []
    for i, nominal in enumerate(stations):
        # Avoid a slice through a mesh vertex; retain the actual station in the result.
        choices = nominal + np.array([0, 32, -32, 64, -64])*tolerance
        valid = [s for s in choices if edges[i] < s < edges[i+1]
                 and np.min(np.abs(depth-s)) > 4*tolerance]
        if not valid:
            raise ValueError('axial.section')
        stations[i] = valid[0]
        area, centroid = _section(mesh.vertices, mesh.faces, axis, origin, stations[i], tolerance)
        if np.linalg.norm(centroid) > 10*tolerance:
            raise ValueError('axial.centered')
        section_areas.append(area)
    face_depth = depth[mesh.faces].mean(axis=1)
    face_sections = np.clip(np.searchsorted(edges, face_depth, side='right')-1, 0, int(count)-1)
    return AxialSections(mesh_id(mesh), origin, axis, edges, stations,
                         np.array(section_areas), float(areas[ids].sum()), ids.copy(),
                         np.flatnonzero(opposite), face_sections)


def section_groups(factors):
    """Reuse numerically equal stresses only; never interpolate probability.

    A 1e-10 relative tolerance absorbs plane-intersection roundoff. Record both
    the exact geometrical factors and the factors actually used for the PDE.
    """
    unique, membership = [], []
    for factor in factors:
        match = next((i for i, value in enumerate(unique)
                      if abs(value-factor) <= 1e-10*abs(factor)), None)
        if match is None:
            match = len(unique); unique.append(float(factor))
        membership.append(match)
    return np.array(unique), np.array(membership, dtype=int)


def run_axial_probability(config, sections, *, nominal_result=None,
                          stop_requested=None, progress=None):
    """Solve each distinct section stress with the unchanged canonical PDE.

    This is one-way force-controlled coupling. Cancellation discards the spatial
    result, so a partially computed specimen is never advertised as complete.
    """
    from solver_v1.kinetic_calibration_workflow import build_time_basis_model
    config.require_material_backend()
    config.validate()
    probability_model = build_time_basis_model(config.energy_model, time_basis=config.time_basis,
                                              calibration=config.time_calibration)
    unique, membership = section_groups(sections.stress_factors)
    histories, floors, residuals, raw_histories = [], [], [], []
    end = config.cycles*config.model_period
    interval = config.model_period/config.steps_per_cycle
    times = np.arange(int(np.floor(config.cycles*config.steps_per_cycle))+1)*interval
    if end-times[-1] > 1e-12*config.model_period:
        times = np.append(times, end)
    else:
        times[-1] = end
    stopped = stop_requested or (lambda: False)
    for i, factor in enumerate(unique):
        if stopped():
            raise InterruptedError('axial.cancelled')
        if progress:
            progress(i+1, len(unique))
        local = replace(config, stress_mean_mpa=config.stress_mean_mpa*factor,
                        stress_amplitude_mpa=config.stress_amplitude_mpa*factor)
        # The caller supplies the just-computed same-config nominal run only.
        if nominal_result is not None and factor == 1.:
            result = nominal_result
        else:
            result = run_ui_analysis(local, stop_requested=stopped, _prepared_model=probability_model)
        if stopped():
            raise InterruptedError('axial.cancelled')
        actual_time = np.asarray(result['model_time'])
        extinct = float(np.asarray(result['intact_probability_mass'])[-1]) == 0.
        if abs(actual_time[-1]-config.cycles*config.model_period) > 1e-9*config.model_period and not extinct:
            raise ValueError('axial.incomplete')
        raw_p = np.asarray(result['local_initiation_probability'])
        raw_floor = np.broadcast_to(result['local_rare_event_floor'], actual_time.shape)
        raw_residual = np.asarray(result['mass_balance_residual'])
        raw_histories.append(dict(model_time=actual_time, local_initiation_probability=raw_p,
            local_rare_event_floor=raw_floor, mass_balance_residual=raw_residual, stress_factor=factor,
            constant_extension_after_extinction=extinct and actual_time[-1] < end))
        # Each adaptive solver retains its own records. Interpolate in TIME only
        # for the common display slider; no interpolation across stress levels.
        histories.append(np.interp(times, actual_time, raw_p))
        floors.append(np.full_like(times, np.max(raw_floor)))
        residuals.append(np.interp(times, actual_time, raw_residual))
    p = np.array(histories).T[:, membership]
    stress = np.asarray(config.stress_mpa(times))[:, None]*sections.stress_factors
    force = np.asarray(config.stress_mpa(times))*sections.loaded_area_mm2
    used_factors = unique[membership]
    return dict(schema='aft.axial-sections/1', mechanics='quasistatic_section_average_N_over_A',
        coupling='one_way_force_controlled_local_reference_PDE',
        boundary_condition='uniform_axial_end_traction_opposite_reaction_lateral_free',
        mesh_sha256=sections.mesh_sha256, origin_mm=sections.origin_mm, axis=sections.axis,
        edges_mm=sections.edges_mm, positions_mm=sections.positions_mm, areas_mm2=sections.areas_mm2,
        loaded_area_mm2=sections.loaded_area_mm2, loaded_faces=sections.loaded_faces,
        opposite_faces=sections.opposite_faces, face_sections=sections.face_sections,
        stress_factors=sections.stress_factors, pde_stress_factors=used_factors,
        distinct_pde_runs=len(unique), stress_group_relative_tolerance=1e-10,
        display_time_sampling='linear_interpolation_of_each_actual_PDE_history',
        raw_histories=raw_histories, section_groups=membership,
        maximum_raw_mass_residual=float(max(np.max(np.abs(r['mass_balance_residual'])) for r in raw_histories)),
        model_time=times, stress_mpa=stress, axial_force_n=force,
        force_balance_residual_n=float(np.max(np.abs(stress*sections.areas_mm2-force[:, None]))),
        local_initiation_probability=p, local_survival_probability=1-p,
        local_rare_event_floor=np.array(floors).T[:, membership],
        mass_balance_residual=np.array(residuals).T[:, membership],
        probability_resolution_certified=False, spatial_resolution_certified=False,
        material_calibrated=False, three_dimensional_stress_concentration=False,
        energy_model=config.energy_model, analysis_quality=config.analysis_quality,
        initialization=config.initialization,
        integration_method=config.integration_method, grid_shape=(config.grid_n_a, config.grid_n_s),
        time_basis=config.time_basis, nominal_stress_mean_mpa=config.stress_mean_mpa,
        nominal_stress_amplitude_mpa=config.stress_amplitude_mpa,
        units=dict(position='mm', area='mm^2', stress='MPa', force='N', probability='1'))


def validate_spatial_result(result, mesh=None):
    """Validate persisted/displayed field dimensions and mesh provenance."""
    if not isinstance(result, dict) or result.get('schema') != 'aft.axial-sections/1':
        raise ValueError('axial.incomplete')
    times = np.asarray(result['model_time'])
    areas = np.asarray(result['areas_mm2'])
    if (times.ndim != 1 or not times.size or areas.ndim != 1 or not 2 <= areas.size <= 64
            or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)
            or not np.isfinite(areas).all() or np.any(areas <= 0)):
        raise ValueError('axial.incomplete')
    for name in ('local_initiation_probability', 'local_survival_probability',
                 'local_rare_event_floor', 'mass_balance_residual', 'stress_mpa'):
        values = np.asarray(result[name])
        if values.shape != (times.size, areas.size) or not np.isfinite(values).all():
            raise ValueError('axial.incomplete')
    p = np.asarray(result['local_initiation_probability'])
    if (np.any((p < 0) | (p > 1)) or np.any(np.diff(p, axis=0) < 0)
            or not np.array_equal(np.asarray(result['local_survival_probability']), 1-p)):
        raise ValueError('axial.incomplete')
    ids = np.asarray(result['face_sections'])
    if ids.ndim != 1 or not ids.size or ids.dtype.kind not in 'iu' or ids.min() < 0 or ids.max() >= areas.size:
        raise ValueError('axial.incomplete')
    if mesh is not None and (result['mesh_sha256'] != mesh_id(mesh) or len(ids) != len(mesh.faces)):
        raise ValueError('axial.stale')
