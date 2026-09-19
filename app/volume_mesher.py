"""Validate TetGen output before FEM; preserve assigned surface facets on retry."""
import numpy as np

MAX_SOLID_NODES = 1_200_000
MAX_SOLID_CELLS = 6_000_000


class SolidMeshLimitError(ValueError):
    def __init__(self, nodes, cells, target_mm, *, max_nodes=None, max_cells=None):
        super().__init__('solid.size')
        self.ui_error_data = dict(key='solid.size_detail', values=dict(
            nodes=int(nodes), cells=int(cells), target=f'{target_mm:.6g}',
            max_nodes=MAX_SOLID_NODES if max_nodes is None else max_nodes,
            max_cells=MAX_SOLID_CELLS if max_cells is None else max_cells))


class SolidMeshError(ValueError):
    def __init__(self, reason, detail=''):
        super().__init__('solid.mesh')
        self.reason = reason
        self.ui_error_data = dict(key='solid.mesh_'+reason, values=dict(detail=detail))


def validate_tetrahedra(mesh, mesher, nodes, cells, volume):
    nodes, cells = np.asarray(nodes), np.asarray(cells)
    if (nodes.ndim != 2 or nodes.shape[1] != 3 or not np.isfinite(nodes).all()
            or cells.ndim != 2 or cells.shape[1] != 4 or not len(cells)
            or cells.dtype.kind not in 'iu' or cells.min() < 0 or cells.max() >= len(nodes)):
        raise SolidMeshError('degenerate', 'invalid node/cell arrays')
    t = nodes[cells]
    edges = t[:, 1:]-t[:, :1]
    determinant = np.linalg.det(edges)
    # Refuse elements whose signed volume is unresolved at float precision.
    arithmetic_scale = np.prod(np.linalg.norm(edges, axis=2), axis=1)
    bad = ~np.isfinite(determinant) | (np.abs(determinant) <= 64*np.finfo(float).eps*arithmetic_scale)
    if np.any(bad):
        raise SolidMeshError('degenerate', f'{np.count_nonzero(bad):,} / {len(cells):,}')
    volumes = np.abs(determinant)/6
    if not np.isclose(volumes.sum(), volume, rtol=1e-8):
        raise SolidMeshError('volume', f'{volumes.sum():.12g} / {volume:.12g} mm³')
    marker = np.asarray(mesher.triface_markers)
    on_surface = marker > 0
    boundary = np.asarray(mesher.trifaces)[on_surface].astype(np.int32)
    parents = marker[on_surface].astype(np.int32)-1
    if (not len(boundary) or boundary.shape[1] != 3 or np.any(parents >= len(mesh.faces))
            or boundary.min() < 0 or boundary.max() >= len(nodes)):
        raise SolidMeshError('boundary', 'invalid facet markers')
    tri = nodes[boundary]
    areas = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)/2
    mapped = np.bincount(parents, weights=areas, minlength=len(mesh.faces))
    match = np.isclose(mapped, mesh.areas, rtol=1e-8, atol=0)
    if not np.all(match):
        raise SolidMeshError('boundary', f'{np.count_nonzero(~match):,} / {len(mesh.faces):,}')
    adjacent = np.asarray(mesher.face2tet)[on_surface]
    if adjacent.shape != (len(boundary), 2) or not np.all(np.sum(adjacent >= 0, axis=1) == 1):
        raise SolidMeshError('boundary', 'invalid boundary adjacency')
    owner = adjacent.max(axis=1)
    for offset in (0, 1):
        candidate = owner-offset
        if np.any(candidate < 0) or np.any(candidate >= len(cells)): continue
        belongs = np.all(np.any(boundary[:, :, None] == cells[candidate][:, None, :], axis=2), axis=1)
        if np.all(belongs):
            return dict(nodes=nodes, cells=cells.astype(np.int32), volumes=volumes,
                boundary_faces=boundary, parents=parents, boundary_areas=areas,
                boundary_cells=candidate.astype(np.int32))
    raise SolidMeshError('boundary', 'facet owner does not contain its vertices')


def tetrahedralize(mesh, target_mm, origin, *, max_nodes=MAX_SOLID_NODES,
                   max_cells=MAX_SOLID_CELLS, progress=None, stop_requested=lambda: False):
    try:
        import tetgen
    except ImportError as exc:
        raise ValueError('solid.dependency') from exc
    volume = mesh.enclosed_volume_mm3
    if len(mesh.vertices) >= max_nodes:
        raise SolidMeshLimitError(len(mesh.vertices), 0, target_mm, max_nodes=max_nodes, max_cells=max_cells)
    attempts = []
    for preserve_surface in (False, True):
        nodes = cells = None
        if stop_requested(): raise InterruptedError('axial.cancelled')
        if preserve_surface and progress: progress('solid.mesh_retry', {})
        mesher = tetgen.TetGen((mesh.vertices-origin).copy(), mesh.faces.astype(np.int32),
                              np.arange(1, len(mesh.faces)+1, dtype=np.int32))
        try:
            nodes, cells, _, _ = mesher.tetrahedralize(plc=True, quality=True,
                fixedvolume=True, maxvolume=target_mm**3/6, facesout=True, neighout=True,
                quiet=True, steinerleft=max_nodes-len(mesh.vertices), nobisect=preserve_surface)
            if len(nodes) >= max_nodes or len(cells) > max_cells:
                raise SolidMeshLimitError(len(nodes), len(cells), target_mm,
                                          max_nodes=max_nodes, max_cells=max_cells)
            result = validate_tetrahedra(mesh, mesher, nodes, cells, volume)
        except (RuntimeError, SolidMeshError, SolidMeshLimitError) as exc:
            attempts.append(dict(preserve_surface=preserve_surface,
                                 reason=getattr(exc, 'reason', str(exc))))
            nodes = cells = None
            del mesher
            if preserve_surface:
                if isinstance(exc, RuntimeError): raise SolidMeshError('engine', str(exc)) from exc
                raise
            continue
        result['meshing_info'] = dict(
            mode='preserved_input_facets' if preserve_surface else 'surface_refinement',
            failed_attempts=attempts, input_surface_faces=len(mesh.faces),
            target_volume_mm3=target_mm**3/6,
            maximum_cell_volume_mm3=float(result['volumes'].max()),
            minimum_cell_volume_mm3=float(result['volumes'].min()),
            input_geometry_and_face_assignments_preserved=True,
            spatial_convergence_certified=False)
        return result
