"""Bounded, mm-based surface geometry; not a volume mechanics solver."""
from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np

MAX_FACES = 20000


@dataclass(frozen=True)
class SurfaceMesh:
    vertices: np.ndarray
    faces: np.ndarray
    source: str

    def __post_init__(self):
        v = np.asarray(self.vertices, dtype=float)
        f = np.asarray(self.faces, dtype=int)
        if (v.ndim != 2 or v.shape[1] != 3 or not np.isfinite(v).all()
                or f.ndim != 2 or f.shape[1] != 3 or not 0 < len(f) <= MAX_FACES
                or f.min() < 0 or f.max() >= len(v)):
            raise ValueError('invalid_mesh')
        v = v.copy(); f = f.copy()
        v.setflags(write=False); f.setflags(write=False)
        object.__setattr__(self, 'vertices', v)
        object.__setattr__(self, 'faces', f)
        if not np.isfinite(self.areas).all() or np.any(self.areas <= 0):
            raise ValueError('invalid_mesh')

    @property
    def areas(self):
        t = self.vertices[self.faces]
        return .5 * np.linalg.norm(np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), axis=1)

    @property
    def closed(self):
        edges = np.sort(np.concatenate([self.faces[:, [0, 1]], self.faces[:, [1, 2]],
                                        self.faces[:, [2, 0]]]), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        return bool(np.all(counts == 2))  # edge incidence, not self-intersection certification

    @property
    def enclosed_volume_mm3(self):
        """Oriented closed-surface volume; no self-intersection certification."""
        if not self.closed: raise ValueError('volume_requires_closed')
        edges = np.concatenate([self.faces[:, [0,1]], self.faces[:, [1,2]], self.faces[:, [2,0]]])
        _, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        signed = np.bincount(inverse, weights=np.where(edges[:, 0] < edges[:, 1], 1., -1.))
        if np.any(signed != 0): raise ValueError('volume_requires_oriented')
        # Translate near the origin to avoid cancellation for translated CAD.
        tri = (self.vertices-self.vertices.mean(axis=0))[self.faces]
        volume = abs(float(np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum()/6))
        if volume <= 0 or not np.isfinite(volume): raise ValueError('volume_requires_closed')
        return volume


def refine_selected(mesh, selected):
    """Split selected edges and triangulate adjacent polygons conformingly.

    Surface facets retain their geometry; this is not curved CAD reconstruction.
    Return child triangles of selected parents for continued interactive editing.
    """
    selected = set(map(int, selected))
    if not selected or min(selected) < 0 or max(selected) >= len(mesh.faces):
        raise ValueError('invalid_selection')
    marked = set()
    for i in selected:
        a,b,c = map(int, mesh.faces[i])
        marked.update(tuple(sorted(edge)) for edge in ((a,b),(b,c),(c,a)))
    vertices = mesh.vertices.tolist()
    midpoints = {}
    for edge in sorted(marked):
        midpoints[edge] = len(vertices)
        vertices.append(mesh.vertices[list(edge)].mean(axis=0).tolist())
    faces, children = [], []
    for i, (a,b,c) in enumerate(mesh.faces):
        polygon = []
        split = False
        for v,w in ((a,b),(b,c),(c,a)):
            polygon.append(int(v))
            mid = midpoints.get(tuple(sorted((int(v),int(w)))))
            if mid is not None: polygon.append(mid); split = True
        start = len(faces)
        if split:
            center = len(vertices)
            vertices.append(mesh.vertices[[a,b,c]].mean(axis=0).tolist())
            faces.extend([[center, v, polygon[(k+1)%len(polygon)]] for k,v in enumerate(polygon)])
        else: faces.append([a,b,c])
        if i in selected: children.extend(range(start,len(faces)))
        if len(faces) > MAX_FACES: raise ValueError('mesh_limit')
    return SurfaceMesh(vertices, faces, mesh.source), np.asarray(children,dtype=int)


def cylinder(radius_mm=5., length_mm=30., target_mm=3.):
    if not all(np.isfinite(x) and x > 0 for x in (radius_mm, length_mm, target_mm)):
        raise ValueError('positive_dimensions')
    if radius_mm/target_mm > MAX_FACES/(2*np.pi) or length_mm/target_mm > MAX_FACES:
        raise ValueError('mesh_limit')
    n = max(12, int(np.ceil(2*np.pi*radius_mm/target_mm)))
    nz = max(1, int(np.ceil(length_mm/target_mm)))
    if 2*n*(nz+1) > MAX_FACES:
        raise ValueError('mesh_limit')
    angle = np.arange(n)*2*np.pi/n
    v = [[radius_mm*np.cos(t), radius_mm*np.sin(t), z]
         for z in np.linspace(0, length_mm, nz+1) for t in angle]
    v.extend([[0, 0, 0], [0, 0, length_mm]])
    bottom, top = len(v)-2, len(v)-1
    faces = []
    for j in range(nz):
        for i in range(n):
            a, b = j*n+i, j*n+(i+1)%n
            faces.extend([[a, b, b+n], [a, b+n, a+n]])
    for i in range(n):
        faces.extend([[bottom, (i+1)%n, i], [top, nz*n+i, nz*n+(i+1)%n]])
    return SurfaceMesh(v, faces, 'cylinder')


def load_surface(path, millimeters_per_unit=1.):
    """Import triangulated STL or OBJ. Source units must be explicitly selected."""
    path = Path(path)
    if not np.isfinite(millimeters_per_unit) or millimeters_per_unit <= 0:
        raise ValueError('positive_dimensions')
    if path.stat().st_size > 10_000_000:
        raise ValueError('mesh_limit')
    data = path.read_bytes()
    if path.suffix.lower() == '.stl':
        count = struct.unpack_from('<I', data, 80)[0] if len(data) >= 84 else -1
        if len(data) == 84+50*count:
            if count > MAX_FACES:
                raise ValueError('mesh_limit')
            vertices = [struct.unpack_from('<9f', data, 84+50*i+12) for i in range(count)]
            triangles = np.asarray(vertices).reshape(-1, 3)
        else:
            rows = [fields[1:] for line in data.decode('ascii').splitlines()
                    if (fields := line.split()) and fields[0] == 'vertex']
            if any(len(row) != 3 for row in rows):
                raise ValueError('invalid_mesh')
            triangles = np.asarray(rows, dtype=float).reshape(-1, 3)
        if not len(triangles) or len(triangles) % 3:
            raise ValueError('invalid_mesh')
        vertices, indices = np.unique(triangles, axis=0, return_inverse=True)
        faces = indices.reshape(-1, 3)
    elif path.suffix.lower() == '.obj':
        vertices, faces = [], []
        for line in data.decode('utf-8-sig').splitlines():
            fields = line.split('#', 1)[0].split()
            if not fields:
                continue
            if fields[0] == 'v':
                vertices.append([float(x) for x in fields[1:4]])
            elif fields[0] == 'f':
                if len(fields) != 4:
                    raise ValueError('triangles_only')
                ids = [int(x.split('/')[0]) for x in fields[1:]]
                if 0 in ids:
                    raise ValueError('invalid_mesh')
                faces.append([x-1 if x > 0 else len(vertices)+x for x in ids])
    else:
        raise ValueError('supported_formats')
    return SurfaceMesh(np.asarray(vertices)*millimeters_per_unit, faces, path.name)


def refine_surface(mesh, target_mm):
    """Conforming midpoint subdivision; no smoothing or geometric repair."""
    if not np.isfinite(target_mm) or target_mm <= 0:
        raise ValueError('positive_dimensions')
    result = mesh
    while True:
        t = result.vertices[result.faces]
        longest = max(np.linalg.norm(t[:, i]-t[:, (i+1)%3], axis=1).max() for i in range(3))
        if longest <= target_mm:
            return result
        if len(result.faces)*4 > MAX_FACES:
            raise ValueError('mesh_limit')
        vertices, faces, cache = result.vertices.tolist(), [], {}
        def midpoint(a, b):
            key = tuple(sorted((a, b)))
            if key not in cache:
                cache[key] = len(vertices)
                vertices.append(((result.vertices[a]+result.vertices[b])/2).tolist())
            return cache[key]
        for a, b, c in result.faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            faces.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]])
        result = SurfaceMesh(vertices, faces, mesh.source)
