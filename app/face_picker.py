"""Visible triangle picking and connected planar CAD-like face selection."""
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d

EDGE_DISPLAY_FACE_LIMIT = 50_000

def ray_triangle_pick(vertices, faces, origin, direction):
    """Closest forward ray intersection, independent of triangle winding."""
    tri = vertices[faces]
    e1, e2 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
    p = np.cross(np.broadcast_to(direction, e2.shape), e2)
    det = np.einsum('ij,ij->i', e1, p)
    valid = np.abs(det) > 1e-12 * np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1)
    inv = np.divide(1., det, out=np.zeros_like(det), where=valid)
    offset = origin - tri[:, 0]
    u = np.einsum('ij,ij->i', offset, p) * inv
    q = np.cross(offset, e1)
    v = q @ direction * inv
    distance = np.einsum('ij,ij->i', e2, q) * inv
    valid &= (u >= -1e-10) & (v >= -1e-10) & (u+v <= 1+1e-10) & (distance >= 0)
    if not valid.any():
        return None
    return int(np.argmin(np.where(valid, distance, np.inf)))


def planar_patch(mesh, index):
    """Connected coplanar triangles; STL/OBJ have no original CAD face IDs."""
    tri = mesh.vertices[mesh.faces]
    normals = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    n = normals[index]
    scale = max(float(np.ptp(mesh.vertices, axis=0).max()), 1e-12)
    mask = (np.abs(normals @ n) > 1-1e-8) & (
        np.max(np.abs((tri-tri[index, 0]) @ n), axis=1) < scale*1e-8)
    edges = {}
    for i in np.flatnonzero(mask):
        f = mesh.faces[i]
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            edges.setdefault(tuple(sorted((int(a), int(b)))), []).append(int(i))
    visited, todo = set(), [int(index)]
    while todo:
        i = todo.pop()
        if i in visited: continue
        visited.add(i)
        f = mesh.faces[i]
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            todo.extend(j for j in edges[tuple(sorted((int(a), int(b))))] if j not in visited)
    return np.array(sorted(visited), dtype=int)


class FacePicker:
    def __init__(self, parent, callback):
        self.callback = callback
        self.mesh = None
        self.figure = Figure(figsize=(5, 5), dpi=90)
        self.ax = self.figure.add_subplot(111, projection='3d')
        self.ax.mouse_init(rotate_btn=3, pan_btn=2, zoom_btn=[])
        self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self.canvas.mpl_connect('button_press_event', self._click)
        self.canvas.mpl_connect('scroll_event', self._zoom)

    def draw(self, mesh, selected):
        changed = mesh is not self.mesh
        self.mesh = mesh
        if changed or mesh is None:
            for artist in list(self.ax.collections): artist.remove()
        if mesh is not None:
            colors = np.tile([.65, .78, .89, 1.], (len(mesh.faces), 1))
            colors[np.asarray(selected, dtype=int)] = [1., .48, .08, 1.]
            if changed or not self.ax.collections:
                # Retain every facet for selection/field display. Dense meshes
                # omit individual edge strokes, not geometry or stress values.
                edges = len(mesh.faces) <= EDGE_DISPLAY_FACE_LIMIT
                self.ax.add_collection3d(Poly3DCollection(mesh.vertices[mesh.faces],
                    facecolors=colors, edgecolors='#38566b' if edges else 'none',
                    linewidth=.3 if edges else 0.))
            else:
                self.ax.collections[0].set_facecolor(colors)
            if changed:
                lo, hi = mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)
                c, r = (lo+hi)/2, max(float((hi-lo).max())*.55, 1e-6)
                self.ax.set(xlim=(c[0]-r,c[0]+r), ylim=(c[1]-r,c[1]+r), zlim=(c[2]-r,c[2]+r))
        self.ax.set_box_aspect((1, 1, 1))
        self.ax.set(xlabel='X [mm]', ylabel='Y [mm]', zlabel='Z [mm]')
        self.canvas.draw_idle()

    def _click(self, event):
        if event.button != 1 or event.inaxes is not self.ax or self.mesh is None: return
        inv = np.linalg.inv(self.ax.get_proj())
        origin = np.asarray(proj3d.inv_transform(event.xdata, event.ydata, -2., inv)).reshape(3)
        far = np.asarray(proj3d.inv_transform(event.xdata, event.ydata, -0.5, inv)).reshape(3)
        direction = far-origin
        direction /= np.linalg.norm(direction)
        index = ray_triangle_pick(self.mesh.vertices, self.mesh.faces, origin, direction)
        if index is not None:
            self.callback(planar_patch(self.mesh, index), event.key in ('control', 'ctrl'))

    def _zoom(self, event):
        if event.inaxes is not self.ax: return
        factor = .85 if event.button == 'up' else 1/.85
        for get, set_ in ((self.ax.get_xlim3d, self.ax.set_xlim3d),
                          (self.ax.get_ylim3d, self.ax.set_ylim3d),
                          (self.ax.get_zlim3d, self.ax.set_zlim3d)):
            low, high = get(); c = (low+high)/2
            set_(c+(low-c)*factor, c+(high-c)*factor)
        self.canvas.draw_idle()
