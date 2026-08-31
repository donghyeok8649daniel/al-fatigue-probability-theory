"""Dependency-light CAD-style viewer for 1D/2D/3D mesh geometry.

Supported imports are mesh/interchange formats that can be parsed safely with
the repository's existing dependencies: OBJ, STL (binary or ASCII), ASCII PLY
and legacy ASCII VTK. STEP/IGES are intentionally not advertised because no
CAD-kernel backend is installed.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import struct
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.widgets import Button, RadioButtons
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import numpy as np


MAX_FILE_BYTES = 256 * 1024 * 1024
MAX_VERTICES = 5_000_000
SUPPORTED_EXTENSIONS = {".obj", ".stl", ".ply", ".vtk"}


@dataclass(frozen=True)
class MeshGeometry:
    vertices: np.ndarray
    faces: tuple[tuple[int, ...], ...] = ()
    lines: tuple[tuple[int, ...], ...] = ()
    source: str = ""

    def __post_init__(self) -> None:
        vertices = np.asarray(self.vertices, dtype=float)
        if vertices.ndim != 2 or vertices.shape[1] not in {1, 2, 3} or vertices.shape[0] == 0:
            raise ValueError("vertices must be a non-empty N-by-1/2/3 array")
        if vertices.shape[0] > MAX_VERTICES or not np.all(np.isfinite(vertices)):
            raise ValueError("mesh contains too many or non-finite vertices")
        if vertices.shape[1] < 3:
            vertices = np.pad(vertices, ((0, 0), (0, 3 - vertices.shape[1])))
        object.__setattr__(self, "vertices", vertices)
        count = vertices.shape[0]
        for cell in self.faces + self.lines:
            if len(cell) < 2 or any(index < 0 or index >= count for index in cell):
                raise ValueError("mesh connectivity references an invalid vertex")

    @property
    def dimension(self) -> int:
        centered = self.vertices - np.mean(self.vertices, axis=0)
        scale = max(float(np.linalg.norm(np.ptp(self.vertices, axis=0))), 1.0)
        rank = int(np.linalg.matrix_rank(centered, tol=scale * 1e-10))
        return max(1, min(rank, 3))


def _check_source(path: Path) -> Path:
    source = Path(path)
    if source.suffix.lower() not in SUPPORTED_EXTENSIONS:
        raise ValueError("supported mesh formats are OBJ, STL, PLY and legacy VTK")
    if not source.is_file() or source.stat().st_size > MAX_FILE_BYTES:
        raise ValueError("mesh file is missing or exceeds the 256 MiB safety limit")
    return source


def _mesh(vertices, faces=(), lines=(), source="") -> MeshGeometry:
    return MeshGeometry(np.asarray(vertices, dtype=float), tuple(map(tuple, faces)),
                        tuple(map(tuple, lines)), source)


def _load_obj(path: Path) -> MeshGeometry:
    vertices, faces, lines = [], [], []
    for raw in path.read_text(encoding="utf-8", errors="strict").splitlines():
        parts = raw.strip().split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "v" and len(parts) >= 4:
            vertices.append(tuple(map(float, parts[1:4])))
        elif parts[0] in {"f", "l"}:
            indices = []
            for item in parts[1:]:
                value = int(item.split("/", 1)[0])
                index = value - 1 if value > 0 else len(vertices) + value
                indices.append(index)
            (faces if parts[0] == "f" else lines).append(tuple(indices))
    return _mesh(vertices, faces, lines, str(path))


def _load_stl(path: Path) -> MeshGeometry:
    data = path.read_bytes()
    binary = False
    if len(data) >= 84:
        triangles = struct.unpack_from("<I", data, 80)[0]
        binary = len(data) == 84 + 50 * triangles
    triangles_xyz = []
    if binary:
        for i in range(triangles):
            offset = 84 + 50 * i + 12
            triangles_xyz.append([struct.unpack_from("<3f", data, offset + 12 * j) for j in range(3)])
    else:
        text = data.decode("utf-8", errors="strict")
        current = []
        for line in text.splitlines():
            parts = line.strip().split()
            if parts and parts[0].lower() == "vertex" and len(parts) == 4:
                current.append(tuple(map(float, parts[1:])))
                if len(current) == 3:
                    triangles_xyz.append(current); current = []
        if current or not triangles_xyz:
            raise ValueError("invalid ASCII STL triangles")
    vertices, lookup, faces = [], {}, []
    for triangle in triangles_xyz:
        face = []
        for point in triangle:
            key = tuple(float(value) for value in point)
            if key not in lookup:
                lookup[key] = len(vertices); vertices.append(key)
            face.append(lookup[key])
        faces.append(tuple(face))
    return _mesh(vertices, faces, source=str(path))


def _load_ply(path: Path) -> MeshGeometry:
    with path.open("r", encoding="utf-8", errors="strict") as stream:
        if stream.readline().strip() != "ply":
            raise ValueError("invalid PLY header")
        vertex_count = face_count = 0
        while True:
            line = stream.readline()
            if not line:
                raise ValueError("unterminated PLY header")
            parts = line.strip().split()
            if parts[:2] == ["format", "ascii"]:
                pass
            elif parts and parts[0] == "format":
                raise ValueError("only ASCII PLY is supported without an external backend")
            elif parts[:2] == ["element", "vertex"]:
                vertex_count = int(parts[2])
            elif parts[:2] == ["element", "face"]:
                face_count = int(parts[2])
            elif parts and parts[0] == "end_header":
                break
        if not 0 < vertex_count <= MAX_VERTICES:
            raise ValueError("invalid PLY vertex count")
        vertices = [tuple(map(float, stream.readline().split()[:3])) for _ in range(vertex_count)]
        faces = []
        for _ in range(face_count):
            values = list(map(int, stream.readline().split()))
            faces.append(tuple(values[1:1 + values[0]]))
    return _mesh(vertices, faces, source=str(path))


def _load_vtk(path: Path) -> MeshGeometry:
    text = path.read_text(encoding="utf-8", errors="strict")
    if "BINARY" in text[:500].upper():
        raise ValueError("only legacy ASCII VTK is supported without an external backend")
    tokens = text.replace("\r", " ").split()
    try:
        point_at = next(i for i, token in enumerate(tokens) if token.upper() == "POINTS")
        count = int(tokens[point_at + 1]); cursor = point_at + 3
        vertices = np.asarray(list(map(float, tokens[cursor:cursor + 3 * count]))).reshape(count, 3)
    except (StopIteration, ValueError, IndexError) as exc:
        raise ValueError("invalid legacy VTK POINTS section") from exc
    faces, lines = [], []
    upper = [token.upper() for token in tokens]
    for keyword, target in (("POLYGONS", faces), ("LINES", lines)):
        if keyword not in upper:
            continue
        at = upper.index(keyword); cell_count = int(tokens[at + 1]); cursor = at + 3
        for _ in range(cell_count):
            width = int(tokens[cursor]); cursor += 1
            target.append(tuple(map(int, tokens[cursor:cursor + width]))); cursor += width
    return _mesh(vertices, faces, lines, str(path))


def load_mesh(path: Path) -> MeshGeometry:
    source = _check_source(path)
    return {".obj": _load_obj, ".stl": _load_stl, ".ply": _load_ply,
            ".vtk": _load_vtk}[source.suffix.lower()](source)


def principal_coordinates(mesh: MeshGeometry, dimension: int) -> tuple[np.ndarray, np.ndarray]:
    """Project geometry to its best-fit 1D/2D basis for display only."""
    if dimension not in {1, 2, 3}:
        raise ValueError("display dimension must be 1, 2 or 3")
    centered = mesh.vertices - np.mean(mesh.vertices, axis=0)
    _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
    basis = np.eye(3) if dimension == 3 else vh[:dimension].T
    return centered @ basis, basis


def _segments(mesh: MeshGeometry) -> tuple[tuple[int, int], ...]:
    result: set[tuple[int, int]] = set()
    cells = mesh.lines + mesh.faces
    for cell in cells:
        pairs = zip(cell[:-1], cell[1:])
        for a, b in pairs:
            result.add(tuple(sorted((a, b))))
        if cell in mesh.faces and len(cell) > 2:
            result.add(tuple(sorted((cell[-1], cell[0]))))
    if not result and mesh.vertices.shape[0] > 1:
        result.update((i, i + 1) for i in range(mesh.vertices.shape[0] - 1))
    return tuple(sorted(result))


class CADNavigation:
    """Wheel zoom, middle-drag pan and reset around one Matplotlib axes."""
    def __init__(self, figure, axes, dimension: int):
        self.figure, self.axes, self.dimension = figure, axes, dimension
        self.home = (axes.get_xlim(), axes.get_ylim(), axes.get_zlim() if dimension == 3 else None)
        self._press = None
        if dimension == 3 and hasattr(axes, "mouse_init"):
            axes.mouse_init(rotate_btn=1, pan_btn=2, zoom_btn=3)
        self.ids = [figure.canvas.mpl_connect("scroll_event", self.on_scroll),
                    figure.canvas.mpl_connect("button_press_event", self.on_press),
                    figure.canvas.mpl_connect("button_release_event", self.on_release),
                    figure.canvas.mpl_connect("motion_notify_event", self.on_motion)]

    @staticmethod
    def _scaled(limits, center, factor):
        return center + (np.asarray(limits, dtype=float) - center) * factor

    def on_scroll(self, event) -> None:
        if event.inaxes is not self.axes:
            return
        factor = 0.82 if event.button == "up" else 1.0 / 0.82
        if self.dimension == 3:
            for getter, setter in ((self.axes.get_xlim3d, self.axes.set_xlim3d),
                                   (self.axes.get_ylim3d, self.axes.set_ylim3d),
                                   (self.axes.get_zlim3d, self.axes.set_zlim3d)):
                limits = getter(); setter(self._scaled(limits, np.mean(limits), factor))
        else:
            xlim, ylim = self.axes.get_xlim(), self.axes.get_ylim()
            xcenter = event.xdata if event.xdata is not None else np.mean(xlim)
            ycenter = event.ydata if event.ydata is not None else np.mean(ylim)
            self.axes.set_xlim(self._scaled(xlim, xcenter, factor))
            self.axes.set_ylim(self._scaled(ylim, ycenter, factor))
        self.figure.canvas.draw_idle()

    def on_press(self, event) -> None:
        if self.dimension < 3 and event.inaxes is self.axes and event.button in {2, 3}:
            self._press = (event.xdata, event.ydata, self.axes.get_xlim(), self.axes.get_ylim())

    def on_release(self, _event) -> None:
        self._press = None

    def on_motion(self, event) -> None:
        if self._press is None or event.inaxes is not self.axes or event.xdata is None or event.ydata is None:
            return
        x0, y0, xlim, ylim = self._press
        self.axes.set_xlim(np.asarray(xlim) + x0 - event.xdata)
        self.axes.set_ylim(np.asarray(ylim) + y0 - event.ydata)
        self.figure.canvas.draw_idle()

    def reset(self, _event=None) -> None:
        self.axes.set_xlim(self.home[0]); self.axes.set_ylim(self.home[1])
        if self.dimension == 3:
            self.axes.set_zlim(self.home[2]); self.axes.view_init(elev=25, azim=-60)
        self.figure.canvas.draw_idle()

    def disconnect(self) -> None:
        for connection in self.ids:
            self.figure.canvas.mpl_disconnect(connection)


class MeshViewport:
    def __init__(self, mesh: MeshGeometry, display_dimension: int | None = None,
                 loading_axis=(1.0, 0.0, 0.0)):
        self.mesh = mesh
        self.dimension = mesh.dimension if display_dimension is None else display_dimension
        self.loading_axis = np.asarray(loading_axis, dtype=float)
        if self.dimension not in {1, 2, 3} or np.linalg.norm(self.loading_axis) == 0:
            raise ValueError("invalid display dimension or loading axis")
        self.figure = plt.figure(figsize=(10, 7))
        self.figure.canvas.manager.set_window_title(f"FTGSim geometry: {Path(mesh.source).name}")
        self.axes = self.figure.add_axes([0.15, 0.12, 0.77, 0.80],
                                         projection="3d" if self.dimension == 3 else None)
        self._draw()
        reset_ax = self.figure.add_axes([0.015, 0.91, 0.08, 0.045])
        self.reset_button = Button(reset_ax, "Reset")
        self.navigation = CADNavigation(self.figure, self.axes, self.dimension)
        self.reset_button.on_clicked(self._on_reset)
        dimension_ax = self.figure.add_axes([0.015, 0.70, 0.09, 0.15])
        self.dimension_radio = RadioButtons(dimension_ax, ("1D", "2D", "3D"),
            active=self.dimension - 1)
        self.dimension_radio.on_clicked(self._on_dimension)
        self.figure.text(0.12, 0.965,
            "Left drag: 3D orbit | Middle drag: pan | Right drag: 3D zoom | Wheel: zoom",
            fontsize=9, va="top")

    def _on_reset(self, _event=None) -> None:
        self.navigation.reset()

    def _on_dimension(self, label: str) -> None:
        dimension = int(label[0])
        if dimension == self.dimension:
            return
        self.navigation.disconnect()
        self.figure.delaxes(self.axes)
        self.dimension = dimension
        self.axes = self.figure.add_axes([0.15, 0.12, 0.77, 0.80],
                                         projection="3d" if dimension == 3 else None)
        self._draw()
        self.navigation = CADNavigation(self.figure, self.axes, dimension)
        self.figure.canvas.draw_idle()

    def _draw(self) -> None:
        coordinates, basis = principal_coordinates(self.mesh, self.dimension)
        edges = _segments(self.mesh)
        if self.dimension == 3:
            if self.mesh.faces:
                polygons = [[coordinates[index] for index in face] for face in self.mesh.faces]
                self.axes.add_collection3d(Poly3DCollection(polygons, facecolor="#78a9d1",
                    edgecolor="#23384d", linewidth=0.45, alpha=0.78))
            if edges:
                self.axes.add_collection3d(Line3DCollection(
                    [[coordinates[a], coordinates[b]] for a, b in edges], colors="#263238", linewidths=0.7))
            self.axes.scatter(*coordinates.T, s=4, color="#17202a", alpha=0.55)
        elif self.dimension == 2:
            if self.mesh.faces:
                self.axes.add_collection(PolyCollection(
                    [[coordinates[index] for index in face] for face in self.mesh.faces],
                    facecolor="#78a9d1", edgecolor="#23384d", linewidth=0.5, alpha=0.8))
            self.axes.add_collection(LineCollection(
                [[coordinates[a], coordinates[b]] for a, b in edges], colors="#263238", linewidths=0.8))
            self.axes.scatter(coordinates[:, 0], coordinates[:, 1], s=5, color="#17202a")
        else:
            x = coordinates[:, 0]; self.axes.plot(x, np.zeros_like(x), "o", ms=3)
            self.axes.add_collection(LineCollection(
                [[(x[a], 0), (x[b], 0)] for a, b in edges], colors="#263238", linewidths=1.2))
        self._fit(coordinates)
        self._draw_loading_axis(coordinates, basis)
        self.axes.set_title(f"{Path(self.mesh.source).name or 'geometry'} | inferred {self.mesh.dimension}D, shown {self.dimension}D")
        self.axes.set_xlabel("x"); self.axes.set_ylabel("y")
        if self.dimension == 3:
            self.axes.set_zlabel("z"); self.axes.set_box_aspect(np.maximum(np.ptp(coordinates, axis=0), 1e-9))
            self.axes.view_init(elev=25, azim=-60)
        else:
            self.axes.set_aspect("equal", adjustable="datalim")

    def _fit(self, coordinates) -> None:
        low, high = np.min(coordinates, axis=0), np.max(coordinates, axis=0)
        span = np.maximum(high - low, max(float(np.max(high - low)), 1.0) * 0.04)
        low -= 0.08 * span; high += 0.08 * span
        self.axes.set_xlim(low[0], high[0])
        if self.dimension >= 2: self.axes.set_ylim(low[1], high[1])
        else: self.axes.set_ylim(-0.1 * span[0], 0.1 * span[0])
        if self.dimension == 3: self.axes.set_zlim(low[2], high[2])

    def _draw_loading_axis(self, coordinates, basis) -> None:
        origin = np.mean(coordinates, axis=0)
        projected = self.loading_axis if self.dimension == 3 else self.loading_axis @ basis
        projected = projected[:self.dimension]
        norm = np.linalg.norm(projected)
        if norm < 1e-12: return
        length = max(float(np.max(np.ptp(coordinates, axis=0))), 1.0) * 0.22
        vector = projected / norm * length
        if self.dimension == 3:
            self.axes.quiver(*origin, *vector, color="crimson", linewidth=2, arrow_length_ratio=0.18)
        else:
            origin_2d = np.asarray([origin[0], origin[1] if len(origin) > 1 else 0.0])
            vector_2d = np.pad(vector, (0, 2-len(vector)))[:2]
            self.axes.annotate("normal load", xy=origin_2d + vector_2d,
                xytext=origin_2d, arrowprops={"arrowstyle": "->", "color": "crimson", "lw": 2})

    def show(self) -> None:
        plt.show()


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="CAD-style viewer for mesh interchange files")
    parser.add_argument("geometry", type=Path)
    parser.add_argument("--dimension", type=int, choices=(1, 2, 3), default=None)
    parser.add_argument("--save-preview", type=Path, default=None)
    args = parser.parse_args(argv)
    viewport = MeshViewport(load_mesh(args.geometry), args.dimension)
    if args.save_preview:
        args.save_preview.parent.mkdir(parents=True, exist_ok=True)
        viewport.figure.savefig(args.save_preview, dpi=170)
    else:
        viewport.show()


if __name__ == "__main__":
    main()
