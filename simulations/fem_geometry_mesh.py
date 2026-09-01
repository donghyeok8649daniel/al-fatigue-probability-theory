# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D normal-only 결과를 표시할 실제 2D/3D mesh를 생성하고 CAD/mesh 파일을 읽는 geometry 계층이다.
# - 주요 클래스: CellBlock, GeometryMesh, AxialProjection
# - 주요 함수/메서드: CellBlock.__post_init__, CellBlock.dimension, GeometryMesh.__post_init__
#   GeometryMesh.topological_dimension, GeometryMesh.cell_count, GeometryMesh.embedding_dimension
#   _point_index_2d, _point_index_3d, structured_rectangle_mesh, structured_box_mesh
#   _deduplicated_triangles, _load_stl, _load_obj, _load_with_meshio, _mesh_cad_with_gmsh
#   load_geometry_mesh, cell_centers, _unit_axis, project_normal_stress, map_axial_element_field
#   boundary_faces, save_mesh_npz, load_mesh_npz
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Geometry and meshing layer for a normal-stress-only fatigue workflow.

Mesh dimension and constitutive-model dimension are deliberately separated.
The mesh may be two- or three-dimensional, while the active probability model
receives only one scalar normal component, ``sigma_nn``.  Mapping a 1D bar
solution onto this mesh is a visualization/post-processing projection; it is
not a 2D/3D elasticity solution.
"""
from __future__ import annotations

from dataclasses import dataclass
import importlib
from pathlib import Path
import struct
from typing import Iterable

import numpy as np


_CELL_DIMENSION = {
    "triangle": 2,
    "quad": 2,
    "tetra": 3,
    "hexahedron": 3,
    "wedge": 3,
    "pyramid": 3,
}

_CELL_NODE_COUNT = {
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "hexahedron": 8,
    "wedge": 6,
    "pyramid": 5,
}

_GMSH_CELL_NAMES = {
    (2, 3): "triangle",
    (2, 4): "quad",
    (3, 4): "tetra",
    (3, 5): "pyramid",
    (3, 6): "wedge",
    (3, 8): "hexahedron",
}


@dataclass(frozen=True)
class CellBlock:
    """One homogeneous block of finite-element connectivity."""

    cell_type: str
    connectivity: np.ndarray

    def __post_init__(self) -> None:
        if self.cell_type not in _CELL_DIMENSION:
            raise ValueError(f"unsupported cell type: {self.cell_type}")
        cells = np.asarray(self.connectivity, dtype=int)
        expected = _CELL_NODE_COUNT[self.cell_type]
        if cells.ndim != 2 or cells.shape[1] != expected:
            raise ValueError(
                f"{self.cell_type} connectivity must have shape (n,{expected})"
            )
        if np.any(cells < 0):
            raise ValueError("cell connectivity must use nonnegative node indices")
        object.__setattr__(self, "connectivity", cells)

    @property
    def dimension(self) -> int:
        return _CELL_DIMENSION[self.cell_type]


@dataclass(frozen=True)
class GeometryMesh:
    """Minimal dependency-free unstructured mesh representation."""

    points_m: np.ndarray
    cell_blocks: tuple[CellBlock, ...]
    source: str = "generated"
    role: str = "geometry/display mesh; not a multidimensional mechanics solution"

    def __post_init__(self) -> None:
        points = np.asarray(self.points_m, dtype=float)
        if points.ndim != 2 or points.shape[1] not in (2, 3) or points.shape[0] == 0:
            raise ValueError("points_m must have shape (n,2) or (n,3)")
        if not np.all(np.isfinite(points)):
            raise ValueError("mesh points must be finite")
        if not self.cell_blocks:
            raise ValueError("mesh must contain at least one cell block")
        blocks = tuple(self.cell_blocks)
        dimensions = {block.dimension for block in blocks}
        if len(dimensions) != 1:
            raise ValueError("all stored cell blocks must have one topological dimension")
        for block in blocks:
            if block.connectivity.size and int(np.max(block.connectivity)) >= points.shape[0]:
                raise ValueError("cell connectivity references a missing point")
        object.__setattr__(self, "points_m", points)
        object.__setattr__(self, "cell_blocks", blocks)

    @property
    def topological_dimension(self) -> int:
        return self.cell_blocks[0].dimension

    @property
    def cell_count(self) -> int:
        return sum(block.connectivity.shape[0] for block in self.cell_blocks)

    @property
    def embedding_dimension(self) -> int:
        return self.points_m.shape[1]


@dataclass(frozen=True)
class AxialProjection:
    """A piecewise-constant 1D axial field copied to mesh cells."""

    values: np.ndarray
    cell_centers_m: np.ndarray
    normalized_axial_coordinate: np.ndarray
    tensile_axis: np.ndarray
    source_role: str = "1D FEM axial-field projection; not 2D/3D elasticity"


def _point_index_2d(i: int, j: int, ny: int) -> int:
    return i * (ny + 1) + j


def _point_index_3d(i: int, j: int, k: int, ny: int, nz: int) -> int:
    return (i * (ny + 1) + j) * (nz + 1) + k


def structured_rectangle_mesh(
    length_m: float,
    width_m: float,
    nx: int,
    ny: int,
) -> GeometryMesh:
    """Generate an actual 2D quadrilateral mesh of a tensile rectangle."""
    if not (np.isfinite(length_m) and length_m > 0.0):
        raise ValueError("length_m must be finite and positive")
    if not (np.isfinite(width_m) and width_m > 0.0):
        raise ValueError("width_m must be finite and positive")
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be at least 1")
    x = np.linspace(0.0, length_m, nx + 1)
    y = np.linspace(-0.5 * width_m, 0.5 * width_m, ny + 1)
    points = np.array([(xi, yj) for xi in x for yj in y], dtype=float)
    cells = []
    for i in range(nx):
        for j in range(ny):
            cells.append(
                [
                    _point_index_2d(i, j, ny),
                    _point_index_2d(i + 1, j, ny),
                    _point_index_2d(i + 1, j + 1, ny),
                    _point_index_2d(i, j + 1, ny),
                ]
            )
    return GeometryMesh(
        points,
        (CellBlock("quad", np.asarray(cells, dtype=int)),),
        source="generated structured tensile rectangle",
    )


def structured_box_mesh(
    length_m: float,
    width_m: float,
    thickness_m: float,
    nx: int,
    ny: int,
    nz: int,
) -> GeometryMesh:
    """Generate an actual 3D hexahedral mesh of a tensile box."""
    dimensions = (length_m, width_m, thickness_m)
    if any(not np.isfinite(value) or value <= 0.0 for value in dimensions):
        raise ValueError("length_m, width_m, and thickness_m must be finite and positive")
    if nx < 1 or ny < 1 or nz < 1:
        raise ValueError("nx, ny, and nz must be at least 1")
    x = np.linspace(0.0, length_m, nx + 1)
    y = np.linspace(-0.5 * width_m, 0.5 * width_m, ny + 1)
    z = np.linspace(-0.5 * thickness_m, 0.5 * thickness_m, nz + 1)
    points = np.array([(xi, yj, zk) for xi in x for yj in y for zk in z], dtype=float)
    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                cells.append(
                    [
                        _point_index_3d(i, j, k, ny, nz),
                        _point_index_3d(i + 1, j, k, ny, nz),
                        _point_index_3d(i + 1, j + 1, k, ny, nz),
                        _point_index_3d(i, j + 1, k, ny, nz),
                        _point_index_3d(i, j, k + 1, ny, nz),
                        _point_index_3d(i + 1, j, k + 1, ny, nz),
                        _point_index_3d(i + 1, j + 1, k + 1, ny, nz),
                        _point_index_3d(i, j + 1, k + 1, ny, nz),
                    ]
                )
    return GeometryMesh(
        points,
        (CellBlock("hexahedron", np.asarray(cells, dtype=int)),),
        source="generated structured tensile box",
    )


def _deduplicated_triangles(vertices: Iterable[tuple[float, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    point_ids: dict[tuple[float, float, float], int] = {}
    points: list[tuple[float, float, float]] = []
    triangles: list[list[int]] = []
    current: list[int] = []
    for vertex in vertices:
        key = tuple(float(value) for value in vertex)
        if key not in point_ids:
            point_ids[key] = len(points)
            points.append(key)
        current.append(point_ids[key])
        if len(current) == 3:
            if len(set(current)) == 3:
                triangles.append(current)
            current = []
    if current:
        raise ValueError("STL vertex stream is not divisible into triangles")
    if not triangles:
        raise ValueError("surface file contains no nondegenerate triangles")
    return np.asarray(points, dtype=float), np.asarray(triangles, dtype=int)


def _load_stl(path: Path, coordinate_scale_to_m: float) -> GeometryMesh:
    raw = path.read_bytes()
    binary = False
    if len(raw) >= 84:
        count = struct.unpack_from("<I", raw, 80)[0]
        binary = 84 + 50 * count == len(raw)
    vertices: list[tuple[float, float, float]] = []
    if binary:
        count = struct.unpack_from("<I", raw, 80)[0]
        for index in range(count):
            offset = 84 + 50 * index + 12
            for local in range(3):
                vertices.append(struct.unpack_from("<fff", raw, offset + 12 * local))
    else:
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError("STL is neither valid binary STL nor UTF-8 ASCII STL") from exc
        for line in text.splitlines():
            fields = line.strip().split()
            if len(fields) == 4 and fields[0].lower() == "vertex":
                vertices.append(tuple(float(value) for value in fields[1:4]))
    points, triangles = _deduplicated_triangles(vertices)
    return GeometryMesh(
        points * coordinate_scale_to_m,
        (CellBlock("triangle", triangles),),
        source=f"surface CAD: {path.name}",
        role="surface geometry/display mesh; volume FEM requires tetrahedralization",
    )


def _load_obj(path: Path, coordinate_scale_to_m: float) -> GeometryMesh:
    points: list[list[float]] = []
    triangles: list[list[int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.strip().split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "v" and len(fields) >= 4:
            points.append([float(fields[1]), float(fields[2]), float(fields[3])])
        elif fields[0] == "f" and len(fields) >= 4:
            face: list[int] = []
            for token in fields[1:]:
                raw_index = int(token.split("/", 1)[0])
                index = raw_index - 1 if raw_index > 0 else len(points) + raw_index
                if index < 0 or index >= len(points):
                    raise ValueError("OBJ face references a missing vertex")
                face.append(index)
            for local in range(1, len(face) - 1):
                triangle = [face[0], face[local], face[local + 1]]
                if len(set(triangle)) == 3:
                    triangles.append(triangle)
    if not points or not triangles:
        raise ValueError("OBJ must contain vertices and polygon faces")
    return GeometryMesh(
        np.asarray(points, dtype=float) * coordinate_scale_to_m,
        (CellBlock("triangle", np.asarray(triangles, dtype=int)),),
        source=f"surface CAD: {path.name}",
        role="surface geometry/display mesh; volume FEM requires tetrahedralization",
    )


def _load_with_meshio(path: Path, coordinate_scale_to_m: float) -> GeometryMesh:
    try:
        meshio = importlib.import_module("meshio")
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            f"{path.suffix} mesh import requires optional dependency 'meshio'; "
            "install requirements-cad.txt"
        ) from exc
    imported = meshio.read(path)
    blocks: list[CellBlock] = []
    aliases = {"hexahedron8": "hexahedron", "tetra4": "tetra", "quad4": "quad"}
    for block in imported.cells:
        cell_type = aliases.get(block.type, block.type)
        if cell_type in _CELL_DIMENSION:
            expected = _CELL_NODE_COUNT[cell_type]
            blocks.append(CellBlock(cell_type, np.asarray(block.data, dtype=int)[:, :expected]))
    if not blocks:
        raise ValueError("mesh file has no supported triangle/quad/tetra/hex/wedge/pyramid cells")
    top_dimension = max(block.dimension for block in blocks)
    blocks = [block for block in blocks if block.dimension == top_dimension]
    points = np.asarray(imported.points, dtype=float)
    if points.shape[1] == 3 and np.allclose(points[:, 2], points[0, 2]):
        points = points[:, :2]
    return GeometryMesh(
        points * coordinate_scale_to_m,
        tuple(blocks),
        source=f"mesh file: {path.name}",
    )


def _mesh_cad_with_gmsh(
    path: Path,
    coordinate_scale_to_m: float,
    target_dimension: int,
    characteristic_length_m: float | None,
) -> GeometryMesh:
    try:
        gmsh = importlib.import_module("gmsh")
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "STEP/IGES/BREP import requires optional dependency 'gmsh'; "
            "install requirements-cad.txt"
        ) from exc
    if target_dimension not in (2, 3):
        raise ValueError("target_dimension must be 2 or 3 for CAD meshing")
    gmsh.initialize()
    try:
        gmsh.model.add("tensile_geometry")
        entities = gmsh.model.occ.importShapes(str(path))
        if not entities:
            raise ValueError(f"Gmsh imported no entities from {path}")
        gmsh.model.occ.synchronize()
        if characteristic_length_m is not None:
            if not np.isfinite(characteristic_length_m) or characteristic_length_m <= 0.0:
                raise ValueError("characteristic_length_m must be finite and positive")
            native_size = characteristic_length_m / coordinate_scale_to_m
            gmsh.option.setNumber("Mesh.MeshSizeMin", native_size)
            gmsh.option.setNumber("Mesh.MeshSizeMax", native_size)
        gmsh.model.mesh.generate(target_dimension)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        node_tags = np.asarray(node_tags, dtype=np.int64)
        points = np.asarray(coordinates, dtype=float).reshape(-1, 3) * coordinate_scale_to_m
        tag_to_index = {int(tag): index for index, tag in enumerate(node_tags)}
        element_types, _, element_nodes = gmsh.model.mesh.getElements(target_dimension)
        blocks: list[CellBlock] = []
        for element_type, flat_nodes in zip(element_types, element_nodes):
            properties = gmsh.model.mesh.getElementProperties(int(element_type))
            element_dimension = int(properties[1])
            node_count = int(properties[3])
            primary_node_count = int(properties[5])
            cell_type = _GMSH_CELL_NAMES.get((element_dimension, primary_node_count))
            if cell_type is None:
                continue
            raw = np.asarray(flat_nodes, dtype=np.int64).reshape(-1, node_count)
            primary = raw[:, :primary_node_count]
            connectivity = np.vectorize(tag_to_index.__getitem__, otypes=[int])(primary)
            blocks.append(CellBlock(cell_type, connectivity))
        if not blocks:
            raise ValueError("Gmsh generated no supported first-order top-dimensional cells")
        return GeometryMesh(
            points,
            tuple(blocks),
            source=f"Gmsh CAD volume/surface mesh: {path.name}",
        )
    finally:
        gmsh.finalize()


def load_geometry_mesh(
    path: Path,
    *,
    coordinate_scale_to_m: float = 1.0,
    target_dimension: int | None = None,
    characteristic_length_m: float | None = None,
) -> GeometryMesh:
    """Read CAD/mesh geometry with explicit units and optional meshing backend.

    STL and OBJ are read without optional dependencies as surface meshes. STEP,
    IGES, and BREP are converted to top-dimensional cells through Gmsh. Common
    solver mesh formats are delegated to meshio.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"geometry file not found: {source}")
    if not np.isfinite(coordinate_scale_to_m) or coordinate_scale_to_m <= 0.0:
        raise ValueError("coordinate_scale_to_m must be finite and positive")
    suffix = source.suffix.lower()
    if suffix == ".stl":
        mesh = _load_stl(source, coordinate_scale_to_m)
    elif suffix == ".obj":
        mesh = _load_obj(source, coordinate_scale_to_m)
    elif suffix in {".step", ".stp", ".iges", ".igs", ".brep"}:
        mesh = _mesh_cad_with_gmsh(
            source,
            coordinate_scale_to_m,
            3 if target_dimension is None else target_dimension,
            characteristic_length_m,
        )
    elif suffix == ".npz":
        mesh = load_mesh_npz(source)
    else:
        mesh = _load_with_meshio(source, coordinate_scale_to_m)
    if target_dimension is not None and mesh.topological_dimension != target_dimension:
        raise ValueError(
            f"loaded topological dimension {mesh.topological_dimension}, "
            f"but target_dimension={target_dimension}; a surface CAD file is not a volume mesh"
        )
    return mesh


def cell_centers(mesh: GeometryMesh) -> np.ndarray:
    """Return centers in the same flattened order as the mesh cell blocks."""
    return np.vstack(
        [np.mean(mesh.points_m[block.connectivity], axis=1) for block in mesh.cell_blocks]
    )


def _unit_axis(axis: np.ndarray | tuple[float, ...], embedding_dimension: int) -> np.ndarray:
    direction = np.asarray(axis, dtype=float)
    if direction.shape != (embedding_dimension,):
        raise ValueError(f"tensile axis must have {embedding_dimension} components")
    norm = float(np.linalg.norm(direction))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("tensile axis must be finite and nonzero")
    return direction / norm


def project_normal_stress(stress_tensor_pa: np.ndarray, tensile_axis: np.ndarray) -> np.ndarray:
    """Project a symmetric or nonsymmetric stress tensor to sigma_nn = n.T sigma n.

    This exact tensor projection does not pass shear components as independent
    variables to the probability model.  The returned scalar is its only stress
    input.
    """
    stress = np.asarray(stress_tensor_pa, dtype=float)
    if stress.ndim < 2 or stress.shape[-1] != stress.shape[-2]:
        raise ValueError("stress_tensor_pa must end in square tensor dimensions")
    axis = _unit_axis(tensile_axis, stress.shape[-1])
    if not np.all(np.isfinite(stress)):
        raise ValueError("stress tensor must be finite")
    return np.einsum("i,...ij,j->...", axis, stress, axis)


def map_axial_element_field(
    mesh: GeometryMesh,
    axial_x_mid_m: np.ndarray,
    axial_values: np.ndarray,
    *,
    tensile_axis: np.ndarray | tuple[float, ...] | None = None,
) -> AxialProjection:
    """Copy a 1D piecewise-constant element field onto 2D/3D mesh cells.

    The mesh axial coordinate is normalized over its geometric extent and mapped
    to the normalized extent of the 1D bar.  This controlled visualization
    approximation preserves axial ordering, but does not solve equilibrium on
    the multidimensional mesh.
    """
    x_mid = np.asarray(axial_x_mid_m, dtype=float)
    values = np.asarray(axial_values, dtype=float)
    if x_mid.ndim != 1 or values.shape != x_mid.shape or x_mid.size == 0:
        raise ValueError("axial_x_mid_m and axial_values must be equal nonempty 1D arrays")
    order = np.argsort(x_mid)
    x_mid = x_mid[order]
    values = values[order]
    if not np.all(np.isfinite(x_mid)) or not np.all(np.isfinite(values)):
        raise ValueError("axial source field must be finite")
    if x_mid.size > 1 and np.any(np.diff(x_mid) <= 0.0):
        raise ValueError("axial element centers must be unique")
    if tensile_axis is None:
        tensile_axis = np.eye(mesh.embedding_dimension, dtype=float)[0]
    axis = _unit_axis(tensile_axis, mesh.embedding_dimension)
    centers = cell_centers(mesh)
    projected = centers @ axis
    extent = float(np.max(projected) - np.min(projected))
    if extent <= 0.0:
        raise ValueError("mesh has zero extent along the tensile axis")
    xi = (projected - float(np.min(projected))) / extent
    if x_mid.size == 1:
        mapped = np.full(mesh.cell_count, values[0], dtype=float)
    else:
        source_xi = (x_mid - x_mid[0]) / (x_mid[-1] - x_mid[0])
        boundaries = 0.5 * (source_xi[:-1] + source_xi[1:])
        indices = np.searchsorted(boundaries, xi, side="right")
        mapped = values[indices]
    return AxialProjection(mapped, centers, xi, axis)


def boundary_faces(mesh: GeometryMesh) -> tuple[list[np.ndarray], np.ndarray]:
    """Extract visible cells/faces and their owner-cell indices."""
    if mesh.topological_dimension == 2:
        faces: list[np.ndarray] = []
        owners: list[int] = []
        offset = 0
        for block in mesh.cell_blocks:
            for local, connectivity in enumerate(block.connectivity):
                faces.append(mesh.points_m[connectivity])
                owners.append(offset + local)
            offset += block.connectivity.shape[0]
        return faces, np.asarray(owners, dtype=int)

    patterns = {
        "tetra": ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
        "hexahedron": (
            (0, 1, 2, 3),
            (4, 7, 6, 5),
            (0, 4, 5, 1),
            (1, 5, 6, 2),
            (2, 6, 7, 3),
            (3, 7, 4, 0),
        ),
        "wedge": ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
        "pyramid": ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
    }
    candidates: dict[tuple[int, ...], list[tuple[np.ndarray, int]]] = {}
    owner_offset = 0
    for block in mesh.cell_blocks:
        for local, connectivity in enumerate(block.connectivity):
            for pattern in patterns[block.cell_type]:
                ids = np.asarray([connectivity[index] for index in pattern], dtype=int)
                key = tuple(sorted(int(value) for value in ids))
                candidates.setdefault(key, []).append((ids, owner_offset + local))
        owner_offset += block.connectivity.shape[0]
    faces = []
    owners = []
    for matches in candidates.values():
        if len(matches) == 1:
            ids, owner = matches[0]
            faces.append(mesh.points_m[ids])
            owners.append(owner)
    return faces, np.asarray(owners, dtype=int)


def save_mesh_npz(path: Path, mesh: GeometryMesh) -> None:
    """Write a dependency-free mesh bundle with explicit cell block metadata."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, np.ndarray] = {
        "points_m": mesh.points_m,
        "cell_types": np.asarray([block.cell_type for block in mesh.cell_blocks]),
        "source": np.asarray(mesh.source),
        "role": np.asarray(mesh.role),
    }
    for index, block in enumerate(mesh.cell_blocks):
        payload[f"cells_{index}"] = block.connectivity
    np.savez_compressed(destination, **payload)


def load_mesh_npz(path: Path) -> GeometryMesh:
    """Read a mesh bundle written by :func:`save_mesh_npz`."""
    with np.load(Path(path), allow_pickle=False) as data:
        cell_types = [str(value) for value in data["cell_types"]]
        blocks = tuple(
            CellBlock(cell_type, np.asarray(data[f"cells_{index}"], dtype=int))
            for index, cell_type in enumerate(cell_types)
        )
        source = str(data["source"].item()) if "source" in data else f"mesh bundle: {path}"
        role = str(data["role"].item()) if "role" in data else "geometry/display mesh"
        return GeometryMesh(np.asarray(data["points_m"], dtype=float), blocks, source, role)
