# === 한국어 파일 안내 시작 ===
# - 파일 역할: 실제 2D/3D cell connectivity와 normal-only scalar field를 mesh 경계 위에 시각화한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: _validate_cell_values, _normalization, plot_geometry_mesh_2d, plot_geometry_mesh_3d
#   save_geometry_mesh_preview
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Visualization of actual 2D/3D mesh connectivity with one axial scalar."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from simulations.fem_geometry_mesh import GeometryMesh, boundary_faces


def _validate_cell_values(mesh: GeometryMesh, cell_values: np.ndarray) -> np.ndarray:
    values = np.asarray(cell_values, dtype=float)
    if values.shape != (mesh.cell_count,):
        raise ValueError(f"cell_values must have shape ({mesh.cell_count},)")
    if not np.all(np.isfinite(values)):
        raise ValueError("cell_values must be finite")
    return values


def _normalization(values: np.ndarray) -> Normalize:
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if np.isclose(vmin, vmax):
        delta = max(abs(vmin) * 0.01, 1.0e-12)
        vmin -= delta
        vmax += delta
    return Normalize(vmin=vmin, vmax=vmax)


def plot_geometry_mesh_2d(
    ax,
    mesh: GeometryMesh,
    cell_values: np.ndarray,
    *,
    norm: Normalize | None = None,
    cmap="viridis",
):
    """Draw every actual triangle/quad cell of a planar two-dimensional mesh."""
    if mesh.topological_dimension != 2 or mesh.embedding_dimension != 2:
        raise ValueError("2D plot requires a planar 2D triangle/quad mesh")
    values = _validate_cell_values(mesh, cell_values)
    faces, owners = boundary_faces(mesh)
    color_norm = _normalization(values) if norm is None else norm
    color_map = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    collection = PolyCollection(
        faces,
        array=values[owners],
        cmap=color_map,
        norm=color_norm,
        edgecolors="black",
        linewidths=0.35,
    )
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    return collection


def plot_geometry_mesh_3d(
    ax,
    mesh: GeometryMesh,
    cell_values: np.ndarray,
    *,
    norm: Normalize | None = None,
    cmap="viridis",
):
    """Draw boundary faces of a volume mesh or an embedded surface mesh."""
    values = _validate_cell_values(mesh, cell_values)
    faces, owners = boundary_faces(mesh)
    color_norm = _normalization(values) if norm is None else norm
    color_map = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    faces_3d = [
        np.column_stack([face, np.zeros(face.shape[0])]) if face.shape[1] == 2 else face
        for face in faces
    ]
    colors = color_map(color_norm(values[owners]))
    collection = Poly3DCollection(
        faces_3d,
        facecolors=colors,
        edgecolors="black",
        linewidths=0.22,
        alpha=0.96,
    )
    ax.add_collection3d(collection)
    points = (
        np.column_stack([mesh.points_m, np.zeros(mesh.points_m.shape[0])])
        if mesh.embedding_dimension == 2
        else mesh.points_m
    )
    mins = np.min(points, axis=0)
    maxs = np.max(points, axis=0)
    spans = np.maximum(maxs - mins, 1.0e-12)
    margin = 0.04 * float(np.max(spans))
    ax.set_xlim(mins[0] - margin, maxs[0] + margin)
    ax.set_ylim(mins[1] - margin, maxs[1] + margin)
    ax.set_zlim(mins[2] - margin, maxs[2] + margin)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    scalar_mappable = cm.ScalarMappable(norm=color_norm, cmap=color_map)
    scalar_mappable.set_array(values)
    return scalar_mappable


def save_geometry_mesh_preview(
    path: Path,
    mesh: GeometryMesh,
    cell_values: np.ndarray,
    *,
    field_label: str,
    title: str,
    dpi: int = 180,
) -> dict[str, float | int | str]:
    """Save one deterministic preview and return exact mesh/value metadata."""
    values = _validate_cell_values(mesh, cell_values)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    planar = mesh.topological_dimension == 2 and mesh.embedding_dimension == 2
    if planar:
        fig, ax = plt.subplots(figsize=(9.0, 4.2))
        artist = plot_geometry_mesh_2d(ax, mesh, values)
    else:
        fig = plt.figure(figsize=(9.0, 5.8))
        ax = fig.add_subplot(111, projection="3d")
        artist = plot_geometry_mesh_3d(ax, mesh, values)
    ax.set_title(title)
    fig.colorbar(artist, ax=ax, shrink=0.78, pad=0.08, label=field_label)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return {
        "points": int(mesh.points_m.shape[0]),
        "cells": int(mesh.cell_count),
        "topological_dimension": int(mesh.topological_dimension),
        "field_min": float(np.min(values)),
        "field_max": float(np.max(values)),
        "source": mesh.source,
        "role": mesh.role,
    }
