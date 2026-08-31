# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: test_obj_lines_faces_and_dimension, test_ascii_and_binary_stl
#   test_ascii_ply_and_legacy_vtk, test_binary_ply_and_unsupported_cad_are_honestly_rejected
#   test_scroll_zoom_reset_and_viewport_dimensions
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
from pathlib import Path
from types import SimpleNamespace
import struct

import matplotlib.pyplot as plt
import numpy as np
import pytest

from simulations.mesh_viewer import CADNavigation, MeshGeometry, MeshViewport, load_mesh


def test_obj_lines_faces_and_dimension(tmp_path: Path):
    path = tmp_path / "plate.obj"
    path.write_text("v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\nl 1 3\n",
                    encoding="utf-8")
    mesh = load_mesh(path)
    assert mesh.dimension == 2
    assert mesh.faces == ((0, 1, 2, 3),)
    assert mesh.lines == ((0, 2),)


def test_ascii_and_binary_stl(tmp_path: Path):
    ascii_path = tmp_path / "one.stl"
    ascii_path.write_text("solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\n"
        "vertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid x\n", encoding="utf-8")
    assert load_mesh(ascii_path).faces == ((0, 1, 2),)

    binary_path = tmp_path / "one_binary.stl"
    header = b"binary triangle".ljust(80, b"\0") + struct.pack("<I", 1)
    triangle = struct.pack("<12fH", 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
    binary_path.write_bytes(header + triangle)
    mesh = load_mesh(binary_path)
    assert mesh.vertices.shape == (3, 3)
    assert mesh.dimension == 2


def test_ascii_ply_and_legacy_vtk(tmp_path: Path):
    ply = tmp_path / "plate.ply"
    ply.write_text("ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\n"
        "property float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_indices\n"
        "end_header\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n", encoding="utf-8")
    assert load_mesh(ply).faces == ((0, 1, 2),)

    vtk = tmp_path / "line.vtk"
    vtk.write_text("# vtk DataFile Version 3.0\nline\nASCII\nDATASET POLYDATA\n"
        "POINTS 3 float\n0 0 0 1 0 0 2 0 0\nLINES 1 4\n3 0 1 2\n", encoding="utf-8")
    mesh = load_mesh(vtk)
    assert mesh.dimension == 1
    assert mesh.lines == ((0, 1, 2),)


def test_binary_ply_and_unsupported_cad_are_honestly_rejected(tmp_path: Path):
    ply = tmp_path / "binary.ply"
    ply.write_bytes(b"ply\nformat binary_little_endian 1.0\nend_header\n")
    with pytest.raises(ValueError, match="only ASCII PLY"):
        load_mesh(ply)
    step = tmp_path / "part.step"; step.write_text("ISO-10303-21", encoding="utf-8")
    with pytest.raises(ValueError, match="OBJ, STL, PLY"):
        load_mesh(step)


def test_scroll_zoom_reset_and_viewport_dimensions():
    mesh = MeshGeometry(np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0]]),
                        faces=((0, 1, 2),), source="memory.obj")
    viewport = MeshViewport(mesh, 2)
    before = np.ptp(viewport.axes.get_xlim())
    viewport.navigation.on_scroll(SimpleNamespace(inaxes=viewport.axes, button="up",
                                                   xdata=0.5, ydata=0.5))
    assert np.ptp(viewport.axes.get_xlim()) < before
    viewport.navigation.reset()
    assert np.isclose(np.ptp(viewport.axes.get_xlim()), before)
    plt.close(viewport.figure)

    one_d = MeshViewport(MeshGeometry(np.array([[0.], [1.], [2.]]), source="line.vtk"), 1)
    assert one_d.dimension == 1
    plt.close(one_d.figure)
