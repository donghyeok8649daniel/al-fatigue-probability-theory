# === 한국어 파일 안내 시작 ===
# - 파일 역할: 경량 FEM mesh UI의 1D/2D/3D 생성, clipping 및 headless rendering helper를 검증한다.
# - 주요 클래스: TestMeshAppConfig, TestMeshVisibility
# - 주요 함수/메서드: TestMeshAppConfig.test_all_three_dimensions_generate_expected_cells
#   TestMeshAppConfig.test_invalid_dimension_and_resolution_are_rejected
#   TestMeshVisibility.test_axial_clip_exposes_subset_and_compacts_nodes
#   TestMeshVisibility.test_headless_smoke_renders_every_dimension
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import tempfile
import unittest
from pathlib import Path

import numpy as np

from simulations.fem_mesh_app import (
    MeshAppConfig,
    axial_visibility_field,
    clipped_mesh_view,
    create_geometry_mesh,
    run_headless_smoke,
    validate_mesh_app_config,
)


class TestMeshAppConfig(unittest.TestCase):
    def test_all_three_dimensions_generate_expected_cells(self) -> None:
        one = create_geometry_mesh(MeshAppConfig(dimension=1, nx=5))
        two = create_geometry_mesh(MeshAppConfig(dimension=2, nx=5, ny=2))
        three = create_geometry_mesh(MeshAppConfig(dimension=3, nx=5, ny=2, nz=2))
        self.assertEqual((one.cell_count, two.cell_count, three.cell_count), (5, 10, 20))
        self.assertEqual(
            (one.topological_dimension, two.topological_dimension, three.topological_dimension),
            (1, 2, 3),
        )

    def test_invalid_dimension_and_resolution_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_mesh_app_config(MeshAppConfig(dimension=4))
        with self.assertRaises(ValueError):
            validate_mesh_app_config(MeshAppConfig(nx=0))


class TestMeshVisibility(unittest.TestCase):
    def test_axial_clip_exposes_subset_and_compacts_nodes(self) -> None:
        mesh = create_geometry_mesh(MeshAppConfig(dimension=3, nx=8, ny=2, nz=2))
        values = axial_visibility_field(mesh)
        clipped, clipped_values = clipped_mesh_view(mesh, values, 0.5)
        self.assertLess(clipped.cell_count, mesh.cell_count)
        self.assertLess(clipped.points_m.shape[0], mesh.points_m.shape[0])
        self.assertEqual(clipped_values.shape, (clipped.cell_count,))
        self.assertTrue(np.all(clipped_values <= 0.5 + 1.0e-12))

    def test_headless_smoke_renders_every_dimension(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            previews = output / "previews"
            counts = run_headless_smoke(output, previews)
            for dimension in (1, 2, 3):
                self.assertTrue((output / f"mesh_{dimension}d" / "geometry_mesh.npz").is_file())
                self.assertTrue((previews / f"mesh_{dimension}d.png").is_file())
            self.assertTrue((previews / "lightweight_mesh_ui.png").is_file())
        self.assertEqual(counts, {"cells_1d": 8, "cells_2d": 24, "cells_3d": 48})


if __name__ == "__main__":
    unittest.main()
