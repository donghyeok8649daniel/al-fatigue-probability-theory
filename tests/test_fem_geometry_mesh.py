# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D/2D/3D mesh 생성, CAD surface 입력, normal projection과 axial field mapping을 회귀검증한다.
# - 주요 클래스: TestStructuredMeshes, TestNormalOnlyProjection, TestGeometryImport
# - 주요 함수/메서드: TestStructuredMeshes.test_line_mesh_is_a_first_class_supported_dimension
#   TestStructuredMeshes.test_rectangle_has_actual_quad_connectivity
#   TestStructuredMeshes.test_box_has_hex_cells_and_only_boundary_faces
#   TestStructuredMeshes.test_dependency_free_npz_round_trip
#   TestNormalOnlyProjection.test_axial_mapping_operates_on_line_mesh
#   TestNormalOnlyProjection.test_tensor_is_reduced_to_sigma_nn_only
#   TestNormalOnlyProjection.test_axial_mapping_is_constant_across_each_cross_section
#   TestGeometryImport.test_ascii_stl_is_read_as_surface_not_volume
#   TestGeometryImport.test_binary_stl_is_detected_by_exact_record_length
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

from simulations.fem_geometry_mesh import (
    boundary_faces,
    cell_centers,
    load_geometry_mesh,
    load_mesh_npz,
    map_axial_element_field,
    project_normal_stress,
    save_mesh_npz,
    structured_box_mesh,
    structured_line_mesh,
    structured_rectangle_mesh,
)


class TestStructuredMeshes(unittest.TestCase):
    def test_line_mesh_is_a_first_class_supported_dimension(self) -> None:
        mesh = structured_line_mesh(0.05, 5)
        self.assertEqual(mesh.topological_dimension, 1)
        self.assertEqual(mesh.embedding_dimension, 1)
        self.assertEqual(mesh.points_m.shape, (6, 1))
        self.assertEqual(mesh.cell_count, 5)
        faces, owners = boundary_faces(mesh)
        self.assertEqual(len(faces), 5)
        np.testing.assert_array_equal(owners, np.arange(5))

    def test_rectangle_has_actual_quad_connectivity(self) -> None:
        mesh = structured_rectangle_mesh(0.05, 0.01, 5, 2)
        self.assertEqual(mesh.topological_dimension, 2)
        self.assertEqual(mesh.points_m.shape, (18, 2))
        self.assertEqual(mesh.cell_count, 10)
        self.assertEqual(mesh.cell_blocks[0].cell_type, "quad")
        self.assertEqual(cell_centers(mesh).shape, (10, 2))
        faces, owners = boundary_faces(mesh)
        self.assertEqual(len(faces), 10)
        np.testing.assert_array_equal(owners, np.arange(10))

    def test_box_has_hex_cells_and_only_boundary_faces(self) -> None:
        mesh = structured_box_mesh(0.05, 0.01, 0.002, 2, 2, 1)
        self.assertEqual(mesh.topological_dimension, 3)
        self.assertEqual(mesh.points_m.shape, (18, 3))
        self.assertEqual(mesh.cell_count, 4)
        faces, owners = boundary_faces(mesh)
        self.assertEqual(len(faces), 16)
        self.assertTrue(np.all((owners >= 0) & (owners < mesh.cell_count)))

    def test_dependency_free_npz_round_trip(self) -> None:
        mesh = structured_box_mesh(0.03, 0.004, 0.001, 3, 1, 1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mesh.npz"
            save_mesh_npz(path, mesh)
            loaded = load_mesh_npz(path)
        np.testing.assert_allclose(loaded.points_m, mesh.points_m)
        np.testing.assert_array_equal(
            loaded.cell_blocks[0].connectivity,
            mesh.cell_blocks[0].connectivity,
        )
        self.assertEqual(loaded.role, mesh.role)


class TestNormalOnlyProjection(unittest.TestCase):
    def test_axial_mapping_operates_on_line_mesh(self) -> None:
        mesh = structured_line_mesh(1.0, 4)
        result = map_axial_element_field(
            mesh,
            np.array([0.125, 0.375, 0.625, 0.875]),
            np.array([1.0, 2.0, 3.0, 4.0]),
        )
        np.testing.assert_allclose(result.values, [1.0, 2.0, 3.0, 4.0])

    def test_tensor_is_reduced_to_sigma_nn_only(self) -> None:
        stress = np.array(
            [
                [[100.0, 25.0, 0.0], [25.0, 40.0, 0.0], [0.0, 0.0, 10.0]],
                [[-20.0, 90.0, 3.0], [90.0, 60.0, 2.0], [3.0, 2.0, 4.0]],
            ]
        )
        projected = project_normal_stress(stress, np.array([1.0, 0.0, 0.0]))
        np.testing.assert_allclose(projected, [100.0, -20.0])

    def test_axial_mapping_is_constant_across_each_cross_section(self) -> None:
        mesh = structured_rectangle_mesh(1.0, 0.2, 4, 3)
        result = map_axial_element_field(
            mesh,
            np.array([0.125, 0.375, 0.625, 0.875]),
            np.array([10.0, 20.0, 30.0, 40.0]),
        )
        mapped = result.values.reshape(4, 3)
        np.testing.assert_allclose(mapped, np.array([[10.0] * 3, [20.0] * 3, [30.0] * 3, [40.0] * 3]))
        self.assertIn("not 2D/3D elasticity", result.source_role)


class TestGeometryImport(unittest.TestCase):
    def test_ascii_stl_is_read_as_surface_not_volume(self) -> None:
        stl = """solid one
facet normal 0 0 1
 outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 0 1 0
 endloop
endfacet
endsolid one
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "triangle.stl"
            path.write_text(stl, encoding="utf-8")
            mesh = load_geometry_mesh(path, coordinate_scale_to_m=1.0e-3)
            with self.assertRaisesRegex(ValueError, "surface CAD file is not a volume mesh"):
                load_geometry_mesh(path, target_dimension=3)
        self.assertEqual(mesh.topological_dimension, 2)
        self.assertEqual(mesh.embedding_dimension, 3)
        self.assertAlmostEqual(float(np.max(mesh.points_m)), 1.0e-3)

    def test_binary_stl_is_detected_by_exact_record_length(self) -> None:
        header = b"binary triangle".ljust(80, b" ")
        record = struct.pack(
            "<12fH",
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "triangle_binary.stl"
            path.write_bytes(header + struct.pack("<I", 1) + record)
            mesh = load_geometry_mesh(path)
        self.assertEqual(mesh.cell_count, 1)
        self.assertEqual(mesh.points_m.shape, (3, 3))


if __name__ == "__main__":
    unittest.main()
