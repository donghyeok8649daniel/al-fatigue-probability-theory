# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D element history를 실제 2D/3D mesh cell field로 내보내는 통합 경로를 검증한다.
# - 주요 클래스: TestTensileMeshProjection
# - 주요 함수/메서드: TestTensileMeshProjection.test_probability_scalar_can_be_mapped_to_real_2d_cells
#   TestTensileMeshProjection.test_stl_surface_is_read_and_colored_with_one_normal_scalar
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import tempfile
import unittest
from pathlib import Path

import numpy as np

from simulations.run_tensile_mesh_projection import run_projection


class TestTensileMeshProjection(unittest.TestCase):
    def test_probability_scalar_can_be_mapped_to_real_2d_cells(self) -> None:
        text = """time_s,step,element,x_mid_m,stress_pa,critical_tail_probability
0.1,4,0,0.125,10,0.01
0.1,4,1,0.375,20,0.02
0.1,4,2,0.625,30,0.03
0.1,4,3,0.875,40,0.04
"""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            history = root / "history.csv"
            history.write_text(text, encoding="utf-8")
            output = root / "out"
            summary = run_projection(
                history_csv=history,
                field="critical_tail_probability",
                step=4,
                output_dir=output,
                dimension=2,
                length_m=1.0,
                width_m=0.2,
                thickness_m=0.01,
                nx=4,
                ny=2,
                nz=1,
            )
            self.assertTrue((output / "geometry_mesh.npz").is_file())
            self.assertTrue((output / "projected_cell_field.csv").is_file())
            self.assertTrue((output / "mesh_field_preview.png").is_file())
            self.assertTrue((output / "summary.json").is_file())
        self.assertEqual(summary["cells"], 8)
        self.assertAlmostEqual(float(summary["field_min"]), 0.01)
        self.assertAlmostEqual(float(summary["field_max"]), 0.04)
        self.assertIn("not 2D/3D elasticity", str(summary["mapping"]))

    def test_stl_surface_is_read_and_colored_with_one_normal_scalar(self) -> None:
        history_text = """time_s,step,element,x_mid_m,stress_pa
0.1,2,0,0.125,10
0.1,2,1,0.375,20
0.1,2,2,0.625,30
0.1,2,3,0.875,40
"""
        stl_text = """solid rectangle
facet normal 0 0 1
 outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 1 1 0
 endloop
endfacet
facet normal 0 0 1
 outer loop
  vertex 0 0 0
  vertex 1 1 0
  vertex 0 1 0
 endloop
endfacet
endsolid rectangle
"""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            history = root / "history.csv"
            geometry = root / "specimen.stl"
            history.write_text(history_text, encoding="utf-8")
            geometry.write_text(stl_text, encoding="utf-8")
            summary = run_projection(
                history_csv=history,
                field="stress_pa",
                step=2,
                output_dir=root / "cad_out",
                dimension=2,
                length_m=1.0,
                width_m=1.0,
                thickness_m=0.1,
                nx=1,
                ny=1,
                nz=1,
                geometry_path=geometry,
                tensile_axis=np.array([1.0, 0.0, 0.0]),
            )
            self.assertTrue((root / "cad_out" / "mesh_field_preview.png").is_file())
        self.assertEqual(summary["cells"], 2)
        self.assertIn("surface CAD", str(summary["source"]))


if __name__ == "__main__":
    unittest.main()
