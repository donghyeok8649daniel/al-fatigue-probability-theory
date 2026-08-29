# === 한국어 파일 안내 시작 ===
# - 파일 역할: 2D/3D tensile UI가 1D 축방향 scalar만 시각적으로 펼치고 transverse mechanics를 만들지 않는지 검증한다.
# - 주요 클래스: TestFemTensionUI
# - 주요 함수/메서드: test_extruded_scalar_is_transversely_constant, test_axial_snapshot_uses_only_requested_scalar, test_invalid_field_rejected
# - 주의: 이 테스트는 시각화 데이터의 구조를 검증하며 C FEM 해석 자체는 fem1d self-test가 담당한다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from simulations.fem_tension_ui import axial_snapshot, extruded_scalar_2d


class TestFemTensionUI(unittest.TestCase):
    @staticmethod
    def _histories():
        node_dtype = [
            ("step", int),
            ("time_s", float),
            ("node", int),
            ("x_m", float),
            ("displacement_m", float),
            ("applied_stress_pa", float),
        ]
        elem_dtype = [
            ("step", int),
            ("time_s", float),
            ("element", int),
            ("x_mid_m", float),
            ("strain", float),
            ("stress_pa", float),
            ("applied_stress_pa", float),
        ]
        nodes = np.array(
            [
                (0, 0.0, 0, 0.0, 0.0, 10.0),
                (0, 0.0, 1, 0.5, 1.0e-4, 10.0),
                (0, 0.0, 2, 1.0, 2.0e-4, 10.0),
            ],
            dtype=node_dtype,
        )
        elements = np.array(
            [
                (0, 0.0, 0, 0.25, 2.0e-4, 20.0e6, 10.0),
                (0, 0.0, 1, 0.75, 4.0e-4, 40.0e6, 10.0),
            ],
            dtype=elem_dtype,
        )
        return nodes, elements

    def test_extruded_scalar_is_transversely_constant(self) -> None:
        x = np.array([0.0, 0.5, 1.0])
        scalar = np.array([2.0, 5.0])
        xg, yg, field = extruded_scalar_2d(x, scalar, 0.1)
        self.assertEqual(xg.shape, (2, 3))
        self.assertEqual(yg.shape, (2, 3))
        self.assertEqual(field.shape, (1, 2))
        np.testing.assert_allclose(field[0], scalar)
        np.testing.assert_allclose(yg[0], -0.1)
        np.testing.assert_allclose(yg[1], 0.1)

    def test_axial_snapshot_uses_only_requested_scalar(self) -> None:
        nodes, elements = self._histories()
        stress = axial_snapshot(nodes, elements, 0, "stress")
        strain = axial_snapshot(nodes, elements, 0, "strain")
        np.testing.assert_allclose(stress["scalar"], [20.0, 40.0])
        np.testing.assert_allclose(strain["scalar"], [2.0e-4, 4.0e-4])
        self.assertNotIn("shear", stress)
        self.assertNotIn("von_mises", stress)

    def test_invalid_field_rejected(self) -> None:
        nodes, elements = self._histories()
        with self.assertRaises(ValueError):
            axial_snapshot(nodes, elements, 0, "von-mises")


if __name__ == "__main__":
    unittest.main()
