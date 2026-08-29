# === 한국어 파일 안내 시작 ===
# - 파일 역할: Python 1D FEM visualizer의 CSV 로딩과 snapshot 선택 규칙을 검증한다.
# - 주요 클래스: TestFem1DVisualizer
# - 주요 함수/메서드: test_select_peak_tension, test_select_peak_absolute, test_select_final
# - 주의: C FEM 자체의 해석 정확도는 fem1d_solver --self-test에서 별도로 검증한다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from simulations.visualize_fem1d import select_snapshot_step


class TestFem1DVisualizer(unittest.TestCase):
    @staticmethod
    def _history():
        dtype = [
            ("step", int),
            ("stress_pa", float),
        ]
        return np.array(
            [
                (0, 0.0),
                (0, 0.0),
                (1, 10.0),
                (1, 10.0),
                (2, -20.0),
                (2, -20.0),
                (3, 5.0),
                (3, 5.0),
            ],
            dtype=dtype,
        )

    def test_select_peak_tension(self) -> None:
        self.assertEqual(select_snapshot_step(self._history(), "peak-tension"), 1)

    def test_select_peak_absolute(self) -> None:
        self.assertEqual(select_snapshot_step(self._history(), "peak-absolute"), 2)

    def test_select_final(self) -> None:
        self.assertEqual(select_snapshot_step(self._history(), "final"), 3)


if __name__ == "__main__":
    unittest.main()
