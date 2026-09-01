# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: TestFem1DVisualizer
# - 주요 함수/메서드: TestFem1DVisualizer._history, TestFem1DVisualizer.test_select_peak_tension
#   TestFem1DVisualizer.test_select_peak_absolute, TestFem1DVisualizer.test_select_final
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
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
