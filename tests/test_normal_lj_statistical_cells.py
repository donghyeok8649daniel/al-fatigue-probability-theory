# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: TestNormalLJStatisticalCells
# - 주요 함수/메서드: TestNormalLJStatisticalCells.test_independent_limit
#   TestNormalLJStatisticalCells.test_fully_identical_limit
#   TestNormalLJStatisticalCells.test_anticorrelation_can_raise_effective_count
#   TestNormalLJStatisticalCells.test_positive_window_estimator
#   TestNormalLJStatisticalCells.test_positive_window_fully_identical_limit
#   TestNormalLJStatisticalCells.test_identical_pair_msd
#   TestNormalLJStatisticalCells.test_independent_any_event_probability
#   TestNormalLJStatisticalCells.test_identical_block_any_event_probability
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.normal_lj_statistical_cells import (
    effective_independent_count,
    finite_correlation_factor,
    identical_block_any_event_probability,
    identical_pair_msd,
    independent_any_event_probability,
    positive_window_empirical_axial_length,
    positive_window_empirical_correlation_factor,
    positive_window_empirical_effective_count,
    variance_equivalent_axial_length,
)


class TestNormalLJStatisticalCells(unittest.TestCase):
    def test_independent_limit(self) -> None:
        m = 16
        rho = np.zeros(m)
        rho[0] = 1.0
        self.assertAlmostEqual(finite_correlation_factor(rho), 1.0, places=14)
        self.assertAlmostEqual(effective_independent_count(rho), float(m), places=14)
        self.assertAlmostEqual(
            variance_equivalent_axial_length(rho, 2.86e-10), 2.86e-10, places=24
        )

    def test_fully_identical_limit(self) -> None:
        m = 16
        rho = np.ones(m)
        self.assertAlmostEqual(finite_correlation_factor(rho), float(m), places=13)
        self.assertAlmostEqual(effective_independent_count(rho), 1.0, places=13)
        self.assertAlmostEqual(
            variance_equivalent_axial_length(rho, 2.86e-10), m * 2.86e-10, places=23
        )

    def test_anticorrelation_can_raise_effective_count(self) -> None:
        rho = np.array([1.0, -0.2, 0.0, 0.0])
        self.assertLess(finite_correlation_factor(rho), 1.0)
        self.assertGreater(effective_independent_count(rho), float(len(rho)))

    def test_positive_window_estimator(self) -> None:
        rho_hat = np.array([1.0, 0.8, 0.4, -0.1, 0.2])
        expected = 1.0 + 2.0 * ((1.0 - 1.0 / 5.0) * 0.8 + (1.0 - 2.0 / 5.0) * 0.4)
        tau_hat = positive_window_empirical_correlation_factor(rho_hat)
        self.assertAlmostEqual(tau_hat, expected, places=14)
        self.assertAlmostEqual(
            positive_window_empirical_effective_count(rho_hat), 5.0 / expected, places=14
        )
        self.assertAlmostEqual(
            positive_window_empirical_axial_length(rho_hat, 2.0), 2.0 * expected, places=14
        )

    def test_positive_window_fully_identical_limit(self) -> None:
        m = 9
        rho_hat = np.ones(m)
        self.assertAlmostEqual(
            positive_window_empirical_correlation_factor(rho_hat), float(m), places=13
        )
        self.assertAlmostEqual(
            positive_window_empirical_effective_count(rho_hat), 1.0, places=13
        )

    def test_identical_pair_msd(self) -> None:
        x = np.array([0.9, 1.0, 1.1])
        self.assertEqual(identical_pair_msd(x, x.copy()), 0.0)
        self.assertGreater(identical_pair_msd(x, x + 1.0e-3), 0.0)

    def test_independent_any_event_probability(self) -> None:
        self.assertAlmostEqual(independent_any_event_probability(0.1, 3), 0.271)
        self.assertEqual(independent_any_event_probability(0.6, 0), 0.0)
        self.assertEqual(independent_any_event_probability(1.0, 5), 1.0)

    def test_identical_block_any_event_probability(self) -> None:
        q = 0.1
        m = 12
        self.assertAlmostEqual(
            identical_block_any_event_probability(q, m, 1),
            independent_any_event_probability(q, m),
        )
        self.assertAlmostEqual(
            identical_block_any_event_probability(q, m, m), q
        )
        self.assertAlmostEqual(
            identical_block_any_event_probability(q, m, 3),
            independent_any_event_probability(q, 4),
        )
        with self.assertRaises(ValueError):
            identical_block_any_event_probability(q, 10, 3)


if __name__ == "__main__":
    unittest.main()
