# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: TestNormalLJProbabilityDynamics
# - 주요 함수/메서드: 
#   TestNormalLJProbabilityDynamics.test_constant_load_preserves_normalization_and_local_equilibrium
#   TestNormalLJProbabilityDynamics.test_finite_relaxation_creates_positive_closed_cycle_area
#   TestNormalLJProbabilityDynamics.test_physical_mean_spacing_is_obtained_by_multiplying_a0
#   TestNormalLJProbabilityDynamics.test_invalid_time_history_is_rejected
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.normal_lj_probability_dynamics import (
    SpacingDynamicsParameters,
    completed_cycle_hysteresis_areas,
    solve_spacing_probability_history,
)


class TestNormalLJProbabilityDynamics(unittest.TestCase):
    def test_constant_load_preserves_normalization_and_local_equilibrium(self) -> None:
        time = np.linspace(0.0, 0.2, 21)
        stress = np.full_like(time, 80.0e6)
        result = solve_spacing_probability_history(
            time,
            stress,
            69.0e9,
            SpacingDynamicsParameters(
                grid_cells=120,
                relaxation_time_s=0.02,
                substeps_per_interval=1,
            ),
        )
        np.testing.assert_allclose(result.normalization, 1.0, rtol=0.0, atol=2.0e-12)
        np.testing.assert_allclose(result.density[0], result.density[-1], rtol=2.0e-5, atol=2.0e-8)
        self.assertLess(abs(result.cumulative_hysteresis_energy_density_j_m3[-1]), 1.0e-5)

    def test_finite_relaxation_creates_positive_closed_cycle_area(self) -> None:
        frequency = 5.0
        steps_per_cycle = 80
        cycles = 3
        time = np.linspace(0.0, cycles / frequency, cycles * steps_per_cycle + 1)
        stress = 120.0e6 + 80.0e6 * np.sin(2.0 * np.pi * frequency * time)
        result = solve_spacing_probability_history(
            time,
            stress,
            69.0e9,
            SpacingDynamicsParameters(
                grid_cells=140,
                relaxation_time_s=0.03,
                substeps_per_interval=1,
            ),
        )
        areas = completed_cycle_hysteresis_areas(result, frequency)
        self.assertEqual(areas.size, cycles)
        self.assertGreater(areas[-1], 0.0)
        self.assertLess(abs(areas[-1] - areas[-2]) / areas[-1], 0.08)

    def test_physical_mean_spacing_is_obtained_by_multiplying_a0(self) -> None:
        time = np.array([0.0, 0.01])
        stress = np.array([50.0e6, 60.0e6])
        result = solve_spacing_probability_history(
            time,
            stress,
            69.0e9,
            SpacingDynamicsParameters(grid_cells=80),
        )
        a0 = 2.86e-10
        mean_spacing = a0 * result.mean_stretch
        self.assertTrue(np.all(mean_spacing > 0.0))
        np.testing.assert_allclose(mean_spacing / a0, result.mean_stretch)

    def test_invalid_time_history_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            solve_spacing_probability_history(
                np.array([0.0, 0.0]),
                np.array([1.0, 1.0]),
                69.0e9,
            )


if __name__ == "__main__":
    unittest.main()
