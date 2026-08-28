# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D layer-LJ chain의 equilibrium, instability, 보존성 및 기본 수치 동작을 회귀검증한다.
# - 주요 클래스: NormalLJChainTests
# - 주요 함수/메서드: NormalLJChainTests.test_normalization, NormalLJChainTests.test_critical_point
#   NormalLJChainTests.test_100mpa_is_reversible_null_case
#   NormalLJChainTests.test_sub_static_critical_dynamic_crossing_exists
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

from theory.normal_lj_chain import (
    NormalLJParameters,
    critical_dimensionless_force,
    critical_stretch,
    normalized_lj_force,
    normalized_lj_stiffness,
    simulate_normal_lj_chain,
    stress_to_dimensionless_force,
)


class NormalLJChainTests(unittest.TestCase):
    def test_normalization(self):
        self.assertAlmostEqual(float(normalized_lj_force(1.0)), 0.0, places=14)
        self.assertAlmostEqual(float(normalized_lj_stiffness(1.0)), 1.0, places=12)

    def test_critical_point(self):
        lam_c = critical_stretch()
        self.assertAlmostEqual(
            float(normalized_lj_stiffness(lam_c)),
            0.0,
            places=12,
        )
        self.assertAlmostEqual(
            critical_dimensionless_force(),
            0.03703426967076833,
            places=12,
        )

    def test_100mpa_is_reversible_null_case(self):
        amplitude = stress_to_dimensionless_force(100.0e6, 69.0e9)
        result = simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=amplitude, omega=0.02),
            atoms=24,
            cycles=4,
        )
        self.assertIsNone(result.first_instability)
        self.assertLess(result.energy_balance_relative_error, 1.0e-7)

    def test_sub_static_critical_dynamic_crossing_exists(self):
        result = simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=0.03, omega=0.02),
            atoms=32,
            cycles=3,
        )
        self.assertIsNotNone(result.first_instability)
        self.assertGreater(result.first_instability.cycle, 2.0)
        self.assertLess(result.first_instability.cycle, 2.5)


if __name__ == "__main__":
    unittest.main()
