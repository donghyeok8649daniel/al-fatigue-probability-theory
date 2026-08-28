# === 한국어 파일 안내 시작 ===
# - 파일 역할: closure-vs-mechanics 비교 지표와 반증시험 계산이 재현되는지 검증한다.
# - 주요 클래스: ClosureFalsificationTests
# - 주요 함수/메서드: ClosureFalsificationTests.setUpClass
#   ClosureFalsificationTests.test_variance_is_close_in_near_equilibrium_case
#   ClosureFalsificationTests.test_full_distribution_is_not_identical
#   ClosureFalsificationTests.test_tail_remains_effectively_zero_before_instability
#   ClosureFalsificationTests.test_skewness_is_not_forced_by_mean_and_energy
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

from theory.normal_lj_chain import NormalLJParameters, simulate_normal_lj_chain
from theory.normal_lj_closure_validation import compare_snapshot_to_closure


class ClosureFalsificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        result = simulate_normal_lj_chain(
            NormalLJParameters(force_amplitude=0.03, omega=0.02),
            cycles=3,
        )
        cls.values = result.cycle_snapshots[2]
        cls.comparison = compare_snapshot_to_closure(
            cls.values,
            closure_quadrature_order=320,
            cdf_quadrature_order=96,
        )

    def test_variance_is_close_in_near_equilibrium_case(self):
        self.assertLess(self.comparison.variance_relative_error, 0.03)

    def test_full_distribution_is_not_identical(self):
        self.assertGreater(self.comparison.kolmogorov_distance, 0.10)

    def test_tail_remains_effectively_zero_before_instability(self):
        self.assertEqual(self.comparison.empirical_critical_tail_probability, 0.0)
        self.assertLess(self.comparison.closure_critical_tail_probability, 1.0e-20)

    def test_skewness_is_not_forced_by_mean_and_energy(self):
        self.assertGreater(
            abs(self.comparison.empirical_skewness - self.comparison.closure_skewness),
            0.05,
        )


if __name__ == "__main__":
    unittest.main()
