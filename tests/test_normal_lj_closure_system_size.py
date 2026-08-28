# === 한국어 파일 안내 시작 ===
# - 파일 역할: system-size sweep의 dynamic-similarity 규칙과 수치 진단 함수를 검증한다.
# - 주요 클래스: NormalLJClosureSystemSizeTests
# - 주요 함수/메서드: NormalLJClosureSystemSizeTests.test_dynamic_similarity_invariant
#   NormalLJClosureSystemSizeTests.test_sharply_concentrated_state_is_resolved
#   NormalLJClosureSystemSizeTests.test_vectorized_cdf_matches_scalar
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from simulations.run_normal_lj_closure_system_size import (
    OMEGA_TIMES_M,
    dynamically_matched_omega,
)
from theory.normal_lj_closure_validation import closure_cdf, closure_cdf_many
from theory.normal_lj_distribution import solve_distribution_closure


class NormalLJClosureSystemSizeTests(unittest.TestCase):
    def test_dynamic_similarity_invariant(self):
        for represented_spacings in (31, 63, 127, 255):
            omega = dynamically_matched_omega(represented_spacings)
            self.assertAlmostEqual(
                omega * represented_spacings,
                OMEGA_TIMES_M,
                places=14,
            )

    def test_sharply_concentrated_state_is_resolved(self):
        target_mu = 0.9885611826976782
        target_energy = 8.704614514402719e-05

        s320 = solve_distribution_closure(
            target_mu,
            target_energy,
            quadrature_order=320,
        )
        s640 = solve_distribution_closure(
            target_mu,
            target_energy,
            quadrature_order=640,
        )

        self.assertAlmostEqual(s640.moments.mean_stretch, target_mu, places=10)
        self.assertAlmostEqual(s640.moments.mean_energy, target_energy, places=12)
        self.assertLess(
            abs(s320.moments.variance_stretch - s640.moments.variance_stretch),
            1.0e-9,
        )

    def test_vectorized_cdf_matches_scalar(self):
        solution = solve_distribution_closure(
            1.0,
            1.0e-3,
            quadrature_order=320,
        )
        points = np.asarray([0.95, 1.0, 1.05, 1.10])
        vectorized = closure_cdf_many(
            points,
            solution,
            quadrature_order=96,
        )
        scalar = np.asarray(
            [
                closure_cdf(x, solution, quadrature_order=96)
                for x in points
            ]
        )
        np.testing.assert_allclose(vectorized, scalar, rtol=0.0, atol=1.0e-13)


if __name__ == "__main__":
    unittest.main()
