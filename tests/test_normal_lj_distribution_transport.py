# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: TestNormalLJDistributionTransport
# - 주요 함수/메서드: TestNormalLJDistributionTransport.test_exact_spacing_acceleration
#   TestNormalLJDistributionTransport.test_moment_rate_identities
#   TestNormalLJDistributionTransport.test_boundary_not_silently_closed
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.normal_lj_chain import normalized_lj_force
from theory.normal_lj_distribution_transport import (
    exact_spacing_acceleration,
    monomial_moment,
    monomial_moment_rate,
    monomial_moment_second_rate,
)


class TestNormalLJDistributionTransport(unittest.TestCase):
    def test_exact_spacing_acceleration(self):
        lam = np.array([0.98, 1.01, 1.03, 0.99, 1.02], dtype=float)
        g = normalized_lj_force(lam)
        acc = exact_spacing_acceleration(lam)
        expected = g[2:] - 2.0 * g[1:-1] + g[:-2]
        np.testing.assert_allclose(acc[1:-1], expected, rtol=0.0, atol=1e-14)

    def test_moment_rate_identities(self):
        lam = np.array([0.97, 1.00, 1.04, 1.02], dtype=float)
        vel = np.array([0.02, -0.01, 0.03, -0.02], dtype=float)
        acc = np.array([0.1, -0.2, 0.05, 0.03], dtype=float)

        self.assertAlmostEqual(monomial_moment(lam, 0), 1.0)
        self.assertAlmostEqual(
            monomial_moment_rate(lam, vel, 2),
            2.0 * float(np.mean(lam * vel)),
        )
        self.assertAlmostEqual(
            monomial_moment_second_rate(lam, vel, acc, 2),
            2.0 * float(np.mean(vel * vel)) + 2.0 * float(np.mean(lam * acc)),
        )

    def test_boundary_not_silently_closed(self):
        lam = np.array([1.0, 1.01, 0.99, 1.0], dtype=float)
        acc = exact_spacing_acceleration(lam)
        self.assertTrue(np.isnan(acc[0]))
        self.assertTrue(np.isnan(acc[-1]))


if __name__ == "__main__":
    unittest.main()
