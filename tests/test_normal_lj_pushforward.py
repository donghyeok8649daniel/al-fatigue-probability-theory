# === 한국어 파일 안내 시작 ===
# - 파일 역할: push-forward 보조 계산의 수학적 identity를 검증한다. harmonic/Taylor 항목은 과거 진단용 테스트로만 유지한다.
# - 주요 클래스: PushForwardDistributionTests
# - 주요 함수/메서드: PushForwardDistributionTests.test_lj_force_taylor_coefficients
#   PushForwardDistributionTests.test_linear_dispersion
#   PushForwardDistributionTests.test_arcsine_normalization_numerically
#   PushForwardDistributionTests.test_arcsine_cdf_endpoints
#   PushForwardDistributionTests.test_single_mode_moments
#   PushForwardDistributionTests.test_two_harmonic_skewness_formula
#   PushForwardDistributionTests.test_two_harmonic_maximum_skewness
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import math
import unittest

import numpy as np

from theory.normal_lj_pushforward import (
    arcsine_cdf,
    arcsine_density,
    lj_force_taylor_coefficients,
    linear_mode_frequency,
    single_mode_moments,
    two_harmonic_max_abs_skewness,
    two_harmonic_moments,
)


class PushForwardDistributionTests(unittest.TestCase):
    def test_lj_force_taylor_coefficients(self):
        c1, c2, c3 = lj_force_taylor_coefficients()
        self.assertAlmostEqual(c1, 1.0, places=14)
        self.assertAlmostEqual(c2, -10.595, places=12)
        self.assertAlmostEqual(c3, 62.97935, places=10)

    def test_linear_dispersion(self):
        q = 0.31
        self.assertAlmostEqual(
            linear_mode_frequency(q) ** 2,
            4.0 * math.sin(0.5 * q) ** 2,
            places=14,
        )

    def test_arcsine_normalization_numerically(self):
        mu = 1.0
        A = 0.03
        theta = (np.arange(500000, dtype=float) + 0.5) * 2.0 * math.pi / 500000
        values = mu + A * np.cos(theta)
        self.assertAlmostEqual(float(np.mean(values)), mu, places=12)
        self.assertAlmostEqual(float(np.var(values)), 0.5 * A * A, places=11)

    def test_arcsine_cdf_endpoints(self):
        mu, A = 1.0, 0.02
        self.assertEqual(float(arcsine_cdf(mu - A, mu, A)), 0.0)
        self.assertEqual(float(arcsine_cdf(mu + A, mu, A)), 1.0)
        self.assertAlmostEqual(float(arcsine_cdf(mu, mu, A)), 0.5, places=14)

    def test_single_mode_moments(self):
        moments = single_mode_moments(0.04)
        self.assertAlmostEqual(moments["variance"], 0.0008, places=14)
        self.assertEqual(moments["skewness"], 0.0)
        self.assertAlmostEqual(moments["kurtosis"], 1.5, places=14)

    def test_two_harmonic_skewness_formula(self):
        A, B = 0.01, 0.003
        moments = two_harmonic_moments(A, B)
        expected = (0.75 * A * A * B) / (0.5 * (A * A + B * B)) ** 1.5
        self.assertAlmostEqual(moments["skewness"], expected, places=14)

    def test_two_harmonic_maximum_skewness(self):
        self.assertAlmostEqual(
            two_harmonic_max_abs_skewness(),
            math.sqrt(2.0 / 3.0),
            places=14,
        )


if __name__ == "__main__":
    unittest.main()
