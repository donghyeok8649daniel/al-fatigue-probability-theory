# === 한국어 파일 안내 시작 ===
# - 파일 역할: historical spacing closure의 normalization, moment recovery, energy relation 및 수치 안정성을 검증한다.
# - 주요 클래스: NormalLJDistributionTests
# - 주요 함수/메서드: NormalLJDistributionTests.test_solver_recovers_mean_and_energy
#   NormalLJDistributionTests.test_energy_decreases_with_beta_at_fixed_mean
#   NormalLJDistributionTests.test_tail_increases_over_reference_energy_sweep
#   NormalLJDistributionTests.test_homogeneous_energy_is_jensen_minimum
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

from theory.normal_lj_distribution import (
    energy_derivative_at_fixed_mean,
    shifted_lj_energy,
    solve_distribution_closure,
)


class NormalLJDistributionTests(unittest.TestCase):
    def test_solver_recovers_mean_and_energy(self):
        target_mu = 1.0
        target_e = 1.0e-3
        sol = solve_distribution_closure(target_mu, target_e, quadrature_order=480)
        self.assertAlmostEqual(sol.moments.mean_stretch, target_mu, places=8)
        self.assertAlmostEqual(sol.moments.mean_energy, target_e, places=8)

    def test_energy_decreases_with_beta_at_fixed_mean(self):
        sol = solve_distribution_closure(1.0, 1.0e-3, quadrature_order=480)
        self.assertLess(energy_derivative_at_fixed_mean(sol.moments), 0.0)

    def test_tail_increases_over_reference_energy_sweep(self):
        energies = [2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3]
        tails = [
            solve_distribution_closure(1.0, e, quadrature_order=480)
            .moments.critical_tail_probability
            for e in energies
        ]
        self.assertTrue(all(b > a for a, b in zip(tails, tails[1:])))

    def test_homogeneous_energy_is_jensen_minimum(self):
        self.assertAlmostEqual(float(shifted_lj_energy(1.0)), 0.0, places=14)


if __name__ == "__main__":
    unittest.main()
