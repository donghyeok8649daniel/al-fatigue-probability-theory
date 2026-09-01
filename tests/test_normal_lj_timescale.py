# === 한국어 파일 안내 시작 ===
# - 파일 역할: 무차원·물리 시간/주파수 변환과 scale 계산을 검증한다.
# - 주요 클래스: NormalLJTimescaleTests
# - 주요 함수/메서드: NormalLJTimescaleTests.test_lowest_mode_decreases_with_chain_length
#   NormalLJTimescaleTests.test_mode_inversion
#   NormalLJTimescaleTests.test_100mpa_stable_branch_has_positive_tangent_stiffness
#   NormalLJTimescaleTests.test_20hz_critical_softening_requires_extreme_proximity
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import math
import unittest

from theory.normal_lj_chain import critical_stretch, normalized_lj_stiffness
from theory.normal_lj_timescale import (
    homogeneous_stretch_for_dimensionless_stress,
    lowest_fixed_free_mode_omega_star,
    moving_atoms_for_target_frequency,
    near_critical_distance_for_target_local_frequency,
)


class NormalLJTimescaleTests(unittest.TestCase):
    def test_lowest_mode_decreases_with_chain_length(self):
        self.assertGreater(
            lowest_fixed_free_mode_omega_star(32),
            lowest_fixed_free_mode_omega_star(320),
        )

    def test_mode_inversion(self):
        t0 = 5.550462614661804e-14
        atoms = moving_atoms_for_target_frequency(20.0, t0)
        omega = lowest_fixed_free_mode_omega_star(atoms)
        recovered = omega / (2.0 * math.pi * t0)
        self.assertAlmostEqual(recovered, 20.0, places=9)

    def test_100mpa_stable_branch_has_positive_tangent_stiffness(self):
        lam = homogeneous_stretch_for_dimensionless_stress(100.0e6 / 69.0e9)
        self.assertLess(lam, critical_stretch())
        self.assertGreater(float(normalized_lj_stiffness(lam)), 0.9)

    def test_20hz_critical_softening_requires_extreme_proximity(self):
        t0 = 5.550462614661804e-14
        delta = near_critical_distance_for_target_local_frequency(20.0, t0)
        self.assertLess(delta, 1.0e-20)


if __name__ == "__main__":
    unittest.main()
