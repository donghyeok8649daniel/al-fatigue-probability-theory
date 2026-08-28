# === 한국어 파일 안내 시작 ===
# - 파일 역할: 비선형 안정점·장벽, 물리 에너지척도, metastable Gibbs P, finite-M fixed-length canonical P의 성질을 검증한다.
# - 주요 클래스: TestNormalLJPhysicalDistribution
# - 주요 함수/메서드: 에너지척도, stationary point, barrier, metastable density, canonical convolution 검증
# - 주의: 이 테스트는 명시된 equilibrium/metastable ensemble의 수학을 검증하며 실제 Al fatigue life를 검증하지 않는다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.normal_lj_chain import (
    critical_dimensionless_force,
    critical_stretch,
    normalized_lj_force,
    normalized_lj_stiffness,
)
from theory.normal_lj_physical_distribution import (
    BOLTZMANN_CONSTANT_J_PER_K,
    fixed_length_canonical_density,
    fixed_length_two_spacing_density,
    inverse_reduced_temperature,
    metastable_barrier_height,
    metastable_gibbs_density,
    metastable_stationary_points,
    metastable_tail_probability,
    physical_energy_scale,
    quasistatic_stable_spacing,
)


class TestNormalLJPhysicalDistribution(unittest.TestCase):
    def test_physical_energy_scale_and_inverse_temperature(self) -> None:
        E = 69.0e9
        A0 = 8.0e-20
        a0 = 2.86e-10
        T = 300.0
        e0 = physical_energy_scale(E, A0, a0)
        self.assertAlmostEqual(e0, E * A0 * a0, places=28)
        chi = inverse_reduced_temperature(E, A0, a0, T)
        self.assertAlmostEqual(
            chi, e0 / (BOLTZMANN_CONSTANT_J_PER_K * T), places=12
        )
        self.assertGreater(chi, 0.0)

    def test_quasistatic_zero_force_is_equilibrium_spacing(self) -> None:
        self.assertEqual(quasistatic_stable_spacing(0.0), 1.0)

    def test_metastable_stationary_points_have_correct_stability(self) -> None:
        fc = critical_dimensionless_force()
        f = 0.5 * fc
        stable, barrier = metastable_stationary_points(f)
        lam_c = critical_stretch()
        self.assertTrue(1.0 < stable < lam_c < barrier)
        self.assertAlmostEqual(float(normalized_lj_force(stable)), f, places=11)
        self.assertAlmostEqual(float(normalized_lj_force(barrier)), f, places=11)
        self.assertGreater(float(normalized_lj_stiffness(stable)), 0.0)
        self.assertLess(float(normalized_lj_stiffness(barrier)), 0.0)

    def test_barrier_is_positive_and_decreases_toward_critical_force(self) -> None:
        fc = critical_dimensionless_force()
        forces = [0.2 * fc, 0.5 * fc, 0.8 * fc, 0.95 * fc]
        barriers = [metastable_barrier_height(f) for f in forces]
        self.assertTrue(all(value > 0.0 for value in barriers))
        self.assertTrue(all(a > b for a, b in zip(barriers, barriers[1:])))

    def test_metastable_density_is_normalized_and_confined_to_basin(self) -> None:
        fc = critical_dimensionless_force()
        f = 0.6 * fc
        _, barrier = metastable_stationary_points(f)
        grid = np.linspace(0.72, barrier * 1.08, 6001)
        density = metastable_gibbs_density(grid, f, 250.0)
        self.assertAlmostEqual(float(np.trapezoid(density, grid)), 1.0, places=9)
        self.assertTrue(np.all(density[grid >= barrier] == 0.0))
        tail = metastable_tail_probability(grid, density)
        self.assertTrue(0.0 <= tail <= 1.0)

    def test_higher_chi_concentrates_metastable_density_near_stable_point(self) -> None:
        fc = critical_dimensionless_force()
        f = 0.5 * fc
        stable, barrier = metastable_stationary_points(f)
        grid = np.linspace(0.78, barrier * 0.9999, 8001)
        p_lo = metastable_gibbs_density(grid, f, 40.0)
        p_hi = metastable_gibbs_density(grid, f, 400.0)
        mean_lo = float(np.trapezoid(grid * p_lo, grid))
        mean_hi = float(np.trapezoid(grid * p_hi, grid))
        var_lo = float(np.trapezoid((grid - mean_lo) ** 2 * p_lo, grid))
        var_hi = float(np.trapezoid((grid - mean_hi) ** 2 * p_hi, grid))
        self.assertLess(var_hi, var_lo)
        self.assertLess(abs(mean_hi - stable), abs(mean_lo - stable))

    def test_fixed_length_two_spacing_density_is_normalized_and_symmetric(self) -> None:
        total = 2.04
        grid = np.linspace(0.30, total - 0.30, 7001)
        density = fixed_length_two_spacing_density(grid, total, 120.0)
        self.assertAlmostEqual(float(np.trapezoid(density, grid)), 1.0, places=9)
        reflected = np.interp(total - grid, grid, density)
        self.assertLess(float(np.max(np.abs(density - reflected))), 2.0e-10)

    def test_convolution_density_matches_exact_M2_formula(self) -> None:
        total = 2.04
        grid = np.linspace(0.0, total, 8193)
        direct = fixed_length_two_spacing_density(grid, total, 120.0)
        convolution = fixed_length_canonical_density(
            grid, total, spacing_count=2, inverse_temperature=120.0
        )
        self.assertLess(float(np.max(np.abs(direct - convolution))), 2.0e-9)

    def test_finite_M_canonical_mean_is_total_length_over_M(self) -> None:
        spacing_count = 5
        total = 5.10
        grid = np.linspace(0.0, total, 16385)
        density = fixed_length_canonical_density(
            grid, total, spacing_count=spacing_count, inverse_temperature=100.0
        )
        self.assertAlmostEqual(float(np.trapezoid(density, grid)), 1.0, places=8)
        mean = float(np.trapezoid(grid * density, grid))
        self.assertAlmostEqual(mean, total / spacing_count, delta=3.0e-4)

    def test_metastable_force_domain_is_enforced(self) -> None:
        fc = critical_dimensionless_force()
        with self.assertRaises(ValueError):
            metastable_stationary_points(0.0)
        with self.assertRaises(ValueError):
            metastable_stationary_points(fc)


if __name__ == "__main__":
    unittest.main()
