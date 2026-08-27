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
