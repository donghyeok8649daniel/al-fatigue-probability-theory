import unittest
import numpy as np

from theory.normal_lj_spatial_correlation import (
    correlation_profile,
    normalized_spatial_correlation,
    random_permutation_expected_rho,
    spatial_covariance,
)


class SpatialCorrelationTests(unittest.TestCase):
    def test_c0_is_empirical_variance(self):
        values = np.asarray([0.98, 1.01, 1.03, 0.99, 1.02])
        self.assertAlmostEqual(spatial_covariance(values, 0), float(np.var(values)), places=15)

    def test_one_point_permutation_invariance_does_not_fix_correlation(self):
        values = np.asarray([0.96, 0.98, 1.00, 1.02, 1.04, 1.06])
        permuted = values[[0, 5, 1, 4, 2, 3]]
        self.assertAlmostEqual(float(np.mean(values)), float(np.mean(permuted)), places=15)
        self.assertAlmostEqual(float(np.mean(values ** 2)), float(np.mean(permuted ** 2)), places=15)
        self.assertNotAlmostEqual(
            normalized_spatial_correlation(values, 1),
            normalized_spatial_correlation(permuted, 1),
            places=6,
        )

    def test_random_permutation_expectation(self):
        self.assertAlmostEqual(random_permutation_expected_rho(31), -1.0 / 30.0)
        self.assertAlmostEqual(random_permutation_expected_rho(255), -1.0 / 254.0)

    def test_profile_starts_at_unity(self):
        values = np.asarray([0.99, 1.00, 1.02, 1.01, 0.98])
        _, _, rho = correlation_profile(values)
        self.assertAlmostEqual(float(rho[0]), 1.0, places=15)


if __name__ == "__main__":
    unittest.main()
