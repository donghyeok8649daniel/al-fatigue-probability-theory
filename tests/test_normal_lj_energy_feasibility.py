import unittest

import numpy as np

from theory.normal_lj_chain import critical_stretch
from theory.normal_lj_energy_feasibility import (
    first_energy_ceiling_crossing_time,
    no_compression_bound_counterexample_energy,
    safe_distribution_exists,
    safe_energy_interval,
    shifted_lj_energy,
)


class NormalLJEnergyFeasibilityTests(unittest.TestCase):
    def test_endpoint_measure_attains_exact_maximum(self):
        lower = 0.95
        mu = 1.0
        interval = safe_energy_interval(mu, lower)
        direct = (
            interval.lower_endpoint_weight * shifted_lj_energy(lower)
            + interval.critical_endpoint_weight
            * shifted_lj_energy(interval.critical_stretch)
        )
        self.assertAlmostEqual(float(direct), interval.maximum_energy, places=14)

    def test_delta_at_mean_attains_exact_minimum(self):
        interval = safe_energy_interval(1.01, 0.95)
        self.assertAlmostEqual(
            interval.minimum_energy,
            float(shifted_lj_energy(1.01)),
            places=14,
        )

    def test_random_safe_measures_stay_below_chord_bound(self):
        rng = np.random.default_rng(20260828)
        lower = 0.95
        upper = critical_stretch()

        for _ in range(200):
            x = np.sort(rng.uniform(lower, upper, size=20))
            w = rng.random(20)
            w /= np.sum(w)
            mu = float(np.dot(w, x))
            energy = float(np.dot(w, shifted_lj_energy(x)))
            interval = safe_energy_interval(mu, lower)
            self.assertLessEqual(energy, interval.maximum_energy + 1.0e-12)
            self.assertGreaterEqual(energy, interval.minimum_energy - 1.0e-12)

    def test_energy_above_ceiling_has_no_safe_distribution(self):
        interval = safe_energy_interval(1.0, 0.95)
        self.assertFalse(
            safe_distribution_exists(
                interval.maximum_energy + 1.0e-6,
                1.0,
                0.95,
            )
        )

    def test_no_lower_compression_bound_allows_unbounded_safe_energy(self):
        e1 = no_compression_bound_counterexample_energy(0.8, 1.0)
        e2 = no_compression_bound_counterexample_energy(0.5, 1.0)
        e3 = no_compression_bound_counterexample_energy(0.3, 1.0)
        self.assertLess(e1, e2)
        self.assertLess(e2, e3)
        self.assertGreater(e3, 100.0)

    def test_first_crossing_is_continuous_time_not_cycle_index(self):
        times = [0.0, 0.1, 0.2, 0.3]
        mu = [1.0] * 4
        lower = [0.95] * 4
        ceiling = safe_energy_interval(1.0, 0.95).maximum_energy
        energies = [0.0, 0.5 * ceiling, 0.99 * ceiling, 1.01 * ceiling]
        self.assertEqual(
            first_energy_ceiling_crossing_time(times, energies, mu, lower),
            0.3,
        )


if __name__ == "__main__":
    unittest.main()
