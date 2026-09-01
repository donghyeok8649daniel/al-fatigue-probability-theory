import unittest

import numpy as np

from theory.fcc_normal_lj import EV_J, FCCNormalLJ, FCCNormalLJParameters


class FCCNormalLJTests(unittest.TestCase):
    def setUp(self):
        self.model = FCCNormalLJ(
            FCCNormalLJParameters(cutoff_lattice_constants=8.0)
        )

    def test_cohesive_energy_calibration(self):
        target = 3.43 * EV_J
        epsilon = self.model.epsilon_for_cohesive_energy(target)
        recovered = self.model.cohesive_energy_j_per_atom(epsilon)
        self.assertAlmostEqual(recovered / target, 1.0, places=12)

    def test_reference_state_has_zero_axial_stress(self):
        epsilon = self.model.epsilon_for_cohesive_energy(3.43 * EV_J)
        stress = self.model.axial_engineering_stress_pa(epsilon, 1.0, 1.0)
        self.assertLess(abs(stress), 2.0e6)

    def test_central_pair_cauchy_relation(self):
        epsilon = self.model.epsilon_for_youngs_modulus(62.702380952380956e9)
        _, C12, C44 = self.model.cubic_elastic_constants_pa(epsilon)
        self.assertLess(abs(C12 - C44) / C12, 2.0e-5)

    def test_normal_fit_exposes_cohesive_energy_gap(self):
        epsilon = self.model.epsilon_for_youngs_modulus(62.702380952380956e9)
        cohesive_ev = self.model.cohesive_energy_j_per_atom(epsilon) / EV_J
        self.assertGreater(cohesive_ev, 0.9)
        self.assertLess(cohesive_ev, 1.1)

    def test_normal_fit_predicts_order_9_gpa_ideal_strength(self):
        epsilon = self.model.epsilon_for_youngs_modulus(62.702380952380956e9)
        stretch = np.linspace(1.0, 1.4, 81)
        _, _, stress = self.model.stress_strain_curve(epsilon, stretch)
        self.assertGreater(np.max(stress), 8.5e9)
        self.assertLess(np.max(stress), 9.6e9)


if __name__ == "__main__":
    unittest.main()
