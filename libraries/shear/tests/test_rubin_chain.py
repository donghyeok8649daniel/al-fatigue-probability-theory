import unittest

import numpy as np

from theory.rubin_chain import (
    RubinParams,
    analytic_response,
    cycle_loop_areas,
    simulate_finite_chain,
)


class RubinChainTests(unittest.TestCase):
    def test_passband_has_positive_hysteresis(self) -> None:
        response = analytic_response(0.5, 0.1, RubinParams())
        self.assertGreater(response["loop_area"], 0.0)
        self.assertGreater(response["z_imag"], 0.0)

    def test_above_band_has_zero_radiation_hysteresis(self) -> None:
        response = analytic_response(2.5, 0.1, RubinParams())
        self.assertAlmostEqual(response["z_imag"], 0.0, places=14)
        self.assertAlmostEqual(response["loop_area"], 0.0, places=14)

    def test_full_newton_chain_matches_analytic_area_and_energy_balance(self) -> None:
        params = RubinParams()
        analytic = analytic_response(0.5, 0.1, params)
        numeric = simulate_finite_chain(
            n_masses=600,
            omega=0.5,
            force_amplitude=0.1,
            dt=0.02,
            n_periods=25,
            ramp_periods=3,
            params=params,
        )
        areas = cycle_loop_areas(
            numeric,
            first_cycle=8,
            last_cycle_exclusive=20,
        )
        mean_area = float(np.mean(areas))

        area_rel_error = abs(mean_area - analytic["loop_area"]) / analytic["loop_area"]
        self.assertLess(area_rel_error, 5.0e-4)

        final_energy = float(np.asarray(numeric["energy"])[-1])
        final_work = float(np.asarray(numeric["work"])[-1])
        energy_rel_error = abs(final_energy - final_work) / abs(final_work)
        self.assertLess(energy_rel_error, 5.0e-4)


if __name__ == "__main__":
    unittest.main()
