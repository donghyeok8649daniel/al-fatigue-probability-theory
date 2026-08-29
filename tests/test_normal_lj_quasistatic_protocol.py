# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: QuasistaticProtocolTests
# - 주요 함수/메서드: QuasistaticProtocolTests.test_zero_force_stable_stretch_is_equilibrium
#   QuasistaticProtocolTests.test_stable_root_recovers_force
#   QuasistaticProtocolTests.test_critical_force_maps_to_critical_stretch
#   QuasistaticProtocolTests.test_static_open_chain_is_exactly_homogeneous
#   QuasistaticProtocolTests.test_zero_mean_sine_has_zero_force_at_integer_cycle
#   QuasistaticProtocolTests.test_nonzero_mean_is_multiplied_by_ramp_envelope
#   QuasistaticProtocolTests.test_zero_variance_snapshot_does_not_create_artificial_cell_length
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.normal_lj_chain import (
    NormalLJParameters,
    critical_dimensionless_force,
    critical_stretch,
    normalized_lj_force,
)
from theory.normal_lj_quasistatic_protocol import (
    cycle_boundary_force,
    quasistatic_open_chain_spacings,
    residual_snapshot_metrics,
    stable_stretch_for_tensile_force,
)


class QuasistaticProtocolTests(unittest.TestCase):
    def test_zero_force_stable_stretch_is_equilibrium(self):
        self.assertEqual(stable_stretch_for_tensile_force(0.0), 1.0)

    def test_stable_root_recovers_force(self):
        for force in (0.001, 0.01, 0.02, 0.03):
            lam = stable_stretch_for_tensile_force(force)
            self.assertGreater(lam, 1.0)
            self.assertLess(lam, critical_stretch())
            self.assertAlmostEqual(float(normalized_lj_force(lam)), force, places=11)

    def test_critical_force_maps_to_critical_stretch(self):
        self.assertAlmostEqual(
            stable_stretch_for_tensile_force(critical_dimensionless_force()),
            critical_stretch(),
            places=12,
        )

    def test_static_open_chain_is_exactly_homogeneous(self):
        values = quasistatic_open_chain_spacings(0.02, 31)
        self.assertEqual(values.shape, (31,))
        self.assertEqual(float(np.var(values)), 0.0)
        self.assertTrue(np.all(values == values[0]))

    def test_zero_mean_sine_has_zero_force_at_integer_cycle(self):
        p = NormalLJParameters(mean_force=0.0, force_amplitude=0.03, omega=0.02, ramp_cycles=2)
        for cycle in range(5):
            self.assertEqual(cycle_boundary_force(p, cycle), 0.0)

    def test_nonzero_mean_is_multiplied_by_ramp_envelope(self):
        p = NormalLJParameters(mean_force=0.01, force_amplitude=0.03, omega=0.02, ramp_cycles=2)
        self.assertEqual(cycle_boundary_force(p, 0), 0.0)
        self.assertAlmostEqual(cycle_boundary_force(p, 1), 0.005, places=14)
        self.assertAlmostEqual(cycle_boundary_force(p, 2), 0.01, places=14)

    def test_zero_variance_snapshot_does_not_create_artificial_cell_length(self):
        values = np.ones(20)
        metrics = residual_snapshot_metrics(values, quasistatic_stretch=1.0)
        self.assertEqual(metrics.variance_c0, 0.0)
        self.assertEqual(metrics.rms_nonuniformity, 0.0)
        self.assertEqual(metrics.rho1, 0.0)
        self.assertEqual(metrics.tau_positive_window, 1.0)
        self.assertEqual(metrics.m_eff_positive_window, 20.0)


if __name__ == "__main__":
    unittest.main()
