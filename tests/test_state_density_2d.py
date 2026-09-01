# === 한국어 파일 안내 시작 ===
# - 파일 역할: theory-core의 (a,s) exact density log-gradient, compatibility curl, 경로복원, G2 평균에너지 적분을 검증한다.
# - 주의: analytic synthetic field를 사용한 수학/수치 검증이며 실제 Al 분포 검증이 아니다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.state_density_2d import (
    compatibility_curl_2d,
    density_log_gradient_2d,
    mean_intrinsic_energy_2d,
    probability_mass_2d,
    reconstruct_density_path_2d,
)


class TestStateDensity2D(unittest.TestCase):
    def setUp(self) -> None:
        self.a = np.linspace(0.86, 1.22, 121)
        self.s = np.linspace(-0.55, 0.55, 151)
        A, S = np.meshgrid(self.a, self.s, indexing="ij")
        self.A, self.S = A, S

        logp = (
            -55.0 * (A - 1.025) ** 2
            -10.0 * (S - 0.28 * (A - 1.0)) ** 2
            -1.5 * S**4
            +1.2 * (A - 1.0) * S
        )
        target = np.exp(logp - float(np.max(logp)))
        self.target = target / probability_mass_2d(self.a, self.s, target)

        self.theta_aa = 0.020 + 0.004 * (A - 1.0) ** 2 + 0.0015 * S**2
        self.theta_ss = 0.035 + 0.003 * S**2 + 0.0010 * (A - 1.0) ** 2
        self.theta_as = 0.0045 + 0.0008 * (A - 1.0) * S

        grad_a = np.gradient(np.log(self.target), self.a, axis=0, edge_order=2)
        grad_s = np.gradient(np.log(self.target), self.s, axis=1, edge_order=2)
        dtu_a = (
            0.003
            * np.sin(2.0 * np.pi * (A - self.a[0]) / (self.a[-1] - self.a[0]))
            * np.cos(np.pi * S / 0.55)
        )
        dtu_s = (
            0.002
            * np.cos(2.0 * np.pi * (A - self.a[0]) / (self.a[-1] - self.a[0]))
            * np.sin(np.pi * S / 0.55)
        )
        div_a = (
            np.gradient(self.theta_aa, self.a, axis=0, edge_order=2)
            + np.gradient(self.theta_as, self.s, axis=1, edge_order=2)
        )
        div_s = (
            np.gradient(self.theta_as, self.a, axis=0, edge_order=2)
            + np.gradient(self.theta_ss, self.s, axis=1, edge_order=2)
        )
        self.dtu_a, self.dtu_s = dtu_a, dtu_s
        self.acc_a = dtu_a + div_a + self.theta_aa * grad_a + self.theta_as * grad_s
        self.acc_s = dtu_s + div_s + self.theta_as * grad_a + self.theta_ss * grad_s

    def test_nonzero_cross_covariance_reconstructs_density(self) -> None:
        grad_a, grad_s = density_log_gradient_2d(
            self.a, self.s,
            self.theta_aa, self.theta_as, self.theta_ss,
            self.acc_a, self.acc_s,
            self.dtu_a, self.dtu_s,
        )
        recovered, mismatch = reconstruct_density_path_2d(
            self.a, self.s, grad_a, grad_s
        )
        l1 = probability_mass_2d(self.a, self.s, np.abs(recovered - self.target))
        self.assertLess(mismatch, 1.0e-8)
        self.assertLess(l1, 2.0e-5)

    def test_exact_fields_have_small_compatibility_curl(self) -> None:
        grad_a, grad_s = density_log_gradient_2d(
            self.a, self.s,
            self.theta_aa, self.theta_as, self.theta_ss,
            self.acc_a, self.acc_s,
            self.dtu_a, self.dtu_s,
        )
        curl = compatibility_curl_2d(self.a, self.s, grad_a, grad_s)
        self.assertLess(float(np.max(np.abs(curl))), 2.0e-9)

    def test_g2_mean_energy_is_preserved_by_reconstruction(self) -> None:
        grad_a, grad_s = density_log_gradient_2d(
            self.a, self.s,
            self.theta_aa, self.theta_as, self.theta_ss,
            self.acc_a, self.acc_s,
            self.dtu_a, self.dtu_s,
        )
        recovered, _ = reconstruct_density_path_2d(self.a, self.s, grad_a, grad_s)
        delta_u0 = (
            (self.A - 1.0) ** 2
            + 0.12 * (1.0 - np.cos(2.0 * np.pi * self.S))
            + 0.08 * (self.A - 1.0) * (1.0 - np.cos(2.0 * np.pi * self.S))
        )
        reference = mean_intrinsic_energy_2d(self.a, self.s, self.target, delta_u0)
        recovered_mean = mean_intrinsic_energy_2d(self.a, self.s, recovered, delta_u0)
        self.assertLess(abs(recovered_mean - reference) / abs(reference), 2.0e-5)

    def test_singular_covariance_is_rejected(self) -> None:
        zeros = np.zeros_like(self.theta_aa)
        with self.assertRaises(ValueError):
            density_log_gradient_2d(
                self.a, self.s,
                self.theta_aa, np.sqrt(self.theta_aa * self.theta_ss), self.theta_ss,
                zeros, zeros, zeros, zeros,
            )


if __name__ == "__main__":
    unittest.main()
