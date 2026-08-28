# === 한국어 파일 안내 시작 ===
# - 파일 역할: exact instantaneous P shape law와 neighbor-joint LJ acceleration 계산을 검증하는 회귀 테스트다.
# - 주요 클래스: TestNormalLJDistributionShape
# - 주요 함수/메서드: test_reconstructs_known_density, test_constant_theta_special_case, test_neighbor_acceleration_independent_case
# - 검증 내용: log-slope 적분으로 알려진 density를 복원하는지, constant-Theta 특수형이 exponential-LJ 구조를 재현하는지, P2 적분식이 조건부 LJ acceleration과 일치하는지 확인한다.
# === 한국어 파일 안내 끝 ===
import unittest

import numpy as np

from theory.normal_lj_chain import normalized_lj_energy, normalized_lj_force
from theory.normal_lj_distribution_shape import (
    density_log_slope,
    lj_conditional_acceleration,
    reconstruct_density_from_shape_fields,
)


class TestNormalLJDistributionShape(unittest.TestCase):
    def test_reconstructs_known_density(self) -> None:
        lam = np.linspace(0.82, 1.28, 2001)
        theta = 0.012 + 0.006 * (lam - 1.0) ** 2
        target = np.exp(-23.0 * (lam - 1.04) ** 2 - 5.0 * (lam - 1.04) ** 4)
        target /= np.trapezoid(target, lam)

        dlogp = np.gradient(np.log(target), lam, edge_order=2)
        dlogtheta = np.gradient(np.log(theta), lam, edge_order=2)

        # Rearranged exact identity:
        # a_bar - D_t u = Theta * (d log P + d log Theta).
        dtu = np.zeros_like(lam)
        abar = theta * (dlogp + dlogtheta)

        recovered = reconstruct_density_from_shape_fields(lam, theta, abar, dtu)
        self.assertLess(float(np.max(np.abs(recovered - target))), 3.0e-4)

    def test_constant_theta_special_case(self) -> None:
        lam = np.linspace(0.86, 1.24, 1801)
        theta_value = 0.02
        theta = np.full_like(lam, theta_value)
        neighbor_force_mean = 0.004
        dtu = np.zeros_like(lam)
        abar = 2.0 * neighbor_force_mean - 2.0 * normalized_lj_force(lam)

        recovered = reconstruct_density_from_shape_fields(lam, theta, abar, dtu)

        # Under the stated special assumptions the exact shape law reduces to
        # P proportional to exp[(2 f_n lambda - 2 phi(lambda))/Theta].
        exponent = (
            2.0 * neighbor_force_mean * lam
            - 2.0 * normalized_lj_energy(lam)
        ) / theta_value
        exponent -= float(np.max(exponent))
        reference = np.exp(exponent)
        reference /= np.trapezoid(reference, lam)

        self.assertLess(float(np.max(np.abs(recovered - reference))), 2.0e-4)

    def test_neighbor_acceleration_independent_case(self) -> None:
        lam = np.linspace(0.9, 1.15, 501)
        nbr = np.linspace(0.9, 1.15, 601)

        p1 = np.exp(-80.0 * (lam - 1.01) ** 2)
        p1 /= np.trapezoid(p1, lam)
        pn = np.exp(-70.0 * (nbr - 1.02) ** 2)
        pn /= np.trapezoid(pn, nbr)

        # Artificial independent-neighbor state used only to test the P2 integral.
        p2 = p1[:, None] * pn[None, :]
        abar = lj_conditional_acceleration(lam, nbr, p2, nbr, p2, p1)

        mean_neighbor_force = float(np.trapezoid(normalized_lj_force(nbr) * pn, nbr))
        expected = 2.0 * mean_neighbor_force - 2.0 * normalized_lj_force(lam)
        self.assertLess(float(np.max(np.abs(abar - expected))), 1.0e-11)

    def test_log_slope_requires_positive_theta(self) -> None:
        lam = np.array([0.9, 1.0, 1.1])
        theta = np.array([0.01, 0.0, 0.01])
        with self.assertRaises(ValueError):
            density_log_slope(lam, theta, np.zeros(3), np.zeros(3))


if __name__ == "__main__":
    unittest.main()
