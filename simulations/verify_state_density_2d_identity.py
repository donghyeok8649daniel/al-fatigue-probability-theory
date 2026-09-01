# === 한국어 파일 안내 시작 ===
# - 파일 역할: theory-core의 (a,s) Theta-tensor density identity를 비대각 covariance를 포함한 analytic field로 검산한다.
# - 주의: 이 스크립트의 target PDF와 energy surface는 수치 검산용 synthetic field이며 물리적 Al 예측이 아니다.
#   실제 (a,s) mechanics 데이터 검증 전 단계의 algebra/numerics verification이다.
# === 한국어 파일 안내 끝 ===
"""Numerically verify the exact two-coordinate density shape identity.

This is deliberately a synthetic self-consistency test.  It verifies the
algebra and grid reconstruction for a non-diagonal conditional velocity
covariance tensor, including non-zero Theta_as and non-zero D_t u.  It does not
claim that the synthetic target density is a physical fatigue distribution.

Run in the integrated tree (or with theory-core available on PYTHONPATH):

    python simulations/verify_state_density_2d_identity.py
"""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np

from theory.state_density_2d import (
    compatibility_curl_2d,
    density_log_gradient_2d,
    mean_intrinsic_energy_2d,
    probability_mass_2d,
    reconstruct_density_path_2d,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "data" / "two_coordinate_density_identity"


def build_case():
    a = np.linspace(0.86, 1.22, 181)
    s = np.linspace(-0.55, 0.55, 221)
    A, S = np.meshgrid(a, s, indexing="ij")

    # Smooth non-Gaussian target used only for an exact identity check.
    logp = (
        -55.0 * (A - 1.025) ** 2
        -10.0 * (S - 0.28 * (A - 1.0)) ** 2
        -1.5 * S**4
        +1.2 * (A - 1.0) * S
    )
    target = np.exp(logp - float(np.max(logp)))
    target /= probability_mass_2d(a, s, target)

    grad_a = np.gradient(np.log(target), a, axis=0, edge_order=2)
    grad_s = np.gradient(np.log(target), s, axis=1, edge_order=2)

    theta_aa = 0.020 + 0.004 * (A - 1.0) ** 2 + 0.0015 * S**2
    theta_ss = 0.035 + 0.003 * S**2 + 0.0010 * (A - 1.0) ** 2
    theta_as = 0.0045 + 0.0008 * (A - 1.0) * S

    # Non-zero material acceleration is included so the test does not collapse
    # to the static special case.
    dtu_a = (
        0.003
        * np.sin(2.0 * np.pi * (A - a[0]) / (a[-1] - a[0]))
        * np.cos(np.pi * S / 0.55)
    )
    dtu_s = (
        0.002
        * np.cos(2.0 * np.pi * (A - a[0]) / (a[-1] - a[0]))
        * np.sin(np.pi * S / 0.55)
    )

    div_theta_a = (
        np.gradient(theta_aa, a, axis=0, edge_order=2)
        + np.gradient(theta_as, s, axis=1, edge_order=2)
    )
    div_theta_s = (
        np.gradient(theta_as, a, axis=0, edge_order=2)
        + np.gradient(theta_ss, s, axis=1, edge_order=2)
    )

    # Rearranged exact identity:
    # A - D_t u = div(Theta) + Theta grad(log P).
    acc_a = dtu_a + div_theta_a + theta_aa * grad_a + theta_as * grad_s
    acc_s = dtu_s + div_theta_s + theta_as * grad_a + theta_ss * grad_s

    # A smooth energy-like surface used only to check G2 quadrature.
    delta_u0 = (
        (A - 1.0) ** 2
        + 0.12 * (1.0 - np.cos(2.0 * np.pi * S))
        + 0.08 * (A - 1.0) * (1.0 - np.cos(2.0 * np.pi * S))
    )
    return (
        a, s, target,
        theta_aa, theta_as, theta_ss,
        acc_a, acc_s, dtu_a, dtu_s,
        delta_u0,
    )


def main() -> None:
    (
        a, s, target,
        theta_aa, theta_as, theta_ss,
        acc_a, acc_s, dtu_a, dtu_s,
        delta_u0,
    ) = build_case()

    grad_a, grad_s = density_log_gradient_2d(
        a, s,
        theta_aa, theta_as, theta_ss,
        acc_a, acc_s,
        dtu_a, dtu_s,
    )
    curl = compatibility_curl_2d(a, s, grad_a, grad_s)
    recovered, path_mismatch = reconstruct_density_path_2d(a, s, grad_a, grad_s)

    l1 = probability_mass_2d(a, s, np.abs(recovered - target))
    mean_u_target = mean_intrinsic_energy_2d(a, s, target, delta_u0)
    mean_u_recovered = mean_intrinsic_energy_2d(a, s, recovered, delta_u0)
    relative_energy_error = abs(mean_u_recovered - mean_u_target) / abs(mean_u_target)

    correlation = theta_as / np.sqrt(theta_aa * theta_ss)
    summary = {
        "classification": "synthetic algebra/numerics verification of exact smooth-moment identity",
        "grid": {"a_points": int(a.size), "s_points": int(s.size)},
        "nonzero_theta_as": True,
        "conditional_velocity_correlation_min": float(np.min(correlation)),
        "conditional_velocity_correlation_max": float(np.max(correlation)),
        "max_abs_compatibility_curl": float(np.max(np.abs(curl))),
        "rms_compatibility_curl": float(np.sqrt(np.mean(curl**2))),
        "path_log_density_mismatch": float(path_mismatch),
        "density_L1_error": float(l1),
        "mean_intrinsic_energy_target": float(mean_u_target),
        "mean_intrinsic_energy_recovered": float(mean_u_recovered),
        "relative_mean_energy_error": float(relative_energy_error),
        "not_a_physical_claim": (
            "Target P and Delta U0 are synthetic verification fields. Physical (a,s) validation "
            "requires a mechanics-generated ensemble supplying A, u and Theta."
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
