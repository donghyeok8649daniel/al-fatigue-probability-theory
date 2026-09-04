"""Analytic/numeric scale audit for the candidate local-traction P0 propagator.

The audit checks whether the dimensionless frequency used in the internal
candidate test corresponds to laboratory fatigue loading under the retained
atomic time calibration. It also quantifies the startup acceleration caused by
using an unrelaxed initial P0 under nonzero mean traction.

No dissipative or empirical fatigue term is introduced.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


M_EXP = 12.19
N_EXP = 6.0
E_PA = 69.0e9
T0_S = 5.55046e-14
STRESS_MEAN_PA = 100.0e6
Q0 = STRESS_MEAN_PA / E_PA
OLD_AUDIT_OMEGA_STAR = 0.02
OLD_AUDIT_SUPPORT = [0.997, 1.0, 1.003]
LAB_FREQUENCIES_HZ = [1.0, 20.0, 100.0, 1000.0]


def dphi(stretch: float) -> float:
    return (
        stretch ** (-N_EXP - 1.0) - stretch ** (-M_EXP - 1.0)
    ) / (M_EXP - N_EXP)


def omega_star_from_hz(f_hz: float) -> float:
    return 2.0 * math.pi * f_hz * T0_S


def hz_from_omega_star(omega_star: float) -> float:
    return omega_star / (2.0 * math.pi * T0_S)


def linear_inertial_correction(omega_star: float, kappa: float = 1.0) -> float:
    return omega_star**2 / (kappa - omega_star**2)


def main() -> None:
    startup = []
    for lam0 in OLD_AUDIT_SUPPORT:
        startup.append(
            {
                "lambda0": lam0,
                "phi_prime_lambda0": dphi(lam0),
                "initial_acceleration_q0_minus_phi_prime": Q0 - dphi(lam0),
            }
        )

    lab_rows = []
    for f_hz in LAB_FREQUENCIES_HZ:
        omega_star = omega_star_from_hz(f_hz)
        lab_rows.append(
            {
                "frequency_hz": f_hz,
                "omega_star": omega_star,
                "linear_relative_inertial_correction_kappa_1": linear_inertial_correction(
                    omega_star
                ),
            }
        )

    result = {
        "classification": "candidate local-traction laboratory-timescale no-go audit",
        "status": "candidate rejected as standalone progressive lab-fatigue mechanism in present conservative atomic-inertia form",
        "retained_time_scale_s": T0_S,
        "old_candidate_audit": {
            "mean_stress_pa": STRESS_MEAN_PA,
            "q_at_initial_time": Q0,
            "omega_star": OLD_AUDIT_OMEGA_STAR,
            "equivalent_frequency_hz_under_retained_t0": hz_from_omega_star(
                OLD_AUDIT_OMEGA_STAR
            ),
            "startup_acceleration_samples": startup,
            "interpretation": "initial support was not mechanically equilibrated under q(0), so the prior strong history signal contains atomic startup transients",
        },
        "laboratory_frequency_mapping": lab_rows,
        "analytic_linearized_result": {
            "equation": "y_ddot + kappa*y = q_a*sin(omega_star*tau)",
            "periodic_particular_amplitude": "q_a/(kappa-omega_star^2)",
            "cycle_work_periodic_particular": 0.0,
            "relative_inertial_correction": "omega_star^2/(kappa-omega_star^2)",
            "interpretation": "away from resonance, the conservative periodic response has no dissipative phase lag and becomes quasistatic as omega_star -> 0",
        },
        "verdict": [
            "The mathematical P0 + stress-history -> P map remains well-defined as a candidate reduced constitutive propagator.",
            "The prior omega*=0.02 audit corresponds to roughly 57.35 GHz under the retained atomic time calibration, not laboratory fatigue frequency.",
            "At 20 Hz the kappa=1 linear inertial correction is approximately 4.865e-23, so the local conservative atomic coordinate is effectively quasistatic.",
            "Without a separately derived slow or irreversible mechanism, the candidate does not explain progressive subcritical laboratory-frequency fatigue accumulation.",
            "Return the main theory path to the exact finite-chain/correlation-hierarchy checkpoint rather than adding arbitrary damping/diffusion/damage."
        ],
    }

    output = Path("results/data/local_traction_lab_timescale/summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
