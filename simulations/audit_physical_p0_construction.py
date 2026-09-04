"""Audit physical scales relevant to candidate P0 construction.

This script does not fit or adopt a physical aluminum P0. It checks the retained
normal calibration and demonstrates why an instantaneous classical harmonic
thermal spacing marginal cannot be silently used as the strict P0-only initial
state of the local-traction candidate.

The preferred strict P0-only interpretation is a slow structural/coarse-grained
spacing distribution obtained by spatial push-forward of measured/computed
residual spacing or microstrain.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


KB = 1.380649e-23
EV_J = 1.602176634e-19

A_REF = 2.8627442948e-10
E_PA = 69.0e9
A0 = 6.0338e-20
M_EXP = 12.19
N_EXP = 6.0

TEMPERATURES_K = [80.0, 293.0, 300.0]
REFERENCE_RESIDUAL_STRAIN_SCALE = 1.0e-4


def normal_curvature_threshold() -> float:
    return ((M_EXP + 1.0) / (N_EXP + 1.0)) ** (1.0 / (M_EXP - N_EXP))


def upper_standard_normal_tail(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def thermal_harmonic_row(temperature_k: float, stiffness_n_per_m: float, lambda_c: float) -> dict:
    sigma_a = math.sqrt(KB * temperature_k / stiffness_n_per_m)
    sigma_lambda = sigma_a / A_REF
    z_to_lambda_c = (lambda_c - 1.0) / sigma_lambda
    return {
        "temperature_k": temperature_k,
        "sigma_a_m": sigma_a,
        "sigma_a_pm": sigma_a * 1.0e12,
        "sigma_lambda": sigma_lambda,
        "z_score_to_lambda_c": z_to_lambda_c,
        "naive_harmonic_upper_tail_above_lambda_c": upper_standard_normal_tail(z_to_lambda_c),
    }


def main() -> None:
    u_ref_j = E_PA * A0 * A_REF
    effective_normal_well_depth_j = u_ref_j / (M_EXP * N_EXP)
    stiffness = E_PA * A0 / A_REF
    lambda_c = normal_curvature_threshold()
    delta_a_c = (lambda_c - 1.0) * A_REF

    thermal_rows = [
        thermal_harmonic_row(t, stiffness, lambda_c) for t in TEMPERATURES_K
    ]

    residual_spacing_scale = REFERENCE_RESIDUAL_STRAIN_SCALE * A_REF

    result = {
        "classification": "candidate physical-P0 calibration and consistency audit",
        "status": "diagnostic only; no physical P0 family adopted",
        "retained_normal_calibration": {
            "a_ref_m": A_REF,
            "E_pa": E_PA,
            "A0_m2": A0,
            "m": M_EXP,
            "n": N_EXP,
            "U_ref_j": u_ref_j,
            "U_ref_eV": u_ref_j / EV_J,
            "effective_1d_normal_well_depth_j": effective_normal_well_depth_j,
            "effective_1d_normal_well_depth_eV": effective_normal_well_depth_j / EV_J,
            "harmonic_spacing_stiffness_N_per_m": stiffness,
            "lambda_c": lambda_c,
            "delta_a_to_lambda_c_m": delta_a_c,
            "delta_a_to_lambda_c_pm": delta_a_c * 1.0e12,
        },
        "structural_P0": {
            "preferred_for_strict_P0_only_candidate": True,
            "definition": "spatial push-forward of a prepared slow structural spacing/microstrain field after fast phonon motion is averaged out",
            "initial_coarse_grained_rate_condition": "c0 = 0",
            "reference_residual_strain_scale_for_resolution_demo_only": REFERENCE_RESIDUAL_STRAIN_SCALE,
            "corresponding_spacing_scale_m": residual_spacing_scale,
            "corresponding_spacing_scale_pm": residual_spacing_scale * 1.0e12,
            "warning": "the 1e-4 strain scale is not adopted as a universal aluminum P0 width",
        },
        "instantaneous_thermal_harmonic_diagnostic": {
            "assumptions": [
                "classical canonical preparation",
                "single harmonic normal spacing coordinate around lambda=1",
                "phi''(1)=1 so K_a=E*A0/a_ref",
                "no quantum phonon correction",
                "no inter-spacing displacement correlation correction",
            ],
            "rows": thermal_rows,
            "warning": "these positional widths are not adopted as the structural P0; a genuine thermal preparation also has nonzero initial rate statistics",
        },
        "consistency_conclusion": [
            "A structural/coarse-grained P0 can be constructed from measured/computed residual spacing or microstrain by spatial push-forward without selecting a named PDF family.",
            "An instantaneous thermal spacing marginal is not sufficient for strict P0-only initialization because the thermal phase-space state also contains nonzero rate statistics.",
            "Naively combining the current classical harmonic thermal width with the current lambda_c gives a non-negligible instantaneous tail at room temperature, so it must not be interpreted as crack probability.",
            "The local-traction P0-only candidate still requires a separate physically justified laboratory-time-scale/cycle-evolution mechanism before promotion.",
        ],
    }

    output = Path("results/data/physical_p0_construction/summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
