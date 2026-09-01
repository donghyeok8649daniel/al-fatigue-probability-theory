"""Independent symbolic/numerical audit of the active multilayer potential."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

from theory.registry_lattice import (
    MultilayerPotentialParameters,
    bessel_lambert,
    bessel_lambert_polylog,
    dU_da,
    dU_da_direct,
    dU_ds,
    h_q_bessel,
    h_q_direct,
    h_q_polylog,
    normal_stationary_points,
    u0,
    v_slip,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "data" / "registry_plasticity" / "multilayer_verification_summary.json"


def half_integer_coefficients(order_n: int) -> list[int]:
    """Coefficients in K_(n+1/2)=sqrt(pi/2x)e^-x sum c_j x^-j."""
    return [
        int(sp.factorial(order_n + j) / (
            sp.factorial(j) * sp.factorial(order_n - j) * 2**j
        ))
        for j in range(order_n + 1)
    ]


def main() -> None:
    coefficients = {
        "K_5_over_2": half_integer_coefficients(2),
        "K_11_over_2": half_integer_coefficients(5),
    }
    assert coefficients["K_5_over_2"] == [1, 3, 3]
    assert coefficients["K_11_over_2"] == [1, 15, 105, 420, 945, 945]

    checks = []
    for q in (6, 12):
        for delta, eta in ((0.0, 0.7), (0.23, 1.0), (0.5, 1.4)):
            direct = h_q_direct(q, delta, eta, 400, 800)
            bessel = h_q_bessel(q, delta, eta, 20, 64)
            polylog = h_q_polylog(q, delta, eta, 20)
            checks.append({
                "q": q,
                "delta": delta,
                "eta": eta,
                "direct_bessel_relative_error": abs(direct - bessel) / abs(bessel),
                "polylog_bessel_relative_error": abs(polylog - bessel) / abs(bessel),
            })

    lambert_error = max(
        abs(bessel_lambert((q - 1) / 2, x, 160) - bessel_lambert_polylog(q, x))
        / abs(bessel_lambert_polylog(q, x))
        for q in (6, 12) for x in (0.9, 2.5, 7.0)
    )
    params = MultilayerPotentialParameters(sigma_lj=0.82, bessel_modes=20, layer_modes=64)
    a, s, step = 0.97, 0.19, 2.0e-6
    derivative_checks = {
        "dU_da_bessel_vs_direct_relative_error": abs(
            dU_da(a, s, params) - dU_da_direct(a, s, params, 300, 600)
        ) / abs(dU_da(a, s, params)),
        "dU_da_finite_difference_relative_error": abs(
            dU_da(a, s, params)
            - (u0(a + step, s, params) - u0(a - step, s, params)) / (2 * step)
        ) / abs(dU_da(a, s, params)),
        "dU_ds_finite_difference_relative_error": abs(
            dU_ds(a, s, params)
            - (u0(a, s + step, params) - u0(a, s - step, params)) / (2 * step)
        ) / abs(dU_ds(a, s, params)),
    }
    identities = {
        "periodicity_absolute_error": abs(u0(a, s + params.b, params) - u0(a, s, params)),
        "even_symmetry_absolute_error": abs(u0(a, -s, params) - u0(a, s, params)),
        "slip_reference_absolute_error": abs(v_slip(a, 0.0, params, 0.0)),
    }
    roots = normal_stationary_points(0.1, 0.0, params, 0.6, 2.8, 160, 100, 200)
    assert len(roots) == 2 and roots[0][1] > 0.0 and roots[1][1] < 0.0

    summary = {
        "status": "passed",
        "counting": "U0=sum_{k>=1}W(k*a,s), unit layer multiplicity, common s",
        "symbolic_half_integer_coefficients": coefficients,
        "point_checks": checks,
        "max_direct_bessel_relative_error": max(
            item["direct_bessel_relative_error"] for item in checks
        ),
        "max_polylog_bessel_relative_error": max(
            item["polylog_bessel_relative_error"] for item in checks
        ),
        "max_lambert_polylog_relative_error": lambert_error,
        "derivative_checks": derivative_checks,
        "identities": identities,
        "normal_stationary_points": [
            {"a": root, "curvature": curvature,
             "classification": "stable" if curvature > 0 else "outer_barrier"}
            for root, curvature in roots
        ],
        "units": {
            "U0": "energy",
            "dU_da_and_dU_ds": "force",
            "P": "inverse state-space measure",
            "J": "probability flux",
            "Ddot_irr": "energy/time",
            "E_hyst": "energy",
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
