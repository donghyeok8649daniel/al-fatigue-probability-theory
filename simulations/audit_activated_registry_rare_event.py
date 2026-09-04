"""Diagnostic audit of a finite-temperature activated registry rare-event route.

This is NOT an active governing law.  It asks whether the already-derived
registry barrier can separate the THz conservative attempt scale from a
many-cycle rare-event time scale without tuning a slow inertia.

The rate law itself is conditional: it assumes a finite-temperature bath,
fast intrawell relaxation, rare activated escape, and a coherent patch of N
reference repeats.  N is scanned only as a sensitivity variable and is not
calibrated here.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

M_EXP = 12.19
N_EXP = 6.0
E_PA = 69.0e9
A0_M2 = 6.0338e-20
A0_PHYS_M = 2.8627442948e-10
T0_S = 5.55046e-14
TEMPERATURE_K = 300.0
LOAD_FREQUENCY_HZ = 20.0
MEAN_STRESS_PA = 100.0e6
AMPLITUDE_STRESS_PA = 100.0e6
KB_EV_PER_K = 8.617333262e-5
EV_J = 1.602176634e-19
K_MAX = 120
P_MAX = 300
THETA_POINTS = 4001

C_MN = M_EXP / (M_EXP - N_EXP) * (M_EXP / N_EXP) ** (
    N_EXP / (M_EXP - N_EXP)
)


def u0(a: float, s: float) -> float:
    """Dimensionless multilayer registry energy, b=sigma=epsilon=1."""
    k = np.arange(1, K_MAX + 1, dtype=float)[:, None]
    p = np.arange(-P_MAX, P_MAX + 1, dtype=float)[None, :]
    r = np.sqrt((k * a) ** 2 + (p + s) ** 2)
    return float(C_MN * np.sum(r ** (-M_EXP) - r ** (-N_EXP)))


def second_derivatives(a: float, s: float, h: float = 1.0e-5) -> tuple[float, float]:
    uc = u0(a, s)
    uaa = (u0(a + h, s) - 2.0 * uc + u0(a - h, s)) / h**2
    uss = (u0(a, s + h) - 2.0 * uc + u0(a, s - h)) / h**2
    return uaa, uss


def phi_prime(stretch: float) -> float:
    return (
        stretch ** (-N_EXP - 1.0) - stretch ** (-M_EXP - 1.0)
    ) / (M_EXP - N_EXP)


def stable_normal_stretch(stress_pa: float) -> float:
    q = stress_pa / E_PA
    if abs(q) < 1.0e-18:
        return 1.0
    if q > 0.0:
        return brentq(lambda x: phi_prime(x) - q, 1.0, 1.1077)
    return brentq(lambda x: phi_prime(x) - q, 0.7, 1.0)


def main() -> None:
    a0_reg = float(
        minimize_scalar(
            lambda a: u0(a, 0.5),
            bounds=(0.8, 1.2),
            method="bounded",
            options={"xatol": 1.0e-12},
        ).x
    )
    uaa, uss = second_derivatives(a0_reg, 0.5)
    curvature_ratio = uss / uaa

    physical_normal_stiffness = E_PA * A0_M2 / A0_PHYS_M
    b_phys = A0_PHYS_M / a0_reg
    energy_scale_j = physical_normal_stiffness * b_phys**2 / uaa
    energy_scale_ev = energy_scale_j / EV_J

    barrier0_dim = u0(a0_reg, 0.0) - u0(a0_reg, 0.5)
    barrier0_ev = barrier0_dim * energy_scale_ev
    attempt_frequency_hz = (
        1.0 / (2.0 * math.pi * T0_S) * math.sqrt(curvature_ratio)
    )

    theta = np.linspace(0.0, 2.0 * math.pi, THETA_POINTS)
    stress = MEAN_STRESS_PA + AMPLITUDE_STRESS_PA * np.sin(theta)
    stretch = np.asarray([stable_normal_stretch(float(s)) for s in stress])
    a_registry = a0_reg * stretch

    a_grid = np.linspace(float(np.min(a_registry)), float(np.max(a_registry)), 501)
    barrier_grid_ev = np.asarray(
        [(u0(a, 0.0) - u0(a, 0.5)) * energy_scale_ev for a in a_grid]
    )
    barrier_ev = np.interp(a_registry, a_grid, barrier_grid_ev)

    kbt_ev = KB_EV_PER_K * TEMPERATURE_K
    rows = []
    for coherent_repeats in range(6, 13):
        rate_hz = attempt_frequency_hz * np.exp(
            -coherent_repeats * barrier_ev / kbt_ev
        )
        hazard_per_cycle = float(
            np.trapezoid(rate_hz, theta)
            / (2.0 * math.pi * LOAD_FREQUENCY_HZ)
        )
        transition_probability_per_cycle = -math.expm1(-hazard_per_cycle)
        local_survival_per_cycle = math.exp(-hazard_per_cycle)
        median_cycles = math.log(2.0) / hazard_per_cycle

        zero_load_rate_hz = attempt_frequency_hz * math.exp(
            -coherent_repeats * barrier0_ev / kbt_ev
        )
        zero_load_hazard_per_20hz_period = zero_load_rate_hz / LOAD_FREQUENCY_HZ

        rows.append(
            {
                "coherent_repeats_N": coherent_repeats,
                "hazard_per_cycle": hazard_per_cycle,
                "transition_probability_per_cycle": transition_probability_per_cycle,
                "local_survival_per_cycle": local_survival_per_cycle,
                "median_cycles_if_rate_model_applied": median_cycles,
                "zero_load_rate_hz": zero_load_rate_hz,
                "zero_load_hazard_per_20hz_period": zero_load_hazard_per_20hz_period,
            }
        )

    barrier_samples = []
    for stress_mpa in (0.0, 50.0, 100.0, 150.0, 200.0):
        lam = stable_normal_stretch(stress_mpa * 1.0e6)
        a_reg = a0_reg * lam
        barrier = (u0(a_reg, 0.0) - u0(a_reg, 0.5)) * energy_scale_ev
        barrier_samples.append(
            {
                "stress_mpa": stress_mpa,
                "stable_normal_stretch": lam,
                "registry_a_over_b": a_reg,
                "barrier_eV_per_repeat": barrier,
            }
        )

    payload = {
        "classification": "conditional activated registry rare-event diagnostic",
        "status": "candidate only; not active law and not characteristic-area calibration",
        "assumptions_required_for_rate": [
            "finite-temperature bath",
            "fast intrawell equilibration compared with interwell escape",
            "rare activated crossing",
            "one coherent registry saddle for a patch of N repeats",
            "harmonic registry frequency used as diagnostic attempt frequency",
            "quasistatic normal response at laboratory loading frequency",
        ],
        "registry_surface": {
            "m": M_EXP,
            "n": N_EXP,
            "a0_registry_over_b": a0_reg,
            "Uaa_dimensionless": uaa,
            "Uss_dimensionless": uss,
            "curvature_ratio": curvature_ratio,
            "energy_scale_eV_per_dimensionless_energy": energy_scale_ev,
            "barrier_at_reference_eV_per_repeat": barrier0_ev,
            "attempt_frequency_hz": attempt_frequency_hz,
            "periodicity_statement": "U0(a,s+b)=U0(a,s); adjacent ideal registry wells are energetically equivalent",
        },
        "loading": {
            "temperature_K": TEMPERATURE_K,
            "frequency_Hz": LOAD_FREQUENCY_HZ,
            "mean_stress_MPa": MEAN_STRESS_PA / 1.0e6,
            "amplitude_stress_MPa": AMPLITUDE_STRESS_PA / 1.0e6,
            "minimum_barrier_over_cycle_eV_per_repeat": float(np.min(barrier_ev)),
            "maximum_barrier_over_cycle_eV_per_repeat": float(np.max(barrier_ev)),
        },
        "barrier_samples": barrier_samples,
        "coherent_patch_sensitivity": rows,
        "verdict": [
            "The barrier factor can separate a THz attempt frequency from many-cycle rare-event probabilities without assigning an unphysical slow inertia.",
            "N is not calibrated or adopted; it is characteristic/coherent patch information deferred to later spatial calibration.",
            "The ideal periodic registry transition changes the unwrapped well index z but does not automatically change the normal interaction branch, so it does not by itself close progressive P_a evolution.",
            "The same activated model predicts zero-load registry hopping; equivalent-well hopping must not be called damage without a physical irreversible post-transition state.",
        ],
    }

    output = Path("results/data/activated_registry_rare_event/summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
