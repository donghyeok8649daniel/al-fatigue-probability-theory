# === 한국어 파일 안내 시작 ===
# - 파일 역할: 최종 P0->survival 축약식의 장벽/특성면적 민감도를 재현한다.
# - 핵심: generalized-LJ stable branch, lambda_c, TST rate를 이용해 A_c/A_0 sweep을 계산한다.
# - 주의: A_c/A_0 값은 calibration이 아니라 sensitivity diagnostic이다.
# === 한국어 파일 안내 끝 ===
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "data" / "final_reduced_closure_sensitivity"

M_EXP = 12.19
N_EXP = 6.0
E = 69e9
A0 = 6.0338e-20
A_REF = 2.8627442948e-10
T0 = 5.55046e-14
KB = 1.380649e-23
TEMP = 300.0
EV = 1.602176634e-19


def phi(lam: float) -> float:
    return lam ** (-M_EXP) / (M_EXP * (M_EXP - N_EXP)) - lam ** (-N_EXP) / (N_EXP * (M_EXP - N_EXP))


def dphi(lam: float) -> float:
    return (lam ** (-N_EXP - 1.0) - lam ** (-M_EXP - 1.0)) / (M_EXP - N_EXP)


def ddphi(lam: float) -> float:
    return ((M_EXP + 1.0) * lam ** (-M_EXP - 2.0) - (N_EXP + 1.0) * lam ** (-N_EXP - 2.0)) / (M_EXP - N_EXP)


LAMBDA_C = ((M_EXP + 1.0) / (N_EXP + 1.0)) ** (1.0 / (M_EXP - N_EXP))
Q_C = dphi(LAMBDA_C)
E0 = E * A0 * A_REF


def stable_lambda(q: float) -> float:
    if q >= Q_C:
        return LAMBDA_C
    return brentq(lambda lam: dphi(lam) - q, 0.8, LAMBDA_C * (1.0 - 1e-12))


def delta_psi(lam: float) -> float:
    q = dphi(lam)
    return (phi(LAMBDA_C) - q * LAMBDA_C) - (phi(lam) - q * lam)


def one_cycle_hazard(area_ratio: float, mean_mpa: float = 100.0, amp_mpa: float = 100.0, freq_hz: float = 20.0) -> float:
    period = 1.0 / freq_hz
    t = np.linspace(0.0, period, 20001)
    sigma = (mean_mpa + amp_mpa * np.sin(2.0 * math.pi * freq_hz * t)) * 1e6
    lam = np.asarray([stable_lambda(s / E) for s in sigma])
    curv = np.asarray([ddphi(x) for x in lam])
    barrier = np.asarray([delta_psi(x) for x in lam])
    prefactor = np.sqrt(np.maximum(curv, 0.0)) / (2.0 * math.pi * T0)
    rate = prefactor * np.exp(-area_ratio * E0 * barrier / (KB * TEMP))
    return float(np.trapezoid(rate, t))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    barrier_rows = []
    for stress_mpa in (0, 50, 100, 150, 200):
        lam = stable_lambda(stress_mpa * 1e6 / E)
        dpsi = delta_psi(lam)
        barrier_rows.append({
            "stress_MPa": stress_mpa,
            "stable_lambda": lam,
            "curvature": ddphi(lam),
            "delta_psi": dpsi,
            "barrier_eV_for_A0": E0 * dpsi / EV,
        })

    sweep = []
    for ratio in (30, 40, 50, 60):
        H = one_cycle_hazard(float(ratio))
        p = -math.expm1(-H)
        n50 = math.log(2.0) / H
        sweep.append({
            "A_c_over_A0": ratio,
            "integrated_hazard_per_cycle": H,
            "escape_probability_per_cycle": p,
            "survival_probability_per_cycle": math.exp(-H),
            "local_median_cycles": n50,
        })

    payload = {
        "classification": "final reduced closure sensitivity audit",
        "status": "sensitivity only; no characteristic-area calibration",
        "m": M_EXP,
        "n": N_EXP,
        "lambda_c": LAMBDA_C,
        "q_c": Q_C,
        "sigma_c_GPa_if_q=sigma/E": Q_C * E / 1e9,
        "E_GPa": E / 1e9,
        "a0_m": A_REF,
        "A0_m2": A0,
        "E0_eV": E0 / EV,
        "temperature_K": TEMP,
        "protocol": {"mean_stress_MPa": 100.0, "amplitude_MPa": 100.0, "frequency_Hz": 20.0},
        "atomic_reference_barrier": barrier_rows,
        "area_ratio_sweep": sweep,
        "verdict": [
            "One atomic reference area gives a barrier of only about 0.02 eV over the tested stress range and is not a rare high-cycle event at 300 K.",
            "The barrier scales linearly with A_c while the TST prefactor is area-independent under coherent mass scaling, creating exponential area sensitivity.",
            "A_c/A0 is deliberately not selected here; the sweep demonstrates why it must be independently calibrated or computed.",
            "For a symbolic A_c, the reduced survivor law is mathematically closed from P0 and sigma(t).",
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
