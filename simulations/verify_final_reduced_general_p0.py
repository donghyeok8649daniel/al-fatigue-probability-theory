# === 한국어 파일 안내 시작 ===
# - 파일 역할: 최종 P0->survival 축약법칙을 일반 이산 structural P0에 대해 검증한다.
# - 검증 항목: 폐주기 기계적 복귀, cycle hazard, survivor selection, 주파수/온도 예측.
# - 주의: 아래 P0 및 A_c/A0=50은 검증용 synthetic input이며 Al 보정값이 아니다.
# === 한국어 파일 안내 끝 ===
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "data" / "final_reduced_general_p0"

M_EXP = 12.19
N_EXP = 6.0
E = 69e9
A0 = 6.0338e-20
A_REF = 2.8627442948e-10
T0 = 5.55046e-14
KB = 1.380649e-23
TEMP = 300.0
AREA_RATIO = 50.0
MEAN_MPA = 100.0
AMP_MPA = 100.0
FREQ_HZ = 20.0


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


def delta_psi_c(lam: float) -> float:
    q = dphi(lam)
    return (phi(LAMBDA_C) - q * LAMBDA_C) - (phi(lam) - q * lam)


def cycle_hazard(lambda0: float, q_ref: float, area_ratio: float, temp: float, freq_hz: float) -> tuple[float, float]:
    period = 1.0 / freq_hz
    t = np.linspace(0.0, period, 20001)
    q = (MEAN_MPA + AMP_MPA * np.sin(2.0 * math.pi * freq_hz * t)) * 1e6 / E
    q_r = dphi(lambda0) - q_ref
    lam = np.asarray([stable_lambda(q_r + qi) for qi in q])
    curv = np.asarray([ddphi(x) for x in lam])
    barrier = np.asarray([delta_psi_c(x) for x in lam])
    prefactor = np.sqrt(np.maximum(curv, 0.0)) / (2.0 * math.pi * T0)
    rate = prefactor * np.exp(-area_ratio * E0 * barrier / (KB * temp))
    H = float(np.trapezoid(rate, t))
    return H, float(abs(lam[-1] - lambda0))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    q_ref = MEAN_MPA * 1e6 / E
    lambda_ref = stable_lambda(q_ref)

    offsets = np.asarray([-0.0012, -0.0008, -0.0004, 0.0, 0.0004, 0.0008, 0.0012])
    lambda0 = lambda_ref + offsets
    weights = np.asarray([0.03, 0.08, 0.16, 0.26, 0.22, 0.16, 0.09], dtype=float)
    weights /= weights.sum()

    rows = []
    hazards = []
    closed_errors = []
    for lam0, w in zip(lambda0, weights):
        H, err = cycle_hazard(float(lam0), q_ref, AREA_RATIO, TEMP, FREQ_HZ)
        hazards.append(H)
        closed_errors.append(err)
        rows.append({
            "lambda0": float(lam0),
            "weight": float(w),
            "integrated_hazard_per_cycle": H,
            "survival_probability_per_cycle": math.exp(-H),
        })

    hazards = np.asarray(hazards)
    S1 = float(np.sum(weights * np.exp(-hazards)))

    survival_rows = []
    for cycles in (1, 10_000, 100_000, 300_000, 1_000_000):
        survivor_weights = weights * np.exp(-cycles * hazards)
        S = float(survivor_weights.sum())
        conditioned = survivor_weights / S
        survival_rows.append({
            "cycles": cycles,
            "survival": S,
            "first_passage_fraction": 1.0 - S,
            "conditional_mean_lambda0": float(np.sum(conditioned * lambda0)),
        })

    frequency_rows = []
    for freq in (1.0, 5.0, 10.0, 20.0, 50.0, 100.0):
        H, _ = cycle_hazard(float(lambda_ref), q_ref, AREA_RATIO, TEMP, freq)
        frequency_rows.append({
            "frequency_Hz": freq,
            "integrated_hazard_per_cycle": H,
            "frequency_times_hazard": freq * H,
            "median_cycles": math.log(2.0) / H,
            "median_time_s": math.log(2.0) / (freq * H),
        })

    temperature_rows = []
    for temp in (250.0, 275.0, 300.0, 325.0, 350.0, 400.0):
        H, _ = cycle_hazard(float(lambda_ref), q_ref, AREA_RATIO, temp, FREQ_HZ)
        temperature_rows.append({
            "temperature_K": temp,
            "integrated_hazard_per_cycle": H,
            "median_cycles": math.log(2.0) / H,
        })

    payload = {
        "classification": "final reduced closure general-P0 verification",
        "status": "synthetic structural P0 and area-ratio diagnostic; not Al life calibration",
        "model": {
            "m": M_EXP,
            "n": N_EXP,
            "lambda_c": LAMBDA_C,
            "q_c": Q_C,
            "E_GPa": E / 1e9,
            "a0_m": A_REF,
            "A0_m2": A0,
            "temperature_K": TEMP,
            "A_c_over_A0": AREA_RATIO,
            "mean_stress_MPa": MEAN_MPA,
            "amplitude_MPa": AMP_MPA,
            "frequency_Hz": FREQ_HZ,
        },
        "reference_lambda": lambda_ref,
        "p0_points": rows,
        "max_closed_cycle_lambda_error": float(max(closed_errors)),
        "mixture_one_cycle_survival": S1,
        "mixture_one_cycle_first_passage": 1.0 - S1,
        "survival_by_cycles": survival_rows,
        "frequency_sweep": frequency_rows,
        "temperature_sweep": temperature_rows,
        "verdict": [
            "The quasistatic characteristic map returns every structural label to its reference spacing after a closed stress cycle to numerical precision.",
            "Despite that reversible mechanical return, the survivor mass decreases because each structural label has a nonzero thermal first-passage hazard.",
            "The normalized survivor population shifts toward lower-hazard initial spacings by selection; no permanent drift of every surviving spacing is imposed.",
            "Within the strict quasistatic/local-equilibrium limit, hazard per cycle scales as 1/f and median life in physical time is frequency-independent; this is a strong falsifiable prediction of the reduced hypothesis.",
            "Temperature sensitivity is exponentially strong and must be checked experimentally before the closure is treated as a material law for Al.",
        ],
    }

    (OUT / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
