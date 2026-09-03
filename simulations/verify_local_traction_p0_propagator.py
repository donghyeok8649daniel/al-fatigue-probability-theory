"""Internal-consistency audit for the candidate local-traction P0 -> P propagator.

This is NOT the active exact finite-chain model. It tests the reduced local
constitutive hypothesis documented in
`docs/CANDIDATE_LOCAL_TRACTION_P0_PROPAGATOR.md` on theory-core.

The diagnostic initial spacing set is deliberately explicit and bounded. It is
not proposed as the physical aluminum P0.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


M_EXP = 12.19
N_EXP = 6.0
E_PA = 69.0e9
STRESS_MEAN_PA = 100.0e6
STRESS_AMPLITUDE_PA = 100.0e6
Q_MEAN = STRESS_MEAN_PA / E_PA
Q_AMPLITUDE = STRESS_AMPLITUDE_PA / E_PA
OMEGA_STAR = 0.02
DT = 0.02
CYCLES = 3
N_CHARACTERISTICS = 5001
LAMBDA_MIN = 0.997
LAMBDA_MAX = 1.003
TAIL_THRESHOLD = 1.004


def phi(stretch: np.ndarray) -> np.ndarray:
    return (
        stretch ** (-M_EXP) / (M_EXP * (M_EXP - N_EXP))
        - stretch ** (-N_EXP) / (N_EXP * (M_EXP - N_EXP))
    )


def dphi(stretch: np.ndarray) -> np.ndarray:
    return (
        stretch ** (-N_EXP - 1.0) - stretch ** (-M_EXP - 1.0)
    ) / (M_EXP - N_EXP)


def q_of_tau(tau: float) -> float:
    return Q_MEAN + Q_AMPLITUDE * np.sin(OMEGA_STAR * tau)


def empirical_ks(sample_a: np.ndarray, sample_b: np.ndarray) -> float:
    grid = np.unique(np.concatenate((sample_a, sample_b)))
    cdf_a = np.searchsorted(np.sort(sample_a), grid, side="right") / sample_a.size
    cdf_b = np.searchsorted(np.sort(sample_b), grid, side="right") / sample_b.size
    return float(np.max(np.abs(cdf_a - cdf_b)))


def main() -> None:
    lambda0 = np.linspace(LAMBDA_MIN, LAMBDA_MAX, N_CHARACTERISTICS)
    lam = lambda0.copy()
    rate = np.zeros_like(lam)
    initial_energy = 0.5 * rate**2 + phi(lam)
    work = np.zeros_like(lam)
    crossed = lam >= TAIL_THRESHOLD

    period = 2.0 * np.pi / OMEGA_STAR
    total_time = CYCLES * period
    n_steps = int(round(total_time / DT))

    target_times = [
        0.0,
        0.5 * period,
        period,
        1.5 * period,
        2.0 * period,
        3.0 * period,
    ]
    target_steps = {int(round(t / DT)): t for t in target_times}
    snapshots: dict[float, dict] = {}

    tau = 0.0
    acceleration = q_of_tau(tau) - dphi(lam)

    for step in range(n_steps + 1):
        crossed |= lam >= TAIL_THRESHOLD

        if step in target_steps:
            label = target_steps[step]
            current_energy = 0.5 * rate**2 + phi(lam)
            snapshots[label] = {
                "actual_tau": tau,
                "q": q_of_tau(tau),
                "mean_lambda": float(np.mean(lam)),
                "variance_lambda": float(np.var(lam)),
                "q99_lambda": float(np.quantile(lam, 0.99)),
                "snapshot_tail_fraction": float(np.mean(lam >= TAIL_THRESHOLD)),
                "first_passage_fraction": float(np.mean(crossed)),
                "mean_energy_change": float(np.mean(current_energy - initial_energy)),
                "mean_external_work": float(np.mean(work)),
                "max_abs_characteristic_work_minus_energy": float(
                    np.max(np.abs((current_energy - initial_energy) - work))
                ),
                "lambda_sample": lam.copy(),
            }

        if step == n_steps:
            break

        rate_half = rate + 0.5 * DT * acceleration
        q_mid = q_of_tau(tau + 0.5 * DT)
        work += q_mid * rate_half * DT
        lam += DT * rate_half

        tau_next = tau + DT
        acceleration_next = q_of_tau(tau_next) - dphi(lam)
        rate = rate_half + 0.5 * DT * acceleration_next
        acceleration = acceleration_next
        tau = tau_next

    same_q_ks_initial_half_cycle = empirical_ks(
        snapshots[0.0]["lambda_sample"], snapshots[0.5 * period]["lambda_sample"]
    )
    same_phase_ks_cycle1_cycle2 = empirical_ks(
        snapshots[period]["lambda_sample"], snapshots[2.0 * period]["lambda_sample"]
    )

    max_energy_work_error = max(
        s["max_abs_characteristic_work_minus_energy"] for s in snapshots.values()
    )

    cycle_rows = []
    for cycle in range(CYCLES + 1):
        target = cycle * period
        snap = snapshots[target]
        cycle_rows.append(
            {
                "cycle": cycle,
                "actual_tau": snap["actual_tau"],
                "snapshot_tail_fraction": snap["snapshot_tail_fraction"],
                "first_passage_fraction": snap["first_passage_fraction"],
                "mean_lambda": snap["mean_lambda"],
                "variance_lambda": snap["variance_lambda"],
                "mean_energy_change": snap["mean_energy_change"],
                "mean_external_work": snap["mean_external_work"],
            }
        )

    result = {
        "classification": "candidate local-traction P0-to-P characteristic audit",
        "status": "candidate reduced model; not exact finite-chain projection",
        "model": {
            "m": M_EXP,
            "n": N_EXP,
            "E_pa": E_PA,
            "stress_mean_pa": STRESS_MEAN_PA,
            "stress_amplitude_pa": STRESS_AMPLITUDE_PA,
            "q_mean": Q_MEAN,
            "q_amplitude": Q_AMPLITUDE,
            "omega_star": OMEGA_STAR,
            "dt": DT,
            "cycles": CYCLES,
            "number_of_characteristics": N_CHARACTERISTICS,
            "diagnostic_P0_support": [LAMBDA_MIN, LAMBDA_MAX],
            "diagnostic_P0_sampling": "uniform quadrature weights over explicit bounded support; not a physical aluminum P0 assumption",
            "initial_rate_condition": "c(lambda0,0)=0",
            "tail_threshold": TAIL_THRESHOLD,
        },
        "checks": {
            "same_applied_q_KS_initial_vs_half_cycle": same_q_ks_initial_half_cycle,
            "same_phase_KS_cycle1_vs_cycle2": same_phase_ks_cycle1_cycle2,
            "max_abs_characteristic_work_minus_energy": max_energy_work_error,
            "first_passage_monotone_at_cycle_ends": all(
                cycle_rows[i + 1]["first_passage_fraction"]
                >= cycle_rows[i]["first_passage_fraction"]
                for i in range(len(cycle_rows) - 1)
            ),
            "snapshot_tail_monotone_at_cycle_ends": all(
                cycle_rows[i + 1]["snapshot_tail_fraction"]
                >= cycle_rows[i]["snapshot_tail_fraction"]
                for i in range(len(cycle_rows) - 1)
            ),
        },
        "cycle_end_results": cycle_rows,
        "interpretation": [
            "Under the candidate local-traction assumption, an initially static F0 is fixed by P0 alone as F0=P0*delta(c).",
            "The supplied stress history and characteristic equations then define a unique push-forward P(t) without ordered atom positions.",
            "The same applied stress at different phases can correspond to different P because characteristic rates retain history internally.",
            "Snapshot tail mass is not cumulative; first-passage mass is cumulative.",
            "The diagnostic bounded P0 is not a proposed physical aluminum initial distribution.",
            "A delta P0 with identical initial rates and identical local stress histories would remain a delta and would not generate a tail by itself.",
        ],
    }

    # Remove arrays retained only for in-memory comparisons.
    for snap in snapshots.values():
        snap.pop("lambda_sample", None)

    output = Path("results/data/local_traction_p0_propagator/summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
