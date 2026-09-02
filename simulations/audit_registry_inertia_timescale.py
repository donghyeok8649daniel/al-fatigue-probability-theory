"""Audit the physical scale of the active reduced registry inertia.

No FCC geometry is used. The calculation uses only:
- the current multiplicity-free row/layer U0(a,s);
- its per-reference-repeat counting convention;
- the existing normal calibration a0, E, A0 and t0.

The direct (k,p) sums below are numerical convergence controls, not physical
cutoffs. The purpose is to evaluate Uaa/Uss at the normalized equilibrium and
map their ratio to a physical registry frequency without freely tuning mu_s.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

M = 12.0
N = 6.0
B = 1.0
SIGMA = 1.0
EPSILON = 1.0
S0 = 0.5
A0_STAR = 0.9919601753795769

# Existing normal calibration retained in theory-core.
A0_PHYS = 2.8627442948e-10  # m, physical equilibrium normal spacing
YOUNGS_MODULUS = 69.0e9     # Pa
REFERENCE_AREA = 6.0338e-20 # m^2
T0 = 5.55046e-14            # s

C_MN = M / (M - N) * (M / N) ** (N / (M - N))


def curvatures(a: float, s: float, kmax: int, pmax: int) -> tuple[float, float]:
    k = np.arange(1, kmax + 1, dtype=float)[:, None]
    p = np.arange(-pmax, pmax + 1, dtype=float)[None, :]
    x = k * a
    y = p * B + s
    r2 = x * x + y * y

    def second(q: float, coordinate: str) -> float:
        if coordinate == "a":
            c2 = k * k
            z2 = x * x
        elif coordinate == "s":
            c2 = 1.0
            z2 = y * y
        else:
            raise ValueError("coordinate must be a or s")
        return float(np.sum(
            -q * c2 * r2 ** (-0.5 * (q + 2.0))
            + q * (q + 2.0) * c2 * z2 * r2 ** (-0.5 * (q + 4.0))
        ))

    uaa = C_MN * EPSILON * (
        SIGMA**M * second(M, "a") - SIGMA**N * second(N, "a")
    )
    uss = C_MN * EPSILON * (
        SIGMA**M * second(M, "s") - SIGMA**N * second(N, "s")
    )
    return uaa, uss


def main() -> None:
    convergence = []
    for kmax, pmax in [(20, 50), (40, 100), (80, 200), (120, 300), (200, 500)]:
        uaa, uss = curvatures(A0_STAR, S0, kmax, pmax)
        convergence.append({
            "kmax": kmax,
            "pmax": pmax,
            "Uaa": uaa,
            "Uss": uss,
            "curvature_ratio": uss / uaa,
        })

    uaa = convergence[-1]["Uaa"]
    uss = convergence[-1]["Uss"]
    ratio = convergence[-1]["curvature_ratio"]

    repeat_mass = T0 * T0 * YOUNGS_MODULUS * REFERENCE_AREA / A0_PHYS
    axial_stiffness = YOUNGS_MODULUS * REFERENCE_AREA / A0_PHYS
    registry_stiffness = axial_stiffness * ratio

    inertia_cases = []
    for rho in (0.5, 1.0):
        fs = math.sqrt(ratio / rho) / (2.0 * math.pi * T0)
        inertia_cases.append({
            "inertia_ratio_mu_over_mrepeat": rho,
            "registry_frequency_hz": fs,
            "principal_loading_frequency_hz": 2.0 * fs,
        })

    resonance = []
    for fload in (1.0, 10.0, 20.0, 100.0, 1000.0):
        rho_req = ratio / (math.pi * fload * T0) ** 2
        resonance.append({
            "loading_frequency_hz": fload,
            "required_inertia_ratio": rho_req,
            "required_registry_inertia_kg": rho_req * repeat_mass,
        })

    out = Path("results/data/registry_inertia_timescale")
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in [
        ("curvature_convergence.csv", convergence),
        ("natural_frequency.csv", inertia_cases),
        ("principal_resonance_required_inertia.csv", resonance),
    ]:
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    print("Uaa", uaa)
    print("Uss", uss)
    print("Uss/Uaa", ratio)
    print("repeat mass [kg]", repeat_mass)
    print("axial stiffness [N/m]", axial_stiffness)
    print("registry stiffness [N/m]", registry_stiffness)
    print("registry f, rho=1 [Hz]", inertia_cases[-1]["registry_frequency_hz"])
    print("20 Hz required rho", resonance[2]["required_inertia_ratio"])


if __name__ == "__main__":
    main()
