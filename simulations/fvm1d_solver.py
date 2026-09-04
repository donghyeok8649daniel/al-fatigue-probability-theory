"""Small-strain 1D finite-volume bar solver with the desktop CSV contract.

The cell balance is written in terms of face tractions.  For the present
uniform bar and prescribed end stress, equilibrium makes the stress constant
over every cell; displacement is reconstructed from the cell strain.  This is
an initial FVM backend, not a fatigue or plasticity constitutive law.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def run(
    *,
    elements: int = 40,
    length_m: float = 0.05,
    area_m2: float = 1.0e-5,
    young_pa: float = 69.0e9,
    stress_mean_mpa: float = 50.0,
    stress_amplitude_mpa: float = 100.0,
    frequency_hz: float = 20.0,
    cycles: int = 2,
    steps_per_cycle: int = 80,
    outdir: Path = Path("fvm1d_output"),
) -> None:
    if elements < 1 or length_m <= 0 or area_m2 <= 0 or young_pa <= 0:
        raise ValueError("elements, length_m, area_m2, and young_pa must be positive")
    if stress_amplitude_mpa < 0 or frequency_hz <= 0 or cycles < 1 or steps_per_cycle < 2:
        raise ValueError("invalid loading controls")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dx = length_m / elements
    total_steps = cycles * steps_per_cycle
    node_path = outdir / "nodes.csv"
    elem_path = outdir / "elements.csv"
    meta_path = outdir / "metadata.csv"

    with node_path.open("w", newline="", encoding="utf-8") as nf, elem_path.open(
        "w", newline="", encoding="utf-8"
    ) as ef:
        nw = csv.writer(nf)
        ew = csv.writer(ef)
        nw.writerow(["time_s", "step", "node", "x_m", "displacement_m", "applied_stress_pa"])
        ew.writerow(["time_s", "step", "element", "x_mid_m", "strain", "stress_pa", "applied_stress_pa"])
        for step in range(total_steps + 1):
            t = step / steps_per_cycle / frequency_hz
            stress = (stress_mean_mpa + stress_amplitude_mpa * math.sin(2.0 * math.pi * frequency_hz * t)) * 1.0e6
            strain = stress / young_pa
            for node in range(elements + 1):
                x = node * dx
                nw.writerow([f"{t:.17g}", step, node, f"{x:.17g}", f"{strain*x:.17g}", f"{stress:.17g}"])
            for element in range(elements):
                ew.writerow([
                    f"{t:.17g}", step, element, f"{(element + 0.5)*dx:.17g}",
                    f"{strain:.17g}", f"{stress:.17g}", f"{stress:.17g}"
                ])

    with meta_path.open("w", newline="", encoding="utf-8") as mf:
        mw = csv.writer(mf)
        mw.writerows([
            ("solver", "quasistatic_linear_bar_fvm"),
            ("elements", elements),
            ("length_m", f"{length_m:.17g}"),
            ("area_m2", f"{area_m2:.17g}"),
            ("young_pa", f"{young_pa:.17g}"),
            ("frequency_hz", f"{frequency_hz:.17g}"),
            ("cycles", cycles),
            ("steps_per_cycle", steps_per_cycle),
        ])


def main() -> None:
    parser = argparse.ArgumentParser(description="Initial 1D finite-volume bar solver")
    parser.add_argument("--elements", type=int, default=40)
    parser.add_argument("--length-m", type=float, default=0.05)
    parser.add_argument("--area-m2", type=float, default=1.0e-5)
    parser.add_argument("--young-pa", type=float, default=69.0e9)
    parser.add_argument("--stress-mean-mpa", type=float, default=50.0)
    parser.add_argument("--stress-amplitude-mpa", type=float, default=100.0)
    parser.add_argument("--frequency-hz", type=float, default=20.0)
    parser.add_argument("--cycles", type=int, default=2)
    parser.add_argument("--steps-per-cycle", type=int, default=80)
    parser.add_argument("--outdir", type=Path, default=Path("fvm1d_output"))
    args = parser.parse_args()
    run(**vars(args))


if __name__ == "__main__":
    main()
