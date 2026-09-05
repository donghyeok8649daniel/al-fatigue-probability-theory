"""Small-strain axial finite-volume solver with 3D cylinder kinematics.

The cell balance is written in terms of face tractions.  For the present
uniform bar and prescribed end stress, equilibrium makes the axial stress
constant over every cell; displacement is reconstructed from axial strain.
Transverse size is a kinematic Poisson post-process with zero applied
transverse stress. This is not a multiaxial constitutive or failure model.
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
    diameter_m: float = 0.006,
    poisson_ratio: float = 0.33,
    tensile_axis: tuple[float, float, float] = (1.0, 0.0, 0.0),
    stress_mean_mpa: float = 50.0,
    stress_amplitude_mpa: float = 100.0,
    frequency_hz: float = 20.0,
    cycles: int = 2,
    steps_per_cycle: int = 80,
    outdir: Path = Path("fvm1d_output"),
) -> None:
    if elements < 1 or length_m <= 0 or area_m2 <= 0 or young_pa <= 0 or diameter_m <= 0:
        raise ValueError("elements, length_m, area_m2, young_pa, and diameter_m must be positive")
    if stress_amplitude_mpa < 0 or frequency_hz <= 0 or cycles < 1 or steps_per_cycle < 2:
        raise ValueError("invalid loading controls")
    if not -1.0 < poisson_ratio < 0.5:
        raise ValueError("poisson_ratio must satisfy -1 < nu < 0.5")
    axis_norm = math.sqrt(sum(float(value) ** 2 for value in tensile_axis))
    if not math.isfinite(axis_norm) or axis_norm <= 0.0:
        raise ValueError("tensile_axis must be a finite nonzero vector")
    axis = tuple(float(value) / axis_norm for value in tensile_axis)

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
        nw.writerow([
            "time_s", "step", "node", "x_m", "displacement_m", "applied_stress_pa",
            "position_x_m", "position_y_m", "position_z_m",
            "displacement_x_m", "displacement_y_m", "displacement_z_m",
        ])
        ew.writerow([
            "time_s", "step", "element", "x_mid_m", "strain", "stress_pa",
            "applied_stress_pa", "transverse_strain", "diameter_m", "transverse_stress_pa",
        ])
        for step in range(total_steps + 1):
            t = step / steps_per_cycle / frequency_hz
            stress = (stress_mean_mpa + stress_amplitude_mpa * math.sin(2.0 * math.pi * frequency_hz * t)) * 1.0e6
            strain = stress / young_pa
            transverse_strain = -poisson_ratio * strain
            current_diameter = diameter_m * (1.0 + transverse_strain)
            for node in range(elements + 1):
                x = node * dx
                displacement = strain * x
                position = tuple(component * x for component in axis)
                displacement_vector = tuple(component * displacement for component in axis)
                nw.writerow([
                    f"{t:.17g}", step, node, f"{x:.17g}", f"{displacement:.17g}",
                    f"{stress:.17g}", *(f"{value:.17g}" for value in position),
                    *(f"{value:.17g}" for value in displacement_vector),
                ])
            for element in range(elements):
                ew.writerow([
                    f"{t:.17g}", step, element, f"{(element + 0.5)*dx:.17g}",
                    f"{strain:.17g}", f"{stress:.17g}", f"{stress:.17g}",
                    f"{transverse_strain:.17g}", f"{current_diameter:.17g}", "0"
                ])

    with meta_path.open("w", newline="", encoding="utf-8") as mf:
        mw = csv.writer(mf)
        mw.writerows([
            ("solver", "quasistatic_linear_bar_fvm"),
            ("elements", elements),
            ("length_m", f"{length_m:.17g}"),
            ("area_m2", f"{area_m2:.17g}"),
            ("young_pa", f"{young_pa:.17g}"),
            ("diameter_m", f"{diameter_m:.17g}"),
            ("poisson_ratio", f"{poisson_ratio:.17g}"),
            ("tensile_axis_x", f"{axis[0]:.17g}"),
            ("tensile_axis_y", f"{axis[1]:.17g}"),
            ("tensile_axis_z", f"{axis[2]:.17g}"),
            ("transverse_applied_stress_pa", "0"),
            ("transverse_kinematics", "epsilon_transverse=-nu*epsilon_axial"),
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
    parser.add_argument("--diameter-m", type=float, default=0.006)
    parser.add_argument("--poisson-ratio", type=float, default=0.33)
    parser.add_argument("--tensile-axis", type=float, nargs=3, default=(1.0, 0.0, 0.0))
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
