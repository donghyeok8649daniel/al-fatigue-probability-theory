"""Generate FCC normal generalized-LJ calibration/validation results."""

from __future__ import annotations

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from theory.fcc_normal_lj import EV_J, FCCNormalLJ, FCCNormalLJParameters

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"

A_LAT_M = 4.0495e-10
COHESIVE_EV = 3.43
C11_REFERENCE_PA = 107.0e9
C12_REFERENCE_PA = 61.0e9
C44_REFERENCE_PA = 29.0e9
DFT_IDEAL_001_GPA = 10.63

E001_REFERENCE_PA = (
    (C11_REFERENCE_PA - C12_REFERENCE_PA)
    * (C11_REFERENCE_PA + 2.0 * C12_REFERENCE_PA)
    / (C11_REFERENCE_PA + C12_REFERENCE_PA)
)
NU001_REFERENCE = C12_REFERENCE_PA / (C11_REFERENCE_PA + C12_REFERENCE_PA)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    parameters = FCCNormalLJParameters(
        lattice_constant_m=A_LAT_M,
        repulsive_exponent=12.19,
        attractive_exponent=6.0,
        cutoff_lattice_constants=15.0,
    )
    model = FCCNormalLJ(parameters)

    epsilon_cohesive = model.epsilon_for_cohesive_energy(COHESIVE_EV * EV_J)
    E_cohesive, nu_cohesive = model.small_strain_properties(epsilon_cohesive)

    epsilon_normal = model.epsilon_for_youngs_modulus(E001_REFERENCE_PA)
    E_normal, nu_normal = model.small_strain_properties(epsilon_normal)
    cohesive_normal_ev = model.cohesive_energy_j_per_atom(epsilon_normal) / EV_J
    C11, C12, C44 = model.cubic_elastic_constants_pa(epsilon_normal)

    stretch = np.linspace(1.0, 1.50, 101)
    _, transverse_cohesive, stress_cohesive = model.stress_strain_curve(
        epsilon_cohesive, stretch
    )
    _, transverse_normal, stress_normal = model.stress_strain_curve(
        epsilon_normal, stretch
    )
    peak = int(np.argmax(stress_normal))

    cutoff_rows = []
    for cutoff in (5.0, 6.0, 8.0, 10.0, 12.0, 15.0):
        cutoff_model = FCCNormalLJ(
            FCCNormalLJParameters(
                lattice_constant_m=A_LAT_M,
                repulsive_exponent=12.19,
                attractive_exponent=6.0,
                cutoff_lattice_constants=cutoff,
            )
        )
        eps_coh = cutoff_model.epsilon_for_cohesive_energy(COHESIVE_EV * EV_J)
        E_coh, nu_coh = cutoff_model.small_strain_properties(eps_coh)
        eps_normal = cutoff_model.epsilon_for_youngs_modulus(E001_REFERENCE_PA)
        coh_normal_ev = cutoff_model.cohesive_energy_j_per_atom(eps_normal) / EV_J
        curve_stretch = np.linspace(1.0, 1.40, 81)
        _, _, curve_stress = cutoff_model.stress_strain_curve(
            eps_normal, curve_stretch
        )
        cutoff_rows.append(
            {
                "cutoff_lattice_constants": cutoff,
                "E001_after_cohesive_fit_GPa": E_coh / 1.0e9,
                "nu001_after_cohesive_fit": nu_coh,
                "cohesive_after_E001_fit_eV_per_atom": coh_normal_ev,
                "ideal_strength_after_E001_fit_GPa": float(np.max(curve_stress) / 1.0e9),
            }
        )

    summary = {
        "status": "normal-only FCC pair-potential validation; not a fatigue-life prediction",
        "fcc_lattice_constant_A": A_LAT_M * 1.0e10,
        "atomic_reference_volume_m3": parameters.atomic_volume_m3,
        "repulsive_exponent_m": parameters.repulsive_exponent,
        "attractive_exponent_n": parameters.attractive_exponent,
        "lattice_sum_cutoff_a_lat": parameters.cutoff_lattice_constants,
        "sigma_lj_A": model.sigma_lj_m * 1.0e10,
        "directional_experimental_reference": {
            "C11_GPa": C11_REFERENCE_PA / 1.0e9,
            "C12_GPa": C12_REFERENCE_PA / 1.0e9,
            "C44_GPa": C44_REFERENCE_PA / 1.0e9,
            "E001_GPa": E001_REFERENCE_PA / 1.0e9,
            "nu001": NU001_REFERENCE,
            "cohesive_energy_eV_per_atom": COHESIVE_EV,
            "DFT_001_ideal_tensile_strength_GPa": DFT_IDEAL_001_GPA,
        },
        "cohesive_energy_calibration": {
            "epsilon_lj_eV": epsilon_cohesive / EV_J,
            "predicted_E001_GPa": E_cohesive / 1.0e9,
            "predicted_nu001": nu_cohesive,
        },
        "normal_E001_calibration": {
            "epsilon_lj_eV": epsilon_normal / EV_J,
            "predicted_E001_GPa": E_normal / 1.0e9,
            "predicted_nu001": nu_normal,
            "predicted_cohesive_energy_eV_per_atom": cohesive_normal_ev,
            "C11_GPa": C11 / 1.0e9,
            "C12_GPa": C12 / 1.0e9,
            "C44_GPa": C44 / 1.0e9,
            "C12_minus_C44_MPa": (C12 - C44) / 1.0e6,
            "ideal_engineering_strength_GPa": stress_normal[peak] / 1.0e9,
            "ideal_engineering_strain": stretch[peak] - 1.0,
            "transverse_stretch_at_peak": transverse_normal[peak],
        },
        "cutoff_convergence": cutoff_rows,
    }

    (DATA / "fcc_normal_lj_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    with (DATA / "fcc_normal_lj_stress_strain.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "axial_stretch",
                "engineering_strain",
                "stress_cohesive_fit_GPa",
                "stress_E001_fit_GPa",
                "transverse_stretch_cohesive_fit",
                "transverse_stretch_E001_fit",
            ]
        )
        for row in zip(
            stretch,
            stretch - 1.0,
            stress_cohesive / 1.0e9,
            stress_normal / 1.0e9,
            transverse_cohesive,
            transverse_normal,
        ):
            writer.writerow(row)

    plt.figure(figsize=(8, 5))
    plt.plot((stretch - 1.0) * 100.0, stress_normal / 1.0e9, label="fit a_lat + E[001]")
    plt.plot((stretch - 1.0) * 100.0, stress_cohesive / 1.0e9, label="fit a_lat + cohesive energy")
    plt.axhline(DFT_IDEAL_001_GPA, linestyle="--", label="DFT [001] ideal-strength reference")
    plt.xlabel("Axial engineering strain (%)")
    plt.ylabel("Axial engineering stress (GPa)")
    plt.title("FCC generalized-LJ [001] normal tension")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "fcc_lj_001_stress_strain.svg")
    plt.close()

    width = 0.35
    names = ["C11", "C12", "C44"]
    x = np.arange(3)
    plt.figure(figsize=(7, 5))
    plt.bar(x - width / 2, [C11 / 1e9, C12 / 1e9, C44 / 1e9], width, label="LJ, E[001]-fit")
    plt.bar(x + width / 2, [107.0, 61.0, 29.0], width, label="experimental reference")
    plt.xticks(x, names)
    plt.ylabel("Elastic constant (GPa)")
    plt.title("Normal elasticity and central-pair Cauchy constraint")
    plt.grid(True, axis="y")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "fcc_lj_elastic_constants.svg")
    plt.close()

    x = np.arange(2)
    plt.figure(figsize=(7, 5))
    plt.bar(x - width / 2, [E_cohesive / 1e9, cohesive_normal_ev], width, label="LJ cross-prediction")
    plt.bar(x + width / 2, [E001_REFERENCE_PA / 1e9, COHESIVE_EV], width, label="reference")
    plt.xticks(x, ["E[001] if cohesion fitted\n(GPa)", "cohesion if E[001] fitted\n(eV/atom)"])
    plt.title("Normal stiffness–cohesion incompatibility")
    plt.grid(True, axis="y")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "fcc_lj_cohesion_conflict.svg")
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(
        [row["cutoff_lattice_constants"] for row in cutoff_rows],
        [row["E001_after_cohesive_fit_GPa"] for row in cutoff_rows],
        marker="o",
    )
    plt.xlabel("Lattice-sum cutoff (a_lat)")
    plt.ylabel("E[001] after cohesive-energy fit (GPa)")
    plt.title("FCC lattice-sum convergence")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIGURES / "fcc_lj_cutoff_convergence.svg")
    plt.close()

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
