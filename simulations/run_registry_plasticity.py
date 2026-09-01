"""Generate verified multilayer-energy and reduced-plasticity demonstrations."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from theory.registry_lattice import (
    MultilayerPotentialParameters,
    RegistryLattice,
    dU_ds,
    h_q_bessel,
    h_q_direct,
    h_q_polylog,
    preferred_registry,
    registry_energy,
    registry_energy_derivative,
    slip_barrier,
    u0,
    v_slip,
)
from theory.registry_plasticity import RegistryTransportConfig, solve_registry


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "registry_plasticity"
FIGURES = ROOT / "results" / "figures" / "registry_plasticity"


def _pulse(time: np.ndarray) -> np.ndarray:
    force = np.zeros_like(time)
    ramp_up = (time >= 2.0) & (time < 4.0)
    force[ramp_up] = 0.55 * (time[ramp_up] - 2.0) / 2.0
    force[(time >= 4.0) & (time < 8.0)] = 0.55
    ramp_down = (time >= 8.0) & (time < 12.0)
    force[ramp_down] = 0.55 * (12.0 - time[ramp_down]) / 4.0
    return force


def _write_history(path: Path, history) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "time_reduced",
                "resolved_force_reduced",
                "mean_registry_s_over_b",
                "mean_well_index",
                "mean_intrawell_registry",
                "variance",
                "mean_lattice_energy_over_epsilon_c",
                "work_over_epsilon_c",
                "entropy_production_reduced",
                "cumulative_hysteresis_energy_reduced",
                "plastic_shear_strain_for_b_over_h_equal_1",
                "plastic_tensile_strain_for_M_equal_1",
                "edge_probability",
            ]
        )
        for row in zip(
            history.time,
            history.generalized_force,
            history.mean_registry,
            history.mean_well_index,
            history.mean_intrawell_registry,
            history.variance,
            history.mean_lattice_energy,
            history.work,
            history.entropy_production,
            history.hysteresis_energy,
            history.plastic_shear_strain,
            history.plastic_tensile_strain,
            history.boundary_probability,
        ):
            writer.writerow([float(value) for value in row])


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    lattice = RegistryLattice(normal_ratio=1.0, bessel_modes=20, layer_modes=48)
    config = RegistryTransportConfig(
        lattice=lattice,
        inverse_temperature=20.0,
        u_min=-6.0,
        u_max=7.0,
        cells=390,
    )
    time = np.linspace(0.0, 30.0, 301)
    force = _pulse(time)
    pulse = solve_registry(time, force, config, max_dt=0.025)
    # Six complete, exactly antisymmetric cycles on the same 0--30 interval.
    cyclic_force = 0.55 * np.sin(2.0 * np.pi * time / 5.0)
    symmetric = solve_registry(time, cyclic_force, config, max_dt=0.025)

    phase = np.linspace(0.0, 1.0, 2001)
    energy = np.asarray(registry_energy(phase, lattice))
    traction = np.asarray(registry_energy_derivative(phase, lattice))
    critical_force = float(np.max(np.abs(traction)))
    _write_history(DATA / "resolved_shear_pulse.csv", pulse)
    _write_history(DATA / "symmetric_cycle.csv", symmetric)
    summary = {
        "status": "dimensionless multiplicity-free multilayer mechanism demonstration",
        "counting_convention": "U0(a,s)=sum_{k>=1} W(k*a,s), no k multiplicity, common s",
        "m": lattice.m,
        "n": lattice.n,
        "a_over_b": lattice.normal_ratio,
        "sigma_LJ_over_b": lattice.sigma_ratio,
        "inverse_temperature_epsilon_c_over_kBT": config.inverse_temperature,
        "bessel_modes": lattice.bessel_modes,
        "bessel_lambert_layer_modes": lattice.layer_modes,
        "preferred_registry_s_over_b": preferred_registry(lattice),
        "ideal_reduced_registry_force": critical_force,
        "pulse_peak_reduced_force": float(np.max(force)),
        "pulse_final_mean_well_index": float(pulse.mean_well_index[-1]),
        "pulse_final_mean_intrawell_registry": float(
            pulse.mean_intrawell_registry[-1]
        ),
        "pulse_final_work_over_epsilon_c": float(pulse.work[-1]),
        "pulse_final_hysteresis_energy_over_epsilon_c": float(
            pulse.hysteresis_energy[-1]
        ),
        "pulse_max_edge_probability": float(np.max(pulse.boundary_probability)),
        "symmetric_final_mean_well_index": float(
            symmetric.mean_well_index[-1]
        ),
        "physical_plastic_shear_mapping": "gamma_p=(b/h_slip)*mean_well_index",
        "physical_axial_mapping": "epsilon_p=M_schmid*gamma_p",
        "limitations": [
            "not calibrated to aluminum",
            "one ideal registry; no dislocation storage or hardening",
            "this transport demonstration prescribes normal separation; the derived 2D PDE is not yet solved here",
            "mobility and representative interface area require atomistic calibration",
        ],
    }
    (DATA / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.6), constrained_layout=True)
    axes[0, 0].plot(phase, energy - np.min(energy), color="black")
    axes[0, 0].plot(phase, -0.55 * phase + energy - np.min(energy), color="tab:red")
    axes[0, 0].set(
        xlabel=r"registry $s/b$",
        ylabel=r"energy $/\varepsilon_c$",
        title="Exact multilayer Bessel--Lambert landscape",
    )
    axes[0, 0].legend(["zero load", "tilted by peak resolved force"])

    mesh = axes[0, 1].pcolormesh(
        pulse.time,
        pulse.registry,
        pulse.density.T,
        shading="auto",
        cmap="magma",
    )
    axes[0, 1].set(
        xlabel="reduced time",
        ylabel=r"unwrapped registry $s/b$",
        title="Probability crosses lattice wells",
    )
    figure.colorbar(mesh, ax=axes[0, 1], label="density")

    axes[1, 0].plot(pulse.time, pulse.generalized_force, label="resolved force")
    axes[1, 0].plot(pulse.time, pulse.mean_well_index, label=r"$\langle z\rangle$")
    axes[1, 0].plot(
        pulse.time,
        pulse.mean_intrawell_registry,
        label="mean intrawell registry",
    )
    axes[1, 0].set(
        xlabel="reduced time",
        title="Unloading recovers intrawell motion but not well population",
    )
    axes[1, 0].legend()

    axes[1, 1].plot(
        pulse.mean_registry, pulse.generalized_force, color="tab:blue", label="biased pulse"
    )
    axes[1, 1].plot(
        symmetric.mean_registry,
        symmetric.generalized_force,
        color="tab:orange",
        alpha=0.8,
        label="symmetric cycles",
    )
    axes[1, 1].set(
        xlabel=r"mean unwrapped registry $\langle s/b\rangle$",
        ylabel="reduced resolved force",
        title="Registry hysteresis and residual translation",
    )
    axes[1, 1].legend()
    figure.savefig(FIGURES / "active_registry_plasticity.png", dpi=180)
    plt.close(figure)

    # Intrinsic U0 surface, fixed-a excess, slip barrier and exact derivative.
    physical = MultilayerPotentialParameters(
        b=1.0, epsilon_lj=1.0, sigma_lj=0.82,
        m=12.0, n=6.0, bessel_modes=20, layer_modes=48,
    )
    eta_values = np.linspace(0.72, 1.55, 90)
    delta_values = np.linspace(0.0, 1.0, 121)
    surface = np.array([
        np.asarray(u0(eta, delta_values, physical)) for eta in eta_values
    ])
    barriers = np.array([slip_barrier(eta, physical, samples=1024) for eta in eta_values])
    fixed_eta = 1.0
    slip_excess = np.asarray(v_slip(fixed_eta, delta_values, physical, s0=0.0))
    traction_u = np.asarray(dU_ds(fixed_eta, delta_values, physical))
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.8), constrained_layout=True)
    mesh = axes[0, 0].pcolormesh(
        delta_values, eta_values, surface, shading="auto", cmap="viridis"
    )
    figure.colorbar(mesh, ax=axes[0, 0], label=r"$U_0/(\mathrm{energy})$")
    axes[0, 0].set(xlabel=r"$s/b$", ylabel=r"$a/b$", title=r"Intrinsic $U_0(a,s)$")
    axes[0, 1].plot(delta_values, slip_excess)
    axes[0, 1].set(xlabel=r"$s/b$", ylabel=r"$V_{\rm slip}$", title=r"Fixed-$a$ slip excess")
    axes[1, 0].plot(eta_values, barriers)
    axes[1, 0].set(xlabel=r"$a/b$", ylabel="barrier energy", title="Slip barrier versus opening")
    axes[1, 1].plot(delta_values, traction_u)
    axes[1, 1].set(xlabel=r"$s/b$", ylabel=r"$\partial_s U_0$", title="Exact intrinsic registry force")
    figure.savefig(FIGURES / "multilayer_energy_landscape.png", dpi=180)
    plt.close(figure)

    # Independent direct/Bessel/polylog and truncation checks for q=6,12.
    verification_rows = []
    for q in (6, 12):
        for delta, eta in ((0.0, 0.7), (0.23, 1.0), (0.5, 1.4)):
            direct = h_q_direct(q, delta, eta, 400, 800)
            bessel = h_q_bessel(q, delta, eta, 20, 48)
            polylog = h_q_polylog(q, delta, eta, 20)
            verification_rows.append((q, delta, eta, direct, bessel, polylog))
    with (DATA / "multilayer_representation_checks.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "q", "delta", "eta", "direct", "bessel_lambert", "polylog",
            "relative_error_direct_bessel", "relative_error_polylog_bessel",
        ])
        for q, delta, eta, direct, bessel, polylog in verification_rows:
            writer.writerow([
                q, delta, eta, direct, bessel, polylog,
                abs(direct - bessel) / abs(bessel),
                abs(polylog - bessel) / abs(bessel),
            ])
    cutoffs = np.array([4, 8, 16, 32, 64, 128, 256, 400])
    reference = h_q_bessel(6, 0.23, 1.0, 24, 64)
    direct_errors = np.array([
        abs(h_q_direct(6, 0.23, 1.0, int(c), int(2 * c)) - reference) / abs(reference)
        for c in cutoffs
    ])
    ell_cutoffs = np.array([1, 2, 3, 4, 6, 8, 12, 20])
    reciprocal_errors = np.array([
        abs(h_q_bessel(6, 0.23, 1.0, int(c), 48) - reference) / abs(reference)
        for c in ell_cutoffs
    ])
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)
    for q in (6, 12):
        rows = [row for row in verification_rows if row[0] == q]
        axes[0].semilogy(
            range(len(rows)),
            [abs(row[3] - row[4]) / abs(row[4]) for row in rows],
            "o-", label=f"q={q}",
        )
    axes[0].set(
        xticks=range(3), xticklabels=["(0,.7)", "(.23,1)", "(.5,1.4)"],
        xlabel=r"$(\delta,\eta)$", ylabel="relative error",
        title="Direct double sum vs exact representation",
    )
    axes[0].legend()
    axes[1].loglog(cutoffs, direct_errors, "o-", label=r"direct $k,p$")
    axes[1].loglog(ell_cutoffs, reciprocal_errors, "s-", label=r"reciprocal $\ell$")
    axes[1].set(xlabel="truncation count", ylabel="relative error", title="Numerical truncation convergence")
    axes[1].legend()
    figure.savefig(FIGURES / "multilayer_representation_verification.png", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
