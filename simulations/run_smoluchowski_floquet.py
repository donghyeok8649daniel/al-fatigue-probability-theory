"""Generate the dimensionless one-cycle survival-spectrum demonstration.

The script does not fit aluminum lifetime.  It evaluates consequences of the
already declared normal-only Smoluchowski model and its tangent-instability
absorbing boundary.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from theory.normal_lj_chain import normalized_lj_energy
from theory.smoluchowski_escape import (
    TransportConfig,
    conditional_equilibrium,
    solve,
    transport_grid,
)
from theory.smoluchowski_floquet import (
    asymptotic_survival_prefactor,
    dense_cycle_spectrum,
    direct_cycle_survival_ratios,
    frozen_principal_escape_rate,
    principal_survival_mode,
)


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures" / "smoluchowski_floquet"
DATA = ROOT / "results" / "data" / "smoluchowski_floquet"


def protocol(period: float, points: int = 101,
             mean_force: float = 0.008, amplitude: float = 0.007):
    time = np.linspace(0.0, period, points)
    force = mean_force + amplitude * np.sin(2.0 * np.pi * time / period)
    return time, force


def absorbing_config(cells: int = 140, inverse_temperature: float = 2000.0):
    return TransportConfig(
        inverse_temperature=inverse_temperature,
        cells=cells,
        boundary="absorbing",
        initiation_definition="tangent_instability",
    )


def phase_observables(result, c):
    x, dx = transport_grid(c)
    p = result.phase_conditional_density
    mean = np.sum(p * x[None, :], axis=1) * dx
    variance = np.sum(p * (x[None, :] - mean[:, None]) ** 2, axis=1) * dx
    shifted = normalized_lj_energy(x, c.m, c.n) - normalized_lj_energy(1.0, c.m, c.n)
    energy = np.sum(p * shifted[None, :], axis=1) * dx
    hazard = result.phase_hazard.copy()
    return x, mean, variance, energy, hazard


def write_rows(path: Path, header, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    period = 12.0
    time, force = protocol(period)
    c = absorbing_config()
    result = principal_survival_mode(
        time, force, c, max_dt=0.03, tolerance=1e-11)
    x, mean, variance, energy, hazard = phase_observables(result, c)

    # Direct iteration starts from the load-conditioned basin density used by
    # the original demonstration, not from the eigenmode.
    _, dx = transport_grid(c)
    generic_initial = conditional_equilibrium(x, dx, float(force[0]), c)
    direct_survival, direct_ratios = direct_cycle_survival_ratios(
        generic_initial, time, force, c, cycles=12, max_dt=0.03)
    eigen_survival, eigen_ratios = direct_cycle_survival_ratios(
        result.start_density, time, force, c, cycles=12, max_dt=0.03)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    axes[0].semilogy(np.arange(1, result.multiplier_history.size + 1),
                     np.abs(result.multiplier_history - result.multiplier) + 1e-18, "o-")
    axes[0].set(xlabel="power iteration", ylabel=r"$|r_j-r|$",
                title="Principal multiplier convergence")
    cycles = np.arange(1, direct_ratios.size + 1)
    axes[1].plot(cycles, direct_ratios, "o-", label="generic initial density")
    axes[1].plot(cycles, eigen_ratios, "s--", label="principal mode")
    axes[1].axhline(result.multiplier, color="black", lw=1, label=r"$r$")
    axes[1].set(xlabel="cycle", ylabel=r"$S_{k+1}/S_k$",
                title="Direct evolution approaches spectrum")
    axes[1].legend(fontsize=8)
    cycle_edges = np.arange(eigen_survival.size)
    axes[2].semilogy(cycle_edges, eigen_survival, "o", label="direct PDE")
    axes[2].semilogy(cycle_edges, result.multiplier ** cycle_edges, "-",
                    label=r"$r^N$")
    axes[2].set(xlabel="cycle", ylabel="survival", title="Eigenmode geometric survival")
    axes[2].legend()
    fig.suptitle("One-cycle operator verification; dimensionless, not an Al life fit")
    fig.savefig(FIG / "operator_verification.png", dpi=180)
    plt.close(fig)

    dense_c = absorbing_config(cells=70)
    dense_spectrum = dense_cycle_spectrum(time, force, dense_c, max_dt=0.03)
    dense_x, dense_dx = transport_grid(dense_c)
    dense_initial = conditional_equilibrium(
        dense_x, dense_dx, float(force[0]), dense_c)
    dense_prefactor = asymptotic_survival_prefactor(
        dense_initial, dense_spectrum, dense_c)
    dense_survival, _ = direct_cycle_survival_ratios(
        dense_initial, time, force, dense_c, cycles=10, max_dt=0.03)
    dense_cycles = np.arange(dense_survival.size)
    scaled_survival = dense_survival / dense_spectrum.multiplier ** dense_cycles
    eigenvalue_moduli = np.sort(np.abs(np.linalg.eigvals(
        dense_spectrum.operator)))[::-1]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    axes[0].plot(dense_x, dense_spectrum.right_density,
                 label="right conditional density")
    left_axis = axes[0].twinx()
    left_axis.plot(dense_x, dense_spectrum.left_survival_weight,
                   color="tab:red", label="left survival weight")
    axes[0].set(xlabel=r"initial $\lambda$", ylabel="right density",
                title="Biorthogonal Perron modes")
    left_axis.set_ylabel("left survival weight", color="tab:red")
    axes[1].semilogy(np.arange(1, 9),
                     eigenvalue_moduli[:8] / eigenvalue_moduli[0], "o-")
    axes[1].set(xlabel="eigenvalue rank", ylabel=r"$|r_j|/r$",
                title="Cycle-spectrum contraction")
    axes[2].plot(dense_cycles, scaled_survival, "o-", label=r"$S_N/r^N$")
    axes[2].axhline(dense_prefactor, color="black", lw=1,
                    label=r"$C(\rho_0)=\langle w_0,\rho_0\rangle$")
    axes[2].set(xlabel="cycle", ylabel="scaled survival",
                title="Initial-state prefactor is predicted")
    axes[2].legend(fontsize=8)
    fig.suptitle("Left mode and transient spectrum; no fitted life prefactor")
    fig.savefig(FIG / "biorthogonal_survival_spectrum.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    image = axes[0, 0].pcolormesh(
        result.phase_time / period, x, result.phase_conditional_density.T,
        shading="auto")
    fig.colorbar(image, ax=axes[0, 0], label=r"$q(\lambda,\theta)$")
    axes[0, 0].set(xlabel="cycle phase", ylabel=r"$\lambda$",
                   title="Periodic conditional intact density")
    ax_load = axes[0, 1]
    ax_load.plot(result.phase_time / period, force, label="reduced force")
    ax_survival = ax_load.twinx()
    ax_survival.plot(result.phase_time / period, result.phase_survival,
                     color="tab:red", label="within-cycle survival")
    ax_load.set(xlabel="cycle phase", ylabel="reduced force",
                title="Sub-Markov mass loss within one cycle")
    ax_survival.set_ylabel("survival from cycle start", color="tab:red")
    axes[1, 0].plot(result.phase_time / period, mean - 1.0,
                    label=r"mean extension $\bar\lambda-1$")
    axes[1, 0].plot(result.phase_time / period, np.sqrt(variance),
                    label=r"standard deviation $\sqrt{\mathrm{Var}}$")
    axes[1, 0].set(xlabel="cycle phase", ylabel="stretch",
                   title="Conditional intact observables")
    energy_axis = axes[1, 0].twinx()
    energy_axis.plot(result.phase_time / period, energy, color="tab:green",
                     label="mean shifted energy")
    energy_axis.set_ylabel("mean shifted energy", color="tab:green")
    axes[1, 0].legend(loc="upper left")
    axes[1, 1].plot(result.phase_time / period, hazard, color="tab:red")
    axes[1, 1].set(xlabel="cycle phase", ylabel="reduced hazard rate",
                   title="First-passage hazard (not an energy source term)")
    fig.suptitle("Principal periodic survival mode")
    fig.savefig(FIG / "phase_resolved_mode.png", dpi=180)
    plt.close(fig)

    # Tail and escape are deliberately evaluated in their different boundary
    # models.  They are plotted together only for comparison, never combined
    # into a constitutive source term.
    reflect_time = np.linspace(0.0, 8.0 * period, 801)
    reflect_force = 0.008 + 0.007 * np.sin(2.0 * np.pi * reflect_time / period)
    reflecting = solve(
        reflect_time, reflect_force,
        TransportConfig(inverse_temperature=2000.0, cells=140,
                        boundary="reflecting", initiation_definition="fixed_coordinate",
                        lambda_max=1.34),
        max_dt=0.03)
    last = reflect_time >= 7.0 * period
    reflect_phase = (reflect_time[last] - 7.0 * period) / period
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    axes[0, 0].plot(reflect_phase, reflect_force[last])
    axes[0, 0].set(xlabel="cycle phase", ylabel="reduced force",
                   title="Common loading protocol")
    axes[0, 1].plot(reflect_phase, reflecting.mean_shifted_energy[last],
                    color="tab:green")
    axes[0, 1].set(xlabel="cycle phase", ylabel="conditional mean shifted energy",
                   title="Reflecting state-function diagnostic")
    axes[1, 0].plot(reflect_phase, reflecting.tail_conditional[last],
                    color="tab:purple", label=r"reflecting tail above $\lambda_c$")
    axes[1, 0].set(xlabel="cycle phase", ylabel="conditional tail probability",
                   title="Instantaneous precursor, not cumulative damage")
    axes[1, 1].plot(result.phase_time / period, hazard, color="tab:red")
    axes[1, 1].set(xlabel="cycle phase", ylabel="absorbing hazard",
                   title="Irreversible first-passage observable")
    fig.suptitle("Energy, reflecting tail, and absorbing escape remain distinct")
    fig.savefig(FIG / "energy_tail_escape_distinction.png", dpi=180)
    plt.close(fig)

    # Frequency response: report both loss per cycle and loss per reduced time.
    periods = np.asarray([0.5, 1.0, 2.0, 4.0, 6.0, 12.0, 30.0, 80.0, 200.0])
    frequency_rows = []
    for sweep_period in periods:
        sweep_time, sweep_force = protocol(float(sweep_period))
        sweep_result = principal_survival_mode(
            sweep_time, sweep_force, absorbing_config(cells=100),
            max_dt=min(0.03, float(sweep_period) / 100.0), tolerance=1e-9)
        n50 = np.log(0.5) / np.log(sweep_result.multiplier)
        frequency_rows.append((
            sweep_period, 1.0 / sweep_period, sweep_result.multiplier,
            sweep_result.escape_per_cycle, sweep_result.integrated_hazard,
            sweep_result.mean_hazard_rate, n50))
    frequency_array = np.asarray(frequency_rows)
    adiabatic_c = absorbing_config(cells=70)
    frozen_rates = np.asarray([
        frozen_principal_escape_rate(value, adiabatic_c) for value in force])
    adiabatic_rate = float(np.trapezoid(
        frozen_rates, result.phase_time / period))

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    frequency = frequency_array[:, 1]
    axes[0].loglog(frequency, frequency_array[:, 3], "o-")
    axes[0].set(xlabel="reduced frequency", ylabel="escape per cycle",
                title=r"$1-r$")
    axes[1].semilogx(frequency, frequency_array[:, 5], "o-")
    axes[1].axhline(adiabatic_rate, color="black", lw=1, ls="--",
                    label="frozen-generator adiabatic limit")
    axes[1].set(xlabel="reduced frequency", ylabel="mean hazard per reduced time",
                title=r"$-\log(r)/T$")
    axes[1].legend(fontsize=8)
    axes[2].loglog(frequency, frequency_array[:, 6], "o-")
    axes[2].set(xlabel="reduced frequency", ylabel="median cycles in principal mode",
                title=r"$N_{50}=\log(0.5)/\log(r)$")
    fig.suptitle("Frequency changes time available for first passage")
    fig.savefig(FIG / "frequency_survival_spectrum.png", dpi=180)
    plt.close(fig)

    grid_cells = np.asarray([50, 70, 100, 140, 200])
    grid_multipliers = []
    for cells in grid_cells:
        grid_multipliers.append(principal_survival_mode(
            time, force, absorbing_config(cells=int(cells)),
            max_dt=0.015, tolerance=1e-9).multiplier)
    timestep_limits = np.asarray([0.08, 0.04, 0.02, 0.01])
    timestep_multipliers = []
    for dt_limit in timestep_limits:
        timestep_multipliers.append(principal_survival_mode(
            time, force, absorbing_config(cells=140),
            max_dt=float(dt_limit), tolerance=1e-9).multiplier)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    axes[0].plot(grid_cells, grid_multipliers, "o-")
    axes[0].set(xlabel="finite-volume cells", ylabel=r"multiplier $r$",
                title="Grid refinement")
    axes[1].semilogx(timestep_limits, timestep_multipliers, "o-")
    axes[1].invert_xaxis()
    axes[1].set(xlabel="maximum reduced timestep", ylabel=r"multiplier $r$",
                title="Timestep refinement")
    fig.suptitle("One-cycle spectrum numerical convergence")
    fig.savefig(FIG / "numerical_refinement.png", dpi=180)
    plt.close(fig)

    # Sensitivity only: beta and loading amplitude are not fitted here.
    inverse_temperatures = np.asarray([1000.0, 1500.0, 2000.0, 3000.0, 5000.0])
    amplitudes = np.asarray([0.003, 0.005, 0.007, 0.009])
    n50_map = np.empty((inverse_temperatures.size, amplitudes.size))
    sensitivity_rows = []
    for i, beta in enumerate(inverse_temperatures):
        for j, amplitude in enumerate(amplitudes):
            local_time, local_force = protocol(period, amplitude=float(amplitude))
            local_result = principal_survival_mode(
                local_time, local_force,
                absorbing_config(cells=90, inverse_temperature=float(beta)),
                max_dt=0.04, tolerance=1e-9)
            n50 = np.log(0.5) / np.log(local_result.multiplier)
            n50_map[i, j] = n50
            sensitivity_rows.append((beta, amplitude, local_result.multiplier,
                                     local_result.escape_per_cycle, n50))
    fig, ax = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    image = ax.imshow(np.log10(n50_map), origin="lower", aspect="auto",
                      extent=[amplitudes[0], amplitudes[-1],
                              inverse_temperatures[0], inverse_temperatures[-1]])
    fig.colorbar(image, ax=ax, label=r"$\log_{10} N_{50}$")
    ax.set(xlabel="reduced force amplitude", ylabel="inverse temperature $E_0/k_BT$",
           title="Uncalibrated sensitivity, not aluminum lifetime")
    fig.savefig(FIG / "parameter_sensitivity.png", dpi=180)
    plt.close(fig)

    write_rows(
        DATA / "phase_mode.csv",
        ["phase", "time_reduced", "force_reduced", "survival_from_cycle_start",
         "outflux", "hazard", "conditional_mean_lambda", "conditional_variance",
         "conditional_mean_shifted_energy"],
        zip(result.phase_time / period, result.phase_time, force,
            result.phase_survival, result.phase_outflux, hazard, mean, variance, energy))
    write_rows(
        DATA / "frequency_spectrum.csv",
        ["period_reduced", "frequency_reduced", "multiplier_r",
         "escape_per_cycle", "integrated_cycle_hazard", "mean_hazard_rate",
         "principal_mode_median_cycles"], frequency_rows)
    write_rows(
        DATA / "frozen_generator_rates.csv",
        ["phase", "force_reduced", "frozen_principal_escape_rate"],
        zip(result.phase_time / period, force, frozen_rates))
    write_rows(
        DATA / "parameter_sensitivity.csv",
        ["inverse_temperature", "force_amplitude", "multiplier_r",
         "escape_per_cycle", "principal_mode_median_cycles"], sensitivity_rows)
    write_rows(
        DATA / "direct_cycle_validation.csv",
        ["cycle", "generic_survival", "principal_mode_survival"],
        zip(np.arange(direct_survival.size), direct_survival, eigen_survival))
    write_rows(
        DATA / "biorthogonal_modes.csv",
        ["lambda", "right_conditional_density", "left_survival_weight"],
        zip(dense_x, dense_spectrum.right_density,
            dense_spectrum.left_survival_weight))
    write_rows(
        DATA / "numerical_refinement.csv",
        ["refinement_kind", "control_value", "multiplier_r"],
        [("grid_cells", int(value), multiplier)
         for value, multiplier in zip(grid_cells, grid_multipliers)]
        + [("max_reduced_timestep", value, multiplier)
           for value, multiplier in zip(timestep_limits, timestep_multipliers)])

    summary = {
        "status": "dimensionless consequence of the declared model; not calibrated aluminum life",
        "active_physics": "one-dimensional normal-tensile Smoluchowski first passage",
        "initiation_boundary": "tangent-instability stretch lambda_c",
        "new_fitted_parameters": 0,
        "m": c.m,
        "n": c.n,
        "inverse_temperature": c.inverse_temperature,
        "period_reduced": period,
        "mean_force_reduced": 0.008,
        "amplitude_reduced": 0.007,
        "principal_multiplier": result.multiplier,
        "escape_per_cycle": result.escape_per_cycle,
        "integrated_cycle_hazard": result.integrated_hazard,
        "mean_hazard_rate_reduced": result.mean_hazard_rate,
        "adiabatic_mean_hazard_rate_reduced_70_cell": adiabatic_rate,
        "principal_mode_median_cycles": float(np.log(0.5) / np.log(result.multiplier)),
        "principal_mode_mean_initiation_cycle": float(
            1.0 / result.escape_per_cycle),
        "principal_mode_initiation_cycle_variance": float(
            result.multiplier / result.escape_per_cycle ** 2),
        "operator_residual_l1": result.residual_l1,
        "power_iterations": result.iterations,
        "mass_balance_residual": float(abs(
            1.0 - result.multiplier
            - np.sum(result.phase_outflux[1:] * np.diff(result.phase_time)))),
        "generic_cycle_ratios": direct_ratios.tolist(),
        "dense_spectrum_cells": dense_c.cells,
        "second_to_principal_spectral_ratio": dense_spectrum.spectral_ratio,
        "generic_initial_asymptotic_prefactor": dense_prefactor,
        "peak_cycle_phases": {
            "force": float(result.phase_time[np.argmax(force)] / period),
            "conditional_mean_stretch": float(
                result.phase_time[np.argmax(mean)] / period),
            "conditional_variance": float(
                result.phase_time[np.argmax(variance)] / period),
            "conditional_mean_shifted_energy": float(
                result.phase_time[np.argmax(energy)] / period),
            "absorbing_hazard": float(
                result.phase_time[np.argmax(hazard)] / period),
        },
        "grid_refinement": dict(zip(map(str, grid_cells), map(float, grid_multipliers))),
        "timestep_refinement": dict(zip(map(str, timestep_limits),
                                         map(float, timestep_multipliers))),
        "frequency_spectrum": [dict(zip(
            ["period", "frequency", "r", "escape", "cycle_hazard",
             "mean_hazard_rate", "N50"], map(float, row))) for row in frequency_rows],
        "interpretation": (
            "The low dimensionless cycle count at this demonstration point is a consequence "
            "of the chosen reduced inverse temperature, force waveform, mobility time scale, "
            "and initiation boundary. It is not a fitted or predicted aluminum lifetime."
        ),
    }
    (DATA / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
