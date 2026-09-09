"""Empirical refinement tables/plots from completed low-stress experiments.

All floors are observed differences, NOT statistical confidence intervals or
proof of physical calibration. Reflecting-box outputs contain no crack sink.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.special import logsumexp

from .low_stress_cyclic_diagnostic import KB_EV_K
from .run_low_stress_cyclic_diagnostic import OUT, write_csv, save_json
from .summarize_low_stress_cycles import rows, load_runs, snapshot


def tag(case, a=61, s=49, dt=128, wells=3, upper=2., period=40.):
    return f"{case}_p{period:g}_a{a}_s{s}_w{wells}_dt{dt}_upper{upper:g}"


METRICS = ("end_outside_mass", "excess_outside_mass", "hold_excess_mean_n",
           "last_cycle_net_transfer", "last_cycle_registry_amplitude",
           "last_cycle_registry_phase_deg", "last_cycle_external_work_ev",
           "cumulative_gross_grid_traffic")


def refinement_tables():
    summaries = {r["run"]: r for r in rows(OUT/"cyclic_summary.csv")}
    comparisons = []; missing = []
    for case in ("axial30", "axial50"):
        pairs = [
            ("a_only", tag(case, 31, 25), tag(case, 61, 25)),
            ("s_only", tag(case, 61, 25), tag(case)),
            ("combined_spatial", tag(case), tag(case, 91, 73)),
            ("dt_64_128", tag(case, dt=64), tag(case)),
            ("dt_128_256", tag(case), tag(case, dt=256)),
            ("dt_256_512", tag(case, dt=256), tag(case, dt=512)),
        ]
        if case == "axial30":
            pairs += [
                ("fine_grid_dt", tag(case, 91, 73), tag(case, 91, 73, dt=256)),
                ("normal_domain", tag(case), tag(case, 85, 49, upper=2.5)),
                ("domain_3_5", tag(case), tag(case, wells=5)),
                ("domain_5_7", tag(case, wells=5), tag(case, wells=7)),
            ]
        for name, low, high in pairs:
            if low not in summaries or high not in summaries:
                missing.append(dict(comparison=name, low=low, high=high)); continue
            x, y = summaries[low], summaries[high]
            for metric in METRICS:
                error = abs(y[metric]-x[metric])
                comparisons.append(dict(case=case, comparison=name, metric=metric,
                    coarse=low, fine=high, coarse_value=x[metric], fine_value=y[metric],
                    absolute_change=error,
                    relative_change=error/abs(y[metric]) if y[metric] else None,
                    status="observed numerical change; not a physical certificate"))
    write_csv(OUT/"refinement_comparison.csv", comparisons)
    # Conservative observed envelope: retain independently changed coordinates,
    # last temporal refinement, and both finite-domain changes. No hardcoded
    # magnitude cutoff or division by a presumed convergence order.
    sources = {"a_only", "s_only", "combined_spatial", "dt_256_512",
               "normal_domain", "domain_3_5", "domain_5_7"}
    floors = []
    for case in ("axial30", "axial50"):
        reference = summaries.get(tag(case, 91, 73, dt=256),
                                  summaries.get(tag(case, 91, 73), summaries[tag(case)]))
        for metric in METRICS:
            selected = [r for r in comparisons if r["case"] == case and r["metric"] == metric
                        and r["comparison"] in sources]
            if not selected: continue
            maximum = max(selected, key=lambda r: r["absolute_change"])
            # Moment/occupation roundoff floor uses actual mass/flux residuals.
            residual = max(reference["mass_residual"], reference["well_balance"],
                           reference["registry_balance"])
            if metric not in ("end_outside_mass", "excess_outside_mass", "hold_excess_mean_n",
                              "last_cycle_net_transfer"):
                residual = 0.  # Do not mix probability units into phase or energy.
            floor = max(maximum["absolute_change"], residual)
            signal = abs(reference[metric])
            floors.append(dict(case=case, metric=metric, reference_run=reference["run"],
                signal=reference[metric], observed_resolution_envelope=floor,
                dominant_comparison=maximum["comparison"], signal_over_envelope=signal/floor if floor else None,
                empirical_numeric_only=True, physical_calibration_accepted=False,
                normal_and_registry_domain_refined=case == "axial30",
                gross_is_not_a_convergent_hop_count=metric == "cumulative_gross_grid_traffic"))
    write_csv(OUT/"observed_resolution_envelopes.csv", floors)
    save_json(OUT/"refinement_execution.json", dict(missing_comparisons=missing,
        completed_comparison_metric_rows=len(comparisons), confidence_interval=False,
        certifies_material_fatigue=False))
    return summaries, comparisons, floors


def static_and_conditional_tables():
    static = rows(OUT/"quasistatic_mpa_cycles.csv"); report = []
    for path in ("direct_110", "shockley_112"):
        for amplitude in (10., 20., 30., 50.):
            r = [r for r in static if r["path"] == path and r["nominal_axial_amplitude_mpa"] == amplitude]
            if not r: continue
            q = np.array([[x["a"], x["s"]] for x in r])
            loads = np.array([[x["normal_traction_GPa"], x["shear_traction_GPa"]] for x in r])
            repeat_error = float(np.max(abs(q[:33]-q[32:])))
            report.append(dict(path=path, nominal_axial_amplitude_mpa=amplitude,
                all_steps_stable=all(x["status"] == "stable_local_equilibrium" for x in r),
                normal_min=min(x["normal_strain"] for x in r), normal_max=max(x["normal_strain"] for x in r),
                min_slip_over_period=min(x["slip_over_period"] for x in r),
                max_slip_over_period=max(x["slip_over_period"] for x in r),
                final_registry_shift=r[-1]["slip_over_period"],
                max_force_residual=max(x["force_residual_ev_coordinate"] for x in r),
                min_hessian_eigenvalue=min(x["min_hessian_eigenvalue"] for x in r),
                repeated_static_cycle_coordinate_error=repeat_error,
                kinetic_or_fatigue_validation=False))
    write_csv(OUT/"quasistatic_summary.csv", report)
    profiles = []; barriers = []; kT = KB_EV_K*293.15
    for file in sorted(OUT.glob("grid_*.npz")):
        meta = json.loads(file.with_suffix(".json").read_text())
        with np.load(file) as data:
            energy = data["energy"]; a = data["a"]; s = data["s"]
            central = abs(s) < .5*float(data["period"])
            e = energy[:, central]; selected_s = s[central]
            logz = logsumexp(-e/kT, axis=0)+np.log(float(data["da"]))
            free = -kT*logz
            conditional = np.exp(-e/kT-logsumexp(-e/kT, axis=0))
            mean_a = np.sum(conditional*a[:, None], axis=0)
            base = int(np.argmin(abs(selected_s)))
            difference = free-free[base]
            for j, value in enumerate(selected_s):
                profiles.append(dict(grid=file.name, s_reduced=float(value),
                    fixed_domain_F_minus_F0_ev=float(difference[j]),
                    conditional_a_over_h=float(mean_a[j]/float(data["h"])),
                    sampled_min_a_over_h=float(a[np.argmin(e[:, j])]/float(data["h"]))))
            barriers.append(dict(grid=file.name, temperature_K=293.15,
                sampled_F_barrier_ev=float(np.max(difference)),
                barrier_over_kT=float(np.max(difference)/kT),
                barrier_s_reduced=float(selected_s[np.argmax(difference)]),
                thermal_normal_strain_at_s0=float(mean_a[base]/float(data["h"])-1),
                interpretation="fixed reflecting-domain one-cell PMF; sampled maximum, not a coupled saddle",
                physical_activation_normalization_validated=False))
    write_csv(OUT/"conditional_profiles.csv", profiles)
    write_csv(OUT/"conditional_barrier_summary.csv", barriers)


def integrator_checks(runs):
    mapped = {r["path"].name: r for r in runs}; comparisons = []
    for steps in (64, 128):
        for case in ("zero", "axial30"):
            name = tag(case, 31, 25, dt=steps, period=.05)+"_c2_h0"
            if name not in mapped or name+"_explicit" not in mapped: continue
            implicit = mapped[name]; explicit = mapped[name+"_explicit"]
            with np.load(implicit["path"]/"snapshots.npz") as x, np.load(explicit["path"]/"snapshots.npz") as y:
                l1 = float(np.sum(abs(snapshot(x, -1)-snapshot(y, -1))))
            comparisons.append(dict(case=case, period_model=.05, driven_cycles=2, steps=steps,
                final_density_L1_explicit_implicit=l1,
                final_registry_difference=abs(implicit["cycles"][-1]["mean_well_index"]-
                                               explicit["cycles"][-1]["mean_well_index"]),
                status="same analytic surface and SG operator, short costly-explicit control"))
    write_csv(OUT/"explicit_implicit_comparison.csv", comparisons)


def plots(runs, summaries):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mapped = {r["path"].name: r for r in runs}
    base = mapped[tag("zero")]; loaded = mapped[tag("axial30")]
    other = mapped[tag("axial50")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for r, label in ((base, "Zero load"), (loaded, "Axial amplitude 30 MPa"),
                     (other, "Axial amplitude 50 MPa")):
        c = r["cycles"]; t = [x["cycle"]*x["period_model"] for x in c]
        axes[0, 0].plot(t, [x["outside_central_mass"] for x in c], label=label)
        control = np.array([x["outside_central_mass"] for x in base["cycles"]])
        axes[0, 1].plot(t, np.array([x["outside_central_mass"] for x in c])-control, label=label)
        axes[1, 0].plot(t, [x["mean_well_index"] for x in c], label=label)
    for ax in axes.flat[:3]:
        ax.axvline(640, color="k", ls="--", lw=.7, label="Unload at 640 model time")
        ax.set_xlabel("Model time (NOT seconds)"); ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[0, 0].set_ylabel("Outside-central-well probability")
    axes[0, 1].set_ylabel("Driven minus matched-zero outside probability")
    axes[1, 0].set_ylabel("Mean registry index <n>")
    axes[0, 0].legend(fontsize=8)
    for r, label in ((loaded, "30 MPa axial"), (other, "50 MPa axial")):
        h = [x for x in r["history"] if x["cycle"] == 16 or x["time_model"] == 600.]
        # Stress at the actual state timestamp, not the midpoint forcing sample.
        sigma = float(r["audit"]["load"]["normal_amplitude_mpa"])*2
        stress = sigma*np.sin(2*np.pi*np.array([x["time_model"] for x in h])/40)
        axes[1, 1].plot([x["registry_shear"] for x in h], stress, label=label)
    axes[1, 1].set_xlabel("Local registry shear <s>/h (unscaled)")
    axes[1, 1].set_ylabel("Example nominal axial stress [MPa]")
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("One-cell reflecting-domain hypothesis test; no crack sink / no physical Hz\n"
                 "Local registry diagnostics: NOT crack probability or macroscopic plastic strain", fontsize=12)
    fig.savefig(OUT/"low_mpa_cycles.png", dpi=150); plt.close(fig)


def main():
    summaries, comparisons, floors = refinement_tables()
    static_and_conditional_tables()
    runs = load_runs(); integrator_checks(runs); plots(runs, summaries)
    manifest = []
    for r in runs:
        a = r["audit"]
        cells = a["grid_s"]//a["wells"]
        upper = a["domain_a_edges"][1]/np.sqrt(2/3)
        command = ("python -m solver_v1.run_low_stress_cyclic_diagnostic --phase cycles"
            f" --case {a['case']} --n-a {a['grid_a']} --cells-per-well {cells}"
            f" --wells {a['wells']} --upper {upper:g} --period {a['period_model']:g}"
            f" --steps {a['steps_per_cycle']} --cycles {len(r['active'])}"
            f" --hold {len(r['cycles'])-len(r['active'])}"
            f" --integrator {a.get('integrator', 'implicit')}")
        if not a["initial_central_only"]: command += " --global-gibbs"
        manifest.append(dict(stored_run=r["path"].name, rerun_command=command,
            grid_prepare_command=("python -m solver_v1.run_low_stress_cyclic_diagnostic --phase grid"
                f" --n-a {a['grid_a']} --cells-per-well {cells} --upper {upper:g}"),
            elapsed_seconds=a["elapsed_seconds"], grid_sha256=a["grid_sha256"],
            parameter_sha256=a["parameter_sha256"], actual_evolution_executed=True))
    save_json(OUT/"execution_manifest.json", manifest)
    print(f"Read {len(runs)} completed runs, {len(comparisons)} refinement metric rows.")
    for r in floors:
        if r["metric"] in ("excess_outside_mass", "hold_excess_mean_n"):
            print(r)


if __name__ == "__main__": main()
