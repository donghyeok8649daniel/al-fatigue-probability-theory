"""Summarize/plot SAVED nonlocal static experiments; no solver or fit replay."""
from __future__ import annotations

import csv
import json
import numpy as np
from .run_nonlocal_interface_reference import OUT
from .run_low_stress_cyclic_diagnostic import save_json


def read(name):
    with (OUT/name).open(encoding="utf-8",newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    meta=json.loads((OUT/"material_metadata.json").read_text(encoding="utf-8"))
    acoustic=read("acoustic_limit_validation.csv")
    dipole=read("dipole_refinement.csv")
    stress=read("small_stress_response.csv")
    material=read("far_field_material_sensitivity.csv")
    fine=[r for r in acoustic if float(r["radius_L0"])==32 and float(r["q_L0"])==.005]
    mesh=sorted((r for r in dipole if float(r["d_over_b"])==128 and float(r["domain_over_d"])==8),
                key=lambda r:float(r["dx_over_b"]),reverse=True)
    domain=sorted((r for r in dipole if float(r["d_over_b"])==128 and float(r["dx_over_b"])==.0625),
                  key=lambda r:float(r["domain_over_d"]))
    shape=read("core_shape_sensitivity.csv")
    finite=read("finite_q_corrected_paths.csv")
    summary=dict(
        preparation_performed=True,static_scenarios_performed=True,
        static_low_stress_case_count=len(stress),variational_dipole_case_count=len(dipole),
        derived_mu_GPa=meta["screw_elastic_factor_Pa"]/1e9,
        max_corrected_acoustic_relative_error=max(float(r["relative_error"]) for r in fine),
        max_wrong_orientation_relative_error=max(float(r["wrong_abc_orientation_relative_error"]) for r in fine),
        min_sampled_finite_q_eigenvalue_eV_L0sq=min(float(r["min_eigenvalue_eV_L0sq"]) for r in finite),
        max_nonlinear_static_residual_MPa=max(float(r["solved_force_residual_MPa"]) for r in stress),
        max_static_unload_slip_over_b=max(float(r["static_unloaded_max_slip_over_b"]) for r in stress),
        last_dipole_mesh_change_MPa=abs(float(mesh[-1]["balance_shear_MPa"])-float(mesh[-2]["balance_shear_MPa"])),
        last_dipole_domain_change_MPa=abs(float(domain[-1]["balance_shear_MPa"])-float(domain[-2]["balance_shear_MPa"])),
        last_periodic_vs_isolated_leading_relative_difference=abs(float(domain[-1]["balance_shear_MPa"])
            /float(domain[-1]["far_field_isolated_MPa"])-1),
        core_shape_far_field_spread_MPa=float(np.ptp([float(r["balance_shear_MPa"]) for r in shape])),
        finest_trial_width_over_b=float(mesh[-1]["width_over_b"]),
        finest_trial_full_euler_residual_MPa=float(mesh[-1]["full_euler_residual_MPa"]),
        core_euler_equation_solved=False,atomic_core_validated=False,
        finite_loop_activation_energy_available=False,
        physical_Al_fatigue_validated=False,physical_mobility_available=False,
        physical_seconds_available=False,physical_Hz_available=False,
        production_PDE_changed=False,static_parameters_refitted=False,
        conclusion="long-wave elastic cost validated; narrow trial core and Al material calibration not certified",
    )
    save_json(OUT/"numerical_and_physical_status.json",summary)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout="constrained")
    a=axes[0,0]
    selected=[r for r in acoustic if r["radius_L0"]=="32.0" and r["direction"]=="[1.0, 0.0, 0.0]"]
    q=np.array([float(r["q_L0"]) for r in selected])
    a.loglog(q,[float(r["wrong_abc_orientation_relative_error"]) for r in selected],"o--",label="Wrong cubic frame")
    a.loglog(q,[float(r["relative_error"]) for r in selected],"o-",label="Correct +ABC frame")
    a.set(xlabel="|q| L0",ylabel="Relative acoustic-matrix error",title="Independent atomistic check")
    a.legend()
    a=axes[0,1]
    selected=[r for r in stress if float(r["shear_amplitude_MPa"])==15]
    a.semilogx([float(r["wavelength_over_b"]) for r in selected],
        [float(r["nonlocal_over_local_amplitude"]) for r in selected],"o-")
    a.set(xlabel="Spatial wavelength / b",ylabel="Nonlocal / local amplitude",title="Intrawell static response")
    a.text(.05,.1,"Short wavelengths are not\natomically certified",transform=a.transAxes,fontsize=9)
    a=axes[1,0]
    for name,label in (("unchanged_research_candidate","Research candidate, isolated far field"),
                       ("Mishin_0K_elastic_reference_only","Source elastic comparator only")):
        rows=[r for r in material if r["material"]==name]
        x=np.array([float(r["d_over_b"]) for r in rows])*meta["length_scale_m"]*1e9
        a.loglog(x,[float(r["isolated_far_field_balance_MPa"]) for r in rows],"o-",label=label)
    a.axhline(4,color="gray",ls="--",label="Specified 4 MPa shear")
    a.set(xlabel="Specified defect separation [nm]",ylabel="Attraction-balance shear [MPa]",
          title="Pre-existing defects, NOT nucleation/yield")
    a.legend(fontsize=8)
    a=axes[1,1]
    a.semilogx([float(r["domain_over_d"]) for r in domain],
        [float(r["balance_shear_MPa"]) for r in domain],"o-",label="Periodized trial, d=36.66 nm")
    a.axhline(float(domain[-1]["far_field_isolated_MPa"]),color="gray",ls="--",label="Isolated leading asymptote")
    a.set(xlabel="Periodic domain / separation",ylabel="Balance shear [MPa]",title="Boundary-size convergence")
    a.legend(fontsize=8)
    for a in axes.flat:
        a.grid(alpha=.25)
    fig.suptitle("Static nonlocal research audit — not fatigue or kinetic validation",fontsize=13)
    fig.savefig(OUT/"nonlocal_static_summary.png",dpi=160)
    svg=OUT/"nonlocal_static_summary.svg"
    fig.savefig(svg)
    # Matplotlib path attributes contain trailing spaces before newlines.
    # Normalize formatting only, so generated vectors pass git diff --check.
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n",encoding="utf-8")
    plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
