"""Execute exact infinite-row harmonic checks and all-profile PN comparisons.

Separate research outputs: the harmonic lattice and continuum core are NOT
silently glued together. Static iterations are never a physical-time path.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import time
import numpy as np
from scipy.optimize import brentq

from .discrete_fcc_screw import FCCScrewRowHessian, row_lj_slip_curvature, row_exponential_slip_derivatives
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor
from .nonlocal_registry_reference import ScrewRegistryFunctional
from .run_low_stress_cyclic_diagnostic import ROOT, build_surface, write_csv, save_json
from .run_nonlocal_interface_reference import load_prepared
from .static_bulk_stability import StaticBulkHessian


OUT = ROOT/"results/fcc111_active_interface/discrete_screw_v6"


def kernel_study():
    started = time.perf_counter()
    surface,units,meta = build_surface(tolerance=2e-11)
    rows = FCCScrewRowHessian.from_surface(surface)
    _,_,elastic_meta = load_prepared()
    tensor = cubic_elastic_tensor(*(elastic_meta["elastic_constants_GPa"][f"C{k}_GPa"]*1e9 for k in (11,12,44)))
    tensor = rotate_elastic_tensor(tensor,surface.interface.bulk.geometry.plane_basis_in_stacked_cubic_axes())
    conversion = units.length_scale_m**3/1.602176634e-19
    tensor *= conversion
    print("row sum",asdict(rows.diagnostics),flush=True)
    rigid = rows.rigid_step_stiffness()
    original = float(surface.packed(surface.h,0.)[5])
    local = rows.converged_jump_stiffness(0.)
    meta.update(row_sum=asdict(rows.diagnostics),rigid_step_hessian=rigid,
        original_interface_Wss=original,rigid_hessian_absolute_discrepancy=abs(original-rigid["total"]),
        discrete_relaxed_zero_wave_stiffness=local["stiffness"],
        rigid_to_relaxed_change_fraction=local["stiffness"]/original-1,
        symbol_units="eV/L0^2 per atom",jump_stiffness_units="eV/L0^2 per atomic interface cell",
        line_direction="e1 / direct_110; invariant along infinite atomic rows",
        row_spacing_L0=rows.d,normal_spacing_L0=rows.h,
        harmonic_lattice_only=True,nonlinear_atomistic_core_solved=False,
        physical_seconds_available=False,physical_hz_available=False,
        local_misfit_was_not_added_to_full_schur_kernel=True)
    # Independent direct 3D evaluator has a known algebraic LJ Hessian tail.
    # Its cutoff is validation ONLY; the new canonical row symbol is infinite.
    checks=[]
    qset=([0,.01,.005],[0,.4,.6],[0,1.2,.3],[0,0,1.5],[0,2.4,2.1])
    for radius in (8.,12.,16.,24.,32.):
        direct=StaticBulkHessian(surface.interface.bulk,cutoff=radius,
            D1=surface.vector_amplitude_ev,D2=surface.quadrupole_amplitude_ev,
            D3=surface.amplitude_ev,angular_decay=surface.angular.kappa)
        for q in qset:
            actual=direct.evaluate(q)
            exact=float(rows.symbol(q[1],q[2]*rows.h+q[1]*rows.d/3))
            value=float(actual["matrix"][0,0])
            checks.append(dict(direct_radius_L0=radius,q_y_L0=q[1],q_z_L0=q[2],
                poisson_symbol=exact,direct_symbol=value,absolute_error=abs(value-exact),
                relative_error=abs(value-exact)/abs(exact),
                scalar_F_second_term=float(actual["scalar_density"][0,0]),
                antiplane_cross_coupling=float(max(abs(actual["matrix"][0,1:])))))
        print(f"independent 3D radius {radius:g} completed",flush=True)
    write_csv(OUT/"direct_3d_validation.csv",checks)
    refinement=[]
    qy,th=np.meshgrid(np.linspace(-np.pi/rows.d,np.pi/rows.d,33),
                      np.linspace(-np.pi,np.pi,34,endpoint=False)+.003,indexing="ij")
    reference=np.concatenate([rows.symbol(qy.ravel()[i:i+256],th.ravel()[i:i+256])
                              for i in range(0,qy.size,256)])
    for ring in (2,3,4,5,6,8,10):
        other=FCCScrewRowHessian.from_surface(surface,fixed_ring=ring)
        values=np.concatenate([other.symbol(qy.ravel()[i:i+256],th.ravel()[i:i+256])
                               for i in range(0,qy.size,256)])
        refinement.append(dict(rings=ring,rows=other.diagnostics.rows_used,
            max_symbol_change=float(max(abs(values-reference))),
            rigid_step_stiffness=other.rigid_step_stiffness()["total"],
            omitted_validation_envelope=other.diagnostics.omitted_validation_ring_envelope))
    write_csv(OUT/"transverse_ring_refinement.csv",refinement)
    meta["minimum_sampled_antiplane_symbol"]=float(min(reference))
    meta["sampled_brillouin_points"]=int(len(reference))
    # Keep the two different finite-q coordinate constraints separate.
    wave_rows=[]; quadrature=[]; slope_rows=[]
    for constraint in ("logical_rows","spectral_same_y"):
        prediction=rows.long_wave_jump_slope(tensor,local["stiffness"],constraint=constraint)
        for q in (0.,.002,.005,.01,.02,.05,.1,.25,.5,1.,2.,np.pi/rows.d):
            solution=rows.converged_jump_stiffness(q,constraint=constraint)
            k=solution["stiffness"]
            wave_rows.append(dict(q_y_L0=q,constraint=constraint,stiffness_eV_L0sq=k,
                rigid_local_stiffness=original,relaxed_local_stiffness=local["stiffness"],
                continuum_additive_comparator=original+elastic_meta["mu_eV_L0cubed"]*rows.bulk.geometry.atomic_cell_area*abs(q)/2,
                observed_quadrature_change=solution["observed_quadrature_change"],
                layer_phase_points=solution["layer_phase_points"],
                derivative_slope=prediction["slope"],
                measured_secant_slope=(k-local["stiffness"])/abs(q) if q else None))
            for stage in solution["refinement_history"]:
                quadrature.append(dict(q_y_L0=q,constraint=constraint,
                    layer_phase_points=stage["layer_phase_points"],stiffness=stage["stiffness"]))
            print(f"{constraint} q={q:.6g}: K={k:.10g}, Ntheta={solution['layer_phase_points']}",flush=True)
        slope_rows.append(prediction)
    write_csv(OUT/"discrete_interface_kernel.csv",wave_rows)
    write_csv(OUT/"layer_phase_refinement.csv",quadrature)
    write_csv(OUT/"acoustic_pole_slopes.csv",slope_rows)
    response=[]
    for row_period in (4,8,16,32,128,512):
        q=2*np.pi/(row_period*rows.d)
        for constraint in ("logical_rows","spectral_same_y"):
            solution=rows.converged_jump_stiffness(q,constraint=constraint)
            for tau_mpa in (4.,15.,25.,50.):
                applied_force=tau_mpa*1e6*conversion*rows.bulk.geometry.atomic_cell_area
                slip=applied_force/solution["stiffness"]
                response.append(dict(atomic_row_period=row_period,wavelength_over_b=row_period*rows.d/rows.b,
                    shear_amplitude_MPa=tau_mpa,constraint=constraint,slip_amplitude_over_b=slip/rows.b,
                    continuum_local_amplitude_over_b=applied_force/original/rows.b,
                    zero_stress_harmonic_slip=0.,opening_dynamics_computed=False,
                    status="exact harmonic discrete static response; not nonlinear plasticity or fatigue"))
    write_csv(OUT/"low_mpa_harmonic_response.csv",response)
    meta["elapsed_seconds"]=time.perf_counter()-started
    save_json(OUT/"kernel_metadata.json",meta)


def profile_study():
    started=time.perf_counter()
    misfit,mu,meta=load_prepared()
    rows=[]; profiles=[]
    # Spatial-grid and domain refinement are distinct. The continuum model
    # may converge numerically while still failing an atomic-core validity test.
    settings=[(d,8.,.0625) for d in (32.,128.,512.,1024.)]
    settings += [(128.,8.,.25),(128.,8.,.125),(128.,8.,.03125),
                 (128.,4.,.0625),(128.,16.,.0625),(128.,32.,.0625)]
    for d,ratio,dx in settings:
        length=d*ratio
        n=2**int(np.ceil(np.log2(length/dx)))
        model=ScrewRegistryFunctional(misfit,elastic_factor=mu,domain_length=length,cells=n)
        trial=model.optimize_trial_width(d,width_bounds=(.02,4.))
        initial=model.dipole_trial(d,trial["width"])[0]
        solved=model.solve_fixed_registry_content(initial,force_tolerance=2e-10)
        if not solved["force_converged"]:
            raise ArithmeticError(f"full profile failed force tolerance at {d,ratio,dx}: {solved['projected_force_residual']}")
        # A midpoint crossing distance is a different observable from fixed
        # mean content; report both instead of silently calling them identical.
        s=solved["slip"]; crossings=[]
        for i in np.flatnonzero((s[:-1]-.5)*(s[1:]-.5)<0):
            crossings.append(model.x[i]+model.dx*(.5-s[i])/(s[i+1]-s[i]))
        crossing_distance=float(crossings[-1]-crossings[0]) if len(crossings)==2 else None
        conv=meta["elastic_conversion_Pa_to_eV_L0cubed"]*1e6
        row=dict(registry_content_over_b2=d,domain_over_content=ratio,cells=n,dx_over_b=model.dx,
            separation_from_half_registry_crossings_over_b=crossing_distance,
            holding_shear_MPa=solved["holding_traction"]/conv,
            projected_force_residual_MPa=solved["projected_force_residual"]/conv,
            mean_constraint_residual=solved["mean_constraint_residual"],
            energy_eV_per_L0_line=solved["energy"],trial_energy_eV_per_L0_line=trial["energy"],
            trial_balance_MPa=trial["balance_traction"]/conv,
            trial_full_force_residual_MPa=trial["full_euler_residual"]/conv,
            spectral_q90_b=solved["spectral_q90_times_period"],
            iterations=solved["iterations"],newton_polish=solved["newton_polish_iterations"],
            force_converged=solved["force_converged"],atomistic_core_validated=False,
            units_warning="energy per dislocation line, not a finite activation energy")
        rows.append(row)
        if d in (32.,128.) and ratio==8. and dx==.0625:
            for i in np.flatnonzero(abs(model.x+d/2)<3.):
                profiles.append(dict(registry_content=d,x_relative_left_content_endpoint=model.x[i]+d/2,
                    relaxed_slip=float(s[i]),arctangent_slip=float(initial[i]),
                    holding_shear_MPa=row["holding_shear_MPa"]))
        print(f"all-profile d={d:g}, L/d={ratio:g}, dx={model.dx:g}: {row['holding_shear_MPa']:.9g} MPa, residual={row['projected_force_residual_MPa']:.3g} MPa",flush=True)
    write_csv(OUT/"all_profile_refinement.csv",rows)
    write_csv(OUT/"all_profile_shapes.csv",profiles)
    stress_rows=[]
    for row in rows:
        if row["domain_over_content"]==8. and row["dx_over_b"]==.0625:
            for tau in (0.,4.,5.,10.,15.,25.,50.):
                stress_rows.append(dict(registry_content_over_b2=row["registry_content_over_b2"],
                    applied_shear_MPa=tau,holding_shear_MPa=row["holding_shear_MPa"],
                    registry_content_growth_drive_MPa=tau-row["holding_shear_MPa"],
                    actual_dynamic_trajectory=False,not_yield_stress=True))
    write_csv(OUT/"all_profile_stress_drives.csv",stress_rows)
    save_json(OUT/"profile_status.json",dict(elapsed_seconds=time.perf_counter()-started,
        all_profiles_force_converged=True,numerically_stationary_continuum_PN=True,
        constraint="fixed mean slip / total registry content; reaction is explicit holding traction",
        atomistic_core_validated=False,kinetic_time_available=False,finite_activation_energy_available=False,
        static_iterations_are_not_time=True,parameter_sha256=meta["parameter_sha256"]))


def matched_separation_domain_study():
    """Compare domains at the SAME actual core crossing separation.

    Fixed total mean content includes an O(L) intrawell far-field offset under
    the holding traction. Holding that content unchanged while increasing L
    changes the core separation. The protocol must correct the constraint,
    not retune the energy or subtract the physical intrawell displacement.
    """
    started=time.perf_counter()
    misfit,mu,meta=load_prepared()
    target=128.; output=[]
    for ratio,dx in ((4.,.0625),(8.,.0625),(16.,.0625),(32.,.0625),(64.,.0625),(8.,.03125)):
        length=ratio*target
        n=2**int(np.ceil(np.log2(length/dx)))
        model=ScrewRegistryFunctional(misfit,elastic_factor=mu,domain_length=length,cells=n)
        cache={}
        def solve(content):
            if content not in cache:
                result=model.solve_fixed_registry_content(model.dipole_trial(content,.278)[0],force_tolerance=2e-10)
                if not result["force_converged"]:
                    raise ArithmeticError("matched separation profile did not reach force tolerance")
                s=result["slip"]
                indices=np.flatnonzero((s[:-1]-.5)*(s[1:]-.5)<0)
                if len(indices)!=2:
                    raise ValueError("dipole branch lost; cannot manufacture a core separation")
                crossings=[model.x[i]+model.dx*(.5-s[i])/(s[i+1]-s[i]) for i in indices]
                cache[content]=(crossings[1]-crossings[0],result)
            return cache[content]
        far_field=mu/(2*length)/np.tan(np.pi*target/length)
        expected_offset=far_field/float(misfit.evaluate(0,2))*length
        lower,upper=target,target+max(2.,2*expected_offset)
        if (solve(lower)[0]-target)*(solve(upper)[0]-target)>=0:
            raise ValueError("no independently bracketed matched-core separation")
        content=brentq(lambda c:solve(c)[0]-target,lower,upper,xtol=2e-8)
        separation,result=solve(content)
        conv=meta["elastic_conversion_Pa_to_eV_L0cubed"]*1e6
        row=dict(target_crossing_separation_over_b=target,actual_crossing_separation_over_b=separation,
            registry_content_over_b2=content,domain_over_separation=ratio,dx_over_b=model.dx,cells=n,
            holding_shear_MPa=result["holding_traction"]/conv,
            projected_force_residual_MPa=result["projected_force_residual"]/conv,
            periodic_far_field_MPa=far_field/conv,
            isolated_far_field_MPa=mu/(2*np.pi*target)/conv,
            distinct_constrained_solves=len(cache),stationary=True,atomistic_certification=False)
        output.append(row)
        print(f"matched actual d={separation:.9g}, L/d={ratio:g}, dx={model.dx:g}: {row['holding_shear_MPa']:.9g} MPa; content={content:.9g}",flush=True)
    write_csv(OUT/"matched_separation_domain_refinement.csv",output)
    save_json(OUT/"matched_separation_status.json",dict(elapsed_seconds=time.perf_counter()-started,
        actual_solutions=sum(r["distinct_constrained_solves"] for r in output),
        target_separation=target,protocol="adjust constrained mean content, NOT material parameters, to match actual core crossings",
        physical_core_or_fatigue_calibrated=False))


def report():
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    def read(name):
        with (OUT/name).open(encoding="utf-8",newline="") as f:
            return list(csv.DictReader(f))
    k=read("discrete_interface_kernel.csv")
    p=read("all_profile_refinement.csv")
    shapes=read("all_profile_shapes.csv")
    fig,axes=plt.subplots(1,3,figsize=(14,4.1),constrained_layout=True)
    for label in ("logical_rows","spectral_same_y"):
        r=[v for v in k if v["constraint"]==label]
        axes[0].plot([float(v["q_y_L0"]) for v in r],[float(v["stiffness_eV_L0sq"]) for v in r],"o-",label=label)
    r=[v for v in k if v["constraint"]=="logical_rows"]
    axes[0].plot([float(v["q_y_L0"]) for v in r],[float(v["continuum_additive_comparator"]) for v in r],"--",label="old continuum + rigid local")
    axes[0].set(xlabel="transverse q L0",ylabel="stiffness [eV/L0^2 per interface cell]",title="Actual discrete constraint matters")
    axes[0].legend(fontsize=8)
    r=[v for v in shapes if float(v["registry_content"])==32.]
    axes[1].plot([float(v["x_relative_left_content_endpoint"]) for v in r],[float(v["arctangent_slip"]) for v in r],"--",label="one-width trial")
    axes[1].plot([float(v["x_relative_left_content_endpoint"]) for v in r],[float(v["relaxed_slip"]) for v in r],label="all-profile stationary")
    axes[1].set(xlabel="(x + registry content / 2) / b",ylabel="s / b",title="Continuum PN only, not atomic core")
    axes[1].legend(fontsize=8)
    r=sorted([v for v in p if float(v["registry_content_over_b2"])==128. and float(v["domain_over_content"])==8.],key=lambda v:float(v["dx_over_b"]))
    axes[2].plot([float(v["dx_over_b"]) for v in r],[float(v["holding_shear_MPa"]) for v in r],"o-",label="all-profile")
    axes[2].plot([float(v["dx_over_b"]) for v in r],[float(v["trial_balance_MPa"]) for v in r],"s--",label="restricted trial")
    axes[2].set(xlabel="numerical dx / b",ylabel="static holding shear [MPa]",title="Mesh convergence is not core validation")
    axes[2].legend(fontsize=8)
    fig.suptitle("Unchanged LJ/Bessel research candidate; no physical clock, no new fit",fontsize=12)
    fig.savefig(OUT/"discrete_screw_audit.png",dpi=160)
    svg=OUT/"discrete_screw_audit.svg"
    fig.savefig(svg)
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n",encoding="utf-8")
    plt.close(fig)
    metadata=json.loads((OUT/"kernel_metadata.json").read_text())
    direct=read("direct_3d_validation.csv")
    final=[v for v in direct if float(v["direct_radius_L0"])==32.]
    summary=dict(rigid_hessian_absolute_error=metadata["rigid_hessian_absolute_discrepancy"],
        direct_3d_radius32_max_absolute_error=max(float(v["absolute_error"]) for v in final),
        direct_3d_radius32_max_relative_error=max(float(v["relative_error"]) for v in final),
        all_profile_max_projected_force_residual_MPa=max(float(v["projected_force_residual_MPa"]) for v in p),
        exact_harmonic_infinite_row_reduction_verified=True,
        full_continuum_profile_force_balance_verified=True,
        full_nonlinear_atomistic_core_verified=False,
        unresolved_physical_items=["nonlinear, vector/normal-relaxed atomic core", "independent Al material fit",
                                  "finite loop / activation-energy normalization", "collective kinetic mobility"],
        physical_hz_available=False,production_solver_changed=False,ui_gate_passed=False)
    if (OUT/"matched_separation_domain_refinement.csv").exists():
        matched=read("matched_separation_domain_refinement.csv")
        m8=[v for v in matched if float(v["domain_over_separation"])==8.]
        final=next(v for v in matched if float(v["domain_over_separation"])==64.)
        summary.update(matched_separation_protocol_performed=True,
            matched_geometry_max_error_b=max(abs(float(v["actual_crossing_separation_over_b"])-128.) for v in matched),
            matched_geometry_last_mesh_traction_change_MPa=abs(float(m8[0]["holding_shear_MPa"])-float(m8[1]["holding_shear_MPa"])),
            matched_domain64_remaining_isolated_difference_fraction=abs(float(final["holding_shear_MPa"])/float(final["isolated_far_field_MPa"])-1))
    save_json(OUT/"physical_and_numerical_status.json",summary)
    print(json.dumps(summary,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("--kernel",action="store_true")
    parser.add_argument("--profiles",action="store_true")
    parser.add_argument("--matched-domains",action="store_true")
    parser.add_argument("--report",action="store_true")
    args=parser.parse_args()
    if not any(vars(args).values()):
        parser.error("choose --kernel, --profiles, --matched-domains and/or --report")
    if args.kernel:
        kernel_study()
    if args.profiles:
        profile_study()
    if args.matched_domains:
        matched_separation_domain_study()
    if args.report:
        report()


if __name__=="__main__":
    main()
