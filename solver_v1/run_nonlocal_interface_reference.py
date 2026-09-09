"""Execute the nonlocal static audit using an unchanged saved LJ/Bessel surface.

No fitting, PDE registration, physical clock, or manufactured activation area.
The arctangent defect family is explicitly VARIATIONAL, not an exact core or
finite-temperature nucleation calculation. Output energies retain /line units.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import brentq

from .aluminum_calibration import EV_J
from .even_moment_calibration import quadrupole_bulk_curvatures
from .full_fcc_calibration_audit import cubic_constants_gpa, relaxed_bulk_observables, independent_bulk_targets
from .nonlocal_interface_elasticity import (
    cubic_elastic_tensor, rotate_elastic_tensor, halfspace_impedance, screw_elastic_factor,
)
from .nonlocal_registry_reference import PeriodicMisfit, ScrewRegistryFunctional
from .run_low_stress_cyclic_diagnostic import build_surface, write_csv, save_json, ROOT
from .static_bulk_stability import StaticBulkHessian


OUT = ROOT/"results/fcc111_active_interface/nonlocal_v5"


def sample_material():
    started = time.perf_counter()
    surface, units, metadata = build_surface(tolerance=2e-11)
    bulk = surface.interface.bulk
    even = quadrupole_bulk_curvatures(surface.quadrupole)
    elastic_rows = []
    for step in (8e-4, 4e-4):
        obs = relaxed_bulk_observables(bulk, strain_step=step)
        obs[2:] += surface.quadrupole_amplitude_ev*even
        constants = cubic_constants_gpa(obs, volume_scale=metadata["lattice_stretch"]**3)
        elastic_rows.append(dict(strain_difference_step=step, **constants))
    constants = elastic_rows[-1]
    c = cubic_elastic_tensor(*(constants[f"C{k}_GPa"]*1e9 for k in (11,12,44)))
    geometry = bulk.geometry
    naive_c = rotate_elastic_tensor(c, [geometry.e1, geometry.e2, geometry.e3])
    c = rotate_elastic_tensor(c, geometry.plane_basis_in_stacked_cubic_axes())
    mu = screw_elastic_factor(c)
    write_csv(OUT/"bulk_elasticity.csv", elastic_rows)
    impedance_rows = []
    for angle in np.linspace(0, np.pi, 13, endpoint=False):
        result = halfspace_impedance(c, [np.cos(angle), np.sin(angle)])
        k = result.jump_per_wave_number/1e9
        impedance_rows.append(dict(angle_radians=angle, min_eigenvalue_GPa=float(np.linalg.eigvalsh(k)[0]),
            normal_GPa=float(k[2,2].real), slip_e1_GPa=float(k[0,0].real),
            normal_slip_GPa=float(k[2,0].real), imag_residual_GPa=float(np.max(abs(k.imag))),
            hermitian_residual=result.hermitian_residual,riccati_residual=result.riccati_residual))
    write_csv(OUT/"halfspace_kernel.csv", impedance_rows)
    print(f"bulk elastic C={constants}; screw factor {mu/1e9:.9f} GPa", flush=True)
    samples = []
    for index in range(128):
        s = index/128*surface.period
        values = surface.packed(surface.h, s)
        samples.append([s,values[0],values[2],values[5]])
        if index % 16 == 0:
            print(f"unchanged analytic registry sample {index}/128; {time.perf_counter()-started:.1f}s",flush=True)
    samples = np.asarray(samples)
    area_reduced = units.atomic_cell_area_m2/units.length_scale_m**2
    fine = PeriodicMisfit.from_samples(samples[:,1]/area_reduced,period=surface.period)
    coarse = PeriodicMisfit.from_samples(samples[::2,1]/area_reduced,period=surface.period)
    checks = []
    for fraction in (.0037,.0231,.077,.149,.249,.333,.411,.497,.623,.783,.951):
        s = fraction*surface.period; value = surface.packed(surface.h,s)
        for order,slot in enumerate((0,2,5)):
            exact = value[slot]/area_reduced
            checks.append(dict(s_over_b=fraction,derivative_order=order,analytic=exact,
                fourier_128=float(fine.evaluate(s,order)),fourier_64=float(coarse.evaluate(s,order)),
                error_128=float(fine.evaluate(s,order)-exact),error_64=float(coarse.evaluate(s,order)-exact)))
    write_csv(OUT/"registry_series_validation.csv", checks)
    write_csv(OUT/"analytic_registry_samples.csv", [dict(s=r[0],W_eV_cell=r[1],W_s=r[2],W_ss=r[3]) for r in samples])
    wave_rows=[]
    conversion=units.length_scale_m**3/EV_J
    volume_reduced=bulk.geometry.atomic_cell_area*bulk.a0
    finite_q_rows=[]
    for radius in (12.,20.,32.):
        direct=StaticBulkHessian(bulk,cutoff=radius,D3=surface.amplitude_ev,
            D1=surface.vector_amplitude_ev,D2=surface.quadrupole_amplitude_ev,
            angular_decay=surface.angular.kappa)
        for magnitude in (.04,.02,.01,.005):
            for direction in ([1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[.6,0.,.8]):
                q=magnitude*np.asarray(direction)
                exact=volume_reduced*np.einsum("ijkl,j,l->ik",c*conversion,q,q)
                actual=direct.evaluate(q)["matrix"]
                wave_rows.append(dict(radius_L0=radius,q_L0=magnitude,direction=str(direction),
                    acoustic_limit_norm=float(np.linalg.norm(exact)),
                    discrepancy=float(np.linalg.norm(actual-exact)),relative_error=float(np.linalg.norm(actual-exact)/np.linalg.norm(exact)),
                    wrong_abc_orientation_relative_error=float(np.linalg.norm(actual-volume_reduced*np.einsum(
                        "ijkl,j,l->ik",naive_c*conversion,q,q))/np.linalg.norm(exact))))
        for label,end in (("GX",[0,1,0]),("GL",[.5,.5,.5]),("GK",[.75,.75,0])):
            for fraction in np.linspace(.05,1.,20):
                q=direct.crystallographic_wavevector(fraction*np.asarray(end))
                eigen=direct.evaluate(q)["eigenvalues"]
                finite_q_rows.append(dict(radius_L0=radius,path=label,path_fraction=fraction,
                    q_plane_frame=str(q.tolist()),min_eigenvalue_eV_L0sq=float(eigen[0]),
                    status="corrected +ABC cubic labels; sampled static stability only"))
        print(f"independent finite-q validation radius {radius:g} L0 complete",flush=True)
    write_csv(OUT/"acoustic_limit_validation.csv",wave_rows)
    write_csv(OUT/"finite_q_corrected_paths.csv",finite_q_rows)
    # Analytic curvature brackets locate all sampled maxima; neither a chosen
    # empirical yield stress nor a coupled a-s spinodal is inferred here.
    locations=[]; line=np.linspace(0,surface.period,257)
    for left,right in zip(line[:-1],line[1:]):
        if fine.evaluate(left,2)*fine.evaluate(right,2)<0:
            locations.append(brentq(lambda s:float(fine.evaluate(s,2)),left,right,xtol=1e-13))
    ideal=max(float(fine.evaluate(s,1)) for s in locations)/conversion/1e6
    metadata.update(elastic_constants_GPa={k:v for k,v in constants.items() if k.startswith("C")},
        screw_elastic_factor_Pa=mu,elastic_conversion_Pa_to_eV_L0cubed=conversion,
        mu_eV_L0cubed=mu*conversion,atomic_area_L0sq=area_reduced,
        registry_curvature_eV_L0fourth=float(fine.evaluate(0,2)),
        fixed_gap_ideal_shear_traction_MPa=ideal,
        max_series_value_error=max(abs(r["error_128"]) for r in checks if r["derivative_order"]==0),
        max_series_force_error=max(abs(r["error_128"]) for r in checks if r["derivative_order"]==1),
        max_series_curvature_error=max(abs(r["error_128"]) for r in checks if r["derivative_order"]==2),
        elapsed_seconds=time.perf_counter()-started,
        surface_role="unchanged angular research candidate; fixed normal gap, direct_110 scalar path",
        continuum_status="long-wavelength harmonic limit, not exact atomic core",
        cubic_frame="+ABC stack: plane basis in actual cubic axes is (-e1,-e2,e3)",
        atom_positions_or_static_parameters_changed=False,
        probability_generator=False,physical_kinetics_available=False)
    np.savez_compressed(OUT/"prepared_misfit.npz",coefficients=fine.coefficients,period=fine.period,
        mu=mu*conversion,tensor_pa=c)
    metadata["prepared_misfit_sha256"]=hashlib.sha256((OUT/"prepared_misfit.npz").read_bytes()).hexdigest()
    save_json(OUT/"material_metadata.json",metadata)
    print(json.dumps(metadata,indent=2),flush=True)


def load_prepared():
    meta=json.loads((OUT/"material_metadata.json").read_text(encoding="utf-8"))
    if hashlib.sha256((ROOT/meta["parameter_source"]).read_bytes()).hexdigest()!=meta["parameter_sha256"]:
        raise ValueError("saved research parameters changed: re-run preparation")
    if hashlib.sha256((OUT/"prepared_misfit.npz").read_bytes()).hexdigest()!=meta["prepared_misfit_sha256"]:
        raise ValueError("prepared misfit checksum mismatch")
    with np.load(OUT/"prepared_misfit.npz") as data:
        misfit=PeriodicMisfit(float(data["period"]),data["coefficients"].copy())
        mu=float(data["mu"])
    return misfit,mu,meta


def static_scenarios():
    started=time.perf_counter(); p,mu,meta=load_prepared()
    b=p.period; conversion=meta["elastic_conversion_Pa_to_eV_L0cubed"]
    wavelength_rows=[]
    for wavelength_b in (4,8,16,32,128,512):
        model=ScrewRegistryFunctional(p,elastic_factor=mu,domain_length=wavelength_b*b,cells=256)
        mode=np.cos(2*np.pi*model.x/model.length)
        for shear_mpa in (4.,15.,25.):
            tau=shear_mpa*1e6*conversion
            s=model.linear_response(tau*mode)
            ratio=tau/(float(p.evaluate(0,2))+mu/2*2*np.pi/model.length)/(tau/float(p.evaluate(0,2)))
            residual=max(abs(model.evaluate(s,tau*mode)["force"]))/conversion/1e6
            solved=model.solve_intrawell_equilibrium(tau*mode)
            unloaded=model.solve_intrawell_equilibrium(0.,initial=solved["slip"])
            wavelength_rows.append(dict(wavelength_over_b=wavelength_b,shear_amplitude_MPa=shear_mpa,
                max_slip_over_b=float(max(abs(s))/b),nonlocal_over_local_amplitude=ratio,
                nonlinear_force_residual_MPa=float(residual),q_b=2*np.pi/wavelength_b,
                solved_max_slip_over_b=float(max(abs(solved["slip"]))/b),
                solved_force_residual_MPa=solved["force_residual"]/conversion/1e6,
                static_unloaded_max_slip_over_b=float(max(abs(unloaded["slip"]))/b),
                newton_iterations=solved["iterations"],
                status="static within-well linearization; short wavelengths not atomically certified"))
    write_csv(OUT/"small_stress_response.csv",wavelength_rows)
    rows=[]; stress_rows=[]; profile_rows=[]
    settings=[(d,8.,.0625) for d in (32.,128.,512.,1024.)]
    settings += [(128.,r,dx) for r,dx in ((4.,.0625),(16.,.0625),(32.,.0625),(64.,.0625),
                                         (8.,.25),(8.,.125),(8.,.03125))]
    for d_b,domain_ratio,dx_b in settings:
        length=domain_ratio*d_b*b
        cells=2**int(np.ceil(np.log2(length/(dx_b*b))))
        model=ScrewRegistryFunctional(p,elastic_factor=mu,domain_length=length,cells=cells)
        fit=model.optimize_trial_width(d_b*b,width_bounds=(.02*b,4*b))
        row=dict(d_over_b=d_b,domain_over_d=domain_ratio,cells=cells,dx_over_b=model.dx/b,
            d_nm=d_b*b*meta["length_scale_m"]*1e9,width_over_b=fit["width"]/b,
            balance_shear_MPa=fit["balance_traction"]/conversion/1e6,
            far_field_periodic_MPa=mu*b/(2*length)*1/np.tan(np.pi*d_b*b/length)/conversion/1e6,
            far_field_isolated_MPa=mu/(2*np.pi*d_b)/conversion/1e6,
            energy_eV_per_L0_line=fit["energy"],elastic_eV_per_L0_line=fit["elastic"],
            misfit_eV_per_L0_line=fit["misfit"],
            width_stationarity=fit["width_force"],full_euler_residual_MPa=fit["full_euler_residual"]/conversion/1e6,
            spectral_q90_b=fit["spectral_q90_times_period"],elastic_series_error=fit["elastic_series_error"],
            optimizer_success=fit["optimizer_success"],width_interior=fit["width_interior"],
            activation_energy_eV=None,physical_core_certified=False,
            status="pre-existing screw dipole; restricted trial-family force, NOT nucleation barrier")
        rows.append(row)
        for shear in (0.,4.,5.,10.,15.,25.,50.):
            drive=shear-row["balance_shear_MPa"]
            stress_rows.append(dict(d_over_b=d_b,domain_over_d=domain_ratio,dx_over_b=model.dx/b,
                applied_shear_MPa=shear,balance_shear_MPa=row["balance_shear_MPa"],
                expansion_drive_MPa=drive,
                energy_derivative_with_separation_eV_L0sq=(row["balance_shear_MPa"]-shear)*1e6*conversion*b,
                tendency="expansion within trial family" if drive>0 else "contraction within trial family",
                dynamic_motion_computed=False))
        if d_b==128. and domain_ratio==8. and dx_b==.0625:
            slip=model.dipole_trial(d_b*b,fit["width"])[0]
            for i in range(0,cells,8):
                profile_rows.append(dict(x_over_b=model.x[i]/b,slip_over_b=slip[i]/b))
        write_csv(OUT/"dipole_refinement.csv",rows)
        print(f"dipole d={row['d_nm']:.4f}nm L/d={domain_ratio:g} dx/b={model.dx/b:g}: "
              f"width/b={row['width_over_b']:.5g}, balance={row['balance_shear_MPa']:.6g} MPa",flush=True)
    write_csv(OUT/"dipole_stress_scenarios.csv",stress_rows)
    write_csv(OUT/"dipole_trial_profile.csv",profile_rows)
    # Independent shape-width sensitivity: separates robust far-field drive
    # from a NOT atomically validated sub-lattice core. No width is adopted.
    sensitivity=[]
    model=ScrewRegistryFunctional(p,elastic_factor=mu,domain_length=8*1024*b,cells=131072)
    for width_b in (.15,.28,.5,1.):
        slip,ds,_=model.dipole_trial(1024*b,width_b*b)
        values=model.evaluate(slip)
        shear=model.dx*np.dot(values["force"],ds)/b/conversion/1e6
        sensitivity.append(dict(d_over_b=1024,width_over_b=width_b,balance_shear_MPa=shear,
            status="specified shape sensitivity, not fitted or physically selected core"))
    write_csv(OUT/"core_shape_sensitivity.csv",sensitivity)
    # Source elastic constants are a separate validation comparator, NOT a
    # replacement of the candidate's elastic tensor to make a desired result.
    source=cubic_constants_gpa(independent_bulk_targets()[0])
    source_mu=np.sqrt(source["C44_GPa"]*(source["C11_GPa"]-source["C12_GPa"])/2)*1e9
    material_rows=[]
    for name,elastic in (("unchanged_research_candidate",mu/conversion),("Mishin_0K_elastic_reference_only",source_mu)):
        for d_b in (128.,512.,1024.):
            material_rows.append(dict(material=name,screw_factor_GPa=elastic/1e9,d_over_b=d_b,
                isolated_far_field_balance_MPa=elastic/(2*np.pi*d_b)/1e6,
                source="repository independent_bulk_targets; Mishin et al. 1999" if name.startswith("Mishin") else meta["parameter_source"],
                not_a_fitted_yield_stress=True,finite_core_certified=False))
    write_csv(OUT/"far_field_material_sensitivity.csv",material_rows)
    save_json(OUT/"execution_status.json",dict(elapsed_seconds=time.perf_counter()-started,
        no_refit=True,no_physical_kinetics=True,no_probability_computed=True,
        exact_core_solved=False,finite_loop_nucleation_solved=False,production_unchanged=True,
        physical_fatigue_validation=False,
        command="python -m solver_v1.run_nonlocal_interface_reference --scenarios"))
    print(f"static scenarios finished in {time.perf_counter()-started:.2f}s",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare",action="store_true")
    parser.add_argument("--scenarios",action="store_true")
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.prepare:
        sample_material()
    if args.scenarios:
        static_scenarios()
    if not (args.prepare or args.scenarios):
        parser.error("select --prepare and/or --scenarios")


if __name__=="__main__":
    main()
