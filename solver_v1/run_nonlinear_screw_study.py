"""Deterministic static defect study; no time, activation rate, or refitting.

python -m solver_v1.run_nonlinear_screw_study --sizes 24 32 48
The transverse grid IS the atomic lattice. Size variation changes periodic
image spacing, not atomic resolution. The supplied dipole is an initial defect,
not a predicted nucleation event or a calibrated Al dislocation population.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np

from .nonlinear_fcc_screw import NonlinearScrewRows
from .run_low_stress_cyclic_diagnostic import ROOT, build_surface, save_json, write_csv


OUT = ROOT/"results/fcc111_active_interface/nonlinear_screw_v7"
EV_J = 1.602176634e-19


def topology_record(cell, field):
    winding = cell.winding(field)
    locations = np.argwhere(winding["charge"] != 0)
    positions = []
    for j, l in locations:
        jc = j+.5-(cell.shape[0]-1)/2
        lc = l+.5-(cell.shape[1]-1)/2
        positions.append(dict(charge=int(winding["charge"][j,l]),
                              y_over_L0=cell.rows.d*(jc+lc/3), z_over_L0=cell.rows.h*lc))
    separation = None
    if len(locations) == 2:
        delta = locations[1]-locations[0]
        # Closest image in the non-orthogonal periodic cross-section.
        separation = min(np.hypot(cell.rows.d*(delta[0]+i*cell.shape[0]+
                         (delta[1]+k*cell.shape[1])/3),
                         cell.rows.h*(delta[1]+k*cell.shape[1]))
                         for i in (-1,0,1) for k in (-1,0,1))
    return dict(positive_cores=int((winding["charge"] > 0).sum()),
                negative_cores=int((winding["charge"] < 0).sum()),
                core_separation_over_L0=None if separation is None else float(separation),
                core_positions=json.dumps(positions, separators=(",", ":")),
                winding_integer_residual=winding["integer_residual"],
                minimum_half_period_margin=winding["minimum_half_period_margin"])


def validation(surface, rows):
    """Independent retained convolution and actual tail changes, not a claim of
    mathematically rigorous infinite-tail certification."""
    shape = (8,9); j,l = np.indices(shape)
    u = .21*np.cos(.71*j+.43*l); gamma=.017
    cell=rows.periodic_cell(shape)
    center=cell.evaluate(u,gamma,direction=np.sin(.3*j+.7*l),gamma_direction=.011,fields=True)
    convolution_error=float(np.max(abs(center['channels']-cell.direct_channels(u,gamma))))
    v=np.sin(.3*j+.7*l); step=2e-6
    plus=cell.evaluate(u+step*v,gamma+step*.011)
    minus=cell.evaluate(u-step*v,gamma-step*.011)
    errors=dict(convolution_channels_max_abs=convolution_error,
        gradient_directional_abs=abs((plus['energy']-minus['energy'])/(2*step)-
                      np.sum(center['gradient']*v)-.011*center['affine_derivative']),
        hessian_vector_max_abs=float(np.max(abs((plus['gradient']-minus['gradient'])/(2*step)-
                                               center['hessian_vector']))),
        affine_hessian_abs=abs((plus['affine_derivative']-minus['affine_derivative'])/(2*step)-
                              center['affine_hessian_vector']))
    results=[]
    for ring in (2,4,6,8):
        refined=NonlinearScrewRows(surface,fixed_ring=ring,tolerance=2e-13)
        out=refined.periodic_cell(shape).evaluate(u,gamma)
        results.append(dict(rings=ring,reciprocal_modes=len(refined.g),
                            energy_ev=out['energy'],
                            energy_difference_reference_ev=abs(out['energy']-center['energy']),
                            force_difference_reference=float(np.max(abs(out['gradient']-center['gradient'])))))
    write_csv(OUT/'tail_refinement.csv',results)
    return errors


def run(sizes, *, check_stability=True):
    started=time.perf_counter()
    surface,units,meta=build_surface(tolerance=2e-11)
    rows=NonlinearScrewRows(surface,tolerance=2e-12)
    mpa_to_internal=1e6*units.length_scale_m**3/EV_J
    histories=[]; stability=[]; snapshots=[]; states={}
    protocol=(0.,4.,15.,25.,50.,25.,0.,-25.,-50.,-25.,0.)
    metadata=dict(**meta,energy_normalization="eV per infinite-row repeat b*L0",
        line_energy_normalization="E_eV*eV_J/(b*L0) [J/m]; NOT finite activation energy",
        state_space="static scalar anti-plane; fixed transverse coordinates; infinite straight e1 lines",
        loading="imposed affine shear stress; minimize E - tau*V_cell*gamma",
        volume_definition="N*b*d*h*L0^3 for the actual periodic atomic cell; no fitted volume",
        kinetic_time=False,production_registered=False,normal_relaxation=False,
        vector_core=False,finite_loop_or_source=False,initial_defect="declared opposite screw pair, seed separation 7.3b",
        row_diagnostics=rows.diagnostics,mpa_to_ev_per_L0_cubed=mpa_to_internal,
        total_affine_shear_is_not_new_plastic_strain=True)
    # Preserve the parameter provenance; this is not a new calibration.
    metadata['energy_unit']='eV per b*L0 infinite-line repeat'
    metadata['verification']=validation(surface,rows)
    save_json(OUT/'metadata.json',metadata)
    for size in sizes:
        cell=rows.periodic_cell((size,size))
        cases=('perfect','dipole') if size==sizes[0] else ('dipole',)
        for case in cases:
            field=np.zeros(cell.shape) if case=='perfect' else cell.dipole_seed(7.3)
            gamma=0.; baseline=None
            for index,tau in enumerate(protocol):
                tick=time.perf_counter()
                result=cell.relax(field,gamma=gamma,shear_stress=tau*mpa_to_internal)
                field=result['displacement'];gamma=result['affine_shear']
                shear=cell.registry_shear(field,gamma)
                if baseline is None:
                    baseline=shear.copy()
                record=dict(size_j=size,size_l=size,case=case,step=index,
                    applied_shear_mpa=tau,internal_shear_mpa=result['internal_shear_stress']/mpa_to_internal,
                    energy_ev_per_line_repeat=result['energy'],
                    line_energy_j_per_m=result['energy']*EV_J/(rows.b*units.length_scale_m),
                    **shear,registry_shear_change_from_initial=shear['registry_shear']-baseline['registry_shear'],
                    intrabond_shear_change_from_initial=shear['intrabond_shear']-baseline['intrabond_shear'],
                    affine_shear_change_from_initial=gamma-baseline['affine_shear'],
                    maximum_force_residual=result['maximum_force_residual'],
                    scaled_stress_residual=result['scaled_stress_residual'],
                    force_converged=result['force_converged'],optimizer_success=result['optimizer_success'],
                    optimizer_message=result['optimizer_message'],iterations=result['iterations'],
                    newton_polish_steps=result['newton_polish_steps'],
                    **topology_record(cell,field),elapsed_seconds=time.perf_counter()-tick)
                histories.append(record)
                write_csv(OUT/'stress_history.csv',histories)
                print(f"{case} {size}x{size} step {index} tau={tau:g} MPa: gamma={gamma:.9g}, "
                      f"d_registry={record['registry_shear_change_from_initial']:.5g}, "
                      f"cores={record['positive_cores']}/{record['negative_cores']}, "
                      f"force={record['maximum_force_residual']:.3g}, {record['elapsed_seconds']:.2f}s",flush=True)
                if not result['force_converged']:
                    raise ArithmeticError('force balance failed; partial output preserved, not classified')
                if index in (0,4,8,10):
                    states[(size,case,index)]=(cell,field.copy(),gamma)
                    for j,l in np.ndindex(cell.shape):
                        snapshots.append(dict(size=size,case=case,step=index,j=j,l=l,u_over_L0=field[j,l],gamma=gamma))
                    write_csv(OUT/'static_states.csv.gz',snapshots)
                if check_stability and case=='dipole' and index in (0,4,8) and size==sizes[0]:
                    print('  evaluating nongauge nonlinear static Hessian ...',flush=True)
                    stable=cell.minimum_curvature(field,gamma,stress_control=True,tolerance=2e-7)
                    stability.append(dict(size=size,case=case,step=index,applied_shear_mpa=tau,**stable))
                    write_csv(OUT/'stability.csv',stability)
                    print(f"  lambda_min={stable['minimum_eigenvalue']:.8g}; residual={stable['eigen_residual']:.3g}",flush=True)
    # Inspect actual relaxed core states with independently extended row neighborhoods.
    tail=[]
    for ring in (4,6,8):
        refined=NonlinearScrewRows(surface,fixed_ring=ring,tolerance=2e-13)
        for (size,case,index),(cell,field,gamma) in states.items():
            if case!='dipole' or index not in (0,4):
                continue
            out=refined.periodic_cell(cell.shape).evaluate(field,gamma)
            original=cell.evaluate(field,gamma)
            tail.append(dict(size=size,step=index,rings=ring,
                energy_change_ev=out['energy']-original['energy'],
                maximum_force_change=float(np.max(abs(out['gradient']-original['gradient']))),
                stress_change_mpa=(out['affine_derivative']-original['affine_derivative'])/
                                  (cell.size*rows.b*rows.d*rows.h*mpa_to_internal)))
    write_csv(OUT/'relaxed_core_tail_refinement.csv',tail)
    summary=dict(elapsed_seconds=time.perf_counter()-started,states=len(histories),
        all_force_balanced=all(r['force_converged'] for r in histories),
        max_force_residual=max(r['maximum_force_residual'] for r in histories),
        max_shear_identity_error=max(abs(r['decomposition_residual']) for r in histories),
        max_absolute_new_registry_shear=max(abs(r['registry_shear_change_from_initial']) for r in histories),
        domain_sizes=list(sizes),stress_protocol_mpa=list(protocol),
        experiment_type='athermal quasistatic local minimization; NOT cyclic kinetics or a hold in physical time',
        physical_al_strength_validated=False,physical_time_available=False,production_gate=False)
    save_json(OUT/'execution_summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes',type=int,nargs='+',default=[24,32,48])
    parser.add_argument('--skip-stability',action='store_true')
    args=parser.parse_args()
    run(args.sizes,check_stability=not args.skip_stability)
