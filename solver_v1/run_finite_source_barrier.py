"""Execute declared finite-source geometry scenarios; no yield/kinetic fit."""
import argparse
from pathlib import Path
import time

import numpy as np

from .aluminum_calibration import EV_J
from .fcc111_geometry import fcc111_geometry_from_b
from .finite_source_barrier import pinned_line_barrier,normal_mode_second_variation
from .finite_source_reference import line_energy_coefficient,pinned_source_branch
from .nonlocal_interface_elasticity import cubic_elastic_tensor,rotate_elastic_tensor
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .vector_material_calibration import LENGTH_M


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve existing source scenarios')
    began=time.perf_counter();b=LENGTH_M
    rotation=fcc111_geometry_from_b(1).plane_basis_in_stacked_cubic_axes()
    tensor=rotate_elastic_tensor(cubic_elastic_tensor(114e9,62e9,32e9),rotation)
    lines={n:line_energy_coefficient(tensor,[b,0.,0.],samples=n) for n in (32,64)}
    cases=[];convergence=[];shapes=[]
    # These are explicit, hypothetical source/core geometries from the existing
    # long-wave research reference, NOT fitted values of an actual Al wire.
    for span in (.2e-6,1e-6,5e-6):
        for cutoff_b in (1.,2.,4.):
            outer_log=np.log(span/(cutoff_b*b))
            critical=pinned_source_branch(lines[64],np.pi/2,span_m=span,
                outer_log_ratio=outer_log)['critical_shear_outer_only_Pa']
            for fraction in (.2,.5,.8,.95,.99,.999):
                value=pinned_line_barrier(lines[64],span_m=span,outer_log_ratio=outer_log,shear_Pa=critical*fraction)
                cases.append(dict(span_um=span*1e6,core_radius_over_b=cutoff_b,
                    outer_radius_convention='R=span, fixed during branch variation',load_fraction=fraction,
                    applied_shear_MPa=critical*fraction/1e6,outer_critical_MPa=critical/1e6,
                    outer_barrier_eV=value['outer_only_barrier_J']/EV_J,
                    additional_swept_area_m2=value['additional_swept_area_m2'],
                    derived_stress_derivative_over_b3=value['outer_only_stress_derivative_m3']/b**3,
                    minimum_angle=value['metastable_endpoint_angle'],saddle_angle=value['saddle_endpoint_angle'],
                    barrier_quadrature_change_eV=value['quadrature_change'][0]/EV_J,
                    geometry_hypothetical=True,core_energy_included=False,actual_yield_prediction=False))
                if cutoff_b==2. and fraction in (.5,.95,.999):
                    coarse=pinned_line_barrier(lines[32],span_m=span,outer_log_ratio=outer_log,
                                               shear_Pa=critical*fraction,quadrature_points=32)
                    finite=[]
                    for delta in (critical*1e-5,critical*5e-6):
                        plus=pinned_line_barrier(lines[64],span_m=span,outer_log_ratio=outer_log,shear_Pa=critical*fraction+delta)
                        minus=pinned_line_barrier(lines[64],span_m=span,outer_log_ratio=outer_log,shear_Pa=critical*fraction-delta)
                        finite.append(-(plus['outer_only_barrier_J']-minus['outer_only_barrier_J'])/(2*delta))
                    convergence.append(dict(span_um=span*1e6,load_fraction=fraction,
                        angular_and_quadrature_barrier_relative_change=coarse['outer_only_barrier_J']/value['outer_only_barrier_J']-1,
                        energy_identity_error_eV=value['total_energy_identity_residual_J']/EV_J,
                        stress_difference_step_Pa=delta,
                        stress_derivative_coarse_relative_error=finite[0]/value['outer_only_stress_derivative_m3']-1,
                        stress_derivative_relative_error=finite[1]/value['outer_only_stress_derivative_m3']-1,
                        stress_derivative_richardson_relative_error=((4*finite[1]-finite[0])/3)/value['outer_only_stress_derivative_m3']-1))
                if span==1e-6 and cutoff_b==2. and fraction in (.5,.95):
                    p=critical*fraction*b;line=lines[64]
                    for branch,key in [('minimum','metastable_endpoint_angle'),('saddle','saddle_endpoint_angle')]:
                        endpoint=value[key];theta=np.linspace(-endpoint,endpoint,129)
                        gamma=line.evaluate(theta)*outer_log;gp=line.evaluate(theta,1)*outer_log
                        Q=gamma*np.sin(theta)+gp*np.cos(theta)
                        r=-gamma*np.cos(theta)+gp*np.sin(theta)
                        re=(-line.evaluate(endpoint)*np.cos(endpoint)+line.evaluate(endpoint,1)*np.sin(endpoint))*outer_log
                        for angle,x,y in zip(theta,Q/p,(re-r)/p):
                            shapes.append(dict(load_fraction=fraction,branch=branch,angle=angle,x_m=x,y_m=y,
                                first_normal_mode_curvature_J_m2=normal_mode_second_variation(
                                    endpoint_angle=endpoint,pressure_J_m2=p,mode=1)))
    write_csv(args.out/'hypothetical_finite_source_barriers.csv',cases)
    write_csv(args.out/'independent_refinement.csv',convergence)
    write_csv(args.out/'metastable_and_saddle_shapes.csv',shapes)
    save_json(args.out/'scope.json',dict(completed=True,cases=len(cases),
        cubic_GPa=[114.,62.,32.],bulk_condition='shared rounded 0K target and exact-bulk analytic calibration',
        underlying_path='same infinite-LJ/Bessel bulk elasticity -> half-space Schur kernel -> outer line energy',
        burgers_m=b,normal=[1,1,1],line_at_zero_bow='direct110 screw',
        source_lengths_measured=False,source_geometry_fitted=False,core_energy_included=False,
        nonlocal_finite_part_included=False,actual_dislocation_source_population=None,
        experimental_continuous_relaxation_matched=False,finite_source_atomistically_validated=False,
        thermally_activated_rate=None,production_probability=None,physical_PDE_time=False,
        elapsed_seconds=time.perf_counter()-began))
    print('finite-source outer-barrier geometry scenarios complete',flush=True)


if __name__=='__main__':main()
