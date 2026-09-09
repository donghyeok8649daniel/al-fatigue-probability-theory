"""Actually solve leading-log pinned-line references; no empirical yield fit."""
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from .fcc111_geometry import fcc111_geometry_from_b
from .full_fcc_calibration_audit import cubic_constants_gpa
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor
from .finite_source_reference import line_energy_coefficient, pinned_source_branch, solve_pinned_graph
from .run_low_stress_cyclic_diagnostic import FIT, write_csv
from .run_vector_registry_audit import save_json
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import VectorCoefficientBasis, observation_matrix, LENGTH_M

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fcc111_active_interface/yield_bridge_v11/finite_source'


def main():
    started=time.perf_counter(); _,obs,_=source_and_targets()
    prior=ROOT/'results/fcc111_active_interface/material_strength_v10/opening_exchange/calibration.json'
    metric=ROOT/'results/fcc111_active_interface/yield_bridge_v11/material_metric_refined/calibration.json'
    fits={'historical_candidate':json.loads(FIT.read_bytes())['angular_monotone_opening'],
          'v10_exchange_rejected':json.loads(prior.read_text())['best']}
    fits.update({f'v11_{k}_not_adopted':v for k,v in json.loads(metric.read_text())['best'].items()})
    materials={}
    for name,p in fits.items():
        matrix=observation_matrix(VectorCoefficientBasis(p['scalar_decay'],p['angular_decay']),obs)
        materials[name]=cubic_constants_gpa(matrix[:5]@p['coefficients'])
    materials['source_elastic_comparator_0K']={f'C{k}_GPa':v for k,v in zip((11,12,44),(114.,62.,32.))}
    metadata=[]; angular=[]; refinement=[]; curves=[]; critical=[]; cases=[]
    for name,c in materials.items():
        tensor=rotate_elastic_tensor(cubic_elastic_tensor(*(c[f'C{k}_GPa']*1e9 for k in (11,12,44))),
            fcc111_geometry_from_b(1.).plane_basis_in_stacked_cubic_axes())
        series={n:line_energy_coefficient(tensor,[LENGTH_M,0.,0.],samples=n) for n in (32,64,128)}
        theta=np.linspace(-np.pi/2,np.pi/2,257)
        fine=series[128]
        for n,line in series.items():
            refinement.append(dict(model=name,kind='angular',samples=n,segments=None,
                coefficient_error_J_m=float(np.max(abs(line.evaluate(theta)-fine.evaluate(theta)))),
                stiffness_error_J_m=float(np.max(abs(line.stiffness(theta)-fine.stiffness(theta)))),
                force_residual=None,shape_error_over_span=None,swept_area_over_span_sq=None))
        for t in theta:
            angular.append(dict(model=name,angle_rad=t,coefficient_J_m=fine.evaluate(t),
                                stiffness_per_log_J_m=fine.stiffness(t)))
        row=dict(model=name,**c,screw_coefficient_J_m=float(fine.evaluate(0)),
                 edge_coefficient_J_m=float(fine.evaluate(np.pi/2)),
                 min_stiffness_per_log_J_m=float(np.min(fine.stiffness(theta))),
                 material_calibrated_for_yield=False)
        metadata.append(row)
        # EXPLICIT hypotheses, not measured source sizes or a fitted core radius.
        # R=L is a leading-log outer-cutoff convention. No finite core/loop claim.
        for span_um in (.5,1.,2.,5.,10.):
            L=span_um*1e-6
            for core_over_b in (.5,1.,2.):
                logarithm=np.log(L/(core_over_b*LENGTH_M))
                try:
                    end=pinned_source_branch(fine,np.pi/2,span_m=L,outer_log_ratio=logarithm)
                except ValueError as error:
                    critical.append(dict(model=name,span_um=span_um,core_radius_over_b=core_over_b,
                        outer_radius_over_span=1.,critical_outer_only_MPa=None,status=str(error))); continue
                stress=end['critical_shear_outer_only_Pa']
                critical.append(dict(model=name,span_um=span_um,core_radius_over_b=core_over_b,
                    outer_radius_over_span=1.,critical_outer_only_MPa=stress/1e6,
                    status='hypothetical geometry; leading-log outer elasticity ONLY; NOT experimental yield'))
                if core_over_b!=1.: continue
                for applied in (2.,4.,10.,25.,50.):
                    if applied*1e6>=stress:
                        cases.append(dict(model=name,span_um=span_um,shear_MPa=applied,
                            endpoint_angle_rad=None,maximum_bow_nm=None,swept_area_m2=None,
                            state='above outer-only local-line fold; no atomistic operation/yield claim'))
                        continue
                    angle=brentq(lambda t:pinned_source_branch(fine,t,span_m=L,outer_log_ratio=logarithm,
                        points=9)['applied_shear_Pa']-applied*1e6,1e-7,np.pi/2)
                    numerical=solve_pinned_graph(fine,span_m=L,outer_log_ratio=logarithm,
                        shear_Pa=applied*1e6,segments=128)
                    if not numerical['converged']: raise RuntimeError('independent pinned graph did not converge')
                    cases.append(dict(model=name,span_um=span_um,shear_MPa=applied,
                        endpoint_angle_rad=angle,maximum_bow_nm=max(numerical['y_m'])*1e9,
                        swept_area_m2=numerical['swept_area_m2'],
                        state='reversible subcritical outer-elastic bow; NOT 0.2% yield'))
        if row['min_stiffness_per_log_J_m']<=0: continue
        L=1e-6; logarithm=np.log(L/LENGTH_M)
        for endpoint in (.35,.8,1.2,1.4):
            exact=pinned_source_branch(fine,endpoint,span_m=L,outer_log_ratio=logarithm,points=513)
            for x,y in zip(exact['x_m'],exact['y_m']):
                curves.append(dict(model=name,endpoint_angle=endpoint,x_over_span=x/L,y_over_span=y/L,
                                   applied_outer_only_MPa=exact['applied_shear_Pa']/1e6))
            for n in (32,64,128,256):
                solve=solve_pinned_graph(fine,span_m=L,outer_log_ratio=logarithm,
                    shear_Pa=exact['applied_shear_Pa'],segments=n)
                if not solve['converged']: raise RuntimeError('refinement graph failed')
                refinement.append(dict(model=name,kind=f'graph_angle_{endpoint}',samples=128,segments=n,
                    coefficient_error_J_m=None,stiffness_error_J_m=None,
                    force_residual=solve['dimensionless_force_residual'],
                    shape_error_over_span=abs(max(solve['y_m'])-max(exact['y_m']))/L,
                    swept_area_over_span_sq=solve['swept_area_m2']/L**2))
        print(f'{name}: line coefficient and independent source branches complete',flush=True)
    for filename,rows in [('elastic_metadata.csv',metadata),('angular_coefficients.csv',angular),
                          ('refinement.csv',refinement),('source_shapes.csv',curves),
                          ('hypothetical_source_thresholds.csv',critical),('stress_scenarios.csv',cases)]:
        write_csv(OUT/filename,rows)
    save_json(OUT/'status.json',dict(completed=True,elapsed_seconds=time.perf_counter()-started,
        energy_origin='bulk Hessian of SAME LJ/Bessel candidate, then derived continuum Schur coefficient',
        reference_only=True,material_promoted=False,empirical_plastic_law_added=False,
        source_lengths_measured=False,core_radius_measured=False,core_energy_included=False,
        finite_part_included=False,atomistic_finite_source_validated=False,
        actual_yield_stress_MPa=None,physical_kinetics=False,production_PDE_changed=False))


if __name__=='__main__': main()
