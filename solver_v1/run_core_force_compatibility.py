"""Actual same-state source-core force validation and fixed-shape lower bound.

Preserve all seven previously exact bulk/initial-interface observations. The
lower bound relaxes sign/stability inequalities and is NOT an adopted refit.
No experimental yield or desired dislocation barrier is a calibration input.
"""
import argparse
import hashlib
import time
from pathlib import Path

import numpy as np

from .core_force_identifiability import force_compatibility_lower_bound
from .current_core_coefficient_basis import current_core_coefficients,NAMES
from .current_material_rows import CurrentMaterialScrewCore
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .report_current_material_core import restore_case
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material,build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh force-audit directory required')
    began=time.perf_counter()
    material,_,binding=load_current_material()
    source,tensor,source_binding=load_source_material()
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    observations=impose_tangents(problem.observations,problem.source_tangents)
    M=problem.matrix(tuple(material.screened_shape))
    selected=[i for i,o in enumerate(observations) if o.role=='exact']
    E=M[selected]; rhs=np.array([observations[i].target for i in selected])
    if len(selected)!=7 or np.max(abs(E@material.coefficients-rhs))>2e-9:
        raise ArithmeticError('the unchanged seven-observation calibration did not replay')
    declared=[('R8_zero','escaped_plus_R8r5',0.),
              ('R10_zero','stable_plus_R10r5',0.),
              ('R8_perturbed_plus','escaped_plus_R8r5',.002),
              ('R8_perturbed_minus','escaped_plus_R8r5',-.002)]
    save_json(args.out/'definition.json',dict(completed=False,shape=material.screened_shape,
        baseline_coefficients=material.coefficients,coefficient_order=NAMES,
        candidate_sha256=binding['parameter_sha256'],source_sha256=source_binding['parameter_sha256'],
        configurations=declared,exact_observations=[observations[i].name for i in selected],
        fit_configuration='R8_zero',heldout=['R10_zero','R8_perturbed_plus','R8_perturbed_minus'],
        force_units='eV/L0 per straight atomic-row repeat', force_metric='equal-weight Euclidean force-component norm',
        source_boundary_frozen=True, same_physical_positions_for_both_energies=True,
        source_potential_is_only_a_target=True, physical_yield_used=False, new_energy_term=False))
    data={}; audit=[]; details=[]
    for name,case,amplitude in declared:
        path=ROOT/'results/current_material_core_v22/source_reference'/case
        source_core,field,meta,_,_=restore_case(path,source,tensor,source_binding,core_builder=build_source_core)
        field=field+amplitude*np.cos(np.arange(field.size)*.37).reshape(field.shape)
        actual_source=source_core.evaluate(field)
        # This is a FROZEN source atomic configuration, not a relaxed candidate
        # with source elasticity spliced into its constitutive response.
        candidate=CurrentMaterialScrewCore(material,source_core.far_field,
            free_radius=source_core.free_radius,ring=7,tolerance=2e-12)
        if not np.array_equal(candidate.indices[candidate.free_ids],source_core.indices[source_core.free_ids]):
            raise ArithmeticError('physical free-site geometry mismatch')
        columns=current_core_coefficients(candidate,field)
        actual=candidate.evaluate(field)
        replay=float(np.max(abs(columns['gradient']@material.coefficients-actual['gradient'])))
        if replay>3e-10:
            raise ArithmeticError('exact current core coefficient identity failed')
        radius=np.linalg.norm(candidate.xyz[candidate.free_ids,1:]-candidate.far_field.center,axis=1)
        data[name]=(columns['gradient'],actual_source['gradient'],radius)
        for window in (1.5,2.,3.):
            mask=radius<window
            difference=actual['gradient'][mask]-actual_source['gradient'][mask]
            audit.append(dict(configuration=name,inner_radius_over_L0=window,rows=int(mask.sum()),
                current_force_rms=float(np.sqrt(np.mean(actual['gradient'][mask]**2))),
                source_force_rms=float(np.sqrt(np.mean(actual_source['gradient'][mask]**2))),
                force_difference_rms=float(np.sqrt(np.mean(difference**2))),
                maximum_force_difference=float(np.max(abs(difference))),coefficient_identity_error=replay))
        for index,(logical,gradient,target,r) in enumerate(zip(candidate.indices[candidate.free_ids],
                columns['gradient'],actual_source['gradient'],radius)):
            if r>=3.: continue
            for component in range(3):
                details.append(dict(configuration=name,j=int(logical[0]),l=int(logical[1]),
                    component=component,radius_over_L0=r,source_force=target[component],
                    **{key:gradient[component,k] for k,key in enumerate(NAMES)}))
        print(f'force columns completed: {name}',flush=True)
    write_csv(args.out/'frozen_source_force_audit.csv',audit)
    write_csv(args.out/'force_coefficient_matrix.csv',details)
    table=[]; coefficients=[]; validation=[]; fits=[]
    F,target,radius=data['R8_zero']
    for window in (1.5,2.,3.):
        mask=radius<window
        fit=force_compatibility_lower_bound(F[mask].reshape(-1,10),target[mask].ravel(),
            E,rhs,material.coefficients,scales=abs(material.coefficients))
        c=fit['coefficients']; violations=[NAMES[i] for i in NONNEGATIVE if c[i]<0]
        fits.append(dict(inner_radius_over_L0=window,**fit,nonnegative_violations=violations,
            strictly_positive_LJ=bool(c[0]>0 and c[1]>0)))
        table.append(dict(inner_radius_over_L0=window,exact_rank=fit['exact_rank'],
            null_dimension=fit['null_dimension'],before_rms=fit['initial_force_rms'],
            lower_bound_rms=fit['unconstrained_lower_bound_force_rms'],
            lower_bound_maximum=fit['maximum_lower_bound_force_error'],
            exact_residual=fit['normalized_exact_residual'],
            normal_residual=fit['least_squares_normal_residual'],violations=';'.join(violations),
            material_accepted=False))
        for key,old,new in zip(NAMES,material.coefficients,c):
            coefficients.append(dict(inner_radius_over_L0=window,parameter=key,baseline=old,
                lower_bound_counterfactual=new,adopted=False))
        for name,(other,target_other,radius_other) in data.items():
            mask=radius_other<2.
            error=other[mask]@c-target_other[mask]
            validation.append(dict(inner_fit_radius_over_L0=window,configuration=name,
                used_in_lower_bound=name=='R8_zero',comparison_radius_over_L0=2.,
                force_rms=float(np.sqrt(np.mean(error**2))),maximum_force_error=float(np.max(abs(error)))))
    write_csv(args.out/'fixed_shape_lower_bounds.csv',table)
    write_csv(args.out/'counterfactual_parameters_NOT_ADOPTED.csv',coefficients)
    write_csv(args.out/'heldout_force_predictions.csv',validation)
    save_json(args.out/'lower_bounds.json',fits)
    save_json(args.out/'completion.json',dict(completed=True,actual_force_calculations=True,
        fixed_shape_only=True,inequalities_relaxed=True,global_family_impossibility_claimed=False,
        material_accepted=False,relaxed_candidate_recomputed=False,production_changed=False,
        physical_yield_calibrated=False,elapsed_seconds=time.perf_counter()-began,
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__':
    main()
