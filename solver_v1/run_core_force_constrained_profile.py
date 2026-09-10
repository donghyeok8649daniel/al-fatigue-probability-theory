"""Fixed-shape, sign/sampled-spectrum-constrained source-core force profile.

All seven static anchors remain exact. Optimize force components ONLY, then
report the other interface observables as tradeoffs, not independent evidence
of an adopted Al material. No new energy term, yield or kinetic parameter.
"""
import argparse
import csv
import hashlib
import time
from pathlib import Path

import numpy as np

from .core_force_identifiability import constrained_force_profile
from .current_core_coefficient_basis import NAMES
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .run_current_material_core import ROOT,load_current_material
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',type=Path,default=ROOT/'results/current_material_core_v22/force_compatibility')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh constrained-force directory required')
    began=time.perf_counter()
    model,_,binding=load_current_material()
    data_path=args.audit/'force_coefficient_matrix.csv'
    with data_path.open(encoding='utf8',newline='') as stream:
        records=list(csv.DictReader(stream))
    selected=[r for r in records if r['configuration']=='R8_zero' and float(r['radius_over_L0'])<2.]
    design=np.array([[float(r[key]) for key in NAMES] for r in selected])
    target=np.array([float(r['source_force']) for r in selected])
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    observations=impose_tangents(problem.observations,problem.source_tangents)
    M=problem.matrix(tuple(model.screened_shape))
    exact=[i for i,o in enumerate(observations) if o.role=='exact']
    E=M[exact]; rhs=np.array([observations[i].target for i in exact])
    operators,tails=problem.operators(tuple(model.screened_shape))
    save_json(args.out/'definition.json',dict(candidate_sha256=binding['parameter_sha256'],
        force_matrix_sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(),
        coefficient_order=NAMES,shape=model.screened_shape,
        exact_observations=[observations[i].name for i in exact],
        fit_configuration='R8_zero',inner_radius_over_L0=2.,
        common_force_component_scale='1 eV/L0 (normalization, not reported uncertainty)',
        force_columns_are_energy_gradients=True,mechanical_force_is_negative_gradient=True,
        all_other_interface_observables_excluded_from_objective=True,
        nonnegative_indices=NONNEGATIVE,spectral_wavepoints=problem.parent_definition['wavepoints_cubic'],
        spectral_stretches=problem.parent_definition['stability_stretches'],
        spectral_radius=problem.parent_definition['radius_over_L0'],
        original_coefficient_shape_unchanged=True,production_changed=False))
    profiles=[]; summary=[]; errors=[]; parameters=[]
    for name,extra in [('sign_only',{}),('sign_and_sampled_spectrum',dict(operators=operators,tails=tails))]:
        fit=constrained_force_profile(design,target,E,rhs,nonnegative=NONNEGATIVE,**extra)
        c=fit['coefficients']; profiles.append(dict(name=name,**fit))
        static_loss=0.
        for o,pred,old in zip(observations,M@c,M@model.coefficients):
            residual=(pred-o.target)/o.scale
            if o.role=='fit': static_loss+=residual**2
            errors.append(dict(profile=name,observable=o.name,previous_role=o.role,
                target=o.target,scale=o.scale,units=o.units,baseline=old,prediction=pred,
                normalized_residual=residual,used_in_core_objective=False,
                retained_as_exact_anchor=o.role=='exact'))
        for parameter,old,new in zip(NAMES,model.coefficients,c):
            parameters.append(dict(profile=name,parameter=parameter,baseline=old,counterfactual=new,adopted=False))
        summary.append(dict(profile=name,force_rms_eV_L0=fit['force_rms'],
            exact_residual=fit['exact_residual'],kkt_residual=fit['kkt_residual'],
            nonnegative_minimum=float(min(c[list(NONNEGATIVE)])),positive_LJ=fit['strictly_positive_LJ'],
            other_interface_loss=static_loss,minimum_robust_margin=fit.get('minimum_robust_margin'),
            material_accepted=False))
        print(f'{name}: core RMS={fit["force_rms"]:.8g}, other interface loss={static_loss:.8g}',flush=True)
    write_csv(args.out/'profile_summary.csv',summary)
    write_csv(args.out/'static_tradeoff.csv',errors)
    write_csv(args.out/'counterfactual_parameters_NOT_ADOPTED.csv',parameters)
    validation=[]
    for fit in profiles:
        for case in sorted({r['configuration'] for r in records}):
            rows=[r for r in records if r['configuration']==case and float(r['radius_over_L0'])<2.]
            F=np.array([[float(r[key]) for key in NAMES] for r in rows])
            t=np.array([float(r['source_force']) for r in rows])
            residual=F@fit['coefficients']-t
            validation.append(dict(profile=fit['name'],configuration=case,
                force_rms_eV_L0=float(np.sqrt(np.mean(residual**2))),
                maximum_force_error_eV_L0=float(np.max(abs(residual))),used_in_loss=case=='R8_zero'))
    write_csv(args.out/'heldout_core_forces.csv',validation)
    save_json(args.out/'profiles.json',profiles)
    save_json(args.out/'completion.json',dict(completed=True,elapsed_seconds=time.perf_counter()-began,
        actual_constrained_optimization=True,fixed_shape_only=True,relaxed_core_recomputed=False,
        material_accepted=False,physical_yield_validated=False,production_changed=False,
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__':
    main()
