"""Predeclared scan of the EXISTING rank1 density screening, not a new law.

Hold the five wide-search shapes fixed; independently profile coefficients
at exponent -1,-.5,0,.5,1. Same sources, seven exact anchors, signs, sampled
bulk spectral inequalities and fractional-error objective. No yield target.
The zero-exponent calculation must replay the previously executed profile.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .coordination_screening import CoordinationScreenedInterface
from .frozen_core_force_target import FrozenCoreForceTarget
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .report_current_material_core import restore_case
from .run_current_material_core import ROOT, load_current_material
from .run_source_core_reference import load_source_material, build_source_core
from .tail_constrained_material import spectral_profile
from .vector_material_calibration import MaterialObservation
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('fresh screening-profile directory required')
    raw=(args.probe/'completion.json').read_bytes();prior=json.loads(raw)
    declared=json.loads((args.probe/'definition.json').read_bytes())
    if not prior['completed'] or not prior['best']['strictly_positive_LJ']:
        raise ValueError('executed positive-LJ source shape profile required')
    model,_,binding=load_current_material()
    source,tensor,source_binding=load_source_material()
    if (declared['candidate_sha256']!=binding['parameter_sha256'] or
            declared['source_sha256']!=source_binding['parameter_sha256']):
        raise ValueError('same material and source bindings required')
    shape=np.asarray(prior['best']['shape'],float)
    if shape.shape!=(6,) or shape[-1]!=0:
        raise ValueError('the declared zero-screening reference is required')
    source_path=ROOT/declared['source_state']
    core,field,*_=restore_case(source_path,source,tensor,source_binding,core_builder=build_source_core)
    target=FrozenCoreForceTarget(core,field,radius=declared['target_radius_over_L0'])
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    original=impose_tangents(problem.observations,problem.source_tangents)
    selected=[i for i,o in enumerate(original) if o.role=='fit']
    static_scale=np.sqrt(declared['static_baseline_squared_loss'])
    core_scale=np.sqrt(declared['core_baseline_squared_loss'])
    observation=[replace(o,scale=o.scale*static_scale) if o.role=='fit' else o for o in original]
    observation += [MaterialObservation(f'frozen_core_gradient_{i}',float(g),core_scale,'eV/L0','fit')
        for i,g in enumerate(target.coefficient_matrix(model)['target'])]
    # Order starts from the mandatory replay, then symmetric fixed probes.
    exponents=(0.,-.5,.5,-1.,1.)
    save_json(args.out/'definition.json',dict(completed=False,exponents=exponents,
        frozen_five_shapes=shape[:5],source_probe_sha256=hashlib.sha256(raw).hexdigest(),
        candidate_sha256=binding['parameter_sha256'],source_sha256=source_binding['parameter_sha256'],
        objective=declared['objective'],existing_family_only=True,new_energy_term=False,
        seven_exact_anchors_unchanged=True,physical_yield_used=False,production_changed=False))
    began=time.perf_counter();profiles=[];summary=[];residuals=[]
    for exponent in exponents:
        shape=shape.copy();shape[-1]=exponent
        tick=time.perf_counter()
        candidate=CoordinationScreenedInterface(shape,np.ones(10),law='power',tolerance=2e-12)
        M=problem.matrix(tuple(shape));F=target.coefficient_matrix(candidate)
        operators,tails=problem.operators(tuple(shape))
        fit=spectral_profile(np.vstack([M,F['design']]),observation,operators,tails,nonnegative=NONNEGATIVE)
        if exponent==0 and abs(fit['squared_loss']-prior['best']['squared_loss'])>2e-8:
            raise ArithmeticError('same-data zero-screening coefficient profile failed replay')
        c=fit['coefficients'];r=(M@c-np.array([o.target for o in original]))/np.array([o.scale for o in original])
        error=F['design']@c-F['target']
        row=dict(exponent=exponent,joint_loss=fit['squared_loss'],
            other_interface_loss=float(r[selected]@r[selected]),
            core_force_RMS_eV_L0=float(np.sqrt(np.mean(error**2))),
            exact_residual=fit['exact_residual'],kkt_residual=fit['kkt_residual'],
            strictly_positive_LJ=fit['strictly_positive_LJ'],
            minimum_sampled_margin=fit['minimum_robust_margin'],seconds=time.perf_counter()-tick,
            fixed_five_shapes=True,relaxed_core_recomputed=False,material_accepted=False)
        for o,pred,ri in zip(original,M@c,r):
            residuals.append(dict(exponent=exponent,observable=o.name,role=o.role,
                target=o.target,prediction=pred,scale=o.scale,units=o.units,normalized_residual=ri))
        summary.append(row);profiles.append(dict(**row,shape=shape,profile=fit))
        save_json(args.out/'checkpoint.json',dict(completed=False,profiles=profiles))
        print(f'screening {exponent:+g}: joint={row["joint_loss"]:.7g}, '
              f'coreRMS={row["core_force_RMS_eV_L0"]:.7g}, other={row["other_interface_loss"]:.7g}',flush=True)
    write_csv(args.out/'profile_summary.csv',summary)
    write_csv(args.out/'observable_residuals.csv',residuals)
    save_json(args.out/'completion.json',dict(completed=True,profiles=profiles,
        actual_profile_count=len(profiles),zero_profile_replayed=True,elapsed_seconds=time.perf_counter()-began,
        relaxed_core_recomputed=False,global_family_fit=False,new_law_required_proved=False,
        material_accepted=False,physical_yield_validated=False,physical_Hz=False))


if __name__=='__main__':main()
