"""Separate row-neighborhood/Bessel errors from same-state material mismatch."""
import argparse
import hashlib
import time
from pathlib import Path

import numpy as np

from .frozen_core_force_target import FrozenCoreForceTarget
from .report_current_material_core import restore_case
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material,build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists(): raise FileExistsError('fresh frozen-force refinement output required')
    began=time.perf_counter(); source,tensor,source_binding=load_source_material()
    source_path=ROOT/'results/current_material_core_v22/source_reference/escaped_plus_R8r5'
    core,field,*_=restore_case(source_path,source,tensor,source_binding,core_builder=build_source_core)
    frozen=FrozenCoreForceTarget(core,field,radius=2.)
    rows=[];differences=[];components=[];bindings={}
    for label,material in [('unchanged',load_current_material()),('joint_trial',load_current_material(args.trial))]:
        model,_,binding=material; bindings[label]=binding['parameter_sha256']; values={}
        for ring,tol in [(5,2e-12),(7,2e-12),(9,2e-12),(7,2e-13)]:
            result=frozen.coefficient_matrix(model,ring=ring,tolerance=tol)
            actual=result['design']@model.coefficients
            values[ring,tol]=actual
            error=actual-result['target']
            rows.append(dict(model=label,ring=ring,reciprocal_tolerance=tol,
                source_force_RMS=float(np.sqrt(np.mean(result['target']**2))),
                material_error_RMS=float(np.sqrt(np.mean(error**2))),
                maximum_material_error=float(np.max(abs(error)))))
            for i,value in enumerate(actual):
                components.append(dict(model=label,ring=ring,tolerance=tol,component=i,
                    candidate_gradient=value,source_gradient=result['target'][i],units='eV/L0'))
        for first,second,kind in [((5,2e-12),(7,2e-12),'row_ring'),
                                  ((7,2e-12),(9,2e-12),'row_ring'),
                                  ((7,2e-12),(7,2e-13),'reciprocal')]:
            change=values[second]-values[first]
            differences.append(dict(model=label,study=kind,first_ring=first[0],second_ring=second[0],
                first_tolerance=first[1],second_tolerance=second[1],
                maximum_force_change=float(np.max(abs(change))),RMS_force_change=float(np.sqrt(np.mean(change**2))),
                units='eV/L0',infinite_domain_certified=False))
        print(label,'frozen-force refinement complete',flush=True)
    write_csv(args.out/'material_signal.csv',rows)
    write_csv(args.out/'numerical_changes.csv',differences)
    write_csv(args.out/'gradient_components.csv',components)
    save_json(args.out/'completion.json',dict(completed=True,actual_represented_state_evaluations=True,
        coefficient_parameter_sha256=bindings,source_sha256=source_binding['parameter_sha256'],
        source_state_sha256=hashlib.sha256((source_path/'state.csv').read_bytes()).hexdigest(),
        source_free_radius=8.,force_observation_radius=2.,source_positions_frozen=True,
        measured_changes_are_not_rigorous_infinite_tail_bounds=True,material_accepted=False,
        elapsed_seconds=time.perf_counter()-began,
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__': main()
