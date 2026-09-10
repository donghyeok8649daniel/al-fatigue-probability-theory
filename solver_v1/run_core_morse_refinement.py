"""Energy/gradient stencil and matrix-free/dense audit at a saved saddle."""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh

from .report_current_material_core import restore_case
from .run_current_material_core import load_current_material,ROOT
from .run_source_core_reference import load_source_material,build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source-reference',action='store_true')
    args=p.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh stencil-audit output required')
    began=time.perf_counter(); loader=load_current_material; extra={}
    if args.source_reference:
        loader=load_source_material; extra['core_builder']=build_source_core
    model,tensor,binding=loader()
    meta=json.loads((args.case/'metadata.json').read_bytes())
    core,_,_,_,_=restore_case(ROOT/meta['baseline'],model,tensor,binding,**extra)
    with (args.case/'state.csv').open(encoding='utf8',newline='') as stream:
        saved={(int(r['j']),int(r['l'])):np.array([float(r[k]) for k in ('ux','uy','uz')])
               for r in csv.DictReader(stream)}
    field=np.array([saved[tuple(i)] for i in core.indices[core.free_ids]])
    value,action=core.linearize(field); matrix=action.explicit_matrix()
    ev,vectors=eigh(matrix,subset_by_index=[0,0]); mode=vectors[:,0].reshape(field.shape)
    analytic=float(ev[0]); eigen_residual=float(np.linalg.norm(action(mode)-analytic*mode))
    checks=[]
    for step in .002*2.**-np.arange(13):
        plus,minus=core.evaluate(field+step*mode),core.evaluate(field-step*mode)
        curvature=(plus['energy']+minus['energy']-2*value['energy'])/step**2
        gradient=float(np.sum(mode*(plus['gradient']-minus['gradient']))/(2*step))
        checks.append(dict(step_over_L0=step,analytic_curvature=analytic,
            energy_curvature=curvature,gradient_curvature=gradient,
            energy_error=abs(curvature-analytic),gradient_error=abs(gradient-analytic),
            plus_energy_change=plus['energy']-value['energy'],minus_energy_change=minus['energy']-value['energy']))
    write_csv(args.out/'stencil_refinement.csv',checks)
    save_json(args.out/'verification.json',dict(completed=True,
        source_reference_only=args.source_reference,material_sha256=binding['parameter_sha256'],
        minimum_eigenvalue=analytic,matrix_free_eigenpair_residual=eigen_residual,
        best_energy_stencil_error=min(r['energy_error'] for r in checks),
        best_gradient_stencil_error=min(r['gradient_error'] for r in checks),
        all_stencils_saved=True,finite_loop_activation_validated=False,
        physical_yield_calibrated=False,elapsed_seconds=time.perf_counter()-began))


if __name__=='__main__':
    main()
