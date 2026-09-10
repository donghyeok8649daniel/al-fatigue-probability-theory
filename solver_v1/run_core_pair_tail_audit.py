"""Execute pair-only infinite-force majorants on actual relaxed v22 cores."""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from .core_pair_tail_bound import pair_force_tail_bound,explicit_pair_gradient
from .current_material_rows import CurrentMaterialScrewCore
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
    if args.out.exists():raise FileExistsError('fresh analytic tail audit output required')
    root=ROOT/'results/current_material_core_v22';rows=[];changes=[];components=[];bindings={}
    old=load_current_material();trial=load_current_material(args.trial)
    source=load_source_material()
    for label,material,path,restore,extra in [
            ('unchanged',old,root/'stable_plus_R8r5',old,{}),
            ('joint_trial',trial,root/'probe_core_plus_R8r5',trial,{}),
            ('unchanged_at_source',old,root/'source_reference/escaped_plus_R8r5',source,dict(core_builder=build_source_core)),
            ('joint_trial_at_source',trial,root/'source_reference/escaped_plus_R8r5',source,dict(core_builder=build_source_core))]:
        model,tensor,binding=material;bindings[label]=binding['parameter_sha256']
        reference,field,*_=restore_case(path,*restore,**extra)
        values={}
        for ring in (3,5,7,9,13):
            core=CurrentMaterialScrewCore(model,reference.far_field,free_radius=reference.free_radius,
                ring=ring,tolerance=2e-13)
            if not np.array_equal(core.indices[core.free_ids],reference.indices[reference.free_ids]):
                raise ValueError('tail check must keep the same relaxed free configuration')
            if np.max(abs(core.boundary[:,1:]))>0:
                raise ValueError('derive a new global displacement bound for a transverse exterior')
            maximum=float(np.max(np.linalg.norm(field[:,1:],axis=1)))
            bound=pair_force_tail_bound(ring=ring,b=core.rows.b,h=core.rows.h,
                u=model.coefficients[0],v=model.coefficients[1],maximum_transverse_displacement=maximum)
            actual=explicit_pair_gradient(core,field,u=model.coefficients[0],v=model.coefficients[1])
            values[ring]=actual
            rows.append(dict(model=label,ring=ring,maximum_free_transverse_displacement=maximum,**bound))
            for logical,gradient in zip(core.indices[core.free_ids],actual):
                components.append(dict(model=label,ring=ring,j=int(logical[0]),l=int(logical[1]),
                    gx=gradient[0],gy=gradient[1],gz=gradient[2],units='eV/L0'))
        for first,second in ((3,5),(5,7),(7,9),(9,13)):
            change=float(np.max(np.linalg.norm(values[second]-values[first],axis=1)))
            bound=next(r['pair_gradient_norm_bound_eV_L0'] for r in rows if r['model']==label and r['ring']==first)
            if change>bound:raise ArithmeticError('actual omitted subset exceeds the analytic majorant')
            changes.append(dict(model=label,first_ring=first,second_ring=second,
                maximum_pair_force_change_eV_L0=change,coarse_infinite_majorant=bound,
                majorant_over_actual_change=bound/change,full_environment_force_certified=False))
    write_csv(args.out/'infinite_pair_force_majorants.csv',rows)
    write_csv(args.out/'actual_pair_force_changes.csv',changes)
    write_csv(args.out/'pair_gradient_components.csv',components)
    save_json(args.out/'completion.json',dict(completed=True,parameter_sha256=bindings,
        source='analytic LJ/Bessel infinite atomic rows and exact Hurwitz-zeta upper series',
        arbitrary_neighbor_cutoff_as_theory=False,force_repaired=False,
        full_environmental_tail_certified=False,finite_free_boundary_certified=False,
        physical_yield_calibrated=False,
        csv_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__':main()
