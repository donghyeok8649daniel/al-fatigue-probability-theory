"""Frozen core path and exact coefficient work; not a fitted/relaxed barrier.

Each model retains its own material/exterior. Interpolate identical interior
atomic coordinates between the original candidate and source stable R8 cores.
Term energies are in the declared gauge and are NOT separately observable
physical energies. Directional contributions sum, component RMS values do not.
"""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from .current_core_coefficient_basis import current_core_coefficients,NAMES
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
    if args.out.exists():raise FileExistsError('fresh mechanism audit output required')
    root=ROOT/'results/current_material_core_v22'
    old=load_current_material();trial=load_current_material(args.trial);source=load_source_material()
    first=root/'stable_plus_R8r5';second=root/'source_reference/escaped_plus_R8r5'
    c0,q0,*_=restore_case(first,*old)
    cs,q1,*_=restore_case(second,*source,core_builder=build_source_core)
    if not np.array_equal(c0.indices[c0.free_ids],cs.indices[cs.free_ids]):
        raise ValueError('same actual interior atomic rows required')
    # Continuous representative near q0; shifting an entire infinite row by
    # b is exact phase periodicity, not multiplying/rescaling displacements.
    delta=q1-q0;delta[:,0]-=c0.rows.b*np.floor(delta[:,0]/c0.rows.b+.5)
    records=[];terms=[];bindings={}
    for name,material,path,extra in [('unchanged',old,first,{}),
            ('joint_trial',trial,root/'probe_core_plus_R8r5',{}),
            ('Mishin_target',source,second,dict(core_builder=build_source_core))]:
        core,_,_,_,_=restore_case(path,*material,**extra)
        bindings[name]=material[2]['parameter_sha256']
        for fraction in (0.,.25,.5,.75,1.):
            field=q0+fraction*delta;value,action=core.linearize(field)
            derivative=float(np.sum(value['gradient']*delta))
            curvature=float(np.sum(delta*action(delta)))
            record=dict(model=name,fraction=fraction,energy_eV_repeat=value['energy'],
                derivative_eV_per_fraction=derivative,curvature_eV_per_fraction_squared=curvature,
                maximum_atomic_gradient_eV_L0=float(np.max(abs(value['gradient']))),
                own_fixed_exterior=True,minimum_energy_path=False,activation_barrier=False)
            if name!='Mishin_target':
                columns=current_core_coefficients(core,field);c=material[0].coefficients
                E=columns['energy']*c
                D=np.einsum('nic,ni->c',columns['gradient'],delta)*c
                if abs(E.sum()-value['energy'])>2e-9 or abs(D.sum()-derivative)>2e-9:
                    raise ArithmeticError('coefficient energy/directional work did not sum')
                record.update(energy_sum_error=float(abs(E.sum()-value['energy'])),
                              work_sum_error=float(abs(D.sum()-derivative)))
                for parameter,energy,work in zip(NAMES,E,D):
                    terms.append(dict(model=name,fraction=fraction,term=parameter,
                        energy_eV_repeat=energy,directional_work_eV=work,
                        term_is_gauge_dependent=True,independent_physical_observable=False))
            records.append(record)
    write_csv(args.out/'frozen_atomic_path.csv',records)
    write_csv(args.out/'gauge_fixed_term_contributions.csv',terms)
    save_json(args.out/'definition.json',dict(completed=True,parameter_sha256=bindings,
        narrow_state=first.relative_to(ROOT).as_posix(),wide_state=second.relative_to(ROOT).as_posix(),
        narrow_state_sha256=hashlib.sha256((first/'state.csv').read_bytes()).hexdigest(),
        wide_state_sha256=hashlib.sha256((second/'state.csv').read_bytes()).hexdigest(),
        interpolation='q(lambda)=q_narrow+lambda*nearest_row_phase(q_wide-q_narrow)',
        radius_over_L0=8.,ring=5,own_material_exterior=True,
        not_minimum_energy_path=True,not_relaxed_barrier=True,no_material_adoption=True,
        physical_yield=False,physical_Hz=False,
        csv_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__':main()
