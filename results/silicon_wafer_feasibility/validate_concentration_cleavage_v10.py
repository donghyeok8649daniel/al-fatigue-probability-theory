"""Independent energy differences, symmetry and fixed-concentration size checks."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_concentration_research import constrained_energy_derivative
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.run_concentration_cleavage_v10 import dump,csv_write


def main(args):
    import torch
    from ase import Atoms
    from ase.io import read
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model identity')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    calls=0;started=time.perf_counter()
    def evaluate(atoms):
        nonlocal calls
        calls+=1
        return force_only(model,atoms)
    all_fd,all_symmetry,all_size,all_backend=[],[],[],[]
    for case in args.cases:
        folder=args.source/case
        result=json.loads((folder/'result.json').read_text())
        base=read(folder/'bulk.extxyz');base.calc=None
        # Direct Cartesian energies at q +/- h; no work/stress formula here.
        for q in [.35,1.1,3.]:
            raw=np.load(folder/f'opening_{q:.3f}.npz')
            true=constrained_energy_derivative(raw['positions'],raw['forces'],raw['cell'],raw['stress'],
                position_tangent=raw['position_tangent'],cell_tangent=raw['cell_tangent'])
            for step in [1e-3,3e-4]:
                energies=[]
                for sign in [-1,1]:
                    cell=raw['cell'].copy();cell[2,2]+=sign*step
                    atom=Atoms(numbers=raw['numbers'],positions=raw['positions'],cell=cell,pbc=True)
                    energy,force=evaluate(atom);energies.append(energy)
                fd=(energies[1]-energies[0])/(2*step)
                all_fd.append(dict(case=case,opening_A=q,step_A=step,
                    energy_minus_eV=energies[0],energy_plus_eV=energies[1],
                    analytical_eV_A=true,finite_difference_eV_A=fd,error_eV_A=fd-true))
        raw=np.load(folder/'opening_0.700.npz')
        atom=Atoms(numbers=raw['numbers'],positions=raw['positions'],cell=raw['cell'],pbc=True)
        e,f=evaluate(atom)
        all_backend.append(dict(case=case,energy_error_eV=e-float(raw['energy']),
            force_error_eV_A=float(np.max(abs(f-raw['forces'])))))
        # Nontrivial proper rotation and reversed atom order independently.
        angle=.431;rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
        ids=np.arange(len(atom))[::-1]
        rotated=Atoms(numbers=atom.numbers[ids],positions=atom.positions[ids]@rotation,
                      cell=atom.cell.array@rotation,pbc=True)
        er,fr=evaluate(rotated)
        all_symmetry.append(dict(case=case,energy_error_eV=er-e,
            force_error_eV_A=float(np.max(abs(fr-f[ids]@rotation)))))
        # Repeat whole chemistry, so reference concentration is EXACTLY fixed.
        # z repeat tests slab thickness; x repeat tests area extensivity.
        for repeats in [(1,1,1),(1,1,2),(2,1,1)]:
            repeated=base.repeat(repeats)
            material_volume=repeated.get_volume()
            count=int(np.count_nonzero(repeated.numbers!=14))
            reference=None
            for q in [0.,.7,3.,10.]:
                opened=repeated.copy();c=opened.cell.array.copy();c[2,2]+=q;opened.set_cell(c,scale_atoms=False)
                energy,force=evaluate(opened)
                if q==0:reference=energy
                area=float(np.linalg.norm(np.cross(c[0],c[1])))
                all_size.append(dict(case=case,repeats='x'.join(map(str,repeats)),atoms=len(opened),
                    opening_A=q,reference_material_volume_A3=material_volume,chemical_concentration_cm3=count*1e24/material_volume,
                    area_A2=area,energy_eV=energy,reference_eV=reference,
                    separation_work_J_m2=(energy-reference)/area*16.02176634))
                np.savez_compressed(args.output/f'{case}_{repeats[0]}{repeats[1]}{repeats[2]}_{q:.1f}.npz',
                    positions=opened.positions,cell=c,numbers=opened.numbers,energy=energy,forces=force)
        print('VALIDATED',case,'calls',calls,flush=True)
    csv_write(args.output/'finite_differences.csv',all_fd)
    csv_write(args.output/'backend.csv',all_backend)
    csv_write(args.output/'symmetry.csv',all_symmetry)
    csv_write(args.output/'size_controls.csv',all_size)
    finest=max(abs(r['error_eV_A']) for r in all_fd if r['step_A']==3e-4)
    size_error=0.
    for row in all_size:
        reference=next(r for r in all_size if r['case']==row['case'] and r['opening_A']==row['opening_A'] and r['repeats']=='1x1x1')
        size_error=max(size_error,abs(row['separation_work_J_m2']-reference['separation_work_J_m2']))
    ok=(finest<2e-4 and max(abs(r['energy_error_eV']) for r in all_symmetry)<1e-8
        and max(r['force_error_eV_A'] for r in all_symmetry)<1e-8
        and max(abs(r['energy_error_eV']) for r in all_backend)<1e-9
        and max(r['force_error_eV_A'] for r in all_backend)<1e-9 and size_error<1e-6)
    dump(args.output/'summary.json',dict(passed=ok,cases=args.cases,force_only_calculator_calls=calls,
        finest_force_difference_max_eV_A=finest,max_size_work_difference_J_m2=size_error,
        elapsed_s=time.perf_counter()-started,model_sha256=MODEL_SHA256,
        scope='independent derivatives and replication; not global stability or material calibration'))
    if not ok:raise AssertionError('independent validation failed; preserve raw data')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+',default=['Si_00_pure','B_08_interface_cluster','P_08_interface_cluster'])
    main(p.parse_args())
