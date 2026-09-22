"""Same LOCAL normal traction on each rigid concentration-dependent path.

Locate actual MACE stationary points of E(q)-A*T*q by bracketed force roots.
This is a constrained 0 K path barrier, not a finite-T activation free energy,
wafer strength, a crack-tip amplification model, or calibrated failure rate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_concentration_research import constrained_energy_derivative
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.run_concentration_cleavage_v10 import dump,csv_write


def main(args):
    import torch
    from ase import Atoms
    from ase.io import read
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model identity')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    if not np.isfinite(args.traction_MPa) or args.traction_MPa<=0:raise ValueError('positive finite local traction')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    cases=args.cases
    if cases is None:
        source_summary=json.loads((args.source/'summary.json').read_text())
        if not source_summary['complete']:raise ValueError('source sweep incomplete')
        cases=[r['case'] for r in source_summary['cases']]
    all_results=[];total_calls=0
    for case in cases:
        if args.follow_source:
            wait_started=time.monotonic()
            while True:
                try:
                    progress=json.loads((args.source/'progress.json').read_text())
                    ready=any(r['case']==case for r in progress['cases'])
                except (FileNotFoundError,json.JSONDecodeError):ready=False
                if ready:break
                if time.monotonic()-wait_started>3600:
                    raise TimeoutError('source case did not finish within one hour: '+case)
                time.sleep(10)
        started=time.perf_counter();folder=args.output/case;folder.mkdir()
        source=args.source/case;metadata=json.loads((source/'result.json').read_text())
        base=read(source/'bulk.extxyz');base.calc=None
        # extxyz keeps eight decimal places for positions; use exact raw npz.
        reference=np.load(source/'opening_0.000.npz')
        r0,c0,numbers=reference['positions'],reference['cell'],reference['numbers']
        area=float(np.linalg.norm(np.cross(c0[0],c0[1])))
        load=area*args.traction_MPa/160217.6634
        cache={};records=[];calls=0
        def evaluate(q):
            nonlocal calls,total_calls
            q=float(q)
            if q in cache:return cache[q]
            cell=c0.copy();cell[2,2]+=q
            if cell[2,2]<=np.max(r0[:,2]):raise ValueError('compression crosses a stored atom; enlarge domain geometry carefully')
            atoms=Atoms(numbers=numbers,positions=r0,cell=cell,pbc=True,calculator=model)
            energy=float(atoms.get_potential_energy());force=atoms.get_forces();stress=atoms.get_stress(voigt=False)
            calls+=1;total_calls+=1
            dc=np.zeros((3,3));dc[2,2]=1
            derivative=constrained_energy_derivative(r0,force,cell,stress,position_tangent=np.zeros_like(r0),cell_tangent=dc)
            row=dict(case=case,opening_A=q,energy_eV=energy,energy_derivative_eV_A=derivative,
                     tilted_energy_eV=energy-load*q,tilted_derivative_eV_A=derivative-load,new_evaluation=True)
            cache[q]=row;records.append(row)
            np.savez_compressed(folder/f'evaluation_{calls:03d}.npz',positions=r0,cell=cell,numbers=numbers,
                energy=energy,forces=force,stress=stress,opening_A=q,energy_derivative=derivative)
            return row
        for file in sorted(source.glob('opening_*.npz')):
            raw=np.load(file);q=float(file.stem.split('_')[1])
            derivative=constrained_energy_derivative(raw['positions'],raw['forces'],raw['cell'],raw['stress'],
                position_tangent=raw['position_tangent'],cell_tangent=raw['cell_tangent'])
            row=dict(case=case,opening_A=q,energy_eV=float(raw['energy']),energy_derivative_eV_A=derivative,
                tilted_energy_eV=float(raw['energy'])-load*q,tilted_derivative_eV_A=derivative-load,new_evaluation=False)
            cache[q]=row;records.append(row)
        # Signed compression extends the path: q=0 is NOT an artificial wall.
        for q in [-.25,-.1,-.03]:evaluate(q)
        grid=sorted(cache)
        if cache[grid[0]]['tilted_derivative_eV_A']>=0:
            evaluate(-.5);grid=sorted(cache)
        roots=[]
        for left,right in zip(grid[:-1],grid[1:]):
            fl,fr=cache[left]['tilted_derivative_eV_A'],cache[right]['tilted_derivative_eV_A']
            if fl*fr<0:
                q=float(brentq(lambda x:evaluate(x)['tilted_derivative_eV_A'],left,right,xtol=2e-7,rtol=1e-12))
                value=evaluate(q)
                h=1e-4;minus,plus=evaluate(q-h),evaluate(q+h)
                curvature=(plus['energy_derivative_eV_A']-minus['energy_derivative_eV_A'])/(2*h)
                fd=(plus['energy_eV']-minus['energy_eV'])/(2*h)
                roots.append(dict(**value,curvature_eV_A2=curvature,
                    finite_difference_eV_A=fd,force_difference_error_eV_A=fd-value['energy_derivative_eV_A'],
                    bracket_A=[left,right],kind='minimum' if fl<fr else 'maximum'))
        minima=[r for r in roots if r['kind']=='minimum' and r['curvature_eV_A2']>0]
        if not minima:
            dump(folder/'failure.json',dict(reason='no metastable minimum bracketed',roots=roots))
            raise RuntimeError('no metastable constrained minimum: '+case)
        well=min(minima,key=lambda r:abs(r['opening_A']))
        saddles=[r for r in roots if r['opening_A']>well['opening_A'] and r['kind']=='maximum' and r['curvature_eV_A2']<0]
        if not saddles:raise RuntimeError('no outward constrained maximum: '+case)
        saddle=min(saddles,key=lambda r:r['opening_A'])
        barrier=saddle['tilted_energy_eV']-well['tilted_energy_eV']
        if barrier<=0:raise RuntimeError('nonpositive constrained barrier')
        error=max(abs(r['force_difference_error_eV_A']) for r in roots)
        if error>5e-4:raise AssertionError('energy/force consistency failure at actual loaded roots')
        result=dict(case=case,species=metadata['species'],dopant_count=metadata['dopant_count'],
            arrangement=metadata['arrangement'],chemical_concentration_cm3=metadata['chemical_concentration_cm3'],
            local_traction_MPa=args.traction_MPa,area_A2=area,conjugate_external_force_eV_A=load,
            minimum_opening_A=well['opening_A'],maximum_opening_A=saddle['opening_A'],
            constrained_path_barrier_eV=barrier,barrier_per_area_J_m2=barrier/area*16.02176634,
            root_force_max_error_eV_A=max(abs(r['tilted_derivative_eV_A']) for r in roots),
            root_energy_difference_max_error_eV_A=error,stationary_points=roots,
            new_standard_calculator_calls=calls,elapsed_s=time.perf_counter()-started,
            specimen_remote_stress_MPa=None,physical_probability=None,physical_rate=None,
            scope='rigid many-atom 0 K path under local normal traction; not a thermal or specimen barrier')
        csv_write(folder/'evaluations.csv',records);dump(folder/'result.json',result);all_results.append(result)
        dump(args.output/'progress.json',dict(completed=len(all_results),planned=len(cases),cases=all_results))
        print('LOADED',case,'barrier',barrier,'q',well['opening_A'],saddle['opening_A'],'calls',calls,flush=True)
    dump(args.output/'summary.json',dict(complete=True,cases=all_results,model_sha256=MODEL_SHA256,
        local_traction_MPa=args.traction_MPa,new_standard_calculator_calls=total_calls,
        physical_probability=None,production_enabled=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+');p.add_argument('--traction-MPa',dest='traction_MPa',type=float,default=100.)
    p.add_argument('--follow-source',action='store_true',help='process explicit cases as their source runs complete')
    main(p.parse_args())
