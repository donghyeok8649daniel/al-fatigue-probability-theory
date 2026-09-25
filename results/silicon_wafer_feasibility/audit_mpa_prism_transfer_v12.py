"""Cross-model residuals on fixed baseline prism states, not re-relaxed candidates.

Model disagreement is not a measured DFT error bar. A source stationary in one
potential need not remain stationary in another; no Hessian is transferred.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.compare_mpa_model_v12 import MODEL_SHA,MODEL_URL
from solver_v1.silicon_initiation_research import grip_observables


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if args.output.exists():raise ValueError('fresh transfer output required')
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA:raise ValueError('MPA checkpoint differs')
    with np.load(args.geometry) as d:
        free=d['free'];lower=d['lower'];upper=d['upper'];area=float(d['area_A2'])
    sources=[(name,args.results/('dense_'+name)/'raw_hessian.npz') for name in ('force100','loading8','return8','loading10')]
    sources.append(('original_zero',args.results.parent/'silicon_initiation_v11/force_controlled_prism_360_tight/state_000/raw.npz'))
    optional=args.results/'return_zero_force/state_000'
    if (optional/'result.json').exists() and json.loads((optional/'result.json').read_text(encoding='utf-8'))['converged']:
        sources.append(('returned_zero',optional/'raw.npz'))
    if any(not path.exists() for _,path in sources):raise ValueError('baseline states unavailable')
    args.output.mkdir(parents=True)
    dump(args.output/'protocol.json',dict(model_sha256=MODEL_SHA,model_url=MODEL_URL,
        baseline_model='MP-0b3; stored predictions on original fixed geometries',
        selection='all four completed Hessian states, original zero force, and returned zero force if converged',
        source_files=[dict(state=name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for name,path in sources],
        geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        new_relaxation=False,truth_reference=None,comparison_is_error_bar=False,physical_clock=None))
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    started=time.perf_counter();rows=[];calls=0;controls=[]
    for name,path in sources:
        if time.perf_counter()-started>args.max_seconds:raise TimeoutError('transfer time budget; prefix retained')
        with np.load(path) as raw:
            positions=raw['positions'];numbers=raw['numbers'];cell=raw['cell'];old_f=raw['forces'];old_e=float(raw['energy'])
        atoms=Atoms(numbers=numbers,positions=positions,cell=cell,pbc=False)
        energy,force=force_only(calc,atoms);calls+=1
        baseline=grip_observables(positions,old_f,lower=lower,upper=upper,area_A2=area)
        candidate=grip_observables(positions,force,lower=lower,upper=upper,area_A2=area)
        np.savez_compressed(args.output/(name+'_raw.npz'),positions=positions,numbers=numbers,cell=cell,
            baseline_energy=old_e,baseline_forces=old_f,candidate_energy=energy,candidate_forces=force,
            free=free,lower=lower,upper=upper,area_A2=area)
        rows.append(dict(state=name,baseline_energy_eV=old_e,candidate_energy_eV=energy,
            baseline_free_force_max_eV_A=baseline['free_force_max_eV_A'],
            candidate_free_force_max_eV_A=candidate['free_force_max_eV_A'],
            baseline_nominal_stress_GPa=baseline['nominal_stress_GPa'],
            candidate_nominal_stress_GPa=candidate['nominal_stress_GPa'],
            free_force_disagreement_component_RMS_eV_A=float(np.sqrt(np.mean((force[free]-old_f[free])**2))),
            candidate_internal_force_norm_eV_A=candidate['total_internal_force_eV_A'],
            candidate_internal_torque_norm_eV=candidate['total_internal_torque_eV']))
        if name=='loading8':
            atoms.calc=calc;standard_e=float(atoms.get_potential_energy());standard_f=atoms.get_forces();calls+=1
            control=dict(state=name,energy_error_eV=standard_e-energy,force_error_eV_A=float(np.max(abs(standard_f-force))))
            if abs(control['energy_error_eV'])>1e-8 or control['force_error_eV_A']>1e-8:raise ValueError('standard ASE model differs')
            controls.append(control)
        dump(args.output/'running.json',dict(completed=len(rows),planned=len(sources),new_model_calls=calls,elapsed_s=time.perf_counter()-started))
    lookup={r['state']:r for r in rows};difference={}
    for model in ('baseline','candidate'):
        difference[model+'_return8_minus_loading8_eV']=lookup['return8'][model+'_energy_eV']-lookup['loading8'][model+'_energy_eV']
    result=dict(complete=True,states=rows,new_model_calls=calls,standard_controls=controls,
        same_grip_energy_comparison=difference,elapsed_s=time.perf_counter()-started,
        baseline_new_calls=0,new_relaxation=False,new_Hessian=False,new_DFT=0,new_MD=0,
        uncertainty_calibrated=False,material_approved=False,physical_clock=None,
        scope='disagreement at identical baseline geometries; no DFT truth or transferred equilibrium/stability claim')
    dump(args.output/'summary.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','geometry','model','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=300.);main(p.parse_args())
