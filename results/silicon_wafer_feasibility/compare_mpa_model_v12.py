"""Predeclared stratified comparison of two pinned neutral foundation models.

New MPA-0 calls on previously audited DFT arrays. Baseline MP-0b3 predictions
are reused exactly. Neither pretraining overlap nor initiation is certified.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.fit_oxide_radial_delta_v12 import energy_differences

MODEL_SHA='75428afe3a1d7d8062e19bcaabd5c433623cabf308242ec9fb493e38604fb638'
MODEL_URL='https://github.com/ACEsuit/mace-foundations/releases/download/mace_mpa_0/mace-mpa-0-medium.model'

def dump(path,x):path.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def select(source,per_group,dataset):
    records=list(csv.DictReader((source/'frames.csv').open(encoding='utf-8')));groups={}
    for row in records:
        group=row.get('observed_group',row.get('group'))
        key=(group,row['composition']) if dataset=='CP2K' else (group,)
        groups.setdefault(key,[]).append(row)
    selected=[]
    for key,rows in sorted(groups.items()):
        ids=np.linspace(0,len(rows)-1,min(per_group,len(rows)),dtype=int)
        for offset in ids:
            row=rows[offset];index=int(row.get('source_frame',row.get('frame')))
            selected.append(dict(dataset=dataset,group=row.get('observed_group',row.get('group')),
                composition=row['composition'],source_frame=index,raw_path=source/'raw'/f'{index:04d}.npz'))
    return selected

def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA:raise ValueError('MPA checkpoint changed')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True);(args.output/'raw').mkdir()
    selected=select(args.cp2k,12,'CP2K')+select(args.qe,16,'QE')
    selection=[{k:v for k,v in row.items() if k!='raw_path'} for row in selected]
    protocol=dict(model_name='MACE-MPA-0-medium',model_sha256=MODEL_SHA,model_url=MODEL_URL,
        baseline='MACE-MP-0b3-medium, stored v11 predictions',
        selection='CP2K: up to12 equispaced source-order frames per group+composition; QE: up to16 per chemistry group; fixed before new model evaluation',
        frames=len(selected),source_selection=selection,source_overlap_status='neither foundation pretraining overlap audited',
        energy_differences='within same dataset and exact composition; no fitted offsets',
        new_DFT=0,new_MD=0,fit_performed=False,production_enabled=False,
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    dump(args.output/'protocol.json',protocol)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    started=time.perf_counter();rows=[];stats={};controls=[];manifest=[]
    control_ids=set(np.linspace(0,len(selected)-1,9,dtype=int).tolist())
    for number,entry in enumerate(selected):
        with np.load(entry['raw_path']) as d:
            z=d['numbers'];r=d['positions'];cell=d['cell'];pbc=d['pbc'] if 'pbc' in d else np.ones(3,bool)
            reference=d['source_forces'];reference_energy=float(d['source_energy'])
            old=d['predicted_forces'];old_energy=float(d['predicted_energy'])
        atoms=Atoms(numbers=z,positions=r,cell=cell,pbc=pbc)
        energy,force=force_only(calc,atoms);b=old-reference;c=force-reference
        row=dict(**selection[number],atoms=len(z),baseline_force_RMSE_eV_A=float(np.sqrt(np.mean(b*b))),
            candidate_force_RMSE_eV_A=float(np.sqrt(np.mean(c*c))),reference_energy_eV=reference_energy,
            baseline_energy_eV=old_energy,candidate_energy_eV=energy)
        rows.append(row)
        np.savez_compressed(args.output/'raw'/f'{number:04d}.npz',numbers=z,positions=r,cell=cell,pbc=pbc,
            source_energy=reference_energy,source_forces=reference,baseline_energy=old_energy,baseline_forces=old,
            predicted_energy=energy,predicted_forces=force)
        manifest.append(dict(dataset=entry['dataset'],source_frame=entry['source_frame'],
            baseline_raw_sha256=hashlib.sha256(entry['raw_path'].read_bytes()).hexdigest()))
        name=entry['dataset']+'__'+entry['group']
        stat=stats.setdefault(name,dict(frames=0,atoms=0,baseline_sq=0.,candidate_sq=0.))
        stat['frames']+=1;stat['atoms']+=len(z);stat['baseline_sq']+=float(np.sum(b*b));stat['candidate_sq']+=float(np.sum(c*c))
        if number in control_ids:
            atoms.calc=calc;e=float(atoms.get_potential_energy());f=atoms.get_forces()
            control=dict(frame=number,energy_difference_eV=e-energy,force_difference_eV_A=float(np.max(abs(f-force))))
            controls.append(control)
            if abs(e-energy)>1e-8 or control['force_difference_eV_A']>1e-8:raise ValueError('standard model calculator differs')
        if number%10==0 or number+1==len(selected):
            dump(args.output/'running.json',dict(completed=number+1,planned=len(selected),elapsed_s=time.perf_counter()-started))
            table(args.output/'frames.csv',rows)
            print('MODEL_FRAME',number+1,len(selected),time.perf_counter()-started,flush=True)
        if time.perf_counter()-started>args.max_seconds:
            dump(args.output/'partial.json',dict(completed=number+1,planned=len(selected),reason='explicit time budget'))
            raise TimeoutError('model comparison budget; prefix retained')
    groups=[dict(group=name,frames=s['frames'],atoms=s['atoms'],
        baseline_force_RMSE_eV_A=float(np.sqrt(s['baseline_sq']/(3*s['atoms']))),
        candidate_force_RMSE_eV_A=float(np.sqrt(s['candidate_sq']/(3*s['atoms'])))) for name,s in stats.items()]
    pairs={dataset:energy_differences([r for r in rows if r['dataset']==dataset]) for dataset in ('CP2K','QE')}
    table(args.output/'source_raw_manifest.csv',manifest);table(args.output/'groups.csv',groups)
    dump(args.output/'summary.json',dict(complete=True,frames=len(rows),new_model_calls=len(rows)+len(controls),
        baseline_calls=0,standard_controls=controls,groups=groups,relative_energy_pairs=pairs,
        elapsed_s=time.perf_counter()-started,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('cp2k','qe','model','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=2400.);main(p.parse_args())
