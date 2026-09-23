"""Neutral MACE versus published CP2K Si/O/H reference configurations.

Evaluate every frame up to an explicit atom-count resource bound. The bound is
chosen before seeing predictions; all excluded frame ids remain in inventory.
Source 'bulk' labels include H-containing slabs; do not silently remove H.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256

SOURCE_SHA='d99170bff4942dcd5a1c38d998fe0a60084a5e132bcb3cf2335f5db4cbe17689'
SOURCE_COMMIT='b61806d3a097e21d7880db39dcc077573e2c3c92'


def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def group_name(numbers,positions,cell):
    n=len(numbers);si=int(np.sum(numbers==14));o=int(np.sum(numbers==8));h=int(np.sum(numbers==1))
    chemical='Si_O_H' if si and o and h else 'Si_H' if si and h else 'Si_O' if si and o else 'Si' if si else 'O' if o else 'H'
    if n<=2:return 'small_'+chemical,[None,None,None]
    inverse=np.linalg.inv(cell);scaled=positions@inverse;gaps=[]
    for j in range(3):
        s=np.sort(scaled[:,j]%1);gaps.append(float(np.diff(np.r_[s,s[0]+1]).max()/np.linalg.norm(inverse[:,j])))
    return chemical+'_vacuum'+str(sum(gap>8 for gap in gaps)),gaps


def main(args):
    import torch
    from ase.io import iread
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.source.read_bytes()).hexdigest()!=SOURCE_SHA:raise ValueError('source changed')
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model changed')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True);(args.output/'raw').mkdir()
    inventory=[];selected=[]
    for index,atoms in enumerate(iread(args.source,index=':',format='extxyz')):
        if not set(atoms.numbers).issubset({1,8,14}):raise ValueError('unexpected chemical species')
        numbers=atoms.numbers.copy();positions=atoms.positions.copy();cell=atoms.cell.array.copy()
        name,gaps=group_name(numbers,positions,cell)
        entry=dict(source_frame=index,atoms=len(atoms),n_Si=int(np.sum(numbers==14)),n_O=int(np.sum(numbers==8)),
            n_H=int(np.sum(numbers==1)),source_config_type=atoms.info.get('config_type'),observed_group=name,
            empty_gap_x_A=gaps[0],empty_gap_y_A=gaps[1],empty_gap_z_A=gaps[2],included=len(atoms)<=args.max_atoms)
        inventory.append(entry)
        if entry['included']:
            if atoms.calc is None or 'energy' not in atoms.calc.results or 'forces' not in atoms.calc.results:
                raise ValueError('DFT energy/force missing')
            energy=float(atoms.calc.results['energy']);forces=np.array(atoms.calc.results['forces'],copy=True)
            if not np.isfinite(energy) or not np.all(np.isfinite(forces)):raise ValueError('nonfinite source')
            atoms.calc=None;selected.append((atoms,energy,forces,entry))
    if len(inventory)!=1466:raise ValueError('source inventory changed')
    with (args.output/'source_inventory.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(inventory[0]));w.writeheader();w.writerows(inventory)
    protocol=dict(source_sha256=SOURCE_SHA,source_commit=SOURCE_COMMIT,
        source_url='https://github.com/lukas-cvitkovich/MLFF-SiOx/tree/'+SOURCE_COMMIT,
        source_doi='10.1063/5.0220091',source_method='CP2K GPW PBE GTH, double-zeta Gaussian basis, 650 Ry cutoff; reported spin protocol in paper',
        total_source_frames=len(inventory),selected_frames=len(selected),max_atoms=args.max_atoms,
        selection='all source frames with N <= max_atoms, before predictions; excluded frames remain in inventory',
        labels='source config_type preserved; auxiliary labels describe composition and observed empty cell directions only',
        vacuum_label_threshold_A=8.,vacuum_label_is_crack_detector=False,
        model_sha256=MODEL_SHA256,model_elements_evaluated=['H','O','Si'],
        force_units='eV/Angstrom',energy_units='eV',actual_new_DFT=0,actual_new_MD=0,
        fit_performed=False,held_out_mace_status='pretraining overlap not audited',production_enabled=False)
    dump(args.output/'protocol.json',protocol)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    started=time.perf_counter();rows=[];references={};groups={};controls=[]
    control_ids={selected[j][3]['source_frame'] for j in np.linspace(0,len(selected)-1,9,dtype=int)}
    with (args.output/'frames.csv').open('w',newline='',encoding='utf-8') as f:
        writer=None
        for number,(atoms,source_energy,source_forces,entry) in enumerate(selected):
            index=entry['source_frame'];energy,forces=force_only(model,atoms)
            error=forces-source_forces;name=entry['observed_group'];composition=atoms.get_chemical_formula()
            if composition not in references:references[composition]=(index,source_energy,energy)
            anchor,dft_ref,ml_ref=references[composition]
            row=dict(source_frame=index,atoms=len(atoms),composition=composition,observed_group=name,
                source_energy_eV=source_energy,mace_energy_eV=energy,
                force_component_RMSE_eV_A=float(np.sqrt(np.mean(error**2))),force_component_MAE_eV_A=float(np.mean(abs(error))),
                force_component_max_error_eV_A=float(np.max(abs(error))),source_force_component_RMS_eV_A=float(np.sqrt(np.mean(source_forces**2))),
                source_total_force_norm_eV_A=float(np.linalg.norm(source_forces.sum(axis=0))),
                mace_total_force_norm_eV_A=float(np.linalg.norm(forces.sum(axis=0))),
                same_composition_reference_frame=anchor,
                relative_energy_error_eV_atom=((energy-ml_ref)-(source_energy-dft_ref))/len(atoms))
            for element in (1,8,14):
                mask=atoms.numbers==element
                row[f'Z{element}_force_RMSE_eV_A']=float(np.sqrt(np.mean(error[mask]**2))) if mask.any() else None
            if writer is None:writer=csv.DictWriter(f,fieldnames=list(row));writer.writeheader()
            writer.writerow(row);f.flush();rows.append(row)
            np.savez_compressed(args.output/'raw'/f'{index:04d}.npz',numbers=atoms.numbers,positions=atoms.positions,
                cell=atoms.cell.array,pbc=atoms.pbc,source_energy=source_energy,source_forces=source_forces,
                predicted_energy=energy,predicted_forces=forces)
            stats=groups.setdefault(name,dict(frames=0,atoms=0,sum_sq=0.,sum_abs=0.))
            stats['frames']+=1;stats['atoms']+=len(atoms);stats['sum_sq']+=float(np.sum(error**2));stats['sum_abs']+=float(np.sum(abs(error)))
            if index in control_ids:
                atoms.calc=model;standard_e=float(atoms.get_potential_energy());standard_f=atoms.get_forces();atoms.calc=None
                control=dict(source_frame=index,energy_difference_eV=standard_e-energy,
                    force_max_difference_eV_A=float(np.max(abs(standard_f-forces))))
                controls.append(control)
                if abs(standard_e-energy)>1e-8 or np.max(abs(standard_f-forces))>1e-8:
                    dump(args.output/'failed_control.json',control);raise RuntimeError('standard MACE comparison failed')
            if number%20==0 or number+1==len(selected):
                dump(args.output/'progress.json',dict(completed=number+1,planned=len(selected),latest=row,
                    elapsed_s=time.perf_counter()-started,groups=groups))
                print('FRAME',number+1,'/',len(selected),'source',index,'seconds',time.perf_counter()-started,flush=True)
            if time.perf_counter()-started>args.max_seconds:
                dump(args.output/'partial.json',dict(completed=number+1,planned=len(selected),reason='explicit wall-time budget'))
                raise RuntimeError('time budget; prefix preserved')
    group_rows=[dict(group=k,frames=s['frames'],atoms=s['atoms'],force_component_RMSE_eV_A=float(np.sqrt(s['sum_sq']/(3*s['atoms']))),
        force_component_MAE_eV_A=s['sum_abs']/(3*s['atoms'])) for k,s in groups.items()]
    with (args.output/'groups.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(group_rows[0]));w.writeheader();w.writerows(group_rows)
    dump(args.output/'summary.json',dict(complete=True,selected_frames=len(rows),excluded_frames=len(inventory)-len(rows),
        atoms=sum(row['atoms'] for row in rows),force_only_calls=len(rows),standard_controls=controls,
        elapsed_s=time.perf_counter()-started,groups=group_rows,protocol=protocol))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--model',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--max-atoms',type=int,default=320)
    p.add_argument('--max-seconds',type=float,default=7200);main(p.parse_args())
