"""Evaluate the public Si/O DFT configurations with the existing neutral MACE.

No fit/new DFT/MD. Source training data for MTP are not held-out MTP data;
overlap with MACE pretraining has not been audited. The source file has no
per-frame morphology labels, so composition groups are not called interfaces.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256

SOURCE_SHA = '949781e648dcf4d0b13fde73520461770232201e1761628e8da61150dcfe5d17'
SOURCE_COMMIT = '4ed842eea9c0b894310fee50e9bbca21404e356e'


def parse_source(path):
    raw = path.read_text(encoding='utf-8')
    if raw.count('BEGIN_CFG') != raw.count('END_CFG'):
        raise ValueError('incomplete CFG input')
    frames = []
    for part in raw.split('BEGIN_CFG')[1:]:
        lines = [s.strip() for s in part.split('END_CFG')[0].splitlines() if s.strip()]
        n = int(lines[lines.index('Size')+1])
        c = lines.index('Supercell')
        cell = np.array([[float(v) for v in s.split()] for s in lines[c+1:c+4]])
        a = next(i for i,s in enumerate(lines) if s.startswith('AtomData:'))
        fields = lines[a].split()[1:]
        atoms = np.array([[float(v) for v in s.split()] for s in lines[a+1:a+1+n]])
        col = {name:atoms[:,i] for i,name in enumerate(fields)}
        types = col['type'].astype(int)
        if (not np.array_equal(col['type'], types) or not set(types).issubset({0,1})
                or not np.array_equal(col['id'], np.arange(1,n+1))):
            raise ValueError('unexpected source type/id convention')
        positions = np.column_stack([col['cartes_'+x] for x in 'xyz'])
        forces = np.column_stack([col['f'+x] for x in 'xyz'])
        energy = float(lines[lines.index('Energy')+1])
        if cell.shape!=(3,3) or np.linalg.det(cell)<=0 or not np.all(np.isfinite(forces)):
            raise ValueError('invalid source cell/force')
        features = [s for s in lines if s.startswith('Feature')]
        frames.append(dict(positions=positions,cell=cell,forces=forces,energy=energy,
            numbers=np.where(types==0,14,8),features=features))
    return frames


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.source.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise ValueError('source identity mismatch')
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:
        raise ValueError('model identity mismatch')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'raw').mkdir()
    frames=parse_source(args.source)
    if len(frames)!=1159:
        raise ValueError('source inventory mismatch')
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    protocol=dict(source_sha256=SOURCE_SHA,source_commit=SOURCE_COMMIT,
        source_url='https://gitlab.com/Kazongogit/MTPu/-/blob/'+SOURCE_COMMIT+'/datasets/Unified_training_set2_1159.cfg',
        source_doi='10.1038/s41524-024-01390-8',source_frames=len(frames),
        source_method='Quantum ESPRESSO PBE PAW DFT; supplied energies eV and forces eV/Angstrom',
        element_mapping={'0':'Si','1':'O'},
        model_sha256=MODEL_SHA256,actual_new_DFT=0,actual_new_MD=0,fit_performed=False,
        held_out_mace_status='training overlap not audited',
        morphology_labels='not supplied per configuration; group only by explicit composition',
        stress_comparison=False,production_enabled=False)
    dump(args.output/'protocol.json',protocol)
    started=time.perf_counter(); rows=[]; references={}; groups={}; controls=[]
    with (args.output/'frames.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=None
        for i,frame in enumerate(frames):
            atoms=Atoms(numbers=frame['numbers'],positions=frame['positions'],cell=frame['cell'],pbc=True)
            energy,force=force_only(model,atoms)
            nsi=int(np.sum(frame['numbers']==14));no=int(np.sum(frame['numbers']==8))
            composition=f'Si{nsi}O{no}'
            group=('Si_only' if no==0 else 'O_only' if nsi==0 else
                   'SiO2_stoichiometry' if no==2*nsi else 'mixed_other_stoichiometry')
            if composition not in references:
                references[composition]=(i,energy,frame['energy'])
            anchor,ref_pred,ref_dft=references[composition]
            delta=force-frame['forces']
            row=dict(frame=i,atoms=len(atoms),n_Si=nsi,n_O=no,composition=composition,group=group,
                dft_energy_eV=frame['energy'],mace_energy_eV=energy,
                force_component_RMSE_eV_A=float(np.sqrt(np.mean(delta**2))),
                force_component_MAE_eV_A=float(np.mean(abs(delta))),
                force_component_max_error_eV_A=float(np.max(abs(delta))),
                source_force_component_RMS_eV_A=float(np.sqrt(np.mean(frame['forces']**2))),
                source_total_force_norm_eV_A=float(np.linalg.norm(frame['forces'].sum(axis=0))),
                mace_total_force_norm_eV_A=float(np.linalg.norm(force.sum(axis=0))),
                same_composition_reference_frame=anchor,
                relative_energy_error_eV_atom=((energy-ref_pred)-(frame['energy']-ref_dft))/len(atoms))
            for element in (8,14):
                mask=frame['numbers']==element
                row[f'Z{element}_force_component_RMSE_eV_A']=(float(np.sqrt(np.mean(delta[mask]**2))) if mask.any() else None)
            if writer is None:
                writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
            writer.writerow(row);stream.flush();rows.append(row)
            np.savez_compressed(args.output/'raw'/f'{i:04d}.npz',positions=frame['positions'],cell=frame['cell'],
                numbers=frame['numbers'],source_energy=frame['energy'],source_forces=frame['forces'],
                predicted_energy=energy,predicted_forces=force)
            g=groups.setdefault(group,dict(frames=0,atoms=0,squared_error=0.,absolute_error=0.,components=0))
            g['frames']+=1;g['atoms']+=len(atoms);g['squared_error']+=float(np.sum(delta**2))
            g['absolute_error']+=float(np.sum(abs(delta)));g['components']+=delta.size
            if i in (0,26,66,163,164,500,800,1158):
                atoms.calc=model; standard_energy=float(atoms.get_potential_energy());standard_force=atoms.get_forces()
                control=dict(frame=i,energy_difference_eV=standard_energy-energy,
                    force_max_difference_eV_A=float(np.max(abs(standard_force-force))))
                controls.append(control)
                if abs(control['energy_difference_eV'])>1e-8 or control['force_max_difference_eV_A']>1e-8:
                    dump(args.output/'failed_control.json',control);raise RuntimeError('force-only mismatch')
            if i%25==0 or i==len(frames)-1:
                dump(args.output/'progress.json',dict(completed=i+1,planned=len(frames),elapsed_s=time.perf_counter()-started,
                    latest=row,groups=groups))
                print('FRAME',i+1,'of',len(frames),'seconds',time.perf_counter()-started,flush=True)
            if time.perf_counter()-started>args.max_seconds:
                dump(args.output/'partial.json',dict(completed=i+1,planned=len(frames),reason='explicit wall-time budget'))
                raise RuntimeError('audit time budget; completed prefix preserved')
    summary=[]
    for name,g in groups.items():
        summary.append(dict(group=name,frames=g['frames'],atoms=g['atoms'],
            force_component_RMSE_eV_A=np.sqrt(g['squared_error']/g['components']),
            force_component_MAE_eV_A=g['absolute_error']/g['components']))
    with (args.output/'groups.csv').open('w',newline='',encoding='utf-8') as stream:
        w=csv.DictWriter(stream,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    dump(args.output/'summary.json',dict(complete=True,frames=len(rows),atoms=sum(x['atoms'] for x in rows),
        elapsed_s=time.perf_counter()-started,force_only_calls=len(rows),standard_controls=controls,
        groups=summary,protocol=protocol))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=5400);main(p.parse_args())
