"""New neutral MACE evaluations on all hash-bound Cambridge Si configurations.

Same 2475 geometries as v5 SW/Tersoff. No fitting or new DFT/MD, and no claim
of held-out status for an independently trained pretrained model.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames, reference_fields, ARCHIVE_SHA256, DATA_SHA256
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256, MODEL_URL
from solver_v1.silicon_atomistic_reference import force_metrics, periodic_basis_audit


def main(args):
    import torch
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:
        raise ValueError('model checksum mismatch')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output folder required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    frames,members=read_frames(args.archive)
    refs=[reference_fields(a) for a in frames]
    anchor=min((i for i,a in enumerate(frames) if a.info['config_type']=='dia' and a.info.get('xc_functional')=='PW91'),key=lambda i:refs[i][0]/len(frames[i]))
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    started=time.perf_counter()
    rows=[];predicted_force=[];predicted_energy=[];transforms=[]
    a=frames[anchor].copy();a.set_cell(periodic_basis_audit(a.cell.array,a.pbc)[0],scale_atoms=False);a.calc=calc
    anchor_prediction=float(a.get_potential_energy()/len(a))
    anchor_dft=refs[anchor][0]/len(frames[anchor])
    with (args.output/'frame_progress.jsonl').open('w',encoding='utf-8') as progress:
        for i,(atoms,(ed,fd,prefix)) in enumerate(zip(frames,refs)):
            a=atoms.copy();cell,signs=periodic_basis_audit(a.cell.array,a.pbc)
            a.set_cell(cell,scale_atoms=False);a.calc=calc
            ep=float(a.get_potential_energy());fp=a.get_forces().copy()
            xc=str(a.info.get('xc_functional','UNSPECIFIED'))
            row=dict(model='MACE-MP-0b3-medium',frame=i,atoms=len(a),config_type=str(a.info['config_type']),
                     declared_xc=xc,dft_field_prefix=prefix,dft_energy_eV=ed,model_energy_eV=ep,
                     basis_signs=' '.join(map(str,signs)),volume_A3=float(abs(np.linalg.det(cell))),
                     **force_metrics(fp,fd))
            row['diamond_referenced_error_eV_atom']=(ep/len(a)-anchor_prediction-(ed/len(a)-anchor_dft) if xc=='PW91' else None)
            rows.append(row);predicted_force.append(fp);predicted_energy.append(ep)
            progress.write(json.dumps(row)+'\n');progress.flush()
            if np.any(signs!=1):transforms.append(dict(frame=i,signs=signs.tolist(),positions_changed=False))
            if (i+1)%100==0:
                print(json.dumps(dict(evaluated=i+1,total=len(frames),elapsed_reported_s=time.perf_counter()-started)),flush=True)
    offsets=np.r_[0,np.cumsum([len(a) for a in frames])]
    np.savez_compressed(args.output/'predictions.npz',offsets=offsets,
        reference_force=np.vstack([f for _,f,_ in refs]),reference_energy=np.array([e for e,_,_ in refs]),
        predicted_force=np.vstack(predicted_force),predicted_energy=np.array(predicted_energy))
    with (args.output/'frame_metrics.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    groups=defaultdict(list)
    for row in rows:groups[(row['config_type'],row['declared_xc'])].append(row)
    summaries=[]
    for (kind,xc),r in sorted(groups.items()):
        components=sum(s['components'] for s in r)
        de=[s['diamond_referenced_error_eV_atom'] for s in r if s['diamond_referenced_error_eV_atom'] is not None]
        summaries.append(dict(config_type=kind,declared_xc=xc,frames=len(r),atoms=sum(s['atoms'] for s in r),
            component_rmse_eV_A=float(np.sqrt(sum(s['squared_error'] for s in r)/components)),
            frame_equal_weight_component_rmse_eV_A=float(np.sqrt(np.mean([s['component_rmse_eV_A']**2 for s in r]))),
            reference_component_rms_eV_A=float(np.sqrt(sum(s['reference_component_rms_eV_A']**2*s['components'] for s in r)/components)),
            diamond_referenced_energy_rmse_meV_atom=float(1000*np.sqrt(np.mean(np.square(de)))) if de else None))
    report=dict(status='completed new pretrained ML energy/force evaluations; material and clock uncalibrated',
        model='MACE-MP-0b3-medium',model_sha256=MODEL_SHA256,model_URL=MODEL_URL,license='MIT',
        model_training='MPTrj PBE+U; target PW91/PBE/unspecified groups kept separate',
        archive_sha256=ARCHIVE_SHA256,member_sha256=DATA_SHA256,source_DOI='10.17863/CAM.65004',
        source_role='GAP training data; MACE training overlap not independently audited',
        frames=len(rows),atoms=int(offsets[-1]),calculator_evaluations=len(rows)+1,
        dtype='float64',device='CPU',threads=2,
        versions={p:version(p) for p in ['mace-torch','torch','e3nn','ase','numpy','scipy','matscipy']},
        anchor=dict(frame=anchor,rule='lowest E/N among explicitly PW91 diamond source frames',
                    relaxed_bulk=False,predicted_eV_atom=anchor_prediction,dft_eV_atom=anchor_dft),
        basis_changes=transforms,group_results=summaries,elapsed_reported_s=time.perf_counter()-started,
        new_DFT_runs=0,new_MD_runs=0,fitting=False,material_calibrated=False,kinetic_calibrated=False)
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('COMPLETE',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
