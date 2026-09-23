"""Independent fixed-format source read and replay of every saved Si/O force."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np


def read_source_fixed_columns(path):
    lines=path.read_text(encoding='utf-8').splitlines();frames=[];i=0
    while i<len(lines):
        if lines[i].strip()!='BEGIN_CFG':i+=1;continue
        if lines[i+1].strip()!='Size' or lines[i+3].strip()!='Supercell':raise ValueError('source format changed')
        n=int(lines[i+2]);cell=np.loadtxt(lines[i+4:i+7])
        expected='AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz'
        if ' '.join(lines[i+7].split())!=' '.join(expected.split()):raise ValueError('source columns changed')
        values=np.loadtxt(lines[i+8:i+8+n],ndmin=2)
        k=i+8+n
        if lines[k].strip()!='Energy':raise ValueError('source energy not found')
        frames.append((cell,values,float(lines[k+1])))
        i=k+2
    return frames


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    frames=read_source_fixed_columns(args.source)
    stored=list(csv.DictReader((args.result/'frames.csv').open(encoding='utf-8')))
    summary=json.loads((args.result/'summary.json').read_text(encoding='utf-8'))
    if len(frames)!=1159 or len(stored)!=len(frames) or not summary['complete']:raise ValueError('incomplete result')
    groups={};worst=[];source_max=0.;metric_max=0.;total_atoms=0
    for index,((cell,values,energy),row) in enumerate(zip(frames,stored)):
        with np.load(args.result/'raw'/f'{index:04d}.npz') as d:
            source_error=max(float(np.max(abs(d['cell']-cell))),
                float(np.max(abs(d['positions']-values[:,2:5]))),
                float(np.max(abs(d['source_forces']-values[:,5:8]))),abs(float(d['source_energy'])-energy))
            source_max=max(source_max,source_error)
            if not np.array_equal(d['numbers'],np.where(values[:,1]==0,14,8)):raise ValueError('element mapping changed')
            error=d['predicted_forces']-values[:,5:8]
            rms=np.sqrt(np.mean(error**2));mae=np.mean(abs(error));total_atoms+=len(values)
            metric_max=max(metric_max,abs(rms-float(row['force_component_RMSE_eV_A'])),
                abs(mae-float(row['force_component_MAE_eV_A'])))
            source_rms=float(np.sqrt(np.mean(values[:,5:8]**2)))
            force_bin='0_to_1' if source_rms<1 else '1_to_5' if source_rms<5 else '5_to_10' if source_rms<10 else '10_or_more'
            for key in (row['group'],row['group']+'__DFT_RMS_'+force_bin):
                g=groups.setdefault(key,dict(frames=0,atoms=0,sum_sq=0.,sum_abs=0.))
                g['frames']+=1;g['atoms']+=len(values);g['sum_sq']+=float(np.sum(error**2));g['sum_abs']+=float(np.sum(abs(error)))
            worst.append(dict(frame=index,composition=row['composition'],atoms=len(values),
                source_force_RMS_eV_A=source_rms,force_RMSE_eV_A=float(rms)))
    group_rows=[]
    for name,g in groups.items():
        group_rows.append(dict(group=name,frames=g['frames'],atoms=g['atoms'],
            force_component_RMSE_eV_A=float(np.sqrt(g['sum_sq']/(3*g['atoms']))),
            force_component_MAE_eV_A=g['sum_abs']/(3*g['atoms'])))
    lookup={r['group']:r for r in group_rows}
    group_max=max(abs(r['force_component_RMSE_eV_A']-lookup[r['group']]['force_component_RMSE_eV_A']) for r in summary['groups'])
    if source_max!=0 or metric_max>1e-13 or group_max>1e-13:raise RuntimeError('source or metric replay failed')
    with (args.output/'force_scale_stratification.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(group_rows[0]));w.writeheader();w.writerows(group_rows)
    result=dict(frames_verified=len(frames),atoms_verified=total_atoms,source_max_absolute_difference=source_max,
        metric_max_absolute_difference=metric_max,group_max_absolute_difference=group_max,
        source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        source_reader='independent fixed-column reader; all original raw arrays compared exactly',
        worst_frames=sorted(worst,key=lambda r:r['force_RMSE_eV_A'],reverse=True)[:10],
        force_bins_status='diagnostic stratification; all source configurations retained, no pass threshold inferred',
        material_approved=False,actual_new_potential_evaluations=0)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
