"""Apply the already fitted conservative delta candidates to larger structures.

No new fit or MACE call. All available frames from the predeclared large-frame
audit are used, including unchanged pure-Si controls. No error-based exclusion.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_oxide_angular_delta_v12 import combined_features
from results.silicon_wafer_feasibility.fit_oxide_radial_delta_v12 import energy_differences


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(args):
    from ase import Atoms
    from ase.neighborlist import neighbor_list
    if args.output.exists():raise ValueError('fresh candidate evaluation output required')
    baseline=json.loads((args.source/'summary.json').read_text(encoding='utf-8'))
    if not baseline['complete']:raise ValueError('large-frame audit incomplete')
    candidates={};pins=[]
    for family in ('radial','angular'):
        path=args.results/('oxide_'+family+'_delta')/'cross_validation.json'
        candidates[family]=np.array(json.loads(path.read_text(encoding='utf-8'))['all_training_fit']['coefficients_eV'])
        pins.append(dict(family=family,fit_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    args.output.mkdir(parents=True);started=time.perf_counter();rows=[];stats={};environment={};sources=[]
    metadata=list(csv.DictReader((args.source/'frames.csv').open(encoding='utf-8')))
    for number,record in enumerate(metadata):
        index=int(record['source_frame']);source=args.source/'raw'/f'{index:04d}.npz'
        with np.load(source) as d:
            z=d['numbers'];positions=d['positions'];cell=d['cell'];pbc=d['pbc']
            truth=d['source_forces'];old=d['predicted_forces'];energy=float(d['predicted_energy']);reference=float(d['source_energy'])
        e_basis,f_basis=combined_features(z,positions,cell,pbc)
        atoms=Atoms(numbers=z,positions=positions,cell=cell,pbc=pbc)
        i,j=neighbor_list('ij',atoms,2.0,self_interaction=False)
        oxygen_count=np.bincount(i[z[j]==8],minlength=len(z))
        labels=np.array(['H' if a==1 else 'O' if a==8 else 'Si_O0' if oxygen_count[k]==0 else 'Si_O1to3' if oxygen_count[k]<=3 else 'Si_O4plus' for k,a in enumerate(z)])
        for family,coefficients in candidates.items():
            count=len(coefficients);delta_e=float(e_basis[:count]@coefficients)
            delta_f=f_basis[:,:,:count]@coefficients
            if not np.any(z==8) and (delta_e!=0 or np.any(delta_f!=0)):raise ValueError('pure Si control was changed')
            before=old-truth;after=old+delta_f-truth;group=record['observed_group']
            row=dict(family=family,source_frame=index,composition=record['composition'],atoms=len(z),group=group,
                reference_energy_eV=reference,baseline_energy_eV=energy,candidate_energy_eV=energy+delta_e,
                baseline_RMSE_eV_A=float(np.sqrt(np.mean(before**2))),candidate_RMSE_eV_A=float(np.sqrt(np.mean(after**2))),
                correction_total_force_norm_eV_A=float(np.linalg.norm(delta_f.sum(axis=0))))
            rows.append(row)
            stats.setdefault((family,group),[]).append((len(z),float(np.sum(before**2)),float(np.sum(after**2))))
            for label in np.unique(labels):
                mask=labels==label
                environment.setdefault((family,group,str(label)),[]).append((int(mask.sum()),float(np.sum(before[mask]**2)),float(np.sum(after[mask]**2))))
        sources.append(dict(source_frame=index,sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
        if number%20==0:dump(args.output/'running.json',dict(completed=number+1,planned=len(metadata),elapsed_s=time.perf_counter()-started))
    def aggregate(values):
        array=np.array(values);atoms=int(array[:,0].sum())
        return dict(frames=len(values),atoms=atoms,baseline_RMSE_eV_A=float(np.sqrt(array[:,1].sum()/(3*atoms))),
            candidate_RMSE_eV_A=float(np.sqrt(array[:,2].sum()/(3*atoms))))
    groups=[dict(family=key[0],group=key[1],**aggregate(values)) for key,values in stats.items()]
    environments=[dict(family=key[0],group=key[1],geometric_label=key[2],**aggregate(values)) for key,values in environment.items()]
    pairs={family:energy_differences([row for row in rows if row['family']==family]) for family in candidates}
    for filename,data in [('frames.csv',rows),('groups.csv',groups),('environment_errors.csv',environments),('source_raw_manifest.csv',sources)]:
        with (args.output/filename).open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    report=dict(complete=True,frames=len(metadata),fits=pins,groups=groups,environment_errors=environments,
        relative_energy_pairs=pairs,elapsed_s=time.perf_counter()-started,
        selection='all completed 320<N<=1200 source frames, excluded from delta coefficient fits',
        caveat='same public DFT collection, correlated configurations, original foundation pretraining overlap unknown; not independent material certification',
        environment_definition='Si counts O neighbors within2.0A; geometric label, not formal oxidation state or charge',
        force_conservation_max=float(max(row['correction_total_force_norm_eV_A'] for row in rows)),
        new_fits=0,new_MACE=0,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None)
    dump(args.output/'summary.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','source','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
