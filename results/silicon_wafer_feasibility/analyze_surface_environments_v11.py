"""Resolve saved force errors by local chemical environment, without refitting.

Coordination is a cutoff-dependent geometric diagnostic, not an oxidation-state
measurement, electronic charge, crack classifier or complete interface label.
Every included atom contributes; large errors are never removed from totals.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np


def main(args):
    from ase import Atoms
    from ase.neighborlist import neighbor_list
    summary=json.loads((args.result/'summary.json').read_text(encoding='utf-8'))
    if not summary['complete']:raise ValueError('completed result required')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    inventory={int(r['source_frame']):r for r in csv.DictReader((args.result/'source_inventory.csv').open(encoding='utf-8'))}
    table={};frame_rows=[];moment_rows=[];atoms_seen=0;sum_sq_total=0.;sum_sq_partition=0.;projection_residual=0.
    # Three diagnostic radius sets are specified before inspecting force errors.
    settings=[(1.8,2.6,1.7),(2.0,2.8,1.9),(2.2,3.0,2.1)]
    for path in sorted((args.result/'raw').glob('*.npz')):
        index=int(path.stem);inv=inventory[index]
        with np.load(path) as data:
            numbers=data['numbers'];source_force=data['source_forces'];predicted_force=data['predicted_forces']
            error=predicted_force-source_force
            atoms=Atoms(numbers=numbers,positions=data['positions'],cell=data['cell'],pbc=data['pbc'])
        n=len(atoms);atoms_seen+=n;sum_sq_total+=float(np.sum(error**2))
        mean_error=error.mean(axis=0);centered_error=error-mean_error
        raw_mse=float(np.mean(error**2));uniform_mse=float(np.mean(mean_error**2));centered_mse=float(np.mean(centered_error**2))
        projection_residual=max(projection_residual,abs(raw_mse-uniform_mse-centered_mse)/max(raw_mse,1e-30))
        source_resultant=float(np.linalg.norm(source_force.sum(axis=0)))
        moment_rows.append(dict(source_frame=index,atoms=n,observed_group=inv['observed_group'],
            source_total_force_eV_A=source_resultant,predicted_total_force_eV_A=float(np.linalg.norm(predicted_force.sum(axis=0))),
            raw_force_component_RMSE_eV_A=float(np.sqrt(raw_mse)),
            centered_error_RMSE_eV_A=float(np.sqrt(centered_mse)),uniform_error_RMS_eV_A=float(np.sqrt(uniform_mse)),
            zero_total_force_model_error_lower_bound_eV_A=source_resultant/(np.sqrt(3.)*n),
            uniform_fraction_of_squared_error=uniform_mse/raw_mse if raw_mse else None))
        first,second,dist=neighbor_list('ijd',atoms,cutoff=3.,self_interaction=False)
        element_name={1:'H',8:'O',14:'Si'}
        for setting,(r_o,r_si,r_h) in enumerate(settings):
            counts={}
            for z,radius in ((8,r_o),(14,r_si),(1,r_h)):
                mask=(numbers[second]==z)&(dist<radius)
                counts[z]=np.bincount(first[mask],minlength=n)
            # Only classify Si chemically; retain all other atoms by element.
            labels=[]
            for atom,z in enumerate(numbers):
                if z==14:
                    label=f'Si_O{counts[8][atom]}_Si{counts[14][atom]}_H{counts[1][atom]}'
                else:label=element_name[int(z)]
                labels.append(label)
            labels=np.array(labels)
            for label in np.unique(labels):
                mask=labels==label;err=error[mask];key=(setting,inv['observed_group'],str(label))
                g=table.setdefault(key,dict(frames=0,atoms=0,sq=0.,absolute=0.,max_error=0.))
                g['frames']+=1;g['atoms']+=int(mask.sum());g['sq']+=float(np.sum(err**2));g['absolute']+=float(np.sum(abs(err)))
                g['max_error']=max(g['max_error'],float(np.max(abs(err))))
                if setting==1:sum_sq_partition+=float(np.sum(err**2))
            if setting==1:
                si=numbers==14
                frame_rows.append(dict(source_frame=index,atoms=n,observed_group=inv['observed_group'],
                    Si_without_O_neighbors=int(np.sum(si&(counts[8]==0))),
                    Si_with_1_to_3_O_neighbors=int(np.sum(si&(counts[8]>=1)&(counts[8]<=3))),
                    Si_with_4_or_more_O_neighbors=int(np.sum(si&(counts[8]>=4))),
                    Si_with_H_neighbors=int(np.sum(si&(counts[1]>0)))))
    rows=[]
    for (setting,group,label),g in table.items():
        rows.append(dict(radius_set=setting,group=group,coordination_label=label,
            Si_O_radius_A=settings[setting][0],Si_Si_radius_A=settings[setting][1],Si_H_radius_A=settings[setting][2],
            frames=g['frames'],atoms=g['atoms'],force_component_RMSE_eV_A=np.sqrt(g['sq']/(3*g['atoms'])),
            force_component_MAE_eV_A=g['absolute']/(3*g['atoms']),force_component_max_error_eV_A=g['max_error']))
    for filename,records in [('environment_force_errors.csv',rows),('frame_environment_inventory.csv',frame_rows),
                              ('net_force_projection_diagnostic.csv',moment_rows)]:
        with (args.output/filename).open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    if atoms_seen!=summary['atoms'] or not np.isclose(sum_sq_partition,sum_sq_total,rtol=1e-12,atol=0) or projection_residual>1e-12:
        raise ValueError('environment partition lost atoms or force errors')
    result=dict(frames=len(frame_rows),atoms=atoms_seen,neighbor_radius_sets_A=settings,
        total_force_squared_error=sum_sq_total,partition_squared_error=sum_sq_partition,
        relative_partition_residual=(sum_sq_partition-sum_sq_total)/sum_sq_total,
        maximum_relative_force_projection_identity_residual=projection_residual,
        largest_source_total_force_frame=max(moment_rows,key=lambda r:r['source_total_force_eV_A']),
        source_net_force_status='measured diagnostic; origin of nonzero source resultant not established; original forces/errors preserved',
        descriptor_status='cutoff-dependent geometric coordination; not electronic oxidation state or unique interface label',
        force_error_selection_performed=False,new_potential_calls=0,new_DFT=0,new_MD=0,
        material_approved=False,crack_initiation_certified=False)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
