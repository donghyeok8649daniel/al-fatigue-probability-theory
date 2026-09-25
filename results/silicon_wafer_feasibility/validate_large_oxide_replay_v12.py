"""Independent text replay of every included CP2K reference state.

No ASE source parser, force model call, fitted offset or omitted outlier.
The original size-based inclusion rule and every excluded id are verified.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,shlex
from pathlib import Path
import numpy as np

SOURCE_SHA='d99170bff4942dcd5a1c38d998fe0a60084a5e132bcb3cf2335f5db4cbe17689'


def frames(path):
    symbols={'H':1,'O':8,'Si':14}
    with path.open(encoding='utf-8') as stream:
        while line:=stream.readline():
            if not line.strip():continue
            n=int(line);meta=dict(word.split('=',1) for word in shlex.split(stream.readline()))
            if meta['Properties']!='species:S:1:pos:R:3:forces:R:3':raise ValueError('source columns changed')
            body=[stream.readline().split() for _ in range(n)]
            if any(len(row)!=7 for row in body):raise ValueError('truncated or extra source columns')
            numbers=np.array([symbols[row[0]] for row in body])
            values=np.array([[float(x) for x in row[1:]] for row in body])
            yield dict(numbers=numbers,positions=values[:,:3],source_forces=values[:,3:],
                cell=np.fromstring(meta['Lattice'],sep=' ').reshape(3,3),
                pbc=np.array([s=='T' for s in meta['pbc'].split()]),
                source_energy=np.array(float(meta['energy'])),config_type=meta['config_type'])


def observed_group(d):
    zs=set(d['numbers']);chemical='_'.join(label for z,label in ((14,'Si'),(8,'O'),(1,'H')) if z in zs)
    if len(d['numbers'])<=2:return 'small_'+chemical
    # Reciprocal-normal cell widths: independent of which direct vectors are skew.
    reciprocal=np.linalg.inv(d['cell']);fractional=d['positions']@reciprocal
    empty=[]
    for axis in range(3):
        coords=np.sort(np.mod(fractional[:,axis],1.))
        gaps=np.r_[np.diff(coords),1+coords[0]-coords[-1]]
        empty.append(float(gaps.max()/np.linalg.norm(reciprocal[:,axis])))
    return chemical+'_vacuum'+str(sum(v>8. for v in empty))


def main(args):
    if hashlib.sha256(args.source.read_bytes()).hexdigest()!=SOURCE_SHA:raise ValueError('source changed')
    summary=json.loads((args.result/'summary.json').read_text(encoding='utf-8'))
    if not summary['complete']:raise ValueError('complete source result required')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    inventory=list(csv.DictReader((args.result/'source_inventory.csv').open(encoding='utf-8')))
    stored={int(r['source_frame']):r for r in csv.DictReader((args.result/'frames.csv').open(encoding='utf-8'))}
    max_atoms=int(summary['protocol']['max_atoms']);min_atoms=int(summary['protocol']['min_atoms_exclusive'])
    seen=[];excluded=[];refs={};groups={};worst=[]
    maximum_source=0.;maximum_metric=0.;total_atoms=0
    def add(name,error):
        g=groups.setdefault(name,dict(frames=0,atoms=0,sq=0.,absolute=0.))
        g['frames']+=1;g['atoms']+=len(error);g['sq']+=float(np.sum(error**2));g['absolute']+=float(np.sum(abs(error)))
    for index,d in enumerate(frames(args.source)):
        inv=inventory[index];n=len(d['numbers'])
        selected=min_atoms<n<=max_atoms and (not summary['protocol'].get('require_oxygen',False) or 8 in d['numbers'])
        if int(inv['source_frame'])!=index or int(inv['atoms'])!=n or (inv['included']=='True')!=selected:
            raise ValueError('source inventory or inclusion mismatch')
        if d['config_type']!=inv['source_config_type'] or observed_group(d)!=inv['observed_group']:
            raise ValueError('source or auxiliary group mismatch')
        if not selected:excluded.append(index);continue
        seen.append(index);total_atoms+=n;row=stored[index]
        with np.load(args.result/'raw'/f'{index:04d}.npz') as raw:
            for key in ('numbers','pbc'):
                if not np.array_equal(d[key],raw[key]):raise ValueError('element or periodicity mismatch')
            for key in ('positions','source_forces','cell','source_energy'):
                maximum_source=max(maximum_source,float(np.max(abs(raw[key]-d[key]))))
            error=raw['predicted_forces']-d['source_forces'];name=inv['observed_group']
            source_rms=float(np.sqrt(np.mean(d['source_forces']**2)))
            expected=dict(force_component_RMSE_eV_A=float(np.sqrt(np.mean(error**2))),
                force_component_MAE_eV_A=float(np.mean(abs(error))),
                force_component_max_error_eV_A=float(np.max(abs(error))),
                source_force_component_RMS_eV_A=source_rms,
                source_total_force_norm_eV_A=float(np.linalg.norm(d['source_forces'].sum(axis=0))),
                mace_total_force_norm_eV_A=float(np.linalg.norm(raw['predicted_forces'].sum(axis=0))))
            composition=tuple(int(np.sum(d['numbers']==z)) for z in (1,8,14))
            refs.setdefault(composition,(index,float(d['source_energy']),float(raw['predicted_energy'])))
            anchor,er,pr=refs[composition]
            expected['relative_energy_error_eV_atom']=((float(raw['predicted_energy'])-pr)-(float(d['source_energy'])-er))/n
            if int(row['same_composition_reference_frame'])!=anchor:raise ValueError('energy reference changed')
            for key,value in expected.items():maximum_metric=max(maximum_metric,abs(value-float(row[key])))
            for z in (1,8,14):
                mask=d['numbers']==z
                if not mask.any():
                    if row[f'Z{z}_force_RMSE_eV_A']!='':raise ValueError('absent element has value')
                    continue
                rms=float(np.sqrt(np.mean(error[mask]**2)))
                maximum_metric=max(maximum_metric,abs(rms-float(row[f'Z{z}_force_RMSE_eV_A'])))
                add(name+f'__Z{z}',error[mask])
            add(name,error)
            force_bin='0_to_1' if source_rms<1 else '1_to_5' if source_rms<5 else '5_to_10' if source_rms<10 else '10_or_more'
            add(name+'__DFT_RMS_'+force_bin,error)
            worst.append(dict(source_frame=index,group=name,atoms=n,source_force_RMS_eV_A=source_rms,
                force_RMSE_eV_A=expected['force_component_RMSE_eV_A']))
    if index+1!=1466 or len(inventory)!=1466 or set(seen)!=set(stored):raise ValueError('source coverage mismatch')
    group_rows=[dict(group=name,frames=g['frames'],atoms=g['atoms'],
        force_component_RMSE_eV_A=float(np.sqrt(g['sq']/(3*g['atoms']))),
        force_component_MAE_eV_A=g['absolute']/(3*g['atoms'])) for name,g in groups.items()]
    lookup={r['group']:r for r in group_rows};maximum_group=0.
    for row in summary['groups']:
        for key in ('frames','atoms','force_component_RMSE_eV_A','force_component_MAE_eV_A'):
            maximum_group=max(maximum_group,abs(row[key]-lookup[row['group']][key]))
    if maximum_source!=0. or maximum_metric>1e-12 or maximum_group>1e-12:raise RuntimeError('independent replay mismatch')
    with (args.output/'group_element_force_scale.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(group_rows[0]));writer.writeheader();writer.writerows(group_rows)
    result=dict(complete=True,source_sha256=SOURCE_SHA,source_frames=index+1,frames_verified=len(seen),
        atoms_verified=total_atoms,excluded_frame_ids=excluded,max_atoms=max_atoms,min_atoms_exclusive=min_atoms,
        require_oxygen=summary['protocol'].get('require_oxygen',False),
        maximum_source_difference=maximum_source,maximum_metric_difference=maximum_metric,
        maximum_group_difference=maximum_group,worst_frames=sorted(worst,key=lambda r:r['force_RMSE_eV_A'],reverse=True)[:12],
        source_reader='independent shlex metadata and fixed atom-column reader; no ASE parser',
        new_potential_calls=0,new_DFT=0,new_MD=0,material_approved=False)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='excluded_frame_ids'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
