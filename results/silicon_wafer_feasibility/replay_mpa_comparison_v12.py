"""Recompute the two-checkpoint comparison from arrays without MACE calls.

The source-selection implementation and pair-energy loop are independent of
the inference runner. This is a serialization/analysis audit, not an independent
DFT or model-accuracy experiment. No statistical independence is assumed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def expected_selection(path,dataset,limit):
    grouped={}
    for row in csv.DictReader((path/'frames.csv').open(encoding='utf-8')):
        label=row.get('observed_group',row.get('group'))
        key=(label,row['composition']) if dataset=='CP2K' else (label,)
        grouped.setdefault(key,[]).append(row)
    result=[]
    for key,rows in sorted(grouped.items()):
        count=min(limit,len(rows))
        offsets=[0] if count==1 else [(len(rows)-1)*i//(count-1) for i in range(count)]
        for offset in offsets:
            row=rows[offset]
            result.append(dict(dataset=dataset,group=row.get('observed_group',row.get('group')),
                composition=row['composition'],source_frame=int(row.get('source_frame',row.get('frame')))))
    return result


def pair_energy(rows):
    result=[]
    for composition in sorted({r['composition'] for r in rows}):
        selected=[r for r in rows if r['composition']==composition]
        if len(selected)<2:continue
        if len({r['atoms'] for r in selected})!=1:raise ValueError('composition atom counts differ')
        squared_before=squared_after=0.;resolved=reverse_before=reverse_after=number=0
        for i,a in enumerate(selected):
            for b in selected[i+1:]:
                n=a['atoms'];truth=(b['source_energy']-a['source_energy'])/n
                before=(b['baseline_energy']-a['baseline_energy'])/n
                after=(b['predicted_energy']-a['predicted_energy'])/n
                squared_before+=(before-truth)**2;squared_after+=(after-truth)**2;number+=1
                if abs(truth)>1e-4:
                    resolved+=1;reverse_before+=truth*before<0;reverse_after+=truth*after<0
        result.append(dict(composition=composition,frames=len(selected),pairs=number,
            baseline_pair_energy_RMSE_eV_atom=float(np.sqrt(squared_before/number)),
            candidate_pair_energy_RMSE_eV_atom=float(np.sqrt(squared_after/number)),
            resolved_pair_count=resolved,baseline_reversed_resolved_pairs=int(reverse_before),
            candidate_reversed_resolved_pairs=int(reverse_after)))
    return result


def main(args):
    if args.output.exists():raise ValueError('fresh audit output required')
    summary=read(args.result/'summary.json');protocol=read(args.result/'protocol.json')
    if not summary['complete']:raise ValueError('inference comparison is incomplete')
    expected=expected_selection(args.cp2k,'CP2K',12)+expected_selection(args.qe,'QE',16)
    if protocol['source_selection']!=expected:raise ValueError('predeclared source selection differs')
    frames=list(csv.DictReader((args.result/'frames.csv').open(encoding='utf-8')))
    manifest=list(csv.DictReader((args.result/'source_raw_manifest.csv').open(encoding='utf-8')))
    if len(frames)!=len(expected) or len(manifest)!=len(expected):raise ValueError('record count differs')
    rows=[];stats={};element_stats={};max_error=0.;source_hashes=[]
    for index,(entry,frame,pin) in enumerate(zip(expected,frames,manifest)):
        for key,value in entry.items():
            if str(frame[key])!=str(value):raise ValueError('frame metadata differs')
        source=(args.cp2k if entry['dataset']=='CP2K' else args.qe)/'raw'/f"{entry['source_frame']:04d}.npz"
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        if (pin['dataset']!=entry['dataset'] or int(pin['source_frame'])!=entry['source_frame']
                or pin['baseline_raw_sha256']!=digest):raise ValueError('source manifest mismatch')
        candidate=args.result/'raw'/f'{index:04d}.npz'
        with np.load(source) as original,np.load(candidate) as current:
            for key in ('numbers','positions','cell','source_energy','source_forces'):
                if not np.array_equal(original[key],current[key]):raise ValueError('source array changed: '+key)
            if not np.array_equal(original['pbc'] if 'pbc' in original else np.ones(3,bool),current['pbc']):
                raise ValueError('periodicity changed')
            for old,new in [('predicted_energy','baseline_energy'),('predicted_forces','baseline_forces')]:
                if not np.array_equal(original[old],current[new]):raise ValueError('baseline array changed')
            before=current['baseline_forces']-current['source_forces']
            after=current['predicted_forces']-current['source_forces']
            if before.shape!=after.shape or not np.all(np.isfinite(after)):raise ValueError('invalid predictions')
            n=len(current['numbers']);b=float(np.sqrt(np.mean(before**2)));a=float(np.sqrt(np.mean(after**2)))
            max_error=max(max_error,abs(b-float(frame['baseline_force_RMSE_eV_A'])),abs(a-float(frame['candidate_force_RMSE_eV_A'])))
            row=dict(**entry,atoms=n,baseline_frame_RMSE=b,candidate_frame_RMSE=a)
            for key,column in [('source_energy','reference_energy_eV'),('baseline_energy','baseline_energy_eV'),('predicted_energy','candidate_energy_eV')]:
                row[key]=float(current[key]);max_error=max(max_error,abs(row[key]-float(frame[column])))
            rows.append(row);label=entry['dataset']+'__'+entry['group']
            stats.setdefault(label,[]).append((n,float(np.sum(before**2)),float(np.sum(after**2)),b,a))
            for z in np.unique(current['numbers']):
                mask=current['numbers']==z
                element_stats.setdefault(label+'__Z'+str(z),[]).append((int(mask.sum()),float(np.sum(before[mask]**2)),float(np.sum(after[mask]**2))))
        source_hashes.append(dict(dataset=entry['dataset'],source_frame=entry['source_frame'],baseline_raw_sha256=digest,
            candidate_raw_sha256=hashlib.sha256(candidate.read_bytes()).hexdigest()))
    groups=[]
    for label,values in stats.items():
        v=np.array(values);atoms=int(v[:,0].sum());b=float(np.sqrt(v[:,1].sum()/(3*atoms)));a=float(np.sqrt(v[:,2].sum()/(3*atoms)))
        group=dict(group=label,frames=len(values),atoms=atoms,baseline_force_RMSE_eV_A=b,candidate_force_RMSE_eV_A=a)
        stored=next(r for r in summary['groups'] if r['group']==label)
        for key in ('frames','atoms','baseline_force_RMSE_eV_A','candidate_force_RMSE_eV_A'):
            max_error=max(max_error,abs(group[key]-stored[key]))
        group.update(frames_improved=int(np.sum(v[:,4]<v[:,3])),frames_worsened=int(np.sum(v[:,4]>v[:,3])),
            median_frame_RMSE_change_eV_A=float(np.median(v[:,4]-v[:,3])),
            worst_frame_RMSE_change_eV_A=float(np.max(v[:,4]-v[:,3])))
        groups.append(group)
    elements=[]
    for label,values in element_stats.items():
        v=np.array(values);atoms=int(v[:,0].sum())
        elements.append(dict(group_element=label,atoms=atoms,baseline_RMSE_eV_A=float(np.sqrt(v[:,1].sum()/(3*atoms))),
            candidate_RMSE_eV_A=float(np.sqrt(v[:,2].sum()/(3*atoms)))))
    pairs={dataset:pair_energy([r for r in rows if r['dataset']==dataset]) for dataset in ('CP2K','QE')}
    for dataset,values in pairs.items():
        if len(values)!=len(summary['relative_energy_pairs'][dataset]):raise ValueError('pair groups differ')
        for row,stored in zip(values,summary['relative_energy_pairs'][dataset]):
            if row['composition']!=stored['composition']:raise ValueError('pair composition order differs')
            for key in row:
                if key!='composition':max_error=max(max_error,abs(row[key]-stored[key]))
    if max_error>1e-10:raise ValueError('stored numerical analysis does not replay: '+str(max_error))
    controls=summary['standard_controls']
    if len(controls)!=9 or max(abs(r['energy_difference_eV']) for r in controls)>1e-8 or max(r['force_difference_eV_A'] for r in controls)>1e-8:
        raise ValueError('standard evaluator controls incomplete or failed')
    args.output.mkdir(parents=True)
    result=dict(complete=True,frames=len(rows),source_array_checks_exact=True,maximum_scalar_replay_error=max_error,
        groups=groups,elements=elements,relative_energy_pairs=pairs,
        caveat='paired descriptive scores; correlated source frames are not independent confidence intervals; foundation pretraining overlap unknown',
        new_model_calls=0,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (args.output/'raw_manifest.json').write_text(json.dumps(source_hashes,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('cp2k','qe','result','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
