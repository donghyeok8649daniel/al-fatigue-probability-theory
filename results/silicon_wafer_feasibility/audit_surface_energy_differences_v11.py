"""Composition-preserving energy-difference audit of completed public snapshots.

Absolute energies from different DFT/model reference conventions are not fitted.
Pair differences cancel one arbitrary constant per exact composition; they are
not transition paths, statistically independent pairs or nucleation barriers.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np


def write_csv(path, records):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)


def main(args):
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    summary=json.loads((args.run/'summary.json').read_text(encoding='utf-8'))
    if not summary['complete']:
        raise ValueError('complete source run required')
    with (args.run/'frames.csv').open(encoding='utf-8',newline='') as stream:
        source_rows=list(csv.DictReader(stream))
    groups={}; relative_replay_error=0.; anchors={}; inventory=[]
    for row in source_rows:
        index=int(row['source_frame']); path=args.run/'raw'/f'{index:04d}.npz'
        with np.load(path) as data:
            composition=tuple(int(np.count_nonzero(data['numbers']==z)) for z in (1,8,14))
            n=len(data['numbers']); dft=float(data['source_energy']); predicted=float(data['predicted_energy'])
        if n!=int(row['atoms']) or abs(dft-float(row['source_energy_eV']))>1e-10 or abs(predicted-float(row['mace_energy_eV']))>1e-10:
            raise ValueError('energy raw replay failed')
        anchor=anchors.setdefault(composition,(index,dft,predicted))
        relative=((predicted-anchor[2])-(dft-anchor[1]))/n
        if anchor[0]!=int(row['same_composition_reference_frame']):
            raise ValueError('reference convention changed')
        relative_replay_error=max(relative_replay_error,abs(relative-float(row['relative_energy_error_eV_atom'])))
        key=(row['observed_group'],composition)
        groups.setdefault(key,[]).append((index,n,dft,predicted))
        inventory.append(dict(source_frame=index,raw_file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    comparisons=[]; group_stats={}; worst=[]; max_identity_error=0.
    for (group,composition),values in sorted(groups.items()):
        indices=np.array([v[0] for v in values]); n=values[0][1]; m=len(values)
        dft=np.array([v[2] for v in values]); predicted=np.array([v[3] for v in values])
        left,right=np.triu_indices(m,1); de=(dft[right]-dft[left])/n; pe=(predicted[right]-predicted[left])/n
        error=pe-de; residual=(predicted-dft)/n
        # Sum_{i<j}(r_i-r_j)^2 = m*sum_i(r_i-mean(r))^2.
        pair_sq=float(error@error); identity=float(m*np.sum((residual-residual.mean())**2))
        identity_error=abs(pair_sq-identity)/max(pair_sq,1.)
        max_identity_error=max(max_identity_error,identity_error)
        row=dict(group=group,n_H=composition[0],n_O=composition[1],n_Si=composition[2],atoms=n,
            frames=m,pairs=len(error),pair_energy_error_RMSE_eV_atom=float(np.sqrt(np.mean(error**2))) if len(error) else None,
            maximum_pair_energy_error_eV_atom=float(np.max(abs(error))) if len(error) else None,
            source_energy_range_eV_atom=float(np.ptp(dft)/n),mace_energy_range_eV_atom=float(np.ptp(predicted)/n),
            pair_variance_identity_relative_error=identity_error)
        comparisons.append(row)
        stats=group_stats.setdefault(group,dict(frames=0,compositions=0,singleton_compositions=0,pairs=0,sq=0.,absolute=0.,
            ordering={str(t):dict(resolved_pairs=0,opposite_sign=0,prediction_tied=0) for t in (0.,1e-6,1e-4)}))
        stats['frames']+=m; stats['compositions']+=1; stats['singleton_compositions']+=int(m==1)
        stats['pairs']+=len(error); stats['sq']+=pair_sq; stats['absolute']+=float(abs(error).sum())
        for threshold in (0.,1e-6,1e-4):
            resolved=abs(de)>threshold; information=stats['ordering'][str(threshold)]
            information['resolved_pairs']+=int(np.sum(resolved))
            information['opposite_sign']+=int(np.sum(resolved&(de*pe<0)))
            information['prediction_tied']+=int(np.sum(resolved&(pe==0)))
        if len(error):
            j=int(np.argmax(abs(error)))
            worst.append(dict(group=group,source_frame_left=int(indices[left[j]]),source_frame_right=int(indices[right[j]]),
                atoms=n,source_delta_eV_atom=float(de[j]),mace_delta_eV_atom=float(pe[j]),error_eV_atom=float(error[j])))
    if relative_replay_error>1e-10 or max_identity_error>1e-10:
        raise ValueError('energy-difference replay/identity failed')
    aggregate=[]
    for group,stats in sorted(group_stats.items()):
        count=stats['pairs']
        aggregate.append(dict(group=group,frames=stats['frames'],compositions=stats['compositions'],
            singleton_compositions=stats['singleton_compositions'],pairs=count,
            pair_energy_error_RMSE_eV_atom=float(np.sqrt(stats['sq']/count)) if count else None,
            pair_energy_error_MAE_eV_atom=stats['absolute']/count if count else None,
            ordering=stats['ordering']))
    write_csv(args.output/'composition_groups.csv',comparisons)
    write_csv(args.output/'maximum_error_pair_per_group.csv',sorted(worst,key=lambda r:abs(r['error_eV_atom']),reverse=True))
    write_csv(args.output/'source_raw_manifest.csv',inventory)
    result=dict(frames=len(source_rows),composition_and_observed_groups=len(comparisons),
        maximum_stored_reference_replay_error_eV_atom=relative_replay_error,
        maximum_pair_variance_identity_relative_error=max_identity_error,groups=aggregate,
        weighting='equal pair weights within each observed group; each pair has identical H/O/Si atom counts; pair counts are not independent samples',
        grouping='exact H/O/Si counts plus observed composition/vacuum group; cell and geometry may vary; not a continuous reaction path',
        energy_offsets_fitted=False,fit_performed=False,new_potential_calls=0,new_DFT=0,new_MD=0,
        initiation_barrier_eV=None,material_approved=False)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True); main(parser.parse_args())
