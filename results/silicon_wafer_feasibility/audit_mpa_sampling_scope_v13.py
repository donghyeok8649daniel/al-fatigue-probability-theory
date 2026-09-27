"""Pre-inference coverage audit of the frozen 199-frame comparison schedule.

Reuses only source CSVs and hashes of existing raw files. Descriptive baseline
RMSE scores are reconstructed from per-frame summaries, not new model calls.
Systematic source-order selection is not a random or independent sample.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.replay_mpa_comparison_v12 import expected_selection
from results.silicon_wafer_feasibility.compare_mpa_model_v12 import select


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def score(rows):
    atoms=sum(int(r['atoms']) for r in rows)
    return dict(frames=len(rows),atoms=atoms,compositions=len({r['composition'] for r in rows}),
        baseline_component_RMSE_eV_A=float(np.sqrt(sum(int(r['atoms'])*float(r['force_component_RMSE_eV_A'])**2 for r in rows)/atoms)))


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    planned=[];groups=[];datasets=[];inputs=[];raw_manifest=[]
    for dataset,folder,limit in [('CP2K',args.cp2k,12),('QE',args.qe,16)]:
        rows=list(csv.DictReader((folder/'frames.csv').open(encoding='utf-8')))
        expected=expected_selection(folder,dataset,limit)
        inference=[{k:v for k,v in r.items() if k!='raw_path'} for r in select(folder,limit,dataset)]
        if inference!=expected:raise ValueError('independent and inference selection disagree')
        ids=[r['source_frame'] for r in expected]
        if len(set(ids))!=len(ids):raise ValueError('duplicate selected source index')
        included=set(ids);kept=[r for r in rows if int(r.get('source_frame',r.get('frame'))) in included]
        if len(kept)!=len(expected):raise ValueError('missing selected source frame')
        before=score(rows);after=score(kept)
        datasets.append(dict(dataset=dataset,full=before,selected=after,
            selected_frame_fraction=after['frames']/before['frames'],selected_atom_fraction=after['atoms']/before['atoms']))
        for group in sorted({r.get('observed_group',r.get('group')) for r in rows}):
            all_group=[r for r in rows if r.get('observed_group',r.get('group'))==group]
            kept_group=[r for r in kept if r.get('observed_group',r.get('group'))==group]
            full_stats=score(all_group);selected_stats=score(kept_group)
            groups.append(dict(dataset=dataset,group=group,
                **{'full_'+k:v for k,v in full_stats.items()},**{'selected_'+k:v for k,v in selected_stats.items()},
                full_dataset_atom_weight=full_stats['atoms']/before['atoms'],
                selected_dataset_atom_weight=selected_stats['atoms']/after['atoms']))
        planned.extend(expected)
        inputs.append(dict(dataset=dataset,frames_csv_sha256=sha(folder/'frames.csv'),source_protocol_sha256=sha(folder/'protocol.json')))
        for r in expected:
            raw_manifest.append(dict(**r,raw_sha256=sha(folder/'raw'/f"{r['source_frame']:04d}.npz")))
    result=dict(complete=True,planned_frames=len(planned),independent_selection_matches_runner=True,
        datasets=datasets,groups=groups,source_tables=inputs,
        selection_sha256=hashlib.sha256(json.dumps(planned,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        selection_fixed_before_candidate_predictions=True,random_sample=False,
        statistical_independence_certified=False,pretraining_overlap_audited=False,
        raw_model_evaluations=0,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None,
        interpretation='CSV baseline scores and sampling coverage only; selected-frame comparisons do not estimate a wafer-population error or an unbiased full-corpus error',
        runner_sha256=sha(Path(__file__)),selector_sources={name:sha(Path(__file__).parent/name) for name in (
            'compare_mpa_model_v12.py','replay_mpa_comparison_v12.py')})
    args.output.mkdir(parents=True)
    for name,value in [('summary.json',result),('planned_selection.json',planned),('source_raw_manifest.json',raw_manifest)]:
        (args.output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (args.output/'groups.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(groups[0]));w.writeheader();w.writerows(groups)
    print(json.dumps(dict(planned_frames=len(planned),datasets=datasets,groups=groups),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('cp2k','qe','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
