"""Composition-held-out diagnostic of a minimal conservative Si-O residual.

Fit only CP2K condensed Si/O(/H) forces using eight fixed radial kernels.
Never fit offsets/relative energies. Test an entire composition at a time and
evaluate the all-training candidate on the separate QE dataset. Differences in
DFT conventions and unknown original MACE training overlap remain explicit.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_oxide_delta_v12 import correction_features,CENTERS_A,WIDTH_A,SUPPORT_A
from solver_v1.silicon_oxide_angular_delta_v12 import combined_features


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def load_frame(path,metadata,basis_model='radial'):
    with np.load(path) as d:
        numbers=d['numbers'];positions=d['positions'];cell=d['cell']
        pbc=d['pbc'] if 'pbc' in d else np.ones(3,bool)
        feature_fn=correction_features if basis_model=='radial' else combined_features
        energy,force=feature_fn(numbers,positions,cell,pbc)
        return dict(metadata=metadata,energy_basis=energy,force_basis=force.reshape(-1,len(energy)),
            reference_force=d['source_forces'].ravel(),baseline_force=d['predicted_forces'].ravel(),
            reference_energy=float(d['source_energy']),baseline_energy=float(d['predicted_energy']),
            atoms=len(numbers),raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def fit(frames):
    # Equal total weight per frame; no hidden hyperparameter search.
    x=np.concatenate([f['force_basis']/np.sqrt(3*f['atoms']) for f in frames])
    y=np.concatenate([(f['reference_force']-f['baseline_force'])/np.sqrt(3*f['atoms']) for f in frames])
    coefficients,residual,rank,singular=np.linalg.lstsq(x,y,rcond=1e-12)
    diagnostics=dict(rank=int(rank),singular_values=singular.tolist(),
        condition_number=float(singular[0]/singular[-1]),coefficients_eV=coefficients.tolist(),
        weighted_force_residual_norm=float(np.linalg.norm(x@coefficients-y)),
        normal_equation_residual_norm=float(np.linalg.norm(x.T@(x@coefficients-y))),
        frame_count=len(frames),force_components=len(y))
    if rank!=x.shape[1]:raise ValueError('correction rank deficient; no regularized candidate accepted')
    return coefficients,diagnostics


def metrics(frames,coefficients):
    baseline=[];candidate=[];per_frame=[]
    for f in frames:
        b=f['baseline_force']-f['reference_force'];c=b+f['force_basis']@coefficients
        baseline.append(b);candidate.append(c)
        per_frame.append(dict(source_frame=int(f['metadata'].get('source_frame',f['metadata'].get('frame'))),
            atoms=f['atoms'],composition=f['metadata']['composition'],
            group=f['metadata'].get('observed_group',f['metadata'].get('group')),
            baseline_force_RMSE_eV_A=float(np.sqrt(np.mean(b*b))),candidate_force_RMSE_eV_A=float(np.sqrt(np.mean(c*c))),
            correction_energy_eV=float(f['energy_basis']@coefficients),
            reference_energy_eV=f['reference_energy'],baseline_energy_eV=f['baseline_energy'],
            candidate_energy_eV=f['baseline_energy']+float(f['energy_basis']@coefficients)))
    b=np.concatenate(baseline);c=np.concatenate(candidate)
    return dict(frames=len(frames),atoms=sum(f['atoms'] for f in frames),baseline_force_RMSE_eV_A=float(np.sqrt(np.mean(b*b))),
        candidate_force_RMSE_eV_A=float(np.sqrt(np.mean(c*c))),
        relative_force_RMSE_change=float(np.sqrt(np.mean(c*c))/np.sqrt(np.mean(b*b))-1) if np.any(b) else None),per_frame


def energy_differences(rows):
    records=[]
    for composition in sorted({r['composition'] for r in rows}):
        group=[r for r in rows if r['composition']==composition]
        if len(group)<2:continue
        if len({r['atoms'] for r in group})!=1:raise ValueError('composition/atom number mismatch')
        n=group[0]['atoms'];ref=np.array([r['reference_energy_eV'] for r in group])/n
        before=np.array([r['baseline_energy_eV'] for r in group])/n
        after=np.array([r['candidate_energy_eV'] for r in group])/n
        i,j=np.triu_indices(len(group),1);truth=ref[j]-ref[i];b=before[j]-before[i];c=after[j]-after[i]
        resolved=abs(truth)>1e-4
        records.append(dict(composition=composition,frames=len(group),pairs=len(i),
            baseline_pair_energy_RMSE_eV_atom=float(np.sqrt(np.mean((b-truth)**2))),
            candidate_pair_energy_RMSE_eV_atom=float(np.sqrt(np.mean((c-truth)**2))),
            resolved_pair_count=int(resolved.sum()),
            baseline_reversed_resolved_pairs=int(np.sum(resolved&(truth*b<0))),
            candidate_reversed_resolved_pairs=int(np.sum(resolved&(truth*c<0)))))
    return records


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True);started=time.perf_counter()
    metadata=list(csv.DictReader((args.cp2k/'frames.csv').open(encoding='utf-8')))
    eligible=[r for r in metadata if r['observed_group'] in ('Si_O_H_vacuum1','Si_O_vacuum0')]
    compositions=sorted({r['composition'] for r in eligible})
    if len(eligible)!=487 or len(compositions)!=5:raise ValueError('predefined CP2K condensed selection changed')
    protocol=dict(selection='all 487 condensed Si/O/H or Si/O CP2K frames already evaluated in v11',
        train_test='leave-one-entire-composition-out, five folds; no random-frame split',
        target='DFT minus baseline MACE forces; equal weight per frame',
        energy_fit=False,energy_offset_fit=False,regularization=False,
        radial_centers_A=CENTERS_A.tolist(),width_A=WIDTH_A,cutoff_A=SUPPORT_A,
        window='quintic rise .6->1.0 A and fall 2.8->3.5 A; energy and first two derivatives vanish at endpoints',
        coefficients=8 if args.basis=='radial' else 17,basis=args.basis,
        angular_definition=None if args.basis=='radial' else 'P0/P1/P2 for O-Si-O, Si-O-Si, Si-Si-O; radial Gaussian 1.8/.45 A and quintic cutoff 2.2->2.8 A',
        architecture_selection='angular extension motivated after seeing radial-only external QE failure; external QE is reused diagnostic, not untouched model-selection test',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        mace_weights_changed=False,production_enabled=False,material_approved=False,
        source_method_CP2K='GPW/PBE/GTH/DZ/650 Ry',source_method_QE='PBE/PAW; distinct basis and cutoff conventions',
        external_QE='not used to fit delta candidate; MACE pretraining overlap remains unaudited',
        complete_source_compounds=compositions,new_DFT=0,new_MD=0,new_MACE_calls=0)
    dump(args.output/'protocol.json',protocol)
    frames=[];manifest=[]
    for row in eligible:
        index=int(row['source_frame']);f=load_frame(args.cp2k/'raw'/f'{index:04d}.npz',row,args.basis)
        frames.append(f);manifest.append(dict(dataset='CP2K',source_frame=index,sha256=f['raw_sha256']))
    folds=[];held_rows=[];held_energy=[]
    for held in compositions:
        training=[f for f in frames if f['metadata']['composition']!=held]
        testing=[f for f in frames if f['metadata']['composition']==held]
        coefficients,diagnostics=fit(training);train_metrics,_=metrics(training,coefficients);test_metrics,rows=metrics(testing,coefficients)
        folds.append(dict(held_composition=held,fit=diagnostics,train=train_metrics,test=test_metrics))
        held_rows.extend(dict(held_composition=held,**r) for r in rows)
        held_energy.extend(energy_differences(rows))
    coefficients,allfit=fit(frames);all_metrics,all_rows=metrics(frames,coefficients)
    dump(args.output/'cross_validation.json',dict(folds=folds,all_training_fit=allfit,all_training_metrics=all_metrics,
        all_training_relative_energy=energy_differences(all_rows),held_composition_relative_energy=held_energy))
    table(args.output/'held_composition_frames.csv',held_rows)
    qe_metadata=list(csv.DictReader((args.qe/'frames.csv').open(encoding='utf-8')))
    qe_rows=[];qe_groups={}
    for row in qe_metadata:
        index=int(row['frame']);f=load_frame(args.qe/'raw'/f'{index:04d}.npz',row,args.basis)
        _,rows=metrics([f],coefficients);qe_rows.extend(rows)
        manifest.append(dict(dataset='QE',source_frame=index,sha256=f['raw_sha256']))
        b=f['baseline_force']-f['reference_force'];c=b+f['force_basis']@coefficients
        stat=qe_groups.setdefault(row['group'],dict(frames=0,atoms=0,baseline_sq=0.,candidate_sq=0.))
        stat['frames']+=1;stat['atoms']+=f['atoms'];stat['baseline_sq']+=float(b@b);stat['candidate_sq']+=float(c@c)
        if index%100==0:
            dump(args.output/'running.json',dict(QE_completed=index+1,QE_planned=len(qe_metadata),elapsed_s=time.perf_counter()-started))
    table(args.output/'external_QE_frames.csv',qe_rows)
    external=[dict(group=name,frames=s['frames'],atoms=s['atoms'],
        baseline_RMSE_eV_A=float(np.sqrt(s['baseline_sq']/(3*s['atoms']))),
        candidate_RMSE_eV_A=float(np.sqrt(s['candidate_sq']/(3*s['atoms'])))) for name,s in qe_groups.items()]
    dump(args.output/'external_QE.json',dict(groups=external,relative_energy=energy_differences(qe_rows),
        DFT_method_mismatch_is_not_resolved=True))
    table(args.output/'source_raw_manifest.csv',manifest)
    # Aggregation reuses held-composition predictions only, never the all-fit errors.
    weights=np.array([3*r['atoms'] for r in held_rows]);before=np.array([r['baseline_force_RMSE_eV_A'] for r in held_rows])
    after=np.array([r['candidate_force_RMSE_eV_A'] for r in held_rows])
    summary=dict(complete=True,training_frames=487,folds=5,external_QE_frames=len(qe_rows),
        held_composition_baseline_RMSE_eV_A=float(np.sqrt(np.sum(weights*before**2)/weights.sum())),
        held_composition_candidate_RMSE_eV_A=float(np.sqrt(np.sum(weights*after**2)/weights.sum())),
        all_fit=allfit,external_QE_groups=external,elapsed_s=time.perf_counter()-started,
        models_fitted=6,new_MACE_calls=0,new_DFT=0,new_MD=0,material_approved=False,
        does_not_determine_initiation_barriers=True,production_enabled=False)
    dump(args.output/'summary.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('cp2k','qe','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--basis',choices=['radial','angular'],default='radial')
    main(p.parse_args())
