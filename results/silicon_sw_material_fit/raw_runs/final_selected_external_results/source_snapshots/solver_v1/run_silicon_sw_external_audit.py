"""Force audit on 28 tagged external GB/di-interstitial/SF/amorphous frames.

Selection is locked before evaluation. All XC labels remain separate; energy
zeros are not transferred from the fitting archive. A cell/wrapped-coordinate
overlap audit is recorded without claiming exhaustive geometry equivalence.
"""
import argparse,hashlib,json,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from ase.io import read
from .silicon_sw_material_fit import force_features
from .run_silicon_sw_material_fit import dump
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames


def fingerprint(atoms):
    frac=np.linalg.solve(atoms.cell.array.T,atoms.positions.T).T
    frac=np.round(frac,6)
    for axis,pbc in enumerate(atoms.pbc):
        if pbc:frac[:,axis]%=1.
    coordinates=np.c_[atoms.numbers,frac];coordinates=coordinates[np.lexsort(coordinates.T[::-1])]
    payload=np.r_[atoms.cell.array.ravel().round(6),atoms.pbc.astype(float),coordinates.ravel()]
    return hashlib.sha256(payload.astype('<f8').tobytes()).hexdigest()


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    sources=['solver_v1/run_silicon_sw_external_audit.py','solver_v1/silicon_sw_material_fit.py',
             'solver_v1/silicon_sw_anharmonic.py','solver_v1/silicon_environment_research.py']
    manifest=[]
    for path in sources:
        raw=(root/path).read_bytes();output=args.output/'source_snapshots'/path;output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(raw)
        manifest.append(dict(path=path,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',manifest)
    raw=args.dataset.read_bytes();selection_raw=args.candidate.read_bytes();selected=json.loads(selection_raw)
    beta=selected.get('angle_beta',selected.get('selected',{}).get('beta'));par=selected['parameters']
    if beta is None:raise ValueError('selected beta required')
    frames=read(args.dataset,':',format='extxyz');archive,_=read_frames(args.archive)
    training=set(json.loads((args.baseline/'fit.json').read_text())['training_frames'])
    archive_fingerprints=defaultdict(list)
    for i,a in enumerate(archive):archive_fingerprints[fingerprint(a)].append(i)
    overlap=[]
    for i,a in enumerate(frames):
        matches=archive_fingerprints.get(fingerprint(a),[])
        overlap.append(dict(frame=i,archive_matches=matches,training_matches=[j for j in matches if j in training]))
    dump(args.output/'overlap_audit.json',dict(
        method='same rounded cell and species/sorted wrapped fractional coordinates at 1e-6; atom ordering and wrapping invariant',
        limitations='global translation/rotation or alternative primitive cell equivalence not exhaustively audited',frames=overlap))
    if any(r['training_matches'] for r in overlap):raise RuntimeError('external audit overlaps a fitting frame')
    dump(args.output/'selection_lock.json',dict(sha256=hashlib.sha256(selection_raw).hexdigest(),heldout_used_for_selection=False,
        declared_before_external_evaluation=True,dataset_sha256=hashlib.sha256(raw).hexdigest()))
    p,_=source_parameters();cases=[('original_SW',p,0.),('selected_plain_control',par,0.),('selected_anharmonic',par,beta)]
    rows=[];predictions=[];references=[];offsets=[0];positions=[];cells=[];calls=0
    for i,a in enumerate(frames):
        if np.any(a.numbers!=14):raise ValueError('external dataset is not pure Si')
        target=np.asarray(a.arrays['dft_force'],float);reference=target.copy();pred=[]
        for name,parameters,b in cases:
            energy,f,s=force_features(a.positions,a.cell.array,a.pbc,parameters,angle_beta=b,strain_derivative=True)
            force=f.sum(axis=2);diff=force-target;calls+=1;pred.append(force)
            rows.append(dict(frame=i,case=name,config_type=str(a.info.get('config_type','UNSPECIFIED')),
                declared_xc=str(a.info.get('xc_functional','UNSPECIFIED')),atoms=len(a),
                force_squared_error=float(np.sum(diff**2)),force_RMSE_eV_A=float(np.sqrt(np.mean(diff**2))),
                maximum_component_error_eV_A=float(np.max(abs(diff))),
                net_model_force_eV_A=float(np.max(abs(force.sum(axis=0)))),energy_eV=float(energy.sum()),
                static_stress_eV_A3=(s.sum(axis=2)/a.get_volume()).tolist()))
        predictions.append(np.stack(pred,axis=2));references.append(reference);positions.append(a.positions.copy());cells.append(a.cell.array)
        offsets.append(offsets[-1]+len(a))
    dump(args.output/'frame_metrics.json',rows)
    np.savez_compressed(args.output/'predictions.npz',offsets=offsets,reference_force=np.vstack(references),
        model_force=np.concatenate(predictions),positions=np.vstack(positions),cells=cells)
    groups=defaultdict(list)
    for row in rows:groups[(row['case'],row['config_type'],row['declared_xc'])].append(row)
    metrics=[]
    for (case,kind,xc),g in sorted(groups.items()):
        metrics.append(dict(case=case,config_type=kind,declared_xc=xc,frames=len(g),atoms=sum(r['atoms'] for r in g),
            force_RMSE_eV_A=float(np.sqrt(sum(r['force_squared_error'] for r in g)/(3*sum(r['atoms'] for r in g)))),
            maximum_component_error_eV_A=max(r['maximum_component_error_eV_A'] for r in g)))
    dump(args.output/'group_metrics.json',metrics)
    dump(args.output/'input_bindings.json',[dict(role=role,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for role,path in [('dataset',args.dataset),('selection',args.candidate),('archive',args.archive)]])
    dump(args.output/'summary.json',dict(complete=True,frames=len(frames),atoms=sum(len(a) for a in frames),cases=len(cases),new_feature_evaluations=calls,
        dataset='tagged silicon-testing-framework v1.0 testing_database.xyz; source manifest binds Git blob',
        fitting_frame_fingerprint_matches=0,not_a_new_DFT_run=True,energy_zero_comparison=False,
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        new_DFT=0,new_MD=0,new_LAMMPS=0,elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['dataset','candidate','archive','baseline','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
