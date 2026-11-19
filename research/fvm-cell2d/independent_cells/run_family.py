"""Run actual independent-cell 2D FVM; never relabel old G4 outputs as new runs."""
import argparse,hashlib,itertools,json,os,time,traceback
from concurrent.futures import ProcessPoolExecutor
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import fvm2d as fvm

ROOT=Path(__file__).resolve().parents[1]
CASES=[(0,0,1),(1,0,1),(0,1,1),(1,1,1),(2,0,1),(0,2,1),(1,1,2),(1,1,4)]
FORCES=(1.,4.5,-4.5,1.)
DURATIONS=(.1,.3,.3,.1)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def save_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')
    os.replace(tmp,path)


def save_npz(path,**values):
    path=Path(path);tmp=path.with_name(path.name+'.tmp')
    with tmp.open('wb') as stream:np.savez_compressed(stream,**values)
    os.replace(tmp,path)


def case_name(a,s,o):return f'a{a}_s{s}_outer{o}'


def gates(results,cycles):
    available={r['case']:r for r in results}
    corner=list(available)[:4]
    space=[]
    def compare(x,y,kind):
        left,right=available[x]['records'],available[y]['records']
        diffs=[]
        for step in range(3,cycles*4,4):
            a,b=left[step]['F'],right[step]['F']
            value=abs(a-b)/min(a,b) if min(a,b)>0 else (0. if a==b else 1e300)
            diffs.append(value)
        return dict(left=x,right=y,kind=kind,cycle_end_raw_relative_differences=diffs,
                    maximum=max(diffs),accepted=max(diffs)<=.01)
    for x,y in itertools.combinations(corner,2):space.append(compare(x,y,'six_corner_pairs'))
    keys=list(available)
    if len(keys)>=6:
        space.extend([compare(keys[4],keys[1],'normal_probe'),compare(keys[5],keys[2],'slip_probe')])
    outer=[]
    if len(keys)==8:
        for x,y in [(keys[3],keys[6]),(keys[6],keys[7])]:
            d=[abs(a['F']-b['F']) for a,b in zip(available[x]['records'],available[y]['records'])]
            outer.append(dict(left=x,right=y,phase_end_absolute_differences=d,maximum=max(d),accepted=max(d)<=1e-12))
    return dict(space=space,outer=outer,space_accepted=all(r['accepted'] for r in space),
                outer_accepted=len(outer)==2 and all(r['accepted'] for r in outer))


def run_case(directory,a,s,outer,target_cycles,source_hashes,refine=1,refine_s=1,resume_revalidated=False):
    case=case_name(a,s,outer);folder=directory/case;folder.mkdir(parents=True,exist_ok=True)
    grid=fvm.Grid2D.original(a,s,outer)
    if refine>1 or refine_s>1:
        def split(edges,count):
            return np.r_[np.concatenate([np.linspace(l,r,count+1)[:-1] for l,r in zip(edges[:-1],edges[1:])]),edges[-1]]
        grid=fvm.Grid2D(split(grid.edges[0],refine),split(grid.edges[1],refine_s))
    path=folder/'checkpoint.npz';record_path=folder/'records.json'
    start=time.perf_counter()
    ops={force:fvm.assemble(grid,force) for force in set(FORCES)}
    if path.exists():
        prior_setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'))
        with np.load(path,allow_pickle=False) as data:
            old_source=str(data['source_id'])
            if old_source!=source_hashes['source_id']:
                if not resume_revalidated:raise RuntimeError('Checkpoint source mismatch')
                if prior_setup['model']!=fvm.model_metadata() or prior_setup['shape']!=list(grid.shape):
                    raise RuntimeError('Physical model/grid changed across checkpoint')
                physics=('fvm2d.py','joint_model.py','mesh_plan_loader.py','cram_action.py')
                for filename in physics:
                    old=next(v for k,v in prior_setup['sources'].items() if k.endswith(filename))
                    new=next(v for k,v in source_hashes.items() if k.endswith(filename))
                    if old!=new:raise RuntimeError(f'Unreviewed physics source change: {filename}')
                check=ROOT/'outputs/fvm-cell2d/FAST_ACTION_CHECK.json'
                if not json.loads(check.read_text())['accepted']:raise RuntimeError('New action not validated')
                save_json(folder/f'resume_{old_source[:12]}_to_{source_hashes["source_id"][:12]}.json',
                    dict(old_source_id=old_source,new_source_id=source_hashes['source_id'],
                         old_checkpoint_sha256=sha(path),physical_sources_identical=True,
                         state_values_modified=False,new_action_check_sha256=sha(check)))
            mass=data['mass'];cdf=data['F'];done=int(data['completed_phase'])
        records=json.loads(record_path.read_text(encoding='utf-8'))[:done]
        for record in records:record.setdefault('source_id',old_source)
        initial=json.loads((folder/'setup.json').read_text(encoding='utf-8'))['initial']
    else:
        p,initial=fvm.initial_mass(grid)
        mass=np.stack([p,p]);cdf=np.zeros(2);done=0;records=[]
        save_npz(folder/'initial.npz',mass=mass,F=cdf,a_edges=grid.edges[0],s_edges=grid.edges[1])
    setup=dict(case=case,axes_levels=[a,s],outer=outer,subdivisions=[refine,refine_s],shape=list(grid.shape),
        total_states_per_cell=grid.size,probability_table_shape=list(mass.shape),
        probability_table_bytes=mass.nbytes,physical_rows='two independent illustrative replicas under identical local load',
        actual_specimen_mesh_solved=False,model=fvm.model_metadata(),initial=initial,
        operators={str(force):op.metadata for force,op in ops.items()},sources=source_hashes)
    save_json(folder/'setup.json',setup)
    for stage in range(done,target_cycles*4):
        phase=stage%4;force=FORCES[phase];duration=DURATIONS[phase]
        mass,cdf,audit=fvm.advance_batch(ops[force],mass,cdf,duration)
        record=dict(phase=stage+1,cycle=stage//4+1,force=force,duration=duration,source_id=source_hashes['source_id'],
            model_time=sum(DURATIONS)*(stage//4)+sum(DURATIONS[:phase+1]),
            F=float(cdf[0]),F_per_cell=cdf.tolist(),U=0.,survivor_mass=mass.sum(axis=1).tolist(),
            two_replica_first_event_probability=float(-np.expm1(np.log1p(-cdf).sum())),
            **audit)
        records.append(record)
        save_json(record_path,records)
        save_npz(path,mass=mass,F=cdf,completed_phase=np.asarray(stage+1),
            source_id=np.asarray(source_hashes['source_id']),a_edges=grid.edges[0],s_edges=grid.edges[1])
        print(json.dumps(dict(case=case,phase=stage+1,F=record['F'],seconds=time.perf_counter()-start)),flush=True)
    result=dict(case=case,setup=setup,records=records,completed_phases=len(records),
        runtime_seconds=time.perf_counter()-start,checkpoint_bytes=path.stat().st_size,
        checkpoint_sha256=sha(path),source_id=source_hashes['source_id'])
    save_json(folder/'result.json',result)
    return result


def initialize_worker(model):
    fvm.set_model(model)
    import fast_propagation
    fvm.advance_batch=fast_propagation.advance


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='outputs/fvm-cell2d')
    parser.add_argument('--model',choices=['coherent_g4','n1_local'],default='coherent_g4')
    parser.add_argument('--base',type=int,default=0)
    parser.add_argument('--base-s',type=int,default=None)
    parser.add_argument('--refine',type=int,default=1)
    parser.add_argument('--refine-s',type=int,default=1)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--resume-revalidated',action='store_true')
    args=parser.parse_args();fvm.set_model(args.model)
    if args.base_s is None:args.base_s=args.base
    out=ROOT/args.out;out.mkdir(parents=True,exist_ok=True)
    sources={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),Path(fvm.__file__),
        Path(fvm.original.__file__),ROOT/'independent_cells/mesh_plan_loader.py',ROOT/'independent_cells/cram_action.py',
        ROOT/'independent_cells/fast_propagation.py']}
    sources['source_id']=hashlib.sha256(json.dumps(sources,sort_keys=True).encode()).hexdigest()
    start=time.perf_counter();manifest=dict(model=fvm.model_metadata(),source_hashes=sources,
        forces=FORCES,durations=DURATIONS,cycles=6,base=[args.base,args.base_s],refine=[args.refine,args.refine_s],
        started_utc=datetime.now(timezone.utc).isoformat(),pid=os.getpid(),workers=args.workers,
        numerical_library_threads_per_worker=1)
    save_json(out/'RUN_MANIFEST.json',manifest)
    results=[]
    try:
        with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize_worker,initargs=(args.model,)) as pool:
            futures=[pool.submit(run_case,out,a+args.base,s+args.base_s,o,1,sources,args.refine,args.refine_s,args.resume_revalidated) for a,s,o in CASES[:4]]
            results=[future.result() for future in futures]
        first=gates(results,1);save_json(out/'FIRST_CYCLE_CHECK.json',first)
        if not first['space_accepted']:
            save_json(out/'FINAL_CHECKS.json',dict(accepted=False,complete=False,
                reason='first-cycle raw six-pair 1% gate failed; no remaining-cycle launch',
                active_dimension=2,model=fvm.model_metadata(),first_cycle=first,source_hashes=sources))
            return 2
        with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize_worker,initargs=(args.model,)) as pool:
            futures=[pool.submit(run_case,out,a+args.base,s+args.base_s,o,6,sources,args.refine,args.refine_s,args.resume_revalidated) for a,s,o in CASES]
            results=[future.result() for future in futures]
        checks=gates(results,6)
        temporal=all(r['maximum_successive_L1']<2e-12 and r['maximum_successive_CDF']<2e-15
            and r['maximum_negative_mass']<2e-12 and r['maximum_mass_residual']<1e-9
            for c in results for r in c['records'])
        final=dict(accepted=checks['space_accepted'] and checks['outer_accepted'] and temporal,
            complete=True,active_dimension=2,model=fvm.model_metadata(),space_and_outer=checks,
            temporal_and_mass_accepted=temporal,cases=8,phase_records=192,
            numerical_unresolved_U=0.,full_state_domain_retained=True,source_hashes=sources,
            seconds=time.perf_counter()-start,probability_output_correction=False,
            physical_fatigue_validation_pass=False,old_4D_runs_reclassified_as_new_2D=False)
        save_json(out/'FINAL_CHECKS.json',final)
        save_json(out/'report_data.json',dict(manifest=manifest,final_checks=final,cases=results))
        return 0 if final['accepted'] else 3
    except Exception:
        save_json(out/'FAILURE.json',dict(error=traceback.format_exc(),source_hashes=sources))
        raise
    finally:
        save_json(out/'PROCESS_EXIT.json',dict(pid=os.getpid(),finished_utc=datetime.now(timezone.utc).isoformat(),
            actual_completed_main_body=True,seconds=time.perf_counter()-start))


if __name__=='__main__':raise SystemExit(main())
