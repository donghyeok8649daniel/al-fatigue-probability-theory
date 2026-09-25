"""Bind v12 source bytes, immutable results, and the unchanged v11 input bundle.

Final manifests must only be made after computations/writers have stopped.
The v11 artifact records pin upstream data, not today's execution of old source
hashes. The intentionally fixed charge validator is pinned separately in v12.
"""
from __future__ import annotations
import argparse,ast,hashlib,json,platform,sys
from datetime import datetime,timezone
from importlib.metadata import version
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.record_initiation_manifest_v11 import digest,source_digest,local_candidates

BASELINE='86360a08abb0a1083dbefb4e4d3024bd959ebd03'


def main(args):
    repo=args.repo.resolve();root=args.results.resolve();upstream=repo/'results/silicon_initiation_v11'
    if not root.is_relative_to(repo):raise ValueError('results outside repository')
    if not args.final:raise ValueError('explicit --final requires stopped/verifiably inactive result writers')
    roots=list((repo/'results/silicon_wafer_feasibility').glob('*v12.py'))
    roots+=list((repo/'solver_v1').glob('*v12.py'))
    roots+=[repo/'solver_v1'/name for name in ('silicon_charge_dynamics.py','test_silicon_charge_dynamics.py',
        'silicon_initiation_probability.py','silicon_initiation_research.py',
        'test_silicon_initiation_probability.py','test_silicon_initiation_research.py')]
    roots+=[repo/'results/silicon_wafer_feasibility'/name for name in
        ('audit_charge_dynamics_v9.py','run_force_controlled_prism_v11.py',
         'analyze_surface_environments_v11.py','summarize_surface_environments_v11.py')]
    pending=[path.resolve() for path in roots];seen=set()
    while pending:
        path=pending.pop()
        if path in seen:continue
        seen.add(path);tree=ast.parse(path.read_text(encoding='utf-8-sig'),filename=path.name)
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):
                for name in node.names:pending.extend(local_candidates(repo,path,name.name))
            elif isinstance(node,ast.ImportFrom):
                pending.extend(local_candidates(repo,path,node.module or '',node.level))
                for name in node.names:
                    if name.name!='*':pending.extend(local_candidates(repo,path,'.'.join(filter(None,(node.module,name.name))),node.level))
    code=[dict(path=path.relative_to(repo).as_posix(),**source_digest(path)) for path in sorted(seen)]
    old_manifest=json.loads((upstream/'artifact_manifest.json').read_text(encoding='utf-8'))
    inputs={}
    for entry in old_manifest['files']:
        path=upstream/entry['path'];actual=digest(path)
        if actual!={key:entry[key] for key in ('bytes','sha256')}:raise ValueError('upstream artifact changed: '+entry['path'])
        inputs[path.relative_to(repo).as_posix()]=actual
    for filename in ('source_manifest.json','artifact_manifest.json'):
        path=upstream/filename;inputs[path.relative_to(repo).as_posix()]=digest(path)
    for filename in ('charge_rate_convergence.csv','spatial_grid_convergence.csv','absorption_convergence.csv',
            'quenched_mixture_counterexample.csv','dynamics_summary.json'):
        path=repo/'results/silicon_doping_v9/dynamics'/filename;inputs[path.relative_to(repo).as_posix()]=digest(path)
    # Do not store private absolute paths in the public provenance bundle.
    external=[]
    expected={
        'baseline_model':('2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012',
            'https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model'),
        'comparison_model':('75428afe3a1d7d8062e19bcaabd5c433623cabf308242ec9fb493e38604fb638',
            'https://github.com/ACEsuit/mace-foundations/releases/download/mace_mpa_0/mace-mpa-0-medium.model'),
        'cp2k_source':('d99170bff4942dcd5a1c38d998fe0a60084a5e132bcb3cf2335f5db4cbe17689',
            'https://github.com/lukas-cvitkovich/MLFF-SiOx/tree/b61806d3a097e21d7880db39dcc077573e2c3c92')}
    for key,(sha,url) in expected.items():
        path=getattr(args,key);record=digest(path)
        if record['sha256']!=sha:raise ValueError('external source changed: '+key)
        external.append(dict(name=key,filename=path.name,url=url,**record))
    source=dict(schema_version=1,created_UTC=datetime.now(timezone.utc).isoformat(),snapshot_label=args.label,
        baseline_git_commit=BASELINE,source_files=code,source_roots=[p.relative_to(repo).as_posix() for p in sorted(roots)],
        upstream_artifact_files=[dict(path=path,**entry) for path,entry in sorted(inputs.items())],
        upstream_scope='v11 immutable data and v9 regression targets; old source manifest is historical evidence, not a claim its validator source is unchanged',
        external_inputs=external,python_version=platform.python_version(),packages={name:version(name) for name in
            ('numpy','scipy','ase','mace-torch','torch','e3nn','matplotlib','pytest')},
        dependency_scope='static local import closure of implemented/prepared v12 tools plus explicit legacy dependencies; source coverage is not an execution ledger or complete external dependency lock',
        material_or_kinetic_validation=False)
    source_path=root/'source_manifest.json'
    source_path.write_text(json.dumps(source,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    artifact_path=root/'artifact_manifest.json';entries=[];excluded=[]
    for path in sorted(root.rglob('*')):
        if not path.is_file():continue
        relative=path.relative_to(root).as_posix()
        if path==artifact_path or '__pycache__' in path.parts or path.suffix in ('.pyc','.log'):
            excluded.append(dict(path=relative,reason='self manifest or runtime log/cache'));continue
        entries.append(dict(path=relative,**digest(path)))
    for entry in entries:
        if digest(root/entry['path'])!={key:entry[key] for key in ('bytes','sha256')}:raise RuntimeError('result changed during final hashing')
    artifact=dict(schema_version=1,created_UTC=datetime.now(timezone.utc).isoformat(),snapshot_label=args.label,
        include_mutable=True,source_manifest_sha256=digest(source_path)['sha256'],files=entries,
        excluded_files=excluded,total_bytes=sum(entry['bytes'] for entry in entries),
        file_hashes_independently_rechecked=True,scope='exact-byte provenance; no scientific adoption')
    artifact_path.write_text(json.dumps(artifact,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(source_files=len(code),artifact_files=len(entries),upstream_files=len(inputs),total_bytes=artifact['total_bytes'])))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('repo','results','baseline-model','comparison-model','cp2k-source'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--label',required=True);p.add_argument('--final',action='store_true');main(p.parse_args())
