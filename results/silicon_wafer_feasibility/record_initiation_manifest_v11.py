"""Hash v11 sources, external input provenance and result snapshots.

Use --include-mutable only after stopping/verifying the associated calculations.
Logs and operating PIDs are always excluded. Raw checkpoints otherwise remain
on disk; omission from an interim manifest is not permission to remove them.
"""
from __future__ import annotations
import argparse,ast,hashlib,json
from datetime import datetime,timezone
from pathlib import Path


def digest(path):
    payload=path.read_bytes()
    return dict(bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())


def source_digest(path):
    payload=path.read_bytes()
    record=digest(path)
    if b'\r\n' in payload:
        # Existing shared Python sources can be CRLF in a Windows worktree and
        # LF in Git. Preserve the actually executed bytes and bind the exact LF
        # representation separately; do not rewrite production/shared sources.
        normalized=payload.replace(b'\r\n',b'\n')
        record['git_lf']=dict(bytes=len(normalized),sha256=hashlib.sha256(normalized).hexdigest())
        record['transport']='CRLF to LF only; original working-file digest retained'
    return record


def local_candidates(repo,current,module,level=0):
    parts=module.split('.') if module else []
    if level:
        bases=[current.parent]
        for _ in range(level-1):bases=[bases[0].parent]
    else:bases=[repo,current.parent]
    for base in bases:
        candidate=base.joinpath(*parts)
        for path in (candidate.with_suffix('.py'),candidate/'__init__.py'):
            if path.is_file() and path.resolve().is_relative_to(repo):yield path.resolve()


def main(args):
    repo=args.repo.resolve();results=args.results.resolve()
    if not results.is_relative_to(repo):raise ValueError('results must belong to the supplied repository')
    roots=list((repo/'results/silicon_wafer_feasibility').glob('*v11.py'))
    roots += [repo/'solver_v1'/name for name in ('silicon_initiation_probability.py','silicon_initiation_research.py',
        'silicon_oxide_research.py','test_silicon_initiation_probability.py','test_silicon_initiation_research.py',
        'test_silicon_charge_dynamics.py')]
    pending=[p.resolve() for p in roots];seen=set()
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
                    module='.'.join(filter(None,(node.module,name.name)))
                    if name.name!='*':pending.extend(local_candidates(repo,path,module,node.level))
    code=[dict(path=p.relative_to(repo).as_posix(),**source_digest(p)) for p in sorted(seen)]
    prepared=json.loads((results/'prepared_inputs_validation.json').read_text(encoding='utf-8'))
    paper=json.loads((results/'thermal_oxide_source.json').read_text(encoding='utf-8'))['source']
    model=json.loads((results/'oxidized_surfaces_max320/protocol.json').read_text(encoding='utf-8'))
    source=dict(schema_version=2,created_UTC=datetime.now(timezone.utc).isoformat(),snapshot_label=args.label,
        external_source_inputs=prepared['sources'],provided_thermal_oxide_manuscript=paper,
        model=dict(name='MACE-MP-0b3-medium',sha256=model['model_sha256'],
            url='https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model',
            mode='neutral, CPU float64',pretraining_overlap_audited=False),
        public_datasets=[dict(name='MTPu',commit='4ed842eea9c0b894310fee50e9bbca21404e356e',
            url='https://gitlab.com/Kazongogit/MTPu',doi='10.1038/s41524-024-01390-8'),
            dict(name='MLFF-SiOx',commit=model['source_commit'],url=model['source_url'],doi=model['source_doi'])],
        runtime=dict(python=prepared['python_version'],packages=prepared['versions']),
        source_roots=[p.relative_to(repo).as_posix() for p in sorted(roots)],source_files=code,
        dependency_scope='v11 runners and named tests with recursively resolved static local imports; external package versions listed, not a complete environment lock')
    source_path=results/'source_manifest.json'
    source_path.write_text(json.dumps(source,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    entries=[];excluded=[]
    manifest_path=results/args.filename
    for path in sorted(results.rglob('*')):
        if not path.is_file():continue
        relative=path.relative_to(results).as_posix()
        if path==manifest_path or (path.suffix=='.json' and 'manifest' in path.stem and path!=source_path):
            excluded.append(dict(path=relative,reason='manifest avoids self/previous-manifest recursion'));continue
        if '__pycache__' in path.parts or path.suffix in ('.pyc','.log') or 'pid' in path.name.lower():
            excluded.append(dict(path=relative,reason='runtime/log metadata, retained locally'));continue
        if not args.include_mutable and (path.name in ('running.json','progress.json') or path.name.startswith('checkpoint')):
            excluded.append(dict(path=relative,reason='mutable calculation state excluded from interim snapshot'));continue
        entries.append(dict(path=relative,**digest(path)))
    # Reject a file changing while the manifest is assembled.
    for entry in entries:
        if digest(results/entry['path'])!={k:entry[k] for k in ('bytes','sha256')}:
            raise RuntimeError('result changed during hashing: '+entry['path'])
    manifest=dict(schema_version=1,created_UTC=datetime.now(timezone.utc).isoformat(),snapshot_label=args.label,
        include_mutable=args.include_mutable,source_manifest_sha256=digest(source_path)['sha256'],
        files=entries,excluded_files=excluded,total_bytes=sum(x['bytes'] for x in entries),
        file_hashes_independently_rechecked=True,scope='artifact integrity, not material or kinetic validation')
    manifest_path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(source_files=len(code),artifact_files=len(entries),excluded=len(excluded),
        total_bytes=manifest['total_bytes'],filename=args.filename,snapshot_label=args.label)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--filename',default='artifact_manifest.json');p.add_argument('--label',required=True)
    p.add_argument('--include-mutable',action='store_true');main(p.parse_args())
