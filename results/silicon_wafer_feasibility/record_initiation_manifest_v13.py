"""Freeze v13 source/results and verify the unchanged v12/v11 input chain.

Run --final only after all scientific jobs and document writers have stopped.
This hashes data; it neither reruns models/tests nor certifies material physics.
Existing verify_initiation_bundle_v12.py accepts this compatible manifest schema.
"""
from __future__ import annotations
import argparse, ast, hashlib, json, platform, sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.record_initiation_manifest_v11 import digest, source_digest, local_candidates

BASELINE = 'c3a93f5d2e4df6e21a68d6080ed99cf6aee840bc'
LEGACY_RUNNERS = (
    'audit_dense_hessians_v12.py', 'analyze_hessian_localization_v12.py',
    'analyze_local_harmonic_v12.py', 'augment_grip_hessian_v12.py',
    'run_large_oxide_audit_v12.py', 'validate_large_oxide_replay_v12.py',
    'compare_mpa_model_v12.py', 'replay_mpa_comparison_v12.py',
    'run_force_controlled_prism_v11.py', 'audit_zero_force_return_v12.py',
    'audit_local_anharmonic_v12.py', 'verify_initiation_bundle_v12.py',
)


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main(args):
    repo = args.repo.resolve(); root = args.results.resolve()
    if root != repo / 'results/silicon_initiation_v13':
        raise ValueError('This recorder is scoped to the v13 result directory')
    if not args.final:
        raise ValueError('--final requires all scientific and record writers to be stopped')
    validation = read_json(root / 'validation.json')
    roots = list((repo / 'results/silicon_wafer_feasibility').glob('*v13.py'))
    roots += list((repo / 'solver_v1').glob('*v13.py'))
    roots += [repo / 'solver_v1' / name for name in ('hessian_resume.py', 'test_hessian_resume.py')]
    roots += [repo / 'results/silicon_wafer_feasibility' / name for name in LEGACY_RUNNERS]
    roots += [repo / item['path'] for item in validation['source_files']]
    pending = [path.resolve() for path in roots]; seen = set()
    while pending:
        path = pending.pop()
        if path in seen: continue
        if not path.is_relative_to(repo): raise ValueError('source outside repository')
        seen.add(path)
        tree = ast.parse(path.read_text(encoding='utf-8-sig'), filename=path.name)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for name in node.names: pending.extend(local_candidates(repo, path, name.name))
            elif isinstance(node, ast.ImportFrom):
                pending.extend(local_candidates(repo, path, node.module or '', node.level))
                for name in node.names:
                    if name.name != '*':
                        pending.extend(local_candidates(repo, path, '.'.join(filter(None, (node.module, name.name))), node.level))
    code = [dict(path=path.relative_to(repo).as_posix(), **source_digest(path)) for path in sorted(seen)]
    # The checkout transport contract is part of reproducibility, not Python code.
    code.append(dict(path='.gitattributes', **source_digest(repo / '.gitattributes')))
    lookup = {item['path']: item['sha256'] for item in code}
    for entry in validation['source_files']:
        if lookup[entry['path']] != entry['sha256']:
            raise ValueError('Tested source changed: ' + entry['path'])
    upstream = repo / 'results/silicon_initiation_v12'
    previous_source = read_json(upstream / 'source_manifest.json')
    previous_artifacts = read_json(upstream / 'artifact_manifest.json')
    if digest(upstream / 'source_manifest.json')['sha256'] != previous_artifacts['source_manifest_sha256']:
        raise ValueError('v12 source/artifact chain changed')
    inputs = {}
    def pin(path, expected=None):
        record = digest(path)
        if expected is not None and record != {key: expected[key] for key in ('bytes', 'sha256')}:
            raise ValueError('Preserved input changed: ' + str(path.relative_to(repo)))
        relative = path.relative_to(repo).as_posix()
        if relative in inputs and inputs[relative] != record:
            raise ValueError('Conflicting upstream provenance: ' + relative)
        inputs[relative] = record
    for entry in previous_source['upstream_artifact_files']:
        pin(repo / entry['path'], entry)
    for entry in previous_artifacts['files']:
        pin(upstream / entry['path'], entry)
    for filename in ('source_manifest.json', 'artifact_manifest.json'):
        pin(upstream / filename)
    external = []
    expected_external = {
        'baseline_model': '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012',
        'comparison_model': '75428afe3a1d7d8062e19bcaabd5c433623cabf308242ec9fb493e38604fb638',
        'cp2k_source': 'd99170bff4942dcd5a1c38d998fe0a60084a5e132bcb3cf2335f5db4cbe17689',
    }
    for entry in previous_source['external_inputs']:
        key = entry['name']; path = getattr(args, key); record = digest(path)
        if key not in expected_external or record['sha256'] != expected_external[key] or record != {k: entry[k] for k in ('bytes', 'sha256')}:
            raise ValueError('External input changed: ' + key)
        external.append(dict(name=key, filename=path.name, url=entry['url'], **record))
    if {item['name'] for item in external} != set(expected_external):
        raise ValueError('External input inventory incomplete')
    source = dict(schema_version=1, created_UTC=datetime.now(timezone.utc).isoformat(), snapshot_label=args.label,
        baseline_git_commit=BASELINE, source_files=code,
        source_roots=sorted({p.relative_to(repo).as_posix() for p in roots}),
        upstream_artifact_files=[dict(path=p, **entry) for p, entry in sorted(inputs.items())],
        upstream_scope='all frozen v12 artifacts and their already pinned upstream v11/v9 data; historical source manifests are evidence, not current execution claims',
        external_inputs=external, python_version=platform.python_version(),
        packages={name: version(name) for name in ('numpy','scipy','ase','mace-torch','torch','e3nn','matplotlib','pytest')},
        dependency_scope='static local imports of v13 tools/tests and explicit reused runners; implemented or prepared source coverage does not imply execution',
        material_or_kinetic_validation=False)
    source_path = root / 'source_manifest.json'
    source_path.write_text(json.dumps(source, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    artifact_path = root / 'artifact_manifest.json'; entries = []; excluded = []
    for path in sorted(root.rglob('*')):
        if not path.is_file(): continue
        relative = path.relative_to(root).as_posix()
        if path == artifact_path or '__pycache__' in path.parts or path.suffix in ('.pyc', '.log'):
            excluded.append(dict(path=relative, reason='self manifest or runtime log/cache')); continue
        entries.append(dict(path=relative, **digest(path)))
    for entry in entries:
        if digest(root / entry['path']) != {key: entry[key] for key in ('bytes', 'sha256')}:
            raise RuntimeError('Result changed during hashing: ' + entry['path'])
    artifact = dict(schema_version=1, created_UTC=datetime.now(timezone.utc).isoformat(), snapshot_label=args.label,
        include_mutable=True, source_manifest_sha256=digest(source_path)['sha256'], files=entries,
        excluded_files=excluded, total_bytes=sum(entry['bytes'] for entry in entries),
        file_hashes_independently_rechecked=True, scope='exact-byte provenance, not scientific adoption')
    artifact_path.write_text(json.dumps(artifact, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(source_files=len(code), artifact_files=len(entries), upstream_files=len(inputs), total_bytes=artifact['total_bytes'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('repo','results','baseline-model','comparison-model','cp2k-source'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--final', action='store_true')
    main(parser.parse_args())
