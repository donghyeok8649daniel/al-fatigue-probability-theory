"""Read-only verification of v11 manifest hashes in working files and Git blobs.

Run after the final manifest is recorded. INDEX checks the staged version;
any other revision must resolve to a commit. This is provenance verification,
not a repeated calculation or a material/kinetic validation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
from pathlib import Path, PurePosixPath


def checked_path(value):
    path = PurePosixPath(value)
    if path.is_absolute() or '..' in path.parts or any(c in value for c in '\r\n\0:'):
        raise ValueError('unsafe manifest path')
    return path.as_posix()


def main(args):
    repo = args.repo.resolve()
    relative = checked_path(args.results)
    root = repo / relative
    artifact = json.loads((root / 'artifact_manifest.json').read_text(encoding='utf-8'))
    source_bytes = (root / 'source_manifest.json').read_bytes()
    source = json.loads(source_bytes)
    if hashlib.sha256(source_bytes).hexdigest() != artifact['source_manifest_sha256']:
        raise ValueError('source manifest does not match artifact manifest')
    if not artifact['include_mutable']:
        raise ValueError('final immutable snapshot required')
    records = {}
    for entry in source['source_files']:
        records[checked_path(entry['path'])] = entry
    for entry in artifact['files']:
        records[relative + '/' + checked_path(entry['path'])] = entry
    # Also compare the manifests themselves with the selected Git revision.
    for name in ('source_manifest.json', 'artifact_manifest.json'):
        payload = (root / name).read_bytes()
        records[relative + '/' + name] = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    worktree_transports = []
    for path, entry in records.items():
        payload = (repo / path).read_bytes()
        actual = dict(bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())
        if actual != {key:entry[key] for key in ('bytes','sha256')}:
            normalized = payload.replace(b'\r\n',b'\n')
            normalized_digest = dict(bytes=len(normalized),sha256=hashlib.sha256(normalized).hexdigest())
            if 'git_lf' not in entry or normalized_digest != entry['git_lf']:
                raise ValueError('working file hash mismatch: ' + path)
            # Existing execution sources include mixed CRLF/LF. A fresh Git
            # checkout can use uniform CRLF or LF. Only the four explicitly
            # bound canonical source representations permit this transport.
            worktree_transports.append(dict(path=path,observed_bytes=actual['bytes'],
                observed_sha256=actual['sha256'],canonical_LF_sha256=entry['git_lf']['sha256']))
        if 'git_lf' in entry:
            normalized = payload.replace(b'\r\n',b'\n')
            if (len(normalized) != entry['git_lf']['bytes']
                    or hashlib.sha256(normalized).hexdigest() != entry['git_lf']['sha256']):
                raise ValueError('recorded CRLF-to-LF transport differs: ' + path)
    validation = json.loads((root / 'validation.json').read_text(encoding='utf-8'))
    for entry in validation['source_files']:
        path = checked_path(entry['path'])
        if path not in records or entry['sha256'] != records[path]['sha256']:
            raise ValueError('tested source differs from final source: ' + path)
    revision = None
    if args.revision:
        if args.revision == 'INDEX':
            revision = 'INDEX'
            prefix = ':'
        else:
            if args.revision.startswith('-') or any(c in args.revision for c in '\r\n\0:'):
                raise ValueError('invalid revision')
            revision = subprocess.check_output(
                ['git', 'rev-parse', '--verify', args.revision + '^{commit}'], cwd=repo, text=True).strip()
            prefix = revision + ':'
        queries = ''.join(prefix + path + '\n' for path in records).encode('utf-8')
        result = subprocess.run(['git', 'cat-file', '--batch'], cwd=repo, input=queries,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        stream = result.stdout
        offset = 0
        for path, entry in records.items():
            end = stream.index(b'\n', offset)
            header = stream[offset:end].split()
            if len(header) != 3 or header[1] != b'blob':
                raise ValueError('Git blob missing or invalid: ' + path)
            size = int(header[2])
            payload = stream[end + 1:end + 1 + size]
            if stream[end + 1 + size:end + 2 + size] != b'\n':
                raise ValueError('invalid Git batch boundary')
            offset = end + 2 + size
            expected_git = entry.get('git_lf',entry)
            if size != expected_git['bytes'] or hashlib.sha256(payload).hexdigest() != expected_git['sha256']:
                raise ValueError('Git stored bytes differ from manifest: ' + path)
        if offset != len(stream):
            raise ValueError('unexpected extra Git batch output')
    print(json.dumps(dict(verified_files=len(records), source_files=len(source['source_files']),
        artifact_files=len(artifact['files']), tested_source_files=len(validation['source_files']),
        recorded_tests_passed=validation['tests_passed'], tests_rerun=False,
        working_tree_hashes_verified=True, git_blob_hashes_verified=revision is not None,
        CRLF_to_LF_source_paths=[path for path,entry in records.items() if 'git_lf' in entry],
        current_worktree_line_ending_transports=worktree_transports,
        result_artifacts_require_exact_bytes=True,
        git_revision=revision, scope='stored byte integrity and test-source binding; no physical validation'), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--results', default='results/silicon_initiation_v11')
    p.add_argument('--revision')
    main(p.parse_args())
