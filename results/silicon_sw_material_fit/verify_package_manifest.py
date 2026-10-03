"""Check every listed package file against its exact Git blob bytes."""
import argparse
import hashlib
import json
from pathlib import Path


def main(root):
    manifest = json.loads((root / 'package_manifest.json').read_text(encoding='utf-8'))
    listed = set()
    for row in manifest['files']:
        name = row['path']
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or name in listed:
            raise ValueError('unsafe or duplicate manifest path')
        listed.add(name)
        raw = (root / relative).read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
        if len(raw) != row['bytes'] or sha != row['sha256'] or blob != row['git_blob_sha1']:
            raise RuntimeError('file differs from declared Git bytes: ' + name)
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    expected = listed | {'package_manifest.json'}
    if actual != expected:
        raise RuntimeError('missing or unlisted package files: ' + str(sorted(actual ^ expected)))
    print(json.dumps(dict(complete=True,files_verified=len(listed),
                         exact_git_bytes=True,manifest_self_hash_excluded=True)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, default=Path(__file__).resolve().parent)
    main(parser.parse_args().package)
