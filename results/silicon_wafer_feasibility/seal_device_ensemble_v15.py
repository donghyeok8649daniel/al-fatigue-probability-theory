"""Bind v15 result and tested source bytes; compare checkout and Git blobs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

def digest(data): return hashlib.sha256(data).hexdigest()

def main(root,write=False,revision=None):
    result=root/'results/silicon_device_ensemble_v15';target=result/'package_manifest.json'
    if write:
        if target.exists():raise ValueError('already sealed')
        paths={p for p in result.rglob('*') if p.is_file()}
        paths.update(root/p for p in json.loads((result/'device_audit/input_manifest.json').read_text())['sha256'])
        paths.update(root/p for p in json.loads((result/'nonlinear_pilot/protocol.json').read_text())['source_sha256'])
        paths.update(root/p for p in (
            'solver_v1/silicon_device_ensemble.py','solver_v1/test_silicon_device_ensemble.py',
            'solver_v1/silicon_thermal_research.py','solver_v1/test_silicon_thermal_research.py',
            'solver_v1/silicon_conditional_research.py',
            'results/silicon_wafer_feasibility/audit_device_ensemble_v15.py',
            'results/silicon_wafer_feasibility/sample_nonlinear_ensemble_v15.py',
            'results/silicon_wafer_feasibility/replay_device_ensemble_v15.py',
            'results/silicon_wafer_feasibility/seal_device_ensemble_v15.py','.gitattributes'))
        files=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=digest(p.read_bytes())) for p in sorted(paths)]
        target.write_text(json.dumps(dict(baseline_commit='b02cb65398ab8ddc6c905c4fca895982da8e46b7',
            scope='v15 results, used upstream data, exact sources; excludes this manifest and external CAD file',files=files),indent=2)+'\n',encoding='utf-8')
    records=json.loads(target.read_text())['files']
    for record in records:
        path=Path(record['path'])
        if path.is_absolute() or '..' in path.parts:raise ValueError('invalid manifest path')
        data=(root/path).read_bytes()
        if len(data)!=record['bytes'] or digest(data)!=record['sha256']:raise ValueError('worktree mismatch '+record['path'])
    if revision:
        prefix=':' if revision=='INDEX' else subprocess.check_output(['git','rev-parse',revision],cwd=root,text=True).strip()+':'
        request=''.join(prefix+r['path']+'\n' for r in records).encode()
        data=subprocess.run(['git','cat-file','--batch'],cwd=root,input=request,stdout=subprocess.PIPE,check=True).stdout
        cursor=0
        for r in records:
            end=data.index(b'\n',cursor);header=data[cursor:end].split()
            if len(header)!=3 or header[1]!=b'blob':raise ValueError('missing Git blob '+r['path'])
            size=int(header[2]);blob=data[end+1:end+1+size];cursor=end+size+2
            if size!=r['bytes'] or digest(blob)!=r['sha256']:raise ValueError('Git mismatch '+r['path'])
        if cursor!=len(data):raise ValueError('unexpected Git bytes')
    print(json.dumps(dict(complete=True,records_verified=len(records),revision=revision,manifest_sha256=digest(target.read_bytes())),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--write',action='store_true');p.add_argument('--revision')
    a=p.parse_args();main(a.root,a.write,a.revision)
