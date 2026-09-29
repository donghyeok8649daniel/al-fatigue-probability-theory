"""Bind the v14 report, raw outputs, tested sources and upstream inputs to bytes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(data): return hashlib.sha256(data).hexdigest()


def main(root, write=False, revision=None):
    result=root/'results/silicon_thermal_validity_v14';target=result/'package_manifest.json'
    if write:
        if target.exists(): raise ValueError('package already sealed')
        # This existing source has CRLF in the Windows worktree and LF in Git.
        # Verify exact logical equality to the original commit; never modify it.
        legacy='solver_v1/probability_pde_2d.py'
        working=(root/legacy).read_bytes()
        baseline=subprocess.check_output(['git','show','4bfe1e196e7be5df242d0d3f7b8d1231c66d84d5:'+legacy],cwd=root)
        if working.replace(b'\r\n',b'\n')!=baseline:
            raise ValueError('existing probability source changed beyond CRLF')
        inputs=json.loads((result/'input_manifest.json').read_text())
        inputs['accepted_line_ending_sha256']={legacy:[sha(baseline)]}
        inputs['line_ending_scope']='only the unchanged legacy probability source; scientific artifacts remain exact-byte checks'
        (result/'input_manifest.json').write_text(json.dumps(inputs,indent=2)+'\n',encoding='utf-8')
        artifacts=json.loads((result/'artifact_manifest.json').read_text())
        artifacts['sha256']['input_manifest.json']=sha((result/'input_manifest.json').read_bytes())
        (result/'artifact_manifest.json').write_text(json.dumps(artifacts,indent=2)+'\n',encoding='utf-8')
        paths={p for p in result.rglob('*') if p.is_file()}
        paths.update(root/p for p in json.loads((result/'input_manifest.json').read_text())['sha256'])
        paths.update(root/p for p in (
            'solver_v1/silicon_thermal_validity.py','solver_v1/test_silicon_thermal_validity.py',
            'results/silicon_wafer_feasibility/audit_thermal_validity_v14.py',
            'results/silicon_wafer_feasibility/probe_thermal_cloud_v14.py',
            'results/silicon_wafer_feasibility/replay_thermal_validity_v14.py',
            'results/silicon_wafer_feasibility/seal_thermal_validity_v14.py',
            'results/silicon_wafer_feasibility/mace_force_only_v9.py',
            'results/silicon_wafer_feasibility/run_mace_boron_v9.py','.gitattributes'))
        records=[]
        for p in sorted(paths):
            data=p.read_bytes();record=dict(path=p.relative_to(root).as_posix(),bytes=len(data),sha256=sha(data))
            if record['path']==legacy:
                record.update(git_bytes=len(baseline),git_sha256=sha(baseline),git_difference='verified CRLF to LF only')
            records.append(record)
        target.write_text(json.dumps(dict(baseline_git_commit='4bfe1e196e7be5df242d0d3f7b8d1231c66d84d5',
                                         scope='v14 package, source code and explicitly used upstream inputs; excludes this manifest and mutable handoffs',files=records),indent=2)+'\n',encoding='utf-8')
    manifest=json.loads(target.read_text());records=manifest['files']
    for r in records:
        p=Path(r['path'])
        if p.is_absolute() or '..' in p.parts: raise ValueError('unsafe manifest path')
        data=(root/p).read_bytes()
        accepted={(r['bytes'],r['sha256']),(r.get('git_bytes',r['bytes']),r.get('git_sha256',r['sha256']))}
        if (len(data),sha(data)) not in accepted: raise ValueError('working bytes differ: '+r['path'])
    verified_revision=None
    if revision:
        verified_revision='INDEX' if revision=='INDEX' else subprocess.check_output(['git','rev-parse',revision],cwd=root,text=True).strip()
        prefix=':' if revision=='INDEX' else verified_revision+':'
        request=''.join(prefix+r['path']+'\n' for r in records).encode()
        stream=subprocess.run(['git','cat-file','--batch'],cwd=root,input=request,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout
        offset=0
        for r in records:
            end=stream.index(b'\n',offset);header=stream[offset:end].split()
            if len(header)!=3 or header[1]!=b'blob': raise ValueError('missing Git blob: '+r['path'])
            size=int(header[2]);data=stream[end+1:end+1+size];offset=end+size+2
            if size!=r.get('git_bytes',r['bytes']) or sha(data)!=r.get('git_sha256',r['sha256']): raise ValueError('Git bytes differ: '+r['path'])
        if offset!=len(stream): raise ValueError('unexpected Git batch bytes')
    print(json.dumps(dict(complete=True,records_verified=len(records),revision=verified_revision,
                         manifest_sha256=sha(target.read_bytes()),new_model_evaluations=0)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--write',action='store_true');p.add_argument('--revision')
    a=p.parse_args();main(a.root,a.write,a.revision)
