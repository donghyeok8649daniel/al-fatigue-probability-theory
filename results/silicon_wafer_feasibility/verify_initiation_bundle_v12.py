"""Verify v12 current/tested bytes plus unchanged upstream data in Git."""
import argparse,contextlib,hashlib,io,json,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.verify_initiation_bundle_v11 import main as verify_current,checked_path


def main(args):
    captured=io.StringIO()
    with contextlib.redirect_stdout(captured):verify_current(args)
    result=json.loads(captured.getvalue())
    source=json.loads((args.repo/args.results/'source_manifest.json').read_text(encoding='utf-8'))
    records={checked_path(entry['path']):entry for entry in source['upstream_artifact_files']}
    for path,entry in records.items():
        payload=(args.repo/path).read_bytes()
        if len(payload)!=entry['bytes'] or hashlib.sha256(payload).hexdigest()!=entry['sha256']:
            raise ValueError('upstream working-file bytes changed: '+path)
    if args.revision:
        prefix=':' if args.revision=='INDEX' else result['git_revision']+':'
        request=''.join(prefix+path+'\n' for path in records).encode('utf-8')
        response=subprocess.run(['git','cat-file','--batch'],cwd=args.repo,input=request,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout
        offset=0
        for path,entry in records.items():
            end=response.index(b'\n',offset);header=response[offset:end].split()
            if len(header)!=3 or header[1]!=b'blob':raise ValueError('upstream Git blob unavailable: '+path)
            size=int(header[2]);payload=response[end+1:end+1+size]
            if response[end+1+size:end+2+size]!=b'\n':raise ValueError('invalid Git batch boundary')
            offset=end+2+size
            if size!=entry['bytes'] or hashlib.sha256(payload).hexdigest()!=entry['sha256']:
                raise ValueError('upstream Git bytes changed: '+path)
        if offset!=len(response):raise ValueError('unexpected extra Git output')
    result.update(upstream_artifact_files_verified=len(records),
        upstream_scope='unchanged stored data; historical v11 source hashes are not asserted against intentionally fixed v12 source',
        total_current_and_upstream_records=result['verified_files']+len(records))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--results',default='results/silicon_initiation_v12');p.add_argument('--revision');main(p.parse_args())
