"""Export exact staged Git blobs to a fresh directory for independent replay."""
import argparse,hashlib,json,subprocess,os
from pathlib import Path


def main(args):
    args.output=args.output.resolve()
    if os.name=='nt' and not str(args.output).startswith('\\\\?\\'):
        args.output=Path('\\\\?\\'+str(args.output))
    if args.output.exists():raise ValueError('fresh export directory required')
    args.output.mkdir(parents=True)
    paths=['solver_v1','results/silicon_wafer_feasibility','results/silicon_sw_material_fit',
           'results/silicon_atomistic_v5/dft_predictions.npz']
    index=subprocess.check_output(['git','ls-files','--stage','-z','--',*paths],cwd=args.repo)
    entries=[]
    for record in index.split(b'\0'):
        if not record:continue
        header,path=record.split(b'\t',1);mode,blob,stage=header.decode().split()
        if stage!='0' or mode not in ['100644','100755']:raise RuntimeError('unsupported Git index entry')
        name=path.decode('utf-8');relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise RuntimeError('unsafe exported path')
        entries.append((name,blob))
    data=subprocess.check_output(['git','cat-file','--batch'],input=''.join(blob+'\n' for _,blob in entries).encode(),cwd=args.repo)
    offset=0;records=[]
    for name,expected in entries:
        end=data.index(b'\n',offset);blob,kind,size=data[offset:end].decode().split();offset=end+1
        raw=data[offset:offset+int(size)];offset+=int(size)+1
        if kind!='blob' or blob!=expected or hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()!=blob:
            raise RuntimeError('Git blob integrity check failed')
        target=args.output/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        records.append(dict(path=name,git_blob_sha1=blob,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)))
    if offset!=len(data):raise RuntimeError('unconsumed Git batch output')
    package=args.output/'results/silicon_sw_material_fit'
    manifest=json.loads((package/'current_source_manifest.json').read_text())
    for r in manifest:
        raw=(args.output/r['path']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=r['lf_sha256']:raise RuntimeError('current source Git byte normalization differs: '+r['path'])
        actual=next(z['git_blob_sha1'] for z in records if z['path']==r['path'])
        if actual!=r['git_blob_sha1']:raise RuntimeError('source Git object identity differs')
    report=dict(complete=True,exported_git_blobs=len(records),current_source_bindings=len(manifest),
                all_files_read_from_index=True,working_tree_bytes_used=False,records=records)
    (args.output/'git_export_manifest.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps({k:v for k,v in report.items() if k!='records'}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['repo','output']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
