"""Archive exact numerical references so raw validation works on a fresh clone."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil

from .run_local_crack_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path('.cache/si-thermal-v4'))
    parser.add_argument('--output',type=Path,default=Path('results/silicon_thermal_v4/references'))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('choose an empty archive output')
    records=[]
    for name in ['front4','front8','multibasin']:
        source=args.source/name;target=args.output/name;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix!='.npz' and path.name!='summary.json':continue
            destination=target/path.name
            shutil.copyfile(path,destination)
            expected=hashlib.sha256(path.read_bytes()).hexdigest()
            if hashlib.sha256(destination.read_bytes()).hexdigest()!=expected:raise RuntimeError('archive copy differs')
            records.append(dict(source=path.as_posix(),archive=destination.as_posix(),sha256=expected,bytes=path.stat().st_size))
    save_json(args.output/'index.json',dict(files=records,total_bytes=sum(r['bytes'] for r in records),
        purpose='exact archival copy, not new physical calculations',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print('archived',len(records),'files',sum(r['bytes'] for r in records),'bytes',flush=True)


if __name__=='__main__':main()
