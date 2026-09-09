"""A predeclared deterministic shape-grid diagnosis before choosing more terms."""
import argparse
from itertools import product
import json
from pathlib import Path
import sys

from .diagnose_interface_calibration import main as diagnose


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve previous research runs')
    args.out.mkdir(parents=True)
    starts=list(product((.9,1.6,2.6,4.5,7.),(2.3,3.2,5.5,8.5,11.5),(2.1,4.,7.,11.9)))
    path=args.out/'predeclared_starts.json'
    path.write_text(json.dumps(starts,indent=2)+'\n',encoding='utf-8')
    old=sys.argv
    try:
        sys.argv=[old[0],'--starts',str(path),'--out',str(args.out/'diagnosis'),
                  '--continue-on-numerical-failure']
        diagnose()
    finally:sys.argv=old


if __name__=='__main__':main()
