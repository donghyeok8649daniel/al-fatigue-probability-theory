"""Independently replay every saved backward vector against the discrete solution."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np


def main(args):
    output=args.result/'vector_replay.json'
    if output.exists():raise ValueError('new replay path required')
    rows=[]
    path=args.result/'selected_backward_solutions.npz'
    with np.load(path) as saved:
        for key in saved.files:
            _,part_rate,part_grid=key.split('_')
            rate=float(part_rate[2:]);n=int(part_grid[1:]);d=.4;p=.2;h=1/n
            theta=np.arccosh(1+rate*h*h/(2*d))
            i=np.arange(n+1);z=n*theta
            ratio=(np.exp(i*theta-z)+np.exp(-i*theta-z))/(-np.expm1(-2*z))
            delta=h*ratio/(p*d*np.sinh(theta))
            average=(1-(i*h)**2)/(2*d)+(1-p)*delta[-1]
            exact=np.r_[average+p*delta,(average-(1-p)*delta)[:-1]]
            numerical=saved[key]
            if exact.shape!=numerical.shape:raise ValueError('backward vector dimension differs')
            error=float(np.max(abs(exact-numerical)));relative=error/float(np.max(abs(exact)))
            if relative>1e-7:raise ValueError('entire backward vector replay differs')
            rows.append(dict(key=key,states=len(exact),maximum_absolute_error_model_time=error,
                maximum_error_over_largest_exact_time=relative))
    result=dict(complete=True,vectors=len(rows),states=sum(r['states'] for r in rows),values=rows,
        maximum_absolute_error_model_time=max(r['maximum_absolute_error_model_time'] for r in rows),
        maximum_error_over_largest_exact_time=max(r['maximum_error_over_largest_exact_time'] for r in rows),
        method='independent arccosh form and full cosh solution; original runner used sparse LU and asinh initial-state formula',
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),new_potential_calls=0,new_DFT=0,new_MD=0,physical_clock=None)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);main(p.parse_args())
