"""Reproduce v11 validation failure from exact Git source and check its repair."""
import argparse,hashlib,json,subprocess,sys,types
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_charge_dynamics import reversible_propagator

def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    baseline='86360a08abb0a1083dbefb4e4d3024bd959ebd03'
    raw=subprocess.check_output(['git','show',baseline+':solver_v1/silicon_charge_dynamics.py'])
    module=types.ModuleType('solver_v1._v11_charge_reference');module.__package__='solver_v1'
    exec(compile(raw,'v11_git_source','exec'),module.__dict__)
    matrix=1e-15*np.array([[-1.,1.],[1.,-1.1]])
    before=module.reversible_propagator(matrix,np.array([.5,.5]),1e15)
    rejection=None
    try:reversible_propagator(matrix,np.array([.5,.5]),1e15)
    except ValueError as exc:rejection=str(exc)
    if rejection is None:raise ValueError('repair did not reject invalid slow generator')
    result=dict(baseline_commit=baseline,baseline_source_git_sha256=hashlib.sha256(raw).hexdigest(),
        leaking_generator=matrix.tolist(),column_residual=matrix.sum(axis=0).tolist(),
        old_propagator=before.tolist(),old_column_masses=before.sum(axis=0).tolist(),
        old_maximum_mass_loss=float(np.max(1-before.sum(axis=0))),new_rejection=rejection,
        scope='synthetic implementation regression; no Si rates or physical time',new_DFT=0,new_MD=0)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
