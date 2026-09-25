import sys,json
from pathlib import Path
import numpy as np
import torch
from ase import Atoms
from mace.calculators import MACECalculator
root=Path(sys.argv[1]);model=Path(sys.argv[2]);out=Path(sys.argv[3])
torch.set_num_threads(2);torch.set_num_interop_threads(1)
with np.load(root/'results/silicon_initiation_v12/hessian_probe/partial.npz') as d:
    b=d['basis'];h=d['hessian_rows']
with np.load(root/'results/silicon_initiation_v11/force_controlled_prism_360_tight/state_001/raw.npz') as d:
    r=d['positions'];n=d['numbers'];c=d['cell']
a=Atoms(numbers=n,positions=r,cell=c,pbc=False)
a.calc=MACECalculator(model_paths=str(model),device='cpu',default_dtype='float64')
v=np.random.default_rng(7713).normal(size=b.shape[1]);v/=np.linalg.norm(v)
rows=[]
for step in (2e-4,1e-4):
    a.positions[:]=r+(b@v).reshape(r.shape)*step;gp=-b.T@a.get_forces().ravel()
    a.positions[:]=r-(b@v).reshape(r.shape)*step;gm=-b.T@a.get_forces().ravel()
    numerical=(gp-gm)/(2*step)
    error=np.linalg.norm(numerical[:len(h)]-h@v)
    rows.append(dict(step_A=step,first_rows_error_norm_eV_A2=float(error),reference_norm=float(np.linalg.norm(h@v))))
    assert error<1e-6,rows
out.write_text(json.dumps(dict(checks=rows,passed=True,standard_calculator_calls=4),indent=2)+'\n',encoding='utf-8')
print(out.read_text())
