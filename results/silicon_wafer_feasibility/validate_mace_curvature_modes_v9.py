"""Independent step refinement along lowest/highest fixed-cell curvature modes."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import eigh,null_space
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def main(args):
    import torch
    from ase.io import read
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model mismatch')
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    args.output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for path in sorted(args.hessians.glob('*_hessian.npz')):
        name=path.name.removesuffix('_hessian.npz')
        atoms=read(args.states/(name+'_final.extxyz'))
        h=np.load(path)['hessian']
        trans=np.tile(np.eye(3),(len(atoms),1))/np.sqrt(len(atoms))
        q=null_space(trans.T)
        values,vectors=eigh(q.T@(.5*(h+h.T))@q)
        for label,index in [('lowest',0),('highest',-1)]:
            direction=q@vectors[:,index]
            for step in [1e-3,3e-4,1e-4,3e-5]:
                plus,minus=atoms.copy(),atoms.copy()
                plus.positions+=step*direction.reshape(-1,3)
                minus.positions-=step*direction.reshape(-1,3)
                hp=-(force_only(calc,plus)[1]-force_only(calc,minus)[1]).ravel()/(2*step)
                curvature=float(direction@hp)
                rows.append(dict(structure=name,mode=label,step_A=step,
                    matrix_curvature_eV_A2=float(values[index]),force_difference_curvature_eV_A2=curvature,
                    absolute_curvature_error=abs(curvature-values[index]),
                    Hv_L2_error=float(np.linalg.norm(hp-values[index]*direction))))
        print('REFINED',name,values[0],flush=True)
    (args.output/'mode_refinement.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    summary=dict(structures=len(set(r['structure'] for r in rows)),directional_checks=len(rows),
        finest_max_curvature_error=max(r['absolute_curvature_error'] for r in rows if r['step_A']==3e-5),
        smallest_finest_checked_curvature=min(r['force_difference_curvature_eV_A2'] for r in rows if r['step_A']==3e-5 and r['mode']=='lowest'),
        scope='fixed-cell Gamma mechanical curvature; no material or kinetic certification')
    (args.output/'mode_refinement_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--states',type=Path,required=True);p.add_argument('--hessians',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
