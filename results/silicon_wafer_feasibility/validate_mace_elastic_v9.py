"""Independent uniaxial, mixed-normal and rotated-shear elastic controls."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.constants import electron_volt
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_calculator


def main(args):
    import torch
    from ase.build import bulk
    from ase.io import write
    from ase.optimize import BFGSLineSearch
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model mismatch')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    calc=force_calculator(model)
    report=json.loads(args.reference.read_text())
    base=bulk('Si','diamond',a=report['a_A'],cubic=True);base.calc=calc
    e0=float(base.get_potential_energy());volume=base.get_volume();step=3e-4
    rows=[]
    def evaluate(strain,label):
        a=base.copy();a.set_cell(base.cell.array@(np.eye(3)+strain).T,scale_atoms=True);a.calc=calc
        opt=BFGSLineSearch(a,logfile=str(args.output/(label+'.log')),maxstep=.03)
        ok=opt.run(fmax=1e-8,steps=100)
        if not ok:raise RuntimeError('force failed')
        energy=float(a.get_potential_energy())
        rows.append(dict(label=label,strain=strain.tolist(),energy_eV=energy,
            max_force_eV_A=float(np.max(np.linalg.norm(a.get_forces(),axis=1))),steps=opt.nsteps))
        write(args.output/(label+'.extxyz'),a)
        return energy
    x=np.diag([1.,0.,0.]);y=np.diag([0.,1.,0.]);shear=np.zeros((3,3));shear[1,2]=1.
    one=[evaluate(s*step*x,f'x_{s:+d}') for s in [-1,1]]
    mixed={(i,j):evaluate(step*(i*x+j*y),f'xy_{i:+d}_{j:+d}') for i in [-1,1] for j in [-1,1]}
    rotated=[evaluate(s*step*shear,f'yz_shear_{s:+d}') for s in [-1,1]]
    factor=electron_volt/1e-21/(volume*step**2)
    values=dict(C11_GPa=(sum(one)-2*e0)*factor,
        C12_GPa=(mixed[1,1]-mixed[1,-1]-mixed[-1,1]+mixed[-1,-1])*factor/4,
        C44_GPa=(sum(rotated)-2*e0)*factor)
    primary=next(r for r in report['finest_constants'] if r['treatment']=='relaxed' and r['method']=='energy')
    difference={k:values[k]-primary[k] for k in values}
    if max(abs(v) for v in difference.values())>.001:raise AssertionError('independent mode mapping failed')
    result=dict(model_sha256=MODEL_SHA256,reference_file=args.reference.name,
        reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        a_A=report['a_A'],strain_step=step,method='uniaxial E second difference, mixed normal E cross difference, rotated engineering yz shear',
        independent_constants=values,differences_from_hydro_tetragonal_xy_GPa=difference,
        states=rows,calculator_evaluations=calc.evaluations,new_DFT_runs=0,new_MD_runs=0,
        material_calibrated=False,kinetic_calibrated=False)
    (args.output/'validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(independent_constants=values,differences_GPa=difference,evaluations=calc.evaluations)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
