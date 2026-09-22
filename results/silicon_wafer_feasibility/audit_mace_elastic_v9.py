"""Independent energy/stress tangents of affine and internally relaxed Si.

Static neutral MACE reference at its own cubic equilibrium. No experimental
temperature matching, fitted modulus, kinetic mobility, or production change.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
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
    class CountedMACE(MACECalculator):
        def __init__(self,**kwargs):
            self.calculator_calls=0
            super().__init__(**kwargs)
        def calculate(self,*call_args,**call_kwargs):
            self.calculator_calls+=1
            return super().calculate(*call_args,**call_kwargs)
    model=CountedMACE(model_paths=str(args.model),device='cpu',default_dtype='float64')
    force=force_calculator(model)
    reference=json.loads(args.bulk.read_text())
    a=float(reference['a_A']);base=bulk('Si','diamond',a=a,cubic=True);base.calc=model
    e0=float(base.get_potential_energy());v0=base.get_volume();s0=base.get_stress(voigt=False)
    conversion=electron_volt/1e-21
    modes={'hydro':np.eye(3),'tetragonal':np.diag([1.,-1.,0.]),
           'engineering_xy_shear':np.array([[0.,1.,0.],[0.,0.,0.],[0.,0.,0.]])}
    records=[];tangents=[];started=time.perf_counter();standard_evaluations=1
    for name,direction in modes.items():
        for step in [3e-3,1e-3,3e-4]:
            pair={'affine':[],'relaxed':[]}
            for sign in [-1,1]:
                deformation=np.eye(3)+sign*step*direction
                atoms=base.copy();atoms.set_cell(base.cell.array@deformation.T,scale_atoms=True)
                for treatment in ['affine','relaxed']:
                    tag=f'{name}_{step:.0e}_{sign:+d}_{treatment}'
                    if treatment=='relaxed':
                        atoms.calc=force
                        opt=BFGSLineSearch(atoms,logfile=str(args.output/(tag+'.log')),maxstep=.03)
                        ok=opt.run(fmax=1e-8,steps=100)
                        if not ok:raise RuntimeError('internal relaxation did not meet force tolerance: '+tag)
                        steps=opt.nsteps
                    else:steps=0
                    atoms.calc=model
                    energy=float(atoms.get_potential_energy());f=atoms.get_forces();sigma=atoms.get_stress(voigt=False)
                    standard_evaluations+=1
                    # P = J sigma F^{-T}; P:D is conjugate to the path parameter.
                    piola=np.linalg.det(deformation)*sigma@np.linalg.inv(deformation).T
                    work_stress=float(np.sum(piola*direction))
                    row=dict(mode=name,strain_step=step,sign=sign,treatment=treatment,
                        energy_eV=energy,conjugate_stress_eV_A3=work_stress,
                        force_max_eV_A=float(np.max(np.linalg.norm(f,axis=1))),optimizer_steps=steps,
                        volume_A3=atoms.get_volume(),energy_reference_eV=e0,reference_volume_A3=v0)
                    pair[treatment].append(row);records.append(row)
                    write(args.output/(tag+'.extxyz'),atoms)
            for treatment,r in pair.items():
                minus,plus=r
                ce=(plus['energy_eV']+minus['energy_eV']-2*e0)/(v0*step**2)*conversion
                cs=(plus['conjugate_stress_eV_A3']-minus['conjugate_stress_eV_A3'])/(2*step)*conversion
                tangents.append(dict(mode=name,strain_step=step,treatment=treatment,
                    energy_tangent_GPa=ce,stress_tangent_GPa=cs,difference_GPa=ce-cs))
            print('TANGENT',name,step,[(r['treatment'],r['energy_tangent_GPa']) for r in tangents[-2:]],flush=True)
    constants=[]
    for step in [3e-3,1e-3,3e-4]:
        for treatment in ['affine','relaxed']:
            for method in ['energy','stress']:
                t={r['mode']:r[method+'_tangent_GPa'] for r in tangents if r['strain_step']==step and r['treatment']==treatment}
                x=t['hydro']/3;delta=t['tetragonal']/2;c44=t['engineering_xy_shear']
                c11=(x+2*delta)/3;c12=(x-delta)/3
                constants.append(dict(strain_step=step,treatment=treatment,method=method,
                    C11_GPa=c11,C12_GPa=c12,C44_GPa=c44,bulk_modulus_GPa=x/3,
                    E100_GPa=(c11-c12)*(c11+2*c12)/(c11+c12),
                    cubic_elastic_positive=bool(x>0 and delta>0 and c44>0)))
    for filename,rows in [('states.csv',records),('tangents.csv',tangents),('elastic_constants.csv',constants)]:
        with (args.output/filename).open('w',newline='',encoding='utf-8') as stream:
            w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    finest=[r for r in constants if r['strain_step']==3e-4]
    error=max(abs(r['difference_GPa']) for r in tangents if r['strain_step']==3e-4)
    if error>.02 or not all(r['cubic_elastic_positive'] for r in constants):
        raise AssertionError('elastic consistency failed; inspect saved raw states')
    summary=dict(status='static neutral ML elasticity consistency checked; no material calibration',
        model_sha256=MODEL_SHA256,bulk_reference_file=args.bulk.name,
        bulk_reference_sha256=hashlib.sha256(args.bulk.read_bytes()).hexdigest(),
        a_A=a,reference_volume_A3=v0,reference_energy_eV=e0,
        reference_stress_max_abs_GPa=float(np.max(abs(s0))*conversion),
        conversion_eV_A3_to_GPa=conversion,states=len(records),tangents=len(tangents),
        force_only_calculator_evaluations=force.evaluations,
        standard_ASE_state_assessments=standard_evaluations,
        standard_ASE_calculator_calls=model.calculator_calls,
        finest_strain_energy_stress_tangent_max_difference_GPa=error,finest_constants=finest,
        elapsed_reported_s=time.perf_counter()-started,new_DFT_runs=0,new_MD_runs=0,
        independent_exp_temperature_target_used=False,material_calibrated=False,kinetic_calibrated=False,
        scope='homogeneous cubic static tangent with periodic internal relaxation; not local cleavage or finite-T free energy')
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--bulk',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
