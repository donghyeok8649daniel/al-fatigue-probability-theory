"""Direct energy curvature for the computed lowest Cartesian Gamma modes."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.run_concentration_cleavage_v10 import dump,csv_write


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model identity')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    raw=np.load(args.modes)
    if len(raw['eigenvalues'])!=2:raise ValueError('two converged modes required')
    atoms=Atoms(numbers=raw['numbers'],positions=raw['positions'],cell=raw['cell'],pbc=True)
    e0,f0=force_only(calc,atoms);rows=[]
    for i,value in enumerate(raw['eigenvalues']):
        direction=raw['eigenvectors'][:,i].reshape(-1,3)
        for step in [1e-2,3e-3,1e-3]:
            energies=[];forces=[]
            for sign in [-1,1]:
                trial=atoms.copy();trial.positions+=sign*step*direction
                e,f=force_only(calc,trial);energies.append(e);forces.append(f)
            curvature=(energies[1]+energies[0]-2*e0)/step**2
            force_curvature=-np.sum(direction*(forces[1]-forces[0]))/(2*step)
            rows.append(dict(mode=i,step_A=step,eigenvalue_eV_A2=float(value),
                energy_reference_eV=e0,energy_minus_eV=energies[0],energy_plus_eV=energies[1],
                energy_curvature_eV_A2=curvature,force_curvature_eV_A2=force_curvature,
                energy_error_eV_A2=curvature-value,force_error_eV_A2=force_curvature-value))
    csv_write(args.output/'mode_curvatures.csv',rows)
    error=float(max(abs(r['energy_error_eV_A2']) for r in rows if r['step_A']==1e-3))
    passed=bool(error<1e-4)
    dump(args.output/'summary.json',dict(passed=passed,force_only_calculator_calls=13,
        finest_energy_curvature_error_eV_A2=error,modes_sha256=hashlib.sha256(args.modes.read_bytes()).hexdigest(),
        model_sha256=MODEL_SHA256,scope='independent energy differences along selected Gamma directions; not full-q stability'))
    if not passed:raise AssertionError('independent mode energy curvature mismatch')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--modes',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
