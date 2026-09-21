"""Prepare distinct-basin Gaussian proposals for exactly the same gap and box.

Each Hessian is a sampling preconditioner; the full unchanged SW energy remains
the Metropolis target. Different preparation basins are deliberately retained.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW
from .run_conditional_dynamics import reduction
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--minima',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('choose an empty output')
    ref=np.load(args.references/'initial.npz',allow_pickle=False)
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    states=np.load(args.minima,allow_pickle=False)
    fixed,bond=ref['fixed'],ref['bond'];free=np.flatnonzero(~fixed)
    coordinates=RelaxedCoordinates(ref['reference_positions'],fixed,bond=bond)
    p,parameter_hash=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    rows=[];positions={}
    for name in ['state00','state04','state01','state02']:
        r=states[name];gap=float(r[bond[1],1]-r[bond[0],1])
        result=model.evaluate(r);gradient,reaction=coordinates.pullback(result.gradient)
        if np.max(abs(gradient))>1e-8:raise RuntimeError('minimum is not stationary')
        h=model.hessian(r,free_atoms=free)
        one,_,_=reduction(h,r,fixed,bond,(1,))
        mean=coordinates.encode(r)
        np.savez_compressed(args.output/(name+'.npz'),positions=r,fixed=fixed,bond=bond,
            reference_positions=ref['reference_positions'],mean=mean,gap_A=gap,
            bath_hessian=one.bath_hessian,response=one.bath_response[:,0],
            cross=one.cross[:,0],bare_curvature=one.bare_curvature,
            relaxed_curvature=one.relaxed_curvature,sites=result.site_energy,
            full_hessian=np.empty((0,0)))
        positions[name]=r
        row=dict(name=name,gap_A=gap,bath_dimension=len(mean),
            energy_difference_eV=float(np.sum(result.site_energy-ref['sites'])),
            static_reaction_eV_A=float(reaction[0]),logdet_bath=one.logdet_bath,
            force_residual_eV_A=float(np.max(abs(gradient))),
            bath_mean_max_coordinate_A=float(abs(mean).max()))
        rows.append(row);print(row,flush=True)
    np.savez_compressed(args.output/'branch_positions.npz',**positions,
        reference=ref['reference_positions'],fixed=fixed,bond=bond)
    save_json(args.output/'summary.json',dict(rows=rows,geometry=meta['geometry'],
        source_parameter_sha256=parameter_hash,
        source_minima_sha256=hashlib.sha256(args.minima.read_bytes()).hexdigest(),
        source_reference_summary_sha256=hashlib.sha256((args.references/'summary.json').read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        purpose='independent preparations for the same finite-box target; no claim of global mixing'))


if __name__=='__main__':main()
