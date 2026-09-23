"""Read a curvature run and distinguish dU from d(U-F*Delta).

The existing finite-difference Hessian of U is also the dead-load enthalpy
Hessian because F*Delta is linear. Gradients and odd energy differences are
different and must include the external force. No Hessian/model call is repeated.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np


def main(args):
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    protocol=json.loads((args.curvature/'protocol.json').read_text(encoding='utf-8'))
    result=json.loads((args.curvature/'summary.json').read_text(encoding='utf-8'))
    state=json.loads((args.state/'result.json').read_text(encoding='utf-8'))
    if protocol['ensemble']!='force' or not state['converged']:
        raise ValueError('force-ensemble curvature and converged source state required')
    if hashlib.sha256((args.state/'raw.npz').read_bytes()).hexdigest()!=protocol['state_sha256']:
        raise ValueError('curvature refers to a different source state')
    with np.load(args.curvature/'modes.npz') as data:
        grad_U=data['gradient']; values=data['values']; vectors=data['vectors']; positions=data['positions']
    with np.load(args.state/'raw.npz') as data:
        expected=data['gradient']; state_positions=data['positions']; extension=float(data['x'][-1])
    if not np.array_equal(positions,state_positions):
        raise ValueError('curvature source coordinates changed')
    load=float(state['target_force_eV_A']); grad_H=grad_U.copy(); grad_H[-1]-=load
    difference=float(np.max(abs(grad_H-expected)))
    if difference>1e-8:
        raise ValueError('fresh curvature gradient differs from force-state gradient')
    rows=[]
    for row in result['modes']:
        i=row['mode']; mode=vectors[:,i]
        projected_U=float(grad_U@mode); projected_H=float(grad_H@mode)
        if abs(projected_U-row['projected_gradient_eV_A'])>1e-10:
            raise ValueError('projected U-gradient replay failed')
        checks=[]
        for check in row['energy_checks']:
            slope_H=check['odd_slope_eV_A']-load*mode[-1]
            checks.append(dict(amplitude_A=check['amplitude_A'],
                even_curvature_U_equals_H_eV_A2=check['even_curvature_eV_A2'],
                odd_slope_U_eV_A=check['odd_slope_eV_A'], odd_slope_enthalpy_eV_A=float(slope_H),
                odd_slope_minus_projected_enthalpy_gradient_eV_A=float(slope_H-projected_H)))
        rows.append(dict(mode=i,ritz_curvature_eV_A2=float(values[i]),
            half_step_rayleigh_eV_A2=row['half_step_rayleigh_eV_A2'],
            residual_norm_eV_A2=row['residual_norm_eV_A2'],
            half_step_residual_eV_A2=row['half_step_residual_eV_A2'],
            projected_U_gradient_eV_A=projected_U,projected_enthalpy_gradient_eV_A=projected_H,
            grip_extension_mode_component=float(mode[-1]), energy_checks=checks))
    audit=dict(target_force_eV_A=load,target_nominal_GPa=state['target_nominal_GPa'],
        enthalpy_eV=result['energy_eV']-load*extension,
        U_reduced_gradient_max_eV_A=float(abs(grad_U).max()),
        enthalpy_reduced_gradient_max_eV_A=float(abs(grad_H).max()),
        maximum_force_state_gradient_replay_difference_eV_A=difference,
        lanczos_converged=result['lanczos_converged'],interruption=result['interruption'],modes=rows,
        saved_eigenmodes=int(len(values)),independently_checked_modes=len(rows),
        force_calls=result['force_calls'],elapsed_s=result['elapsed_s'],
        krylov_dimension=protocol['krylov_dimension'],coordinates=protocol['coordinates'],
        source_force_tolerance_note='source state identified by raw SHA; do not transfer stability to another relaxed state',
        linear_dead_load_has_zero_Hessian=True,
        original_curvature_result_preserved=True,
        distinction='original reduced_gradient and odd_slope refer to U; force-control stationarity uses H=U-F*Delta',
        new_potential_calls=0,new_DFT=0,new_MD=0,global_minimum_certified=False,initiation_probability=None,
        source_manifest=[dict(file=name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for name,path in
            [('curvature/protocol.json',args.curvature/'protocol.json'),('curvature/summary.json',args.curvature/'summary.json'),
             ('curvature/modes.npz',args.curvature/'modes.npz'),('state/raw.npz',args.state/'raw.npz'),
             ('state/result.json',args.state/'result.json')]])
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--curvature',type=Path,required=True)
    parser.add_argument('--state',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
