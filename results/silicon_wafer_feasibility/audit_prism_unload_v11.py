"""Compare converged loading/unloading structures at the same imposed grips.

This measures a static preparation dependence. It supplies no physical time,
plastic-flow label, dissipated work, endurance limit, or crack-initiation event.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from audit_prism_continuation_v11 import all_pairs


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    protocol=json.loads((args.unload/'protocol.json').read_text(encoding='utf-8'))
    unloaded=json.loads((args.unload/'summary.json').read_text(encoding='utf-8'))
    baseline=json.loads((args.baseline/'result.json').read_text(encoding='utf-8'))
    if not baseline['converged']:raise ValueError('force-converged loading baseline required')
    if abs(unloaded['strain']-baseline['strain'])>1e-12:raise ValueError('same target strain required')
    if hashlib.sha256(args.geometry.read_bytes()).hexdigest()!=protocol['geometry_sha256']:
        raise ValueError('geometry source changed')
    with np.load(args.geometry) as d:
        reference=d['reference'];lower=d['lower'];upper=d['upper'];free=d['free'];area=float(d['area_A2'])
    with np.load(args.baseline/'raw.npz') as d:r0=d['positions'];e0=float(d['energy']);n0=d['numbers']
    with np.load(args.unload/'raw.npz') as d:r1=d['positions'];e1=float(d['energy']);force=d['forces'];n1=d['numbers']
    if not np.array_equal(n0,n1) or not np.all(n0==14):raise ValueError('same pure-Si body required')
    fixed=lower|upper
    grip_error=float(abs(r1[fixed]-r0[fixed]).max())
    if grip_error>1e-12:raise ValueError('loading and unloading boundary coordinates differ')
    target=reference.copy();target[lower,2]-=unloaded['extension_A']/2;target[upper,2]+=unloaded['extension_A']/2
    if np.max(abs(r1[fixed]-target[fixed]))>1e-12:raise ValueError('unload target grips differ')
    axial=float((force[lower,2].sum()-force[upper,2].sum())/2)
    support=-force[fixed].sum(axis=0);torque=-np.cross(r1[fixed],force[fixed]).sum(axis=0)
    expected=dict(energy_eV=e1,conjugate_force_eV_A=axial,nominal_stress_GPa=axial/area*160.2176634,
        free_force_max_eV_A=float(np.linalg.norm(force[free],axis=1).max()),
        reaction_force_residual_eV_A=float(np.linalg.norm(support)),reaction_torque_eV=float(np.linalg.norm(torque)),
        total_internal_force_eV_A=float(np.linalg.norm(force.sum(axis=0))),
        total_internal_torque_eV=float(np.linalg.norm(np.cross(r1,force).sum(axis=0))))
    replay_error=max(abs(expected[k]-unloaded[k]) for k in expected)
    if replay_error>1e-11 or abs(e0-baseline['energy_eV'])>1e-11:raise ValueError('raw observable replay failed')
    # Compare labelled positions and geometry modulo permutations of identical
    # free Si atoms. The latter is a descriptor, not a dynamical return test.
    differences=np.linalg.norm(r1[free]-r0[free],axis=1)
    cost=cdist(r0[free],r1[free],metric='sqeuclidean')
    i,j=linear_sum_assignment(cost)
    matched=np.linalg.norm(r0[free][i]-r1[free][j],axis=1)
    if not np.isclose(float(cost[i,j].sum()),float(matched@matched),rtol=1e-13,atol=1e-13):
        raise ValueError('permutation-matching cost replay failed')
    graphs=[]
    for row in unloaded['connectivity']:
        cutoff=row['cutoff_A'];old=all_pairs(r0,cutoff);new=all_pairs(r1,cutoff)
        if len(new)!=row['pair_count']:raise ValueError('direct pair count differs from stored result')
        graphs.append(dict(cutoff_A=cutoff,loading_pairs=len(old),unloading_pairs=len(new),
            lost_relative_to_loading=len(old-new),formed_relative_to_loading=len(new-old),
            retained_fraction=len(old&new)/len(old)))
    result=dict(loading_converged=True,unloading_converged=unloaded['converged'],
        interruption=unloaded['interruption'],target_strain=unloaded['strain'],
        same_grip_coordinates_error_A=grip_error,maximum_observable_replay_difference=replay_error,
        loading_energy_eV=e0,unloading_energy_eV=e1,energy_difference_at_same_grips_eV=e1-e0,
        loading_stress_GPa=baseline['nominal_stress_GPa'],unloading_stress_GPa=unloaded['nominal_stress_GPa'],
        loading_free_force_max_eV_A=baseline['free_force_max_eV_A'],
        unloading_free_force_max_eV_A=unloaded['free_force_max_eV_A'],
        labelled_free_atom_RMS_difference_A=float(np.sqrt(np.mean(differences**2))),
        labelled_free_atom_max_difference_A=float(differences.max()),
        minimum_matching_RMS_difference_A=float(np.sqrt(np.mean(matched**2))),
        minimum_matching_max_difference_A=float(matched.max()),
        minimum_matching_permuted_atoms=int(np.count_nonzero(i!=j)),
        permutation_invariant_matching='minimum total squared distance, identical free Si atoms only; fixed grips excluded',
        force_identity_error_eV_A=float(np.linalg.norm(support-force[free].sum(axis=0))),
        torque_identity_error_eV=float(np.linalg.norm(torque-np.cross(r1[free],force[free]).sum(axis=0))),
        pair_diagnostics=graphs,new_potential_calls=0,new_DFT=0,new_MD=0,
        physical_dissipated_work=None,physical_clock=None,initiation_probability=None,
        scope='static structural preparation dependence at identical boundary displacement; not a committor, crack label, or physical fatigue cycle',
        sources=[dict(role=role,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for role,path in
            [('geometry',args.geometry),('baseline_raw',args.baseline/'raw.npz'),('baseline_result',args.baseline/'result.json'),
             ('unload_raw',args.unload/'raw.npz'),('unload_summary',args.unload/'summary.json'),('unload_protocol',args.unload/'protocol.json')]])
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('baseline','unload','geometry','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
