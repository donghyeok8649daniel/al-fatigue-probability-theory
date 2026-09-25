"""Independent raw comparison of two same-load, different-history states.

Zero axial force does not force equal grip extensions. This is not a physical
fatigue cycle, residual plasticity claim, or free-energy/kinetic calibration.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_structure_diagnostics_v12 import pair_set,edge_bottleneck

def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    with np.load(args.geometry) as d:
        ref=d['reference'];free=d['free'];lower=d['lower'];upper=d['upper'];area=float(d['area_A2']);gauge=float(d['gauge_length_A'])
    snapshots=[];observations=[];manifest=[]
    for name,folder in [('original',args.baseline),('after_rearrangement',args.returned)]:
        stored=json.loads((folder/'result.json').read_text(encoding='utf-8'))
        if abs(stored['target_nominal_GPa'])>1e-14:raise ValueError('zero external axial load required')
        if not stored['converged'] or stored['active_guard_bound']:
            raise ValueError('unconverged or guard-limited state is not a zero-force equilibrium')
        protocol=json.loads((folder.parent/'protocol.json').read_text(encoding='utf-8'))
        with np.load(folder/'raw.npz') as d:
            r=d['positions'].copy();f=d['forces'].copy();energy=float(d['energy']);numbers=d['numbers'].copy()
        if not np.all(numbers==14):raise ValueError('same pure Si specimen required')
        extension=float((r[upper,2]-ref[upper,2]).mean()-(r[lower,2]-ref[lower,2]).mean())
        expected=ref.copy();expected[lower,2]-=extension/2;expected[upper,2]+=extension/2
        grip_error=float(np.max(abs(r[~free]-expected[~free])))
        conjugate=float((f[lower,2].sum()-f[upper,2].sum())/2)
        values=dict(energy_eV=energy,extension_A=extension,nominal_stress_GPa=conjugate/area*160.2176634,
            free_force_max_eV_A=float(np.linalg.norm(f[free],axis=1).max()),
            reaction_force_residual_eV_A=float(np.linalg.norm(f[~free].sum(axis=0))),
            reaction_torque_eV=float(np.linalg.norm(np.cross(r[~free],f[~free]).sum(axis=0))))
        error=max(abs(v-stored[k]) for k,v in values.items())
        if error>1e-10 or grip_error>1e-12:raise ValueError('stored observations or rigid constraints do not replay')
        if values['free_force_max_eV_A']>protocol['fmax_eV_A'] or abs(conjugate)>protocol['fmax_eV_A']:
            raise ValueError('independent forces fail the declared zero-force convergence tolerance')
        observations.append(dict(role=name,converged=stored['converged'],**values,raw_replay_error=error,grip_error_A=grip_error))
        snapshots.append(r);manifest.append(dict(role=name,raw_sha256=hashlib.sha256((folder/'raw.npz').read_bytes()).hexdigest(),
            result_sha256=hashlib.sha256((folder/'result.json').read_bytes()).hexdigest()))
    a,b=snapshots;design=np.column_stack((a[free],np.ones(free.sum())))
    affine=np.linalg.lstsq(design,b[free],rcond=None)[0];residual=b[free]-design@affine
    i,j=linear_sum_assignment(cdist(a[free],b[free],metric='sqeuclidean'))
    matched=np.linalg.norm(a[free][i]-b[free][j],axis=1)
    graphs=[]
    for cutoff in (2.7,2.8,3.1,3.4,3.7):
        old=pair_set(a,cutoff);new=pair_set(b,cutoff)
        graph=edge_bottleneck(len(a),new,lower,upper)
        graphs.append(dict(cutoff_A=cutoff,original_pairs=len(old),returned_pairs=len(new),
            lost_pairs=len(old-new),new_pairs=len(new-old),returned_edge_disjoint_paths=graph['edge_disjoint_paths']))
    result=dict(complete=True,states=observations,geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),source_states=manifest,
        same_external_load=True,same_boundary_displacement=False,
        grip_extension_difference_A=observations[1]['extension_A']-observations[0]['extension_A'],
        extension_difference_over_original_force_free_gauge=(observations[1]['extension_A']-observations[0]['extension_A'])/(gauge+observations[0]['extension_A']),
        energy_difference_at_zero_external_force_eV=observations[1]['energy_eV']-observations[0]['energy_eV'],
        labelled_free_RMS_A=float(np.sqrt(np.mean(np.sum((b[free]-a[free])**2,axis=1)))),
        permutation_matching_RMS_A=float(np.sqrt(np.mean(matched**2))),permuted_atoms=int(np.sum(i!=j)),
        best_affine_removed_RMS_A=float(np.sqrt(np.mean(np.sum(residual**2,axis=1)))),
        graph_diagnostics=graphs,new_model_calls=0,new_DFT=0,new_MD=0,
        actual_plasticity_certified=False,first_crack_certified=False,physical_clock=None,
        scope='static history dependence after zero axial load; no thermal hold, physical dissipation or initiation classification')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('geometry','baseline','returned','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
