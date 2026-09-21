"""Analyze actual finite-box conditional ensembles; do not invent convergence."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline, CubicHermiteSpline
from scipy.optimize import brentq

from solver_v1.silicon_thermal_research import split_rhat, block_statistics
from .run_local_crack_audit import save_json


def profile_from_forces(gaps, forces, anchor):
    interpolator = PchipInterpolator(gaps, forces, extrapolate=False)
    primitive = interpolator.antiderivative()
    roots = interpolator.roots(extrapolate=False)
    roots = roots[np.isfinite(roots)]
    roots = roots[(roots >= min(gaps)) & (roots <= max(gaps))]
    stationary = [dict(gap_A=float(q), free_energy_eV=float(primitive(q)-primitive(anchor)),
        mean_force_derivative_eV_A2=float(interpolator.derivative()(q)),
        kind='minimum' if interpolator.derivative()(q)>0 else 'maximum') for q in roots]
    return dict(stationary=stationary,
        free_energy_at_samples_eV=(primitive(gaps)-primitive(anchor)).tolist()), interpolator


def profile_with_static_control(gaps,forces,static_energies,static_forces,anchor):
    """Integrate only thermal excess force; retain known static energies/jets.

    The exact computed static branch is interpolated with its own derivative.
    This removes the avoidable error from integrating a large known 0 K force
    on a sparse thermal grid. Interpolation is still checked independently.
    """
    gaps,forces,static_energies,static_forces=map(lambda a:np.asarray(a,float),
        (gaps,forces,static_energies,static_forces))
    mechanical=CubicHermiteSpline(gaps,static_energies,static_forces,extrapolate=False)
    excess=PchipInterpolator(gaps,forces-static_forces,extrapolate=False)
    integral=excess.antiderivative()
    free=lambda q:mechanical(q)-mechanical(anchor)+integral(q)-integral(anchor)
    force=lambda q:mechanical.derivative()(q)+excess(q)
    curvature=lambda q:mechanical.derivative(2)(q)+excess.derivative()(q)
    dense=np.linspace(gaps[0],gaps[-1],2001);sampled=force(dense)
    roots=list(dense[sampled == 0])
    for i in range(len(dense)-1):
        if sampled[i]*sampled[i+1]<0:
            roots.append(brentq(force,dense[i],dense[i+1],xtol=1e-13))
    roots=sorted(roots)
    roots=[q for i,q in enumerate(roots) if i==0 or abs(q-roots[i-1])>1e-10]
    stationary=[dict(gap_A=float(q),free_energy_eV=float(free(q)),
        mean_force_derivative_eV_A2=float(curvature(q)),kind='minimum' if curvature(q)>0 else 'maximum') for q in roots]
    cubic=CubicSpline(gaps,forces-static_forces).antiderivative()
    difference=(integral(dense)-integral(anchor))-(cubic(dense)-cubic(anchor))
    report=dict(stationary=stationary,free_energy_at_samples_eV=free(gaps).tolist(),
        excess_force_pchip_cubic_max_difference_eV=float(np.max(abs(difference))))
    return report,free,force


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, required=True)
    parser.add_argument('--allow-incomplete',action='store_true',help='analyze only saved complete replica groups; retain incomplete status')
    args = parser.parse_args()
    out = args.calculation
    metadata_path=out/'summary.json'
    if not metadata_path.exists() and args.allow_incomplete:metadata_path=out/'interrupted_summary.json'
    meta = json.loads(metadata_path.read_text(encoding='utf-8'))
    records = meta['records']
    names = sorted(set(row['reference'] for row in records),
        key=lambda name:next(row['gap_A'] for row in records if row['reference']==name))
    temperatures = sorted(set(row['temperature_K'] for row in records))
    all_rows, profiles = [], {}
    for temperature in temperatures:
        temperature_rows, forces_by_seed = [], {}
        for name in names:
            runs = sorted([row for row in records if row['temperature_K']==temperature and row['reference']==name],
                          key=lambda row:row['seed'])
            chains = np.stack([np.load(out/(row['label']+'.npz'),allow_pickle=False)['observations'] for row in runs])
            rhat = split_rhat(chains)
            means = chains.mean(axis=1)
            sem = means.std(axis=0,ddof=1)/np.sqrt(len(runs))
            row = dict(temperature_K=temperature,name=name,gap_A=runs[0]['gap_A'],
                replicas=len(runs),draws_per_replica=chains.shape[1],
                raw_force_eV_A=float(means[:,0].mean()), controlled_force_eV_A=float(means[:,1].mean()),
                raw_replica_sem_eV_A=float(sem[0]),controlled_replica_sem_eV_A=float(sem[1]),
                raw_split_rhat=float(rhat[0]),controlled_split_rhat=float(rhat[1]),
                anharmonic_split_rhat=float(rhat[3]),
                mean_anharmonic_residual_kBT=float(means[:,3].mean()),
                maximum_bath_coordinate_A=float(chains[:,:,4].max()),
                raw_control_difference_eV_A=float((means[:,0]-means[:,1]).mean()),
                minimum_acceptance=float(min(row['acceptance_fraction'] for row in runs)),
                control_variance_reduction=float(np.var(chains[:,:,0])/np.var(chains[:,:,1])),
                block_sem_16=[np.asarray(block_statistics(chain,blocks=16)['block_standard_error']).tolist() for chain in chains],
                block_sem_32=[np.asarray(block_statistics(chain,blocks=32)['block_standard_error']).tolist() for chain in chains],
                replica_means=means.tolist())
            temperature_rows.append(row); all_rows.append(row)
            for run, mean in zip(runs,means):
                forces_by_seed.setdefault(str(run['seed']),[]).append(float(mean[1]))
        gaps = np.array([row['gap_A'] for row in temperature_rows])
        if len(gaps) >= 5:
            anchor = next(row['gap_A'] for row in records if row['reference']=='initial')
            forces = np.array([row['controlled_force_eV_A'] for row in temperature_rows])
            profile, interpolator = profile_from_forces(gaps,forces,anchor)
            profile['replicas'] = {seed:profile_from_forces(gaps,f,anchor)[0] for seed,f in forces_by_seed.items()}
            cube = CubicSpline(gaps,forces).antiderivative()
            dense = np.linspace(gaps[0],gaps[-1],501)
            primitive = interpolator.antiderivative()
            profile['pchip_cubic_max_integral_difference_eV'] = float(np.max(abs(
                (primitive(dense)-primitive(anchor))-(cube(dense)-cube(anchor)))))
            reference_rows = {r['name']:r for r in meta['reference_summary']['rows']}
            static_forces = np.array([reference_rows[name]['static_reaction_eV_A'] for name in names])
            static, static_interpolator = profile_from_forces(gaps,static_forces,anchor)
            exact_energy = np.array([reference_rows[name]['energy_difference_eV'] for name in names])
            profile['static_integral_max_error_at_samples_eV'] = float(np.max(abs(
                np.array(static['free_energy_at_samples_eV'])-exact_energy)))
            corrected,_,_=profile_with_static_control(gaps,forces,exact_energy,static_forces,anchor)
            corrected['replicas']={seed:profile_with_static_control(gaps,f,exact_energy,static_forces,anchor)[0]
                                  for seed,f in forces_by_seed.items()}
            profile['static_controlled_integration']=corrected
            profiles[str(round(temperature))] = profile
    result = dict(rows=all_rows,profiles=profiles,
        interpretation='classical finite fixed-bath-box SW ensemble estimates; incomplete branch exploration and quantum validation',
        uncertainty='replica spread and batch diagnostics; split-Rhat is classical, not rank-normalized; no automatic material gate',
        finite_T_global_PMF_certified=False,physical_mobility=None,production_t0_seconds=None,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        computation_completed=meta.get('computation_completed',False),
        input_summary_file=metadata_path.name,
        input_summary_sha256=hashlib.sha256(metadata_path.read_bytes()).hexdigest())
    save_json(out/'analysis.json',result)
    scalar_keys=[key for key in all_rows[0] if key not in ('block_sem_16','block_sem_32','replica_means')]
    with (out/'mean_forces.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=scalar_keys); writer.writeheader()
        writer.writerows({key:row[key] for key in scalar_keys} for row in all_rows)
    for temperature in temperatures:
        rows=[r for r in all_rows if r['temperature_K']==temperature]
        print('T',temperature,'max force Rhat',max(r['controlled_split_rhat'] for r in rows),
            'max residual Rhat',max(r['anharmonic_split_rhat'] for r in rows),
            'max force SEM',max(r['controlled_replica_sem_eV_A'] for r in rows),flush=True)
    print(json.dumps(profiles,indent=2),flush=True)


if __name__=='__main__':
    main()
