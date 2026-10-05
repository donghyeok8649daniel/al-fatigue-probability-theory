"""Independently replay the bounded atomistic pCN pilot and thermal identity.

No new potential calls or fracture counts. Diagnostics of short correlated
chains are not convergence certification or a measurement of physical time.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import Boltzmann, electron_volt

from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_thermal_identity import canonical_temperature_terms
from solver_v1.silicon_thermal_research import block_statistics, split_rhat


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_arrays(path):
    with np.load(path) as data:
        return {key:data[key].copy() for key in data.files}


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def replay_chain(root, pilot, chain, source, target, protocol, common, direction):
    folder = pilot/f'chain_{chain}'
    files = sorted(folder.glob('evaluation_*.npz'))
    expected_count = 1+protocol['warmup']+protocol['draws']
    if len(files) != expected_count:
        raise ValueError('incomplete chain evaluation record')
    evaluations = [read_arrays(path) for path in files]
    for index, path in enumerate(files):
        if path.name != f'evaluation_{index:03d}.npz':
            raise ValueError('missing or duplicate proposal evaluation')
    rng = np.random.default_rng(protocol['seeds'][chain])
    current = target.reference.whiten(np.zeros(target.basis.shape[1]))
    retained, accepted, observation, potential, retained_indices = [], [], [], [], []
    current_index, accepted_warmup = 0, 0
    persistence = np.sqrt(1-protocol['proposal_scale']**2)
    max_position_error, max_correction_error, max_proposal_error = 0., 0., 0.
    for index, value in enumerate(evaluations):
        u = current if index == 0 else persistence*current+protocol['proposal_scale']*rng.standard_normal(len(current))
        max_proposal_error = max(max_proposal_error, float(np.max(abs(u-value['u']))))
        positions = target.positions+(target.basis@(target.reference.mean+np.sqrt(Boltzmann/electron_volt*300.)
                        *np.linalg.solve(target.reference.lower.T,u))).reshape(target.positions.shape)
        max_position_error = max(max_position_error, float(np.max(abs(positions-value['positions']))))
        inside = target.inside(positions)
        if inside != bool(value['inside']):
            raise ValueError('domain membership changed')
        expected_phi = ((float(value['energy'])-float(source['energy']))
                        /(Boltzmann/electron_volt*300.)-.5*u@u) if inside else np.inf
        if inside:
            max_correction_error = max(max_correction_error, abs(expected_phi-float(value['correction'])))
        elif not np.isposinf(value['correction']):
            raise ValueError('outside point has a finite correction')
        if index == 0:
            if (abs(float(value['energy'])-float(source['energy'])) > 1e-8
                    or np.max(abs(value['forces']-source['forces'])) > 1e-8):
                raise ValueError('initial model replay mismatch')
            continue
        accept = bool(np.log(rng.random()) < min(0., float(evaluations[current_index]['correction'])-expected_phi))
        if accept:
            current, current_index = u.copy(), index
        if index <= protocol['warmup']:
            accepted_warmup += int(accept)
        else:
            retained.append(current.copy()); accepted.append(accept)
            observation.append(evaluations[current_index]['observation'])
            potential.append(float(evaluations[current_index]['correction']))
            retained_indices.append(current_index)
    saved = read_arrays(folder/'chain.npz')
    for key, calculated in [('samples',retained), ('observations',observation), ('potentials',potential),
                             ('accepted_at_saved_step',accepted)]:
        np.testing.assert_allclose(saved[key], calculated, rtol=0, atol=1e-10)
    metadata = json.loads((folder/'summary.json').read_text())
    if (abs(metadata['acceptance_fraction']-np.mean(accepted)) > 1e-15
            or abs(metadata['warmup_acceptance_fraction']-accepted_warmup/protocol['warmup']) > 1e-15
            or metadata['final_rng_state'] != rng.bit_generator.state):
        raise ValueError('acceptance or final RNG replay mismatch')
    if max_proposal_error > 1e-12 or max_position_error > 1e-12 or max_correction_error > 1e-10:
        raise ValueError('independent proposal/position/energy replay error')
    thermal, features, terms, own_rms = [], [], [], []
    for index in retained_indices:
        value = evaluations[index]
        r, f = value['positions'], value['forces']
        row = canonical_temperature_terms(r, f, target.free, common, source['positions'],
            temperature_K=300., halfwidth_A=target.halfwidth,
            minimum_pair_A=target.minimum_pair, taper_width_A=.4)
        thermal.append([row['numerator_eV'],row['denominator'],row['residual_eV']])
        terms.append(dict(proposal_evaluation=index, pair_taper=row['pair_taper'],
                          active_pair_tapers=row['active_pair_tapers']))
        projection = float(((r-common).ravel()@direction))
        delta = r[target.free]-source['positions'][target.free]
        own_rms.append(float(np.sqrt(np.mean(np.sum(delta*delta,axis=1)))))
        common_delta = r[target.free]-common[target.free]
        features.append([float(value['energy']), projection,
                         float(np.sqrt(np.mean(np.sum(common_delta*common_delta,axis=1))))])
    thermal, features = np.asarray(thermal), np.asarray(features)
    blocks = block_statistics(thermal, blocks=8)
    numerator, denominator, residual = thermal.mean(axis=0)
    result = dict(chain=chain, evaluation_count=len(files), saved_correlated_draws=len(retained),
        unique_retained_evaluations=len(set(retained_indices)), acceptance_fraction=float(np.mean(accepted)),
        warmup_acceptance_fraction=accepted_warmup/protocol['warmup'],
        maximum_proposal_error=max_proposal_error, maximum_position_error_A=max_position_error,
        maximum_correction_error=max_correction_error,
        mean_energy_eV=float(features[:,0].mean()), mean_projection_A=float(features[:,1].mean()),
        mean_rms_from_own_source_A=float(np.mean(own_rms)),
        mean_rms_from_common_center_A=float(features[:,2].mean()),
        temperature_identity_mean_numerator_eV=float(numerator),
        temperature_identity_mean_denominator=float(denominator),
        temperature_ratio_diagnostic_K=float(numerator/denominator/(Boltzmann/electron_volt)) if denominator > 0 else None,
        mean_identity_residual_eV=float(residual),
        block_residual_standard_error_diagnostic_eV=float(blocks['block_standard_error'][2]),
        first_half_residual_eV=float(blocks['first_half_mean'][2]),
        second_half_residual_eV=float(blocks['second_half_mean'][2]),
        no_confidence_interval_claim=True, actual_physical_temperature_inferred=False,
        thermal_equilibrium_certified=False, physical_time_ps=None,
        retained_evaluation_indices=retained_indices, pair_boundary_controls=terms)
    return result, thermal, features, files


def main(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    protocol = json.loads((args.pilot/'protocol.json').read_text())
    completion = json.loads((args.pilot/'summary.json').read_text())
    if not completion['complete']:
        raise ValueError('pilot incomplete; preserve its partial record')
    if protocol['temperature_K'] != 300. or protocol['planned_proposals'] != 128:
        raise ValueError('this replay targets the declared two-chain 300 K pilot')
    source_paths = [args.root/relative for relative in protocol['source_sha256']]
    for path in source_paths:
        relative = str(path.relative_to(args.root)).replace('\\','/')
        if sha(path) != protocol['source_sha256'][relative]:
            raise ValueError('source hash changed')
    for key, relative in [('runner_sha256','results/silicon_wafer_feasibility/sample_nonlinear_ensemble_v15.py'),
                          ('target_module_sha256','solver_v1/silicon_device_ensemble.py'),
                          ('sampler_module_sha256','solver_v1/silicon_thermal_research.py')]:
        if sha(args.root/relative) != protocol[key]:
            raise ValueError('source module hash changed')
    sources = [read_arrays(path) for path in source_paths]
    common = (sources[0]['positions']+sources[1]['positions'])/2
    direction = (sources[1]['positions']-sources[0]['positions']).ravel()
    endpoint_distance = float(np.linalg.norm(direction)); direction /= endpoint_distance
    args.output.mkdir(parents=True)
    results, thermals, features, evaluations = [], [], [], []
    for chain, source in enumerate(sources):
        target = FixedGripTarget(source['positions'],source['basis'],source['hessian'],source['gradient'],300.,common)
        result, thermal, feature, files = replay_chain(args.root,args.pilot,chain,source,target,protocol,common,direction)
        results.append(result); thermals.append(thermal); features.append(feature); evaluations.extend(files)
    feature_array = np.asarray(features)
    rhat = split_rhat(feature_array)
    summary = dict(completed_utc=datetime.now(timezone.utc).isoformat(), chains=results,
        source_endpoint_distance_A=endpoint_distance, feature_names=['total_energy_eV','endpoint_projection_A','rms_from_common_center_A'],
        classical_split_rhat_diagnostic=rhat.tolist(),
        rhat_is_not_a_confidence_interval=True, full_covariance_rank_at_most=protocol['draws']-1,
        new_potential_calls=completion['new_model_calls'], new_calls_in_this_replay=0,
        new_MD=0, new_DFT=0, equilibrium_certified=False, free_energy_estimated=False,
        material_approved=False, first_crack_states_verified=False, crack_probability_estimated=False,
        physical_clock=None, production_changed=False, analysis_runner_sha256=sha(Path(__file__)),
        interpretation='short two-start canonical-model pilot; no calibration of kBT or real Si')
    np.savez_compressed(args.output/'chain_diagnostics.npz',features=feature_array,thermal_terms=np.asarray(thermals))
    write_json(args.output/'summary.json',summary)
    raw_files = sorted(args.pilot.rglob('*'))
    write_json(args.output/'input_manifest.json',{
        'pilot/'+str(path.relative_to(args.pilot)).replace('\\','/'):sha(path) for path in raw_files if path.is_file()})
    write_json(args.output/'source_manifest.json',{
        str(path.relative_to(args.root)).replace('\\','/'):sha(path) for path in source_paths+[args.root/'solver_v1/silicon_thermal_identity.py',
        args.root/'solver_v1/silicon_device_ensemble.py',args.root/'solver_v1/silicon_thermal_research.py']})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,3,figsize=(12,3.4),layout='constrained')
    for chain in range(2):
        label = ['initial 8% start','return 8% start'][chain]
        axes[0].plot(feature_array[chain,:,0],label=label)
        axes[1].plot(feature_array[chain,:,1],label=label)
        axes[2].plot(np.asarray(thermals)[chain,:,2],label=label)
    for axis, title, ylabel in zip(axes,['Total atomistic energy','Projection between initial states','Boundary-safe thermal identity'],
                                   ['eV','Angstrom','Residual (eV)']):
        axis.set(title=title,xlabel='Retained MC index (NOT time)',ylabel=ylabel);axis.grid(alpha=.2)
    axes[2].axhline(0,color='black',lw=.7);axes[0].legend(fontsize=8)
    fig.savefig(args.output/'thermal_chain_diagnostics.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--pilot',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
