"""Compare separately replayed force-control states at two force tolerances.

This is a numerical convergence check in one small bare specimen, not material
validation. Source hashes and the reference geometry are checked before pairing.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    paths = [args.baseline, args.refined]
    replays = [args.baseline_replay, args.refined_replay]
    protocols = [json.loads((p/'protocol.json').read_text()) for p in paths]
    audits = [json.loads((p/'summary.json').read_text()) for p in replays]
    for key in ('model_sha256', 'geometry_sha256', 'area_A2',
                'bulk_reference_gauge_A', 'nominal_stresses_GPa', 'energy'):
        if protocols[0][key] != protocols[1][key]:
            raise ValueError(f'comparison protocols differ: {key}')
    if not protocols[1]['fmax_eV_A'] < protocols[0]['fmax_eV_A']:
        raise ValueError('refined run must have a smaller force tolerance')
    for folder, audit in zip(paths, audits):
        for source in audit['source_raw_files']:
            if sha(folder/source['file']) != source['sha256']:
                raise ValueError('source raw changed since independent replay')
    a, b = audits
    if a['states_verified'] != b['states_verified']:
        raise ValueError('different state counts')
    rows = []
    for old, new in zip(a['records'], b['records']):
        if old['target_nominal_GPa'] != new['target_nominal_GPa']:
            raise ValueError('different loads')
        if not old['converged'] or not new['converged']:
            raise ValueError('both compared states must be converged')
        index = old['state']
        with np.load(paths[0]/f'state_{index:03d}/raw.npz') as data:
            r0, g0 = data['positions'].copy(), data['gradient'].copy()
        with np.load(paths[1]/f'state_{index:03d}/raw.npz') as data:
            r1, g1 = data['positions'].copy(), data['gradient'].copy()
        difference = np.linalg.norm(r1-r0, axis=1)
        rows.append(dict(target_MPa=1000*old['target_nominal_GPa'],
            old_fmax_eV_A=old['free_force_max_eV_A'],
            new_fmax_eV_A=new['free_force_max_eV_A'],
            old_reduced_gradient_max_eV_A=float(abs(g0).max()),
            new_reduced_gradient_max_eV_A=float(abs(g1).max()),
            old_axial_error_MPa=old['axial_target_residual_MPa'],
            new_axial_error_MPa=new['axial_target_residual_MPa'],
            old_support_residual_MPa=old['support_resultant_equivalent_MPa'],
            new_support_residual_MPa=new['support_resultant_equivalent_MPa'],
            old_support_torque_eV=old['reaction_torque_eV'],
            new_support_torque_eV=new['reaction_torque_eV'],
            extension_change_A=new['extension_A']-old['extension_A'],
            enthalpy_change_eV=new['enthalpy_eV']-old['enthalpy_eV'],
            atom_displacement_RMS_A=float(np.sqrt(np.mean(difference**2))),
            atom_displacement_max_A=float(difference.max())))
    result = dict(records=rows, old_force_tolerance_eV_A=protocols[0]['fmax_eV_A'],
        new_force_tolerance_eV_A=protocols[1]['fmax_eV_A'],
        old_secant_GPa=a['comparison']['secant_stress_strain_ratio_GPa'],
        new_secant_GPa=b['comparison']['secant_stress_strain_ratio_GPa'],
        scope='numerical force-tolerance control of the same finite bare prism; not a bulk modulus or initiation validation',
        exact_support_equilibrium_certified=False,
        old_curvature_not_transferred_to_refined_states=True,
        new_potential_calls=0, new_DFT=0, new_MD=0, initiation_probability=None,
        inputs=[dict(role=label, file=p.name, sha256=sha(p)) for label, p in
            [('baseline_protocol', paths[0]/'protocol.json'),
             ('refined_protocol', paths[1]/'protocol.json'),
             ('baseline_replay', replays[0]/'summary.json'),
             ('refined_replay', replays[1]/'summary.json')]])
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with (args.output/'comparison.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('baseline', 'refined', 'baseline-replay', 'refined-replay', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())
