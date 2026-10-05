"""Verify exact artifact/source bytes and core numerical scope without new calls."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def main(root):
    result = root/'results/silicon_crack_thermal_validation'
    manifest = json.loads((result/'package_manifest.json').read_text())
    checks = 0
    for relative, digest in manifest['files_sha256'].items():
        path = root/relative
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('hash mismatch: '+relative)
        checks += 1
    pilot = json.loads((result/'nonlinear_pilot/summary.json').read_text())
    analysis = json.loads((result/'analysis/summary.json').read_text())
    geometry = json.loads((result/'geometry/summary.json').read_text())
    if not pilot['complete'] or pilot['new_model_calls'] != 130 or pilot['evaluations_recorded'] != 130:
        raise ValueError('pilot completion/call count changed')
    if any(pilot[name] for name in ('equilibrium_certified','free_energy_estimated','crack_probability_estimated')):
        raise ValueError('unverified physical claim in pilot')
    if analysis['physical_clock'] is not None or analysis['new_calls_in_this_replay'] != 0:
        raise ValueError('clock or replay scope changed')
    for chain in analysis['chains']:
        if chain['evaluation_count'] != 65 or chain['saved_correlated_draws'] != 48:
            raise ValueError('chain record changed')
        if max(chain['maximum_proposal_error'],chain['maximum_position_error_A'],chain['maximum_correction_error']) > 1e-10:
            raise ValueError('unresolved independent replay')
    if geometry['first_crack_label'] is not None or geometry['new_potential_calls'] != 0:
        raise ValueError('unverified crack label or geometry scope changed')
    for relative, digest in geometry['source_sha256'].items():
        if hashlib.sha256((root/relative).read_bytes()).hexdigest() != digest:
            raise ValueError('geometry source changed')
        checks += 1
    with np.load(result/'analysis/chain_diagnostics.npz') as data:
        if data['features'].shape != (2,48,3) or data['thermal_terms'].shape != (2,48,3):
            raise ValueError('retained diagnostic dimensions changed')
        means = data['thermal_terms'].mean(axis=1)
    for chain, mean in zip(analysis['chains'],means):
        expected = [chain['temperature_identity_mean_numerator_eV'],chain['temperature_identity_mean_denominator'],
                    chain['mean_identity_residual_eV']]
        np.testing.assert_allclose(mean,expected,rtol=0,atol=1e-12)
    checks += 9
    print(json.dumps(dict(verified=True,checks=checks,new_potential_calls=0,
                         material_approved=False,physical_clock=None)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('.'))
    main(parser.parse_args().root)
