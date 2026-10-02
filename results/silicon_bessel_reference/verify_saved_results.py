"""Audit saved raw arrays and source hashes without repeating the lattice run."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np


def main(directory):
    directory = Path(directory).resolve()
    root = Path(__file__).resolve().parents[2]
    rows = json.loads((directory / 'states.json').read_text())
    summary = json.loads((directory / 'summary.json').read_text())
    maxima = np.zeros(3)
    fields = ['energy_eV', 'gradient_eV_per_A', 'hessian_eV_per_A2']
    assert len(rows) == 12
    assert {r['cut'] for r in rows} == {'shuffle', 'glide'}
    for row in rows:
        errors = []
        for field in fields:
            a = np.asarray(row['bessel'][field], float)
            b = np.asarray(row['direct'][field], float)
            assert np.all(np.isfinite(a)) and np.all(np.isfinite(b))
            error = float(np.max(abs(a - b)))
            assert error == row['absolute_error'][field]
            errors.append(error)
        h = np.asarray(row['bessel']['hessian_eV_per_A2'])
        assert np.max(abs(h - h.T)) < 1e-12
        maxima = np.maximum(maxima, errors)
    for key, value in zip(['max_energy_error_eV', 'max_gradient_error_eV_per_A',
                           'max_hessian_error_eV_per_A2'], maxima):
        assert summary[key] == value
    refinement = json.loads((directory / 'refinement.json').read_text())
    assert len(refinement) == 7
    for row in refinement:
        for field, error_key in zip(fields, ['energy_error_eV',
                    'gradient_error_eV_per_A', 'hessian_error_eV_per_A2']):
            a, b = [np.asarray(row[method][field]) for method in ['bessel', 'direct']]
            assert float(np.max(abs(a - b))) == row[error_key]
    angular = json.loads((directory / 'angular_cross_check.json').read_text())
    assert abs(angular['combined_moments_energy_eV'] -
               angular['explicit_three_body_energy_eV']) < 3e-12
    assert abs(angular['missing_cross_term_eV']) > 1e-4
    manifest = json.loads((directory / 'source_manifest.json').read_text())
    for item in manifest:
        assert hashlib.sha256((root / item['path']).read_bytes()).hexdigest() == item['sha256']
    assert summary['parameter_sha256'] == manifest[-1]['sha256']
    assert not any(summary[key] for key in ['new_DFT', 'new_MD', 'material_fit',
                     'first_initiation_validated', 'probability_evaluations',
                     'physical_clock_validated', 'refinement_is_certified_tail_bound'])
    result = dict(saved_array_reaggregation='passed', source_hashes='passed',
                  interface_states=len(rows), refinement_controls=len(refinement),
                  max_errors=maxima.tolist(), source_files=len(manifest),
                  new_lattice_evaluations=0, new_MD=0, new_DFT=0)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent)
