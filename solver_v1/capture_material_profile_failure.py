"""Capture a real profiled-QP boundary case as a compact independent fixture."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from .interface_even_development_targets import even_development_observations
from .isotropic_bulk_validation import IsotropicBulkBasis
from .material_calibration_controls import append_fixed_pair
from .quadrupole_saturation import SaturatedQuadrupoleCache
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from . import tail_constrained_material as qp
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve the original numerical regression fixture')
    raw = (args.directory/'calibration.json').read_bytes()
    data = json.loads(raw); definition = json.loads((args.directory/'definition.json').read_bytes())
    source, obs, states = source_and_targets()
    if not definition.get('quadrupole_saturation_extension'):
        raise ValueError('this reproduction is for the actual v16 saturated-family profile')
    obs, _ = even_development_observations(source, obs, states)
    shape = data['best']['decays']
    raw_matrix = SaturatedQuadrupoleCache(obs).matrix(shape)
    matrix, obs = cubic_metric_problem(raw_matrix[:, :8], obs)
    extra = raw_matrix[:, 8:].copy()
    extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
    matrix = np.column_stack([matrix, extra])
    obs = [replace(o, role='exact') if i < 5 else o for i, o in enumerate(obs)]
    matrix, obs = append_fixed_pair(matrix, obs, definition['fixed_pair_control'])
    values = []
    for stretch in definition['stability_stretches']:
        basis = IsotropicBulkBasis(shape[:3], stretch=stretch, radius=definition['radius_over_L0'])
        values.extend(basis.evaluate(q) for q in definition['wavepoints_cubic'])
    columns, tails = map(np.asarray, zip(*values))
    original = qp.convex_profile
    last = {}

    def capture(M, observations, inequalities=None, *, nonnegative):
        result = original(M, observations, inequalities, nonnegative=nonnegative)
        last.update(matrix=M, target=np.array([o.target for o in observations]),
            scale=np.array([o.scale for o in observations]),
            role=np.array([o.role for o in observations]),
            inequalities=np.asarray(inequalities), nonnegative=np.asarray(nonnegative),
            returned_coefficients=result['coefficients'])
        return result

    qp.convex_profile = capture
    try:
        result = qp.spectral_profile(matrix, obs, columns, tails, nonnegative=(0, 1, 2, 4, 5, 8, 9))
    finally:
        qp.convex_profile = original
    args.out.mkdir(parents=True)
    np.savez_compressed(args.out/'actual_qp.npz', **last)
    save_json(args.out/'reproduction.json', dict(shape=shape,
        original_calibration_sha256=hashlib.sha256(raw).hexdigest(), source_sha256=source.reference.sha256,
        coefficients=result['coefficients'], original_saved_coefficients=data['best']['coefficients'],
        maximum_replay_change=float(np.max(abs(result['coefficients']-data['best']['coefficients']))),
        loss=result['squared_loss'], kkt=result['kkt_residual'],
        physical_nonnegative_min=float(np.min(result['coefficients'][[0, 1, 2, 4, 5, 8, 9]])),
        interpretation='numerical coefficient/QP reproduction, not a material fit or clipped parameter'))
    print('actual numerical fixture captured', result['coefficients'], flush=True)


if __name__ == '__main__':
    main()
