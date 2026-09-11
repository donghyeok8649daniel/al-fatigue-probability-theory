"""One existing rational-shape ablation after a quantified jet incompatibility.

Fix all six current shapes. Compare beta=alpha3*I_reference, where I_reference
is the actual maximum site Q3:Q3 at the pre-existing source saddle state.
This is a dimensionless moment normalization, not an area/length or yield fit.
No source-core or production implementation of this ablation is presumed.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .core_interface_compatibility import minimax_compatibility
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .quartic_jet_ablation import RationalQuarticJetColumn
from .run_current_material_core import ROOT, load_current_material
from .tail_constrained_material import spectral_profile
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--model', type=Path, default=ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh explicit quartic-ablation output required')
    began = time.perf_counter()
    model, _, binding = load_current_material(args.model)
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    obs = impose_tangents(problem.observations, problem.source_tangents)
    cache = RationalQuarticJetColumn(model, obs)
    saddle = next(o.state for o in obs if o.name == 'saddle_Hxx')
    reference = cache.reference_invariant(saddle)
    if not reference > 0:
        raise ValueError('positive resolved saddle invariant required for numerical normalization')
    betas = (0., .1, 1., 10., 100.)
    save_json(args.out/'definition.json', dict(candidate_sha256=binding['parameter_sha256'],
        source_sha256=problem.source.reference.sha256, frozen_shape=model.screened_shape,
        beta_values=betas, reference_saddle=saddle, reference_maximum_site_Q3_norm_squared=reference,
        energy_law='K3 * I3^2 / (1 + alpha3*I3)', alpha3='beta/reference_invariant',
        explicit_new_combination_of_previously_separate_ablation=True,
        previous_v14_ablation_was_not_adopted=True, shapes_refitted=False,
        core_force_fit=False, new_production_energy=False, physical_yield_used=False,
        core_bridge_supported=False, material_accepted=False))
    M = problem.matrix(tuple(model.screened_shape))
    zero = cache.column(0.)
    replay = float(np.max(abs(zero-M[:, 9])))
    if replay > 2e-10:
        raise ArithmeticError(f'zero-cap column failed exact old-family replay: {replay}')
    ops, tails = problem.operators(tuple(model.screened_shape))
    if np.max(abs(ops[:, 9])) > 1e-14 or np.max(abs(tails[:, 9])) > 1e-14:
        raise ArithmeticError('centrosymmetric-bulk quartic harmonic column is not zero')
    fit_rows = [i for i, o in enumerate(obs) if o.role == 'fit']
    jets = [i for i, o in enumerate(obs) if o.role != 'exact' and o.units == 'eV/L0^2']
    summaries, profiles, residuals = [], [], []
    for beta in betas:
        tick = time.perf_counter()
        alpha = beta/reference
        matrix = M.copy()
        matrix[:, 9] = cache.column(alpha)
        lp = minimax_compatibility(matrix, obs, jets, nonnegative=NONNEGATIVE)
        fit = spectral_profile(matrix, obs, ops, tails, nonnegative=NONNEGATIVE)
        # Sensitivity of the new column against D3 at this fixed shape.
        columns = matrix[fit_rows][:, [5, 9]]/np.array([obs[i].scale for i in fit_rows])[:, None]
        columns /= np.linalg.norm(columns, axis=0)
        singular = np.linalg.svd(columns, compute_uv=False)
        row = dict(beta=beta, alpha3=alpha, static_fit_loss=fit['squared_loss'],
                   inspected_curvature_lower_bound=lp.get('minimax_normalized_error'),
                   positive_LJ=fit['strictly_positive_LJ'], exact_residual=fit['exact_residual'],
                   kkt=fit['kkt_residual'], bulk_margin=fit['minimum_robust_margin'],
                   D3_K3_column_correlation=float(columns[:, 0]@columns[:, 1]),
                   D3_K3_singular_ratio=float(singular[-1]/singular[0]),
                   elapsed_seconds=time.perf_counter()-tick, material_accepted=False)
        for o, pred, ri in zip(obs, fit['predictions'], fit['residuals']):
            residuals.append(dict(beta=beta, observable=o.name, previous_role=o.role,
                target=o.target, scale=o.scale, prediction=pred, units=o.units,
                normalized_error=ri, used_for_LS=o.role == 'fit',
                used_for_diagnostic_LP=o.role != 'exact' and o.units == 'eV/L0^2',
                blind_validation=False))
        summaries.append(row)
        profiles.append(dict(beta=beta, alpha3=alpha, fixed_shape=model.screened_shape,
                             minimax=lp, coefficient_fit=fit, material_accepted=False))
        save_json(args.out/'checkpoint.json', dict(completed=False, profiles=profiles))
        print(f'quartic beta={beta:g}: static={row["static_fit_loss"]:.7g}, '
              f'curvature lower={row["inspected_curvature_lower_bound"]:.7g}', flush=True)
    write_csv(args.out/'profile_summary.csv', summaries)
    write_csv(args.out/'observable_residuals.csv', residuals)
    save_json(args.out/'completion.json', dict(completed=True, profiles=profiles,
        column_zero_replay_error=replay, actual_profiles=len(profiles),
        elapsed_seconds=time.perf_counter()-began, material_accepted=False,
        core_implementation_validated=False, whole_family_failure_proved=False,
        physical_Hz=False, production_changed=False,
        file_hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__ == '__main__':
    main()
