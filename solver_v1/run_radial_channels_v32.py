"""Minimal extra per-site radial invariants: necessary compatibility audit only.

The pair and existing ten columns stay intact. Each extra column is extracted
from the verified infinite-sum jet, NOT an independently fitted force curve.
Excluded q points never enter selection. No candidate is production-approved.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import null_space

from .coordination_screening import CoordinationScreenedBulk
from .core_interface_compatibility import minimax_compatibility
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT
from .run_finite_q_compatibility_v31 import points
from .run_source_core_reference import load_source_material
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_low_frequency_forcing_v29 import sha
from .vector_material_calibration import LENGTH_M, MaterialObservation


def channel_shape(shape, kind, decay):
    """One invariant only; other altered columns must never be included."""
    result = np.array(shape, dtype=float, copy=True)
    if result.shape != (6,) or np.any(~np.isfinite(result)) or not np.isfinite(decay) or decay <= 0:
        raise ValueError('six shape parameters and positive finite decay required')
    if kind == 'screened_rank1':
        result[4], column = decay, 6
    elif kind == 'rank3':
        result[1], column = decay, 5
    elif kind == 'quadratic_rank2':
        result[2], result[3], column = decay, 0., 7
    elif kind == 'density_linear':
        result[0], column = decay, 3
    elif kind == 'density_quadratic':
        result[0], column = decay, 4
    else:
        raise ValueError('unknown analytic invariant')
    return tuple(result), column


def conditional_svd(matrix, obs, rows):
    scales = np.array([o.scale for o in obs])
    exact = [i for i, o in enumerate(obs) if o.role == 'exact']
    norm = np.linalg.norm(matrix[exact+rows]/scales[exact+rows, None], axis=0)
    D = np.divide(1., norm, out=np.ones_like(norm), where=norm > 0)
    N = null_space(matrix[exact]*D/scales[exact, None])
    J = matrix[rows]*D/scales[rows, None] @ N
    return dict(exact_rank=len(D)-N.shape[1], null_dimension=N.shape[1],
                singular_values=np.linalg.svd(J, compute_uv=False),
                conditional_rank=int(np.linalg.matrix_rank(J)))


def run(study, out, *, powers_only=False, density_only=False, cross_only=False, vector_mix=False, vector_density_shape=False, spatial_gradient=False, rank_four=False, rank_four_decay=None):
    if rank_four_decay is not None and (not rank_four or not np.isfinite(rank_four_decay) or rank_four_decay<=0):
        raise ValueError('independent rank-four decay requires its explicit study and positive finite value')
    if sum((powers_only, density_only, cross_only, vector_mix, vector_density_shape, spatial_gradient, rank_four)) > 1:
        raise ValueError('alternative studies must be separate')
    started = time.perf_counter()
    study, out = Path(study), Path(out)
    if out.exists():
        raise FileExistsError('fresh output required; preserve previous studies')
    parent = json.loads((study/'summary.json').read_bytes())
    if not parent['completed']:
        raise ValueError('completed parent required')
    shape = tuple(parent['best']['shape'])
    choices = [] if powers_only or cross_only or vector_mix or vector_density_shape or spatial_gradient or rank_four else [(kind, decay) for kind in ('screened_rank1', 'rank3', 'quadratic_rank2')
                                    for decay in (3., 5., 8., 11.)]
    if density_only:
        choices = [(kind, decay) for kind in ('density_linear', 'density_quadratic')
                   for decay in (3., 5., 8., 11.)]
    out.mkdir(parents=True)
    save_json(out/'definition.json', dict(parent_sha256=sha(study/'summary.json'),
        shape=shape, channels=choices, extra_amplitude_nonnegative=not density_only,
        alternative_positive_exponents=[0., 1., 1.5, 2., 3., 4., 6., 8.] if powers_only else [],
        positive_square_cross_shapes=[-100., -10., -1., -.1, .1, 1., 10., 100.] if cross_only else [],
        coherent_vector_shapes=[(k, t) for k in (3., 5.) for t in (-3., -1., -.3, 1.)] if vector_mix else [],
        positive_density_shape_weights=[-.5, .5, 1., 2., 4., 8.] if vector_density_shape else [],
        signed_extra_channels=['density_linear'],
        density_linear_gauge='per-site linear density is equivalent to an exponential pair redistribution; not independently identifiable EAM physics',
        selection='same 115 inspected interface curvature rows plus four fit-q matrices',
        excluded_q_used_for_selection=False, production_changed=False,
        fresh_validation_states=[(.856,.213,.057),(.962,.348,.173),(1.183,.274,.143)] if rank_four else
                                [(.827, .071, -.033), (.879, .267, .113), (.941, .361, .207), (1.127, .419, .139)],
        fresh_states_used_for_selection=False,
        spatial_gradient_energy='D |grad x_i|^2/x_i, D>=0, SAME scalar density' if spatial_gradient else None,
        rank_four_energy='D4 sum_i ||STF(sum_R R^4 f(R))||^2; NONZERO bulk background' if rank_four else None,
        rank_four_decay=(shape[0] if rank_four_decay is None else rank_four_decay) if rank_four else None,
        interpretation='single-column trials first; dictionary is capacity diagnostic only'))
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    obs = list(impose_tangents(problem.observations, problem.source_tangents))
    base = problem.matrix(shape)
    source, _, binding = load_source_material()
    operator = SourceEAMBlochHessian(source.source)
    components = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    weights = np.array([1., 1., 1., np.sqrt(2), np.sqrt(2), np.sqrt(2)])

    def bulk_rows(sh):
        bulk = CoordinationScreenedBulk(sh, radius=16., law='power')
        result = []
        for _, q in points():
            columns, _ = bulk.evaluate(q)
            result.append(weights[:, None]*np.array([columns[:, i, j] for i, j in components]))
        return np.vstack(result)

    matrix = np.vstack([base, bulk_rows(shape)])
    for k, (role, q) in enumerate(points()):
        H = operator.evaluate(np.array(q)*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
        for w, (i, j) in zip(weights, components):
            obs.append(MaterialObservation(f'v32_q{k}_H{i}{j}', float(w*H[i, j]),
                float(.05*np.linalg.norm(H)), 'eV/L0^2', 'heldout' if role == 'excluded' else role))
    rows = [i for i, o in enumerate(obs[:len(base)]) if o.role != 'exact' and o.units == 'eV/L0^2']
    rows += [i for i in range(len(base), len(obs)) if obs[i].role == 'fit']
    reports, summaries, extras = [], [], []

    def solve(label, design, extra_signs=None):
        fit = minimax_compatibility(design, obs, rows,
            nonnegative=NONNEGATIVE+(tuple(range(10, design.shape[1])) if extra_signs is None else tuple(extra_signs)))
        if not fit['completed']:
            raise ArithmeticError('compatibility LP failed')
        excluded = []
        for k, (role, _) in enumerate(points()):
            if role == 'excluded':
                ix = np.arange(len(base)+6*k, len(base)+6*(k+1))
                target = np.array([obs[i].target for i in ix])
                excluded.append(float(np.linalg.norm(fit['predictions'][ix]-target)/np.linalg.norm(target)))
        report = dict(label=label, **fit, identifiability=conditional_svd(design, obs, rows),
                      excluded_q_relative_errors=excluded)
        reports.append(report)
        summaries.append(dict(label=label, eta=fit['minimax_normalized_error'],
            max_excluded_q_error=max(excluded), positive_LJ=fit['strictly_positive_LJ'],
            exact_residual=fit['equality_residual'], duality_gap=fit['duality_gap'],
            conditional_rank=report['identifiability']['conditional_rank'],
            coefficient_count=design.shape[1], material_accepted=False))
        save_json(out/'profiles.json', reports)
        print(label, summaries[-1], flush=True)
        return fit

    baseline = solve('baseline', matrix)
    if abs(baseline['minimax_normalized_error']-parent['best']['value']) > 2e-5:
        raise ArithmeticError('baseline replay mismatch')
    all_signs = []
    for kind, decay in choices:
        sh, column = channel_shape(shape, kind, decay)
        extra = np.r_[problem.matrix(sh)[:, column], bulk_rows(sh)[:, column]]
        extras.append(extra)
        positive = kind != 'density_linear'
        if positive:
            all_signs.append(9+len(extras))
        solve(f'{kind}_k{decay:g}', np.column_stack([matrix, extra]), [10] if positive else [])
    if extras:
        solve('dictionary_capacity_NOT_model_selection', np.column_stack([matrix, *extras]), all_signs)
    if powers_only:
        from .coordination_power_research import positive_power_site_jet
        site_model = problem.cache.site_model(shape[0], shape[4])
        for exponent in (0., 1., 1.5, 2., 3., 4., 6., 8.):
            design = matrix.copy()
            for i, o in enumerate(problem.raw_observations):
                if o.bulk_index is None:
                    density, moments, _ = site_model.site_inputs(o.state)
                    design[i, 6] = 2*np.asarray(o.jet_weights)@positive_power_site_jet(density, moments, exponent)
            # At pristine cubic bulk x=1,Q1=0: this column is EXACTLY p-independent.
            solve(f'alternative_power_{exponent:g}', design)
    if cross_only:
        from .density_angular_cross_material import CrossDensityAngularCache
        raw = CrossDensityAngularCache(problem.raw_observations).cross(shape[0], shape[1])
        cross = np.r_[raw, np.zeros(len(obs)-len(raw))]
        # Lambda * [(x-1)+t I3]^2 per atom, Lambda>=0. Its C/K/cross
        # 2x2 coefficient block is PSD by construction, not an unchecked
        # signed Taylor term. At fixed t only ONE added amplitude is fitted.
        squares = []
        for t in (-100., -10., -1., -.1, .1, 1., 10., 100.):
            column = matrix[:, 4]+2*t*cross+t*t*matrix[:, 9]
            squares.append(column)
            solve(f'positive_cross_square_t{t:g}', np.column_stack([matrix, column]))
        solve('cross_square_cone_capacity_NOT_model_selection', np.column_stack([matrix, *squares]))
    if vector_mix:
        from .mixed_vector_channel import mix_site_moments, mixed_vector_bulk_column
        from .coordination_screening import screened_site_jet
        first_site = problem.cache.site_model(shape[0], shape[4])
        first_bulk = CoordinationScreenedBulk(shape, radius=16., law='power')
        mixed_columns = []
        tails = []
        for k in (3., 5.):
            sh, _ = channel_shape(shape, 'screened_rank1', k)
            second_site = problem.cache.site_model(shape[0], k)
            second_bulk = CoordinationScreenedBulk(sh, radius=16., law='power')
            for t in (-3., -1., -.3, 1.):
                column = np.zeros(len(obs))
                for i, o in enumerate(problem.raw_observations):
                    if o.bulk_index is None:
                        density, a, _ = first_site.site_inputs(o.state)
                        _, b, _ = second_site.site_inputs(o.state)
                        jet = 2*screened_site_jet(density, mix_site_moments(a, b, t), shape[5], law='power')
                        column[i] = np.asarray(o.jet_weights)@jet
                for j, (_, q) in enumerate(points()):
                    H, tail = mixed_vector_bulk_column(first_bulk, second_bulk, q, t)
                    column[len(base)+6*j:len(base)+6*(j+1)] = [w*H[i, h] for w, (i, h) in zip(weights, components)]
                    tails.append(dict(decay=k, weight=t, point_index=j, unit_coefficient_tail=tail))
                mixed_columns.append(column)
                solve(f'coherent_vector_k{k:g}_t{t:g}', np.column_stack([matrix, column]))
        solve('coherent_vector_capacity_NOT_model_selection', np.column_stack([matrix, *mixed_columns]))
        write_csv(out/'mixed_vector_tails.csv', tails)
    if vector_density_shape:
        from .coordination_power_research import positive_power_site_jet
        site_model = problem.cache.site_model(shape[0], shape[4])
        powers = np.zeros((len(obs), 3))
        powers[len(base):] = matrix[len(base):, 6, None]
        for i, o in enumerate(problem.raw_observations):
            if o.bulk_index is None:
                density, moments, _ = site_model.site_inputs(o.state)
                powers[i] = [2*np.asarray(o.jet_weights)@positive_power_site_jet(density, moments, p)
                             for p in (0., 1., 2.)]
        shape_columns = []
        for t in (-.5, .5, 1., 2., 4., 8.):
            # lambda*[1+t(x-1)]²*|Q1|² per site, lambda>=0.
            # Polynomial in the full local density, never in applied stress.
            column = powers@np.array([(1-t)**2, 2*t*(1-t), t*t])
            shape_columns.append(column)
            solve(f'positive_density_vector_shape_t{t:g}', np.column_stack([matrix, column]))
        solve('density_vector_cone_capacity_NOT_model_selection', np.column_stack([matrix, *shape_columns]))
    if spatial_gradient:
        from .spatial_density_gradient import SpatialDensityGradient, gradient_bulk_column
        from .coordination_screening import screened_site_jet
        environment = SpatialDensityGradient(shape[0])
        site = problem.cache.site_model(shape[0], shape[4])
        bulk = CoordinationScreenedBulk(shape, radius=16., law='power')
        column = np.zeros(len(obs)); tails = []
        for i, o in enumerate(problem.raw_observations):
            if o.bulk_index is None:
                density, _, _ = site.site_inputs(o.state)
                gradients, diagnostics = environment.site_jets(o.state)
                jet = 2*screened_site_jet(density, gradients, -1., law='power')
                column[i] = np.asarray(o.jet_weights)@jet
                tails.append(dict(state=o.state, **diagnostics))
        for j, (_, q) in enumerate(points()):
            H, tail = gradient_bulk_column(bulk, q)
            column[len(base)+6*j:len(base)+6*(j+1)] = [w*H[i, h] for w, (i, h) in zip(weights, components)]
            tails.append(dict(q=q, unit_coefficient_bloch_tail=tail))
        solve('true_spatial_density_gradient_square', np.column_stack([matrix, column]))
        save_json(out/'gradient_diagnostics.json', tails)
    if rank_four:
        from .rank_four_environment import RankFourEnvironment, rank_four_bulk_column
        from .yield_elastic_metric import cubic_to_mode_matrix
        environment = RankFourEnvironment(shape[0] if rank_four_decay is None else rank_four_decay)
        bulk = CoordinationScreenedBulk(shape, radius=16., law='power')
        column = np.zeros(len(obs));column[:5] = environment.bulk_column()
        column[2:5] = np.linalg.solve(cubic_to_mode_matrix(),column[2:5])
        tails=[]
        for i,o in enumerate(problem.raw_observations):
            if o.bulk_index is None:
                jet,diagnostics=environment.interface_jet(o.state)
                column[i]=np.asarray(o.jet_weights)@jet
                tails.append(dict(state=o.state,**diagnostics))
        for j,(_,q) in enumerate(points()):
            H,tail=rank_four_bulk_column(environment,bulk,q)
            column[len(base)+6*j:len(base)+6*(j+1)]=[w*H[i,h] for w,(i,h) in zip(weights,components)]
            tails.append(dict(q=q,unit_coefficient_bloch_tail=tail))
        design=np.column_stack([matrix,column])
        solve('rank_four_nonnegative',design)
        solve('rank_four_signed_DIAGNOSTIC',design,[])
        save_json(out/'rank_four_diagnostics.json',tails)
    save_json(out/'observations.json', [asdict(o) for o in obs])
    write_csv(out/'comparison.csv', summaries)
    save_json(out/'summary.json', dict(completed=True, source=binding,
        elapsed_seconds=time.perf_counter()-started, profiles=summaries,
        material_accepted=False, physical_time_calibrated=False, actual_yield_validated=False,
        production_changed=False, interpretation='necessary fixed-shape compatibility only; no full spectral or fresh validation approval'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--powers-only', action='store_true')
    parser.add_argument('--density-only', action='store_true')
    parser.add_argument('--cross-only', action='store_true')
    parser.add_argument('--vector-mix', action='store_true')
    parser.add_argument('--vector-density-shape', action='store_true')
    parser.add_argument('--spatial-gradient', action='store_true')
    parser.add_argument('--rank-four', action='store_true')
    parser.add_argument('--rank-four-decay', type=float)
    args = parser.parse_args()
    if sum((args.powers_only, args.density_only, args.cross_only, args.vector_mix, args.vector_density_shape, args.spatial_gradient, args.rank_four)) > 1:
        parser.error('alternative studies must be separate')
    run(args.study, args.out, powers_only=args.powers_only, density_only=args.density_only,
        cross_only=args.cross_only, vector_mix=args.vector_mix, vector_density_shape=args.vector_density_shape,
        spatial_gradient=args.spatial_gradient,rank_four=args.rank_four,rank_four_decay=args.rank_four_decay)
