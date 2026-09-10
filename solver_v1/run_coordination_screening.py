"""Actual deterministic fixed-range calibration of the v19 screening hypothesis.

Fits the SAME v17 objective with/without inherited pair CONTROL. Preserves old
results. New off-grid observations are declared before optimization, excluded
from selection, and saved separately. No production parameter is changed.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar

from .coordination_screening import CoordinationScreenedCache, screening_factor
from .rank_one_range_material import RankOneRangeBulk
from .run_vector_material_calibration import source_and_targets
from .interface_normal_development_targets import normal_development_observations
from .material_calibration_controls import append_fixed_pair
from .tail_constrained_material import spectral_profile
from .vector_material_calibration import IDEAL_H, UNITS, MaterialObservation
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .profiled_identifiability import equality_tangent_sensitivity


def declared_new_observations(source, *, followup=False):
    # Deliberately different from v17/v16 scalar locations. Do not choose these
    # after seeing predictions or fit to these states in the current experiment.
    states = [(r*IDEAL_H, 0., 0.) for r in (.992, 1.028, 1.145, 1.325, 1.705, 2.175)]
    states += [(1.058*IDEAL_H, .173, -.043), (1.19*IDEAL_H, .367, .081),
               (IDEAL_H, .317, 0.), (1.015*IDEAL_H, .413, .219)]
    if followup:
        states = [(r*IDEAL_H, 0., 0.) for r in (.982, 1.052, 1.168, 1.385, 1.765, 2.285)]
        states += [(1.042*IDEAL_H, .193, -.063), (1.175*IDEAL_H, .347, .061),
                   (IDEAL_H, .337, 0.), (1.022*IDEAL_H, .393, .229)]
    rows = []
    for number, q in enumerate(states):
        result = source.evaluate(q)
        jet = np.r_[result.energy, result.gradient, result.hessian[0], result.hessian[1, 1:], result.hessian[2, 2]]
        for index, field, unit in ((0, 'energy', 'eV/cell'), (1, 'force', 'eV/L0'),
                                  (4, 'Haa', 'eV/L0^2'), (7, 'Hxx', 'eV/L0^2')):
            value = float(jet[index])
            scale = (float(UNITS.traction_mpa_to_force(250.)) if index == 1 else
                     max(.1*abs(value), .01/float(UNITS.energy_to_surface(1.))) if index == 0 else
                     max(.1*abs(value), .1))
            prefix = 'v19_power_new' if followup else 'v19_new'
            rows.append(MaterialObservation(f'{prefix}_{number}_{field}', value, scale, unit,
                                           'heldout', tuple(q), tuple(np.eye(10)[index])))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--maxiter', type=int, default=40)
    parser.add_argument('--law', choices=('rational', 'power'), default='rational')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('use a fresh result directory, never overwrite old fits')
    began = time.perf_counter()
    parent = json.loads((args.parent/'calibration.json').read_bytes())
    definition = json.loads((args.parent/'definition.json').read_bytes())
    if not parent['completed'] or not definition.get('rank_one_range_extension'):
        raise ValueError('completed v17 independent-range parent required')
    shape = np.array(parent['best']['decays'])
    source, previous, states = source_and_targets()
    if definition['source_sha256'] != source.reference.sha256:
        raise ValueError('source binding mismatch')
    observations, provenance = normal_development_observations(source, previous, states)
    new = declared_new_observations(source)
    new_prefix = 'v19_new_'
    if args.law == 'power':
        # The rational probe's states have already been inspected. Keep them
        # excluded but classify them as retrospective, and declare new ones.
        new += declared_new_observations(source, followup=True)
        new_prefix = 'v19_power_new_'
    lower, upper = (-1., 1.) if args.law == 'power' else (0., 1.)
    grid = [-1., -.5, 0., .5, 1.] if args.law == 'power' else [0., .25, .5, .75, .9, .99, 1.]
    existing = {(o.state, o.jet_weights) for o in observations if o.state is not None}
    if any((o.state, o.jet_weights) in existing for o in new):
        raise ValueError('new validation overlaps old observations')
    observations += new
    save_json(args.out/'definition.json', dict(
        family='coordination_screened_rank1_v19', source_sha256=source.reference.sha256,
        parent_calibration_sha256=hashlib.sha256((args.parent/'calibration.json').read_bytes()).hexdigest(),
        inherited_shape=shape, fixed_shape_control=True, actual_radial_optimization=False,
        screening_law=args.law, screening_definition='x_site**p' if args.law == 'power' else '1/(1-z+z*x_site)',
        screening_bounds=[lower, upper], screening_grid=grid, optimizer_maxiter=args.maxiter,
        observations=[asdict(o) for o in observations],
        retrospective_validation='all previously inspected excluded states', new_validation_prefix=new_prefix,
        fixed_pair_control=definition.get('fixed_pair_control'),
        source_temperature_K=0., length_scale_m=UNITS.length_scale_m,
        energy_scale_J=1.602176634e-19,
        stability_stretches=definition['stability_stretches'],
        wavepoints_cubic=definition['wavepoints_cubic'], radius_over_L0=definition['radius_over_L0'],
        physical_kinetics_calibrated=False, production_changed=False))
    cache = CoordinationScreenedCache(observations, law=args.law)
    inherited_spectra = []
    for stretch in definition['stability_stretches']:
        basis = RankOneRangeBulk(shape, radius=definition['radius_over_L0'], stretch=stretch)
        matrices, tails = map(np.asarray, zip(*(basis.evaluate(q) for q in definition['wavepoints_cubic'])))
        inherited_spectra.append((basis.x, basis.density_series_tail, matrices, tails))
    print('source/new-validation defined; inherited finite-q bases ready', flush=True)
    results = []
    for mode in ('fixed_pair_control', 'free_pair'):
        profiles = []

        def problem(z):
            raw = cache.matrix(np.r_[shape, z])
            matrix, obs = cubic_metric_problem(raw[:, :8], observations)
            extra = raw[:, 8:].copy()
            extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
            matrix = np.column_stack([matrix, extra])
            obs = [replace(o, role='exact') if i < 5 else o for i, o in enumerate(obs)]
            if mode == 'fixed_pair_control':
                matrix, obs = append_fixed_pair(matrix, obs, definition['fixed_pair_control'])
            return matrix, obs

        def objective(z):
            M, obs = problem(z); operators, errors = [], []
            for x, dx, oldH, oldtail in inherited_spectra:
                H, tail = oldH.copy(), oldtail.copy()
                factor = float(screening_factor(x, z, law=args.law)[0])
                dg = abs(float(screening_factor(x-dx, z, law=args.law)[1]))
                tail[:, 6] = factor*tail[:, 6]+dg*dx*(np.linalg.norm(H[:, 6], ord=2, axis=(1, 2))+tail[:, 6])
                H[:, 6] *= factor
                operators.extend(H); errors.extend(tail)
            fit = spectral_profile(M, obs, operators, errors, nonnegative=(0, 1, 2, 4, 5, 8, 9))
            row = dict(screening=float(z), **fit)
            profiles.append(row)
            save_json(args.out/mode/'checkpoint.json', dict(completed=False, profiles=profiles))
            print(mode, len(profiles), 'z', z, 'loss', fit['squared_loss'], flush=True)
            return fit['squared_loss']

        for z in grid:
            objective(z)
        optimization = minimize_scalar(objective, bounds=(lower, upper), method='bounded',
            options=dict(maxiter=args.maxiter, xatol=1e-6))
        best = min(profiles, key=lambda p: p['squared_loss'])
        z, c = best['screening'], np.asarray(best['coefficients'])
        M, obs = problem(z)
        target, scales = np.array([o.target for o in obs]), np.array([o.scale for o in obs])
        predictions = M@c; norm = (predictions-target)/scales
        rows = [dict(observable=o.name, role=o.role,
            validation_status='new_predeclared' if o.name.startswith(new_prefix) else
                              'retrospective' if o.role == 'heldout' else 'development_or_control',
            target=o.target, prediction=float(p), units=o.units, scale=o.scale,
            normalized_residual=float(r)) for o, p, r in zip(obs, predictions, norm)]
        write_csv(args.out/mode/'residuals.csv', rows)
        # Screening tangent is taken at fixed inherited ranges. Boundary endpoints
        # use second-order one-sided differences; no fictitious negative screening.
        sensitivities = []
        for step in (2e-4, 1e-4):
            if z < lower+step:
                derivative = (-3*M+4*problem(z+step)[0]-problem(z+2*step)[0])/(2*step)
            elif z > upper-step:
                derivative = (3*M-4*problem(z-step)[0]+problem(z-2*step)[0])/(2*step)
            else:
                derivative = (problem(z+step)[0]-problem(z-step)[0])/(2*step)
            sensitivities.append(equality_tangent_sensitivity(M, obs, c, derivative[None],
                exact_rows=[i for i, o in enumerate(obs) if o.role == 'exact']))
        identity = sensitivities[-1]
        identity.update(fixed_inherited_ranges=True, screening_coordinate='linear shape', screening_law=args.law,
            jacobian_step_change=float(np.max(abs(identity['jacobian']-sensitivities[0]['jacobian']))),
            statistical_confidence_claimed=False)
        save_json(args.out/mode/'identifiability.json', identity)
        save_json(args.out/mode/'calibration.json', dict(completed=True, best=best,
            shape=np.r_[shape, z], screening_law=args.law, profiles=profiles,
            scalar_optimizer=dict(success=bool(optimization.success), message=str(optimization.message),
                                  nfev=optimization.nfev, x=optimization.x, fun=optimization.fun),
            material_accepted=False, actual_radial_optimization=False,
            physical_kinetics_calibrated=False))
        summary = dict(mode=mode, screening=z, squared_loss=best['squared_loss'],
            positive_LJ=best['strictly_positive_LJ'], scalar_optimizer_success=bool(optimization.success),
            profiles=len(profiles), D1=float(c[6]))
        for label, choose in (
            ('retrospective', [i for i, o in enumerate(obs) if o.role == 'heldout' and not o.name.startswith(new_prefix)]),
            ('new', [i for i, o in enumerate(obs) if o.name.startswith(new_prefix)])):
            summary[label+'_rms'] = float(np.sqrt(np.mean(norm[choose]**2)))
        for name in ('perfect_Haa', 'perfect_Hxx'):
            i = next(i for i, o in enumerate(obs) if o.name == name)
            summary[name] = float(predictions[i]); summary[name+'_source'] = float(target[i])
        results.append(summary)
        write_csv(args.out/'summary.csv', results)
    save_json(args.out/'completion.json', dict(completed=True, summary=results,
        entire_family_impossibility_proved=False, material_accepted=False,
        elapsed_seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
