"""Execute one provenance-bound current-material Bessel screw-core study.

Static research only. Existing positive-LJ parameters are never optimized.
Save checkpoints, distinguish iteration/budget failure, preserve all old runs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .coordination_screening import CoordinationScreenedInterface
from .core_continuation import continue_core_field
from .core_stability import lowest_core_mode
from .current_material_rows import CurrentMaterialScrewCore
from .current_core_diagnostics import (phase_winding_cells, adjacent_registry_profile,
    translated_core_seed,reference_configuration_seed)
from .isolated_screw_core import ScrewFarField
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .vector_material_calibration import LENGTH_M


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATERIAL = ROOT/'results/fcc111_active_interface/tangent_calibration_v20/shape_refinement/old_family/calibration.json'
EV_J = 1.602176634e-19


def load_current_material(path=DEFAULT_MATERIAL, *, tolerance=2e-12):
    path = Path(path).resolve()
    raw = path.read_bytes(); data = json.loads(raw)
    if not data['completed'] or not data['best']['strictly_positive_LJ']:
        raise ValueError('completed positive-LJ research candidate required')
    fit = data['best']
    shape = np.asarray(data.get('shape', fit.get('shape')), float)
    model = CoordinationScreenedInterface(shape, fit['coefficients'],
        law=data['screening_law'], tolerance=tolerance)
    elastic = np.asarray(fit['predictions'][2:5], float)
    tensor = rotate_elastic_tensor(cubic_elastic_tensor(*elastic),
        model.geometry.plane_basis_in_stacked_cubic_axes())*1e9*LENGTH_M**3/EV_J
    metadata = dict(parameter_source=path.relative_to(ROOT).as_posix(),
        parameter_sha256=hashlib.sha256(raw).hexdigest(), material_shape=shape,
        coefficients=fit['coefficients'], screening_law=data['screening_law'],
        material_accepted=False, energy_terms_omitted=False, length_scale_m=LENGTH_M,
        energy_unit='eV per infinite straight-row repeat b*L0',
        tensor_unit='eV/L0^3', elastic_tensor_same_candidate=True,
        C11_C12_C44_GPa=elastic, physical_time=False, physical_Hz=False,
        actual_yield_calibrated=False, production_changed=False)
    return model, tensor, metadata


class CalculationBudgetReached(RuntimeError):
    pass


def run(args, *, material_loader=load_current_material, core_builder=CurrentMaterialScrewCore):
    for key in ('radius', 'max_seconds', 'force_tolerance', 'tolerance', 'escape_amplitude'):
        value = getattr(args, key)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'{key} must be finite and positive')
    if args.iterations < 1 or args.ring < 1:
        raise ValueError('positive iteration and ring counts required')
    if not np.isfinite(args.shear_mpa) or not np.isfinite(args.center_y):
        raise ValueError('finite loading and center coordinates required')
    out = Path(args.out)
    if out.exists():
        raise FileExistsError('fresh case directory required; prior data preserved')
    began = time.perf_counter()
    model, tensor, metadata = material_loader(args.material, tolerance=args.tolerance)
    center = (np.sqrt(3)*model.geometry.b/12+args.center_y, model.h/2)
    far = ScrewFarField(tensor, model.geometry.b, center)
    traction = args.shear_mpa*1e6*LENGTH_M**3/EV_J
    core = core_builder(model, far, free_radius=args.radius, ring=args.ring,
        tolerance=args.tolerance, shear_traction=traction, burgers_sign=args.burgers_sign)
    metadata.update(radius_over_L0=args.radius, ring=args.ring, center_over_L0=center,
        b_reduced=model.geometry.b, h_reduced=model.h,
        shear_traction_MPa=args.shear_mpa, burgers_sign=args.burgers_sign,
        force_tolerance_eV_L0=args.force_tolerance, reciprocal_tolerance=args.tolerance,
        free_sites=len(core.free_ids), energy_sites=len(core.energy_ids),
        explicit_environment_rows=len(core.offsets), channels=core.channel_count,
        farfield_log_coefficient_eV_per_repeat=far.log_energy_coefficient,
        boundary='fixed same-material Volterra plus prescribed affine shear',
        boundary_is_not_measured_pinning=True, domain_converged=False,
        max_iterations=args.iterations, wall_budget_seconds=args.max_seconds,
        stability_method=args.stability_method)
    seed = core.initial
    if getattr(args,'resume_checkpoint',None):
        if args.from_case or getattr(args,'configuration_seed',None) or args.escape_sign or args.translate_j:
            raise ValueError('checkpoint recovery is a separate explicitly incomplete-state calculation')
        origin=Path(args.resume_checkpoint).resolve()
        if not (origin/'failure.json').exists():
            raise ValueError('checkpoint recovery requires a stopped/failed run, not a live or certified case')
        prior=json.loads((origin/'metadata.json').read_bytes())
        failure=json.loads((origin/'failure.json').read_bytes())
        if prior['parameter_sha256']!=metadata['parameter_sha256'] or failure.get('completed',True):
            raise ValueError('same bound material and explicitly uncompleted source required')
        if prior['ring']!=args.ring or prior['reciprocal_tolerance']!=args.tolerance:
            raise ValueError('checkpoint recovery cannot silently change the represented operator')
        with (origin/'checkpoint.csv').open(newline='',encoding='utf8') as stream:
            rows=list(csv.DictReader(stream))
        saved={(int(r['j']),int(r['l'])):np.array([float(r[k]) for k in ('ux','uy','uz')]) for r in rows}
        if len(saved)!=len(rows):raise ValueError('duplicate checkpoint rows')
        seed=reference_configuration_seed(core,prior,saved,length_scale_m=LENGTH_M)
        metadata['checkpoint_recovery']=dict(source=origin.relative_to(ROOT).as_posix(),
            checkpoint_sha256=hashlib.sha256((origin/'checkpoint.csv').read_bytes()).hexdigest(),
            source_was_unconverged=True,source_failure=failure,
            independent_final_force_and_Morse_checks_required=True)
    if getattr(args,'configuration_seed',None):
        if args.from_case or args.escape_sign or args.translate_j:
            raise ValueError('reference configuration is a separate multistart, not continuation')
        origin=Path(args.configuration_seed).resolve()
        prior=json.loads((origin/'metadata.json').read_bytes())
        previous=json.loads((origin/'summary.json').read_bytes())
        if not previous['completed'] or not previous['final_stability_probe']['stable_on_tested_fixed_boundary']:
            raise ValueError('reference seed requires a separately verified stable source')
        with (origin/'state.csv').open(newline='',encoding='utf8') as stream:
            rows=list(csv.DictReader(stream))
        saved={(int(r['j']),int(r['l'])):np.array([float(r[k]) for k in ('ux','uy','uz')]) for r in rows}
        if len(saved)!=len(rows):raise ValueError('duplicate logical source rows')
        seed=reference_configuration_seed(core,prior,saved,length_scale_m=LENGTH_M)
        metadata['configuration_seed']=dict(source=origin.relative_to(ROOT).as_posix(),
            source_parameter_sha256=prior['parameter_sha256'],
            source_state_sha256=hashlib.sha256((origin/'state.csv').read_bytes()).hexdigest(),
            coordinates_only=True,source_energy_adopted=False,source_elasticity_adopted=False,
            target_exterior_unchanged=True,final_position_constrained=False)
    if args.from_case:
        origin = Path(args.from_case).resolve()
        prior = json.loads((origin/'metadata.json').read_bytes())
        previous = json.loads((origin/'summary.json').read_bytes())
        if not previous.get('completed'):
            raise ValueError('continuation source must be a completed calculation')
        if prior.get('burgers_sign', 1) != args.burgers_sign:
            raise ValueError('continuation cannot change the imposed Burgers topology')
        with (origin/'state.csv').open(newline='', encoding='utf-8') as stream:
            rows = list(csv.DictReader(stream))
        saved = {(int(r['j']), int(r['l'])):np.array([float(r[k]) for k in ('ux', 'uy', 'uz')]) for r in rows}
        if len(saved) != len(rows):
            raise ValueError('duplicate logical source rows')
        seed, transfer = continue_core_field(core, prior, saved,
            parameter_sha256=metadata['parameter_sha256'], length_scale_m=LENGTH_M,
            shear_mpa=args.shear_mpa, allow_expansion=args.extend_domain,
            require_stable=not args.allow_unstable_source, prior_summary=previous)
        metadata.update(continued_from=origin.relative_to(ROOT).as_posix(), continuation=transfer,
            source_state_sha256=hashlib.sha256((origin/'state.csv').read_bytes()).hexdigest())
    if args.escape_sign:
        if not args.from_case or not args.allow_unstable_source:
            raise ValueError('mode descent requires explicitly selected unstable saved source')
        mode = lowest_core_mode(core, seed, force_tolerance=args.force_tolerance, method=args.stability_method)
        if not mode['minimum_eigenvalue'] < -max(mode['eigenpair_residual'], 1e-9):
            raise ValueError('no resolved negative direction')
        if not all(r['curvature'] < 0 for r in mode['energy_curvature_checks']):
            raise ValueError('energy differences do not verify the negative mode')
        seed += args.escape_sign*args.escape_amplitude*mode['eigenvector']
        metadata['mode_descent'] = dict(sign=args.escape_sign, amplitude=args.escape_amplitude,
            minimum_eigenvalue=mode['minimum_eigenvalue'], energy_checks=mode['energy_curvature_checks'])
    if args.translate_j:
        if not args.from_case or args.allow_unstable_source or args.escape_sign:
            raise ValueError('translation requires a verified stable source and a separate experiment')
        seed = translated_core_seed(core, seed, args.translate_j)
        metadata['core_translation_initial_guess'] = dict(j_steps=args.translate_j,
            glide_y_over_L0=args.translate_j*core.rows.d, imposed_final_translation=False,
            exterior_boundary_translated=False)
    save_json(out/'metadata.json', metadata)
    latest = seed.copy()
    def save_field(path, field, gradient=None):
        rows = []
        for index, (logical, u) in enumerate(zip(core.indices[core.free_ids], field)):
            row = dict(j=int(logical[0]), l=int(logical[1]), ux=u[0], uy=u[1], uz=u[2])
            if gradient is not None:
                row.update(gx=gradient[index, 0], gy=gradient[index, 1], gz=gradient[index, 2])
            rows.append(row)
        write_csv(path, rows)
    def checkpoint(iteration, field):
        nonlocal latest
        latest = field.copy()
        elapsed = time.perf_counter()-began
        if iteration == 1 or iteration % 10 == 0:
            save_field(out/'checkpoint.csv', field)
            save_json(out/'progress.json', dict(completed=False, iteration=iteration, elapsed_seconds=elapsed))
            print(f'{out.name}: iteration {iteration}, elapsed={elapsed:.2f}s', flush=True)
        if elapsed >= args.max_seconds:
            raise CalculationBudgetReached('declared wall budget reached; checkpoint is not a converged result')
    try:
        initial = core.evaluate(seed)
        if args.minimizer == 'stationary_root':
            # Independent full-matrix root for a near-stationary DOMAIN check.
            # This is not energy-descent dynamics or a load-history integrator.
            # A saddle can result; the separate final Morse gate is mandatory.
            from .core_stationary_point import stationary_core
            root = stationary_core(core,seed,max_iterations=args.iterations,
                force_tolerance=args.force_tolerance,
                callback=lambda i,q,row:checkpoint(i,q))
            value = core.evaluate(root['field'])
            result = dict(field=root['field'],energy=root['energy'],
                maximum_free_force=root['maximum_force'],force_converged=root['force_converged'],
                optimizer_success=root['force_converged'],
                optimizer_message='stationary Newton root; final Morse gate separately decides minimum',
                iterations=len(root['history'])-1,evaluations=None,polish_steps=0,
                minimizer='stationary_root',energy_descent_path=False,physical_dynamics=False,
                stationary_history=root['history'],boundary_winding=core.boundary_winding(root['field']),
                minimum_density=value['minimum_density'],minimum_transverse_radius=value['minimum_transverse_radius'])
        else:
            result = core.relax(seed, max_iterations=args.iterations,
                force_tolerance=args.force_tolerance, callback=checkpoint, method=args.minimizer)
        final = core.evaluate(result['field'])
        mode = lowest_core_mode(core, result['field'], force_tolerance=args.force_tolerance, method=args.stability_method)
        result['final_stability_probe'] = {k:v for k,v in mode.items() if k != 'eigenvector'}
        save_field(out/'state.csv', result['field'], final['gradient'])
        write_csv(out/'radial_energy.csv', core.radial_energy(result['field'], np.arange(.5, args.radius+.01, .5)))
        save_json(out/'projected_circulation.json', phase_winding_cells(core, result['field']))
        write_csv(out/'adjacent_registry.csv', adjacent_registry_profile(core, result['field']))
        save_json(out/'summary.json', dict(**{k:v for k,v in result.items() if k != 'field'},
            completed=True, initial_energy=initial['energy'],
            initial_maximum_free_force=float(np.max(abs(initial['gradient']))),
            reciprocal_modes=final['reciprocal_modes'], last_mode_envelope=final['last_mode_envelope'],
            analytic_affine_hessian=core.affine_antiplane_hessian(), continuum_antiplane_hessian=far.matrix,
            actual_yield_calibrated=False, domain_converged=False, material_accepted=False,
            elapsed_seconds=time.perf_counter()-began))
        save_json(out/'progress.json', dict(completed=True, elapsed_seconds=time.perf_counter()-began))
        print(f'{out.name}: completed force={result["maximum_free_force"]:.7g}, '
              f'minH={mode["minimum_eigenvalue"]:.7g}, winding={result["boundary_winding"]:.5g}', flush=True)
    except Exception as exc:
        save_field(out/'checkpoint.csv', latest)
        save_json(out/'failure.json', dict(completed=False, error_type=type(exc).__name__, error=str(exc),
            elapsed_seconds=time.perf_counter()-began, actual_yield_calibrated=False))
        raise


def argument_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--material', type=Path, default=DEFAULT_MATERIAL)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--radius', type=float, default=2.)
    p.add_argument('--ring', type=int, default=3)
    p.add_argument('--shear-mpa', type=float, default=0.)
    p.add_argument('--center-y', type=float, default=0.)
    p.add_argument('--burgers-sign', type=int, choices=(-1, 0, 1), default=1)
    p.add_argument('--from-case', type=Path)
    p.add_argument('--configuration-seed', type=Path)
    p.add_argument('--resume-checkpoint', type=Path)
    p.add_argument('--extend-domain', action='store_true')
    p.add_argument('--allow-unstable-source', action='store_true')
    p.add_argument('--escape-sign', type=int, choices=(-1, 0, 1), default=0)
    p.add_argument('--escape-amplitude', type=float, default=.02)
    p.add_argument('--translate-j', type=int, default=0)
    p.add_argument('--iterations', type=int, default=150)
    p.add_argument('--max-seconds', type=float, default=1200.)
    p.add_argument('--force-tolerance', type=float, default=5e-7)
    p.add_argument('--tolerance', type=float, default=2e-12)
    p.add_argument('--minimizer', choices=('lbfgs', 'newton_cg', 'stationary_root'), default='lbfgs')
    p.add_argument('--stability-method', choices=('operator', 'assembled'), default='operator')
    return p


def main():
    run(argument_parser().parse_args())


if __name__ == '__main__':
    main()
