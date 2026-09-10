"""Replay completed v22 states under their bound full per-atom energy.

This report does not rerun relaxation or calibrate material/kinetics. Compare
stable branches, neighborhoods, inner fields and stresses without turning a
fixed-boundary energy into a finite-loop activation or specimen yield.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .core_matching import finite_part_samples
from .current_core_diagnostics import phase_winding_cells, adjacent_registry_profile
from .current_material_rows import CurrentMaterialScrewCore
from .isolated_screw_core import ScrewFarField
from .run_current_material_core import load_current_material, EV_J, ROOT
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def restore_case(path, model, tensor, metadata, *, core_builder=CurrentMaterialScrewCore):
    meta = json.loads((path/'metadata.json').read_bytes())
    summary = json.loads((path/'summary.json').read_bytes())
    if not summary['completed'] or meta['parameter_sha256'] != metadata['parameter_sha256']:
        raise ValueError('completed state with the same bound material required')
    with (path/'state.csv').open(encoding='utf8', newline='') as stream:
        data = list(csv.DictReader(stream))
    state = {(int(r['j']), int(r['l'])): np.array([float(r[k]) for k in ('ux', 'uy', 'uz')]) for r in data}
    far = ScrewFarField(tensor, model.geometry.b, tuple(meta['center_over_L0']))
    traction = meta['shear_traction_MPa']*1e6*metadata['length_scale_m']**3/EV_J
    core = core_builder(model, far, free_radius=meta['radius_over_L0'], ring=meta['ring'],
        tolerance=meta['reciprocal_tolerance'], shear_traction=traction, burgers_sign=meta['burgers_sign'])
    if set(state) != {tuple(i) for i in core.indices[core.free_ids]}:
        raise ValueError('state indices differ from the declared free disk')
    field = np.array([state[tuple(i)] for i in core.indices[core.free_ids]])
    return core, field, meta, summary, state


def compare_states(first, second, *, radius, b, h):
    ma, a = first; mb, other = second
    if ma['parameter_sha256'] != mb['parameter_sha256']:
        raise ValueError('different materials cannot certify domain convergence')
    if not np.allclose(ma['center_over_L0'], mb['center_over_L0'], atol=1e-13, rtol=0):
        raise ValueError('different prescribed centers are a boundary study, not identical-domain refinement')
    center = np.asarray(ma['center_over_L0'])
    records = []
    for index, u in a.items():
        j, layer = index
        yz = np.array([np.sqrt(3)*b/2*(j+layer/3), h*layer])
        if index not in other or np.linalg.norm(yz-center) >= radius:
            continue
        diff = other[index]-u
        diff[0] -= b*np.floor(diff[0]/b+.5)
        records.append(diff)
    if not records:
        raise ValueError('no common interior rows')
    records = np.asarray(records)
    return dict(common_rows=len(records), comparison_radius_over_L0=radius,
        maximum_vector_change_over_L0=float(np.max(np.linalg.norm(records, axis=1))),
        rms_vector_change_over_L0=float(np.sqrt(np.mean(np.sum(records**2, axis=1)))),
        maximum_screw_change_over_L0=float(np.max(abs(records[:, 0]))),
        maximum_transverse_change_over_L0=float(np.max(np.linalg.norm(records[:, 1:], axis=1))))


def report(root, out, *, pairs=None, material_loader=load_current_material,
           core_builder=CurrentMaterialScrewCore):
    if out.exists():
        raise FileExistsError('fresh report directory required')
    began = time.perf_counter()
    model, tensor, metadata = material_loader()
    summaries, checks, profiles, circulation, radial = [], [], [], [], []
    stored = {}; other_materials = []
    for path in sorted(root.iterdir()):
        if not (path/'metadata.json').exists():
            continue
        if not (path/'summary.json').exists():
            continue  # incomplete cases are listed separately below, never PASS
        declared = json.loads((path/'metadata.json').read_bytes())
        if declared['parameter_sha256'] != metadata['parameter_sha256']:
            other_materials.append(dict(case=path.name,parameter_sha256=declared['parameter_sha256']))
            continue  # explicitly listed, never mixed into convergence pairs
        core, field, meta, saved, state = restore_case(path, model, tensor, metadata, core_builder=core_builder)
        actual = core.evaluate(field)
        error = abs(actual['energy']-saved['energy'])
        if error > 2e-9:
            raise ArithmeticError(f'saved energy failed replay: {path.name}')
        stored[path.name] = meta, state
        probe = saved['final_stability_probe']
        summaries.append(dict(case=path.name, radius=meta['radius_over_L0'], ring=meta['ring'],
            shear_MPa=meta['shear_traction_MPa'], free_rows=len(field),
            energy_eV_per_row_repeat=saved['energy'], maximum_force_eV_L0=saved['maximum_free_force'],
            minimum_H_eV_L0sq=probe['minimum_eigenvalue'], eigenpair_residual=probe['eigenpair_residual'],
            stable_fixed_boundary=probe['stable_on_tested_fixed_boundary'],
            method=saved['minimizer'], iterations=saved['iterations'], elapsed_seconds=saved['elapsed_seconds'],
            winding=core.boundary_winding(field), energy_replay_error_eV=error,
            affine_bulk_error_eV_L0cubed=float(np.max(abs(core.affine_antiplane_hessian()-core.far_field.matrix))),
            reciprocal_modes=saved['reciprocal_modes'], last_mode_envelope=saved['last_mode_envelope'],
            material_accepted=False, actual_yield_calibrated=False))
        for check in probe['energy_curvature_checks']:
            checks.append(dict(case=path.name, minimum_H_eV_L0sq=probe['minimum_eigenvalue'], **check))
        for row in adjacent_registry_profile(core, field):
            profiles.append(dict(case=path.name, **row))
        projected = phase_winding_cells(core, field)
        for row in projected['cells']:
            circulation.append(dict(case=path.name, **row))
        if meta['shear_traction_MPa'] == 0 and core.free_radius >= 2:
            distance = np.linalg.norm(core.xyz[core.energy_ids, 1:]-core.far_field.center, axis=1)
            for row in finite_part_samples(distance, actual['quadratic_remainder_site_energy'],
                    np.arange(2., core.free_radius+.01, .5), log_coefficient=core.far_field.log_energy_coefficient,
                    b=core.rows.b, free_radius=core.free_radius):
                radial.append(dict(case=path.name, **row))
        print(f'replayed {path.name}: E error={error:.3g}', flush=True)
    comparisons = []
    for first, second, kind in (pairs or []):
        if first not in stored or second not in stored:
            raise ValueError(f'declared comparison missing: {first}/{second}')
        for radius in (1., 1.5, 2.):
            if radius <= min(stored[first][0]['radius_over_L0'], stored[second][0]['radius_over_L0']):
                comparisons.append(dict(first=first, second=second, study=kind,
                    **compare_states(stored[first], stored[second], radius=radius, b=model.geometry.b, h=model.h)))
    write_csv(out/'core_summary.csv', summaries)
    write_csv(out/'curvature_checks.csv', checks)
    write_csv(out/'registry_profiles.csv', profiles)
    write_csv(out/'projected_circulation.csv', circulation)
    write_csv(out/'inner_field_changes.csv', comparisons)
    write_csv(out/'finite_part_samples.csv', radial)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    selected = [r for r in summaries if r['stable_fixed_boundary']]
    for row in selected:
        name = row['case']
        p = [r for r in profiles if r['case'] == name and r['both_rows_free']]
        axes[0].plot([r['y_over_L0'] for r in p], [r['affine_subtracted_delta_x'] for r in p],
                     '.-', label=name)
        if row['shear_MPa'] == 0:
            p = [r for r in radial if r['case'] == name and r['inner_fraction'] == .75]
            axes[1].plot([r['radius_over_L0'] for r in p],
                [r['finite_part_at_reference_b_eV'] for r in p], '.-', label=name)
    axes[0].set(xlabel='Adjacent-row midpoint y/L0', ylabel='Affine-subtracted screw displacement / L0',
                title='Actual three-component relaxed row state')
    axes[1].set(xlabel='Integration radius R/L0', ylabel='Finite part [eV/straight-row repeat]',
                title='Same-energy logarithm; no fitted core radius')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(fontsize=5)
    fig.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'; plt.rcParams['svg.hashsalt'] = 'current-core-v22'
    stream = io.StringIO(); fig.savefig(stream, format='svg', metadata={'Date': None})
    (out/'current_core.svg').write_text('\n'.join(s.rstrip() for s in stream.getvalue().splitlines())+'\n', encoding='utf8')
    cache = ROOT/'.cache/current_material_core_v22'; cache.mkdir(parents=True, exist_ok=True)
    fig.savefig(cache/(out.name+'.png'), dpi=140); plt.close(fig)
    incomplete = [path.name for path in root.iterdir()
                  if (path/'metadata.json').exists() and not (path/'summary.json').exists()
                  and not (path/'saddle_summary.json').exists()]
    save_json(out/'decision.json', dict(completed=True, completed_cases=len(summaries),
        incomplete_cases=incomplete, other_material_cases_not_replayed=other_materials,
        actual_relaxation_reexecuted=False, energies_replayed=True,
        all_energy_terms_retained=True, parameter_sha256=metadata['parameter_sha256'],
        actual_yield_calibrated=False, material_accepted=False,
        infinite_transverse_domain_certified=False, finite_loop_activation_validated=False,
        physical_a_s_mobility_available=False, physical_seconds=False, physical_Hz=False,
        production_changed=False, elapsed_seconds=time.perf_counter()-began,
        files_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.csv')}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT/'results/current_material_core_v22')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--pairs', type=Path)
    p.add_argument('--source-reference', action='store_true')
    p.add_argument('--material', type=Path)
    args = p.parse_args()
    extra = {}
    if args.source_reference:
        if args.material:
            raise ValueError('source and candidate material selectors are distinct')
        from .run_source_core_reference import load_source_material, build_source_core
        extra.update(material_loader=load_source_material, core_builder=build_source_core)
    elif args.material:
        from functools import partial
        extra['material_loader'] = partial(load_current_material,args.material)
    report(args.root, args.out, pairs=json.loads(args.pairs.read_bytes()) if args.pairs else None, **extra)


if __name__ == '__main__':
    main()
