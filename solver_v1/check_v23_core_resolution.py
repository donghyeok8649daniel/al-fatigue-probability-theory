"""New-shape neighborhood/tolerance checks without another material fit.

Frozen-source force comparison uses the identical source atomic coordinates.
The candidate stationary-field checks use its own far field and energy. A
fixed-field check is not re-relaxation or infinite free-domain convergence.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .current_material_rows import CurrentMaterialScrewCore
from .frozen_core_force_target import FrozenCoreForceTarget
from .report_current_material_core import restore_case
from .report_core_research_v22 import read_csv, sampled_span
from .run_current_material_core import ROOT, load_current_material
from .run_source_core_reference import load_source_material, build_source_core
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh resolution study required')
    began = time.perf_counter()
    root = ROOT/'results/core_interface_compatibility_v23'
    material_path = root/'radial_full_validation/joint_LS/research_candidate_snapshot.json'
    model, tensor, binding = load_current_material(material_path)
    source, source_tensor, source_binding = load_source_material()
    source_path = ROOT/'results/current_material_core_v22/source_reference/escaped_plus_R8r5'
    source_core, source_field, *_ = restore_case(source_path, source, source_tensor,
        source_binding, core_builder=build_source_core)
    frozen = FrozenCoreForceTarget(source_core, source_field, radius=2.)
    candidate_path = root/'radial_joint_escape_plus_newton_R8r5'
    core, field, meta, summary, saved = restore_case(candidate_path, model, tensor, binding)
    if not summary['final_stability_probe']['stable_on_tested_fixed_boundary']:
        raise ValueError('separately verified reference core required')
    save_json(args.out/'definition.json', dict(parameter_sha256=binding['parameter_sha256'],
        source_sha256=source_binding['parameter_sha256'],
        source_state_sha256=hashlib.sha256((source_path/'state.csv').read_bytes()).hexdigest(),
        candidate_state_sha256=hashlib.sha256((candidate_path/'state.csv').read_bytes()).hexdigest(),
        frozen_source_rings=[5, 7, 9], frozen_source_tolerances=[2e-12, 2e-14],
        candidate_fixed_field_rings=[5, 7], material_fit_performed=False,
        static_research_only=True, infinite_domain_certified=False))
    force_rows, last = [], None
    for ring, tol in ((5, 2e-12), (7, 2e-12), (9, 2e-12), (9, 2e-14)):
        design = frozen.coefficient_matrix(model, ring=ring, tolerance=tol)
        value = design['design']@model.coefficients
        error = value-design['target']
        row = dict(ring=ring, reciprocal_tolerance=tol,
            force_RMS_eV_L0=float(np.sqrt(np.mean(error*error))),
            maximum_force_error_eV_L0=float(np.max(abs(error))),
            change_RMS_eV_L0=None if last is None else float(np.sqrt(np.mean((value-last)**2))),
            change_maximum_eV_L0=None if last is None else float(np.max(abs(value-last))),
            interpretation='frozen source configuration, no new fit')
        force_rows.append(row); last = value
        print('frozen source', row, flush=True)
    state_rows, previous = [], None
    for ring in (5, 7):
        trial = CurrentMaterialScrewCore(model, core.far_field, free_radius=8., ring=ring,
            tolerance=2e-12, shear_traction=core.shear_traction, burgers_sign=meta['burgers_sign'])
        q = np.array([saved[tuple(i)] for i in trial.indices[trial.free_ids]])
        value = trial.evaluate(q)
        actual = {tuple(i): g for i, g in zip(trial.indices[trial.free_ids], value['gradient'])}
        ordered = np.array([actual[tuple(i)] for i in core.indices[core.free_ids]])
        state_rows.append(dict(ring=ring, energy_eV_repeat=value['energy'],
            force_maximum_eV_L0=float(np.max(abs(ordered))),
            force_change_maximum_eV_L0=None if previous is None else float(np.max(abs(ordered-previous))),
            force_converged_on_this_operator=bool(np.max(abs(ordered)) <= meta['force_tolerance_eV_L0']),
            actual_re_relaxation=False, material_accepted=False))
        previous = ordered
    spans = []
    for name, path, b in (
        ('source_R8', source_path, source_core.rows.b),
        ('v22_wider_R8', ROOT/'results/current_material_core_v22/wider_core_near_R8r5', model.geometry.b),
        ('v23_radial_R8', candidate_path, model.geometry.b)):
        spans.append(dict(model=name, **sampled_span(read_csv(path/'adjacent_registry.csv'), b),
            used_for_fitting=False, continuum_core_radius=False, physical_yield_inferred=False))
    write_csv(args.out/'frozen_source_force_refinement.csv', force_rows)
    write_csv(args.out/'fixed_candidate_field_refinement.csv', state_rows)
    write_csv(args.out/'sampled_core_spans.csv', spans)
    save_json(args.out/'completion.json', dict(completed=True, elapsed_seconds=time.perf_counter()-began,
        fit_force_error_much_larger_than_last_resolution_change=bool(
            force_rows[2]['force_RMS_eV_L0'] > 100*force_rows[2]['change_RMS_eV_L0']),
        hundredfold_is_a_reported_ratio_check_not_a_physical_probability_threshold=True,
        infinite_domain_certified=False, material_accepted=False,
        actual_yield_validated=False, physical_Hz=False, production_changed=False))


if __name__ == '__main__':
    main()
