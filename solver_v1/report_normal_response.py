"""Compare completed material fits on one declared v17 jet dataset; no refit.

Old and new fitted losses used different data. This report explicitly replays
ALL named materials against the SAME new data before comparing residuals.
It does not turn a new validation failure into an optimization observation.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .interface_normal_development_targets import normal_development_observations
from .report_public_calibration_evidence import save_figure
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .validate_tail_calibration import load_material
from .vector_material_calibration import UNITS
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix


def grouped_validation(rows):
    """Return predeclared heldout errors, never a selected best subset."""
    result = []
    held = [r for r in rows if r['role'] == 'heldout']
    for field in ('energy', 'normal_force', 'Haa', 'Hxx'):
        selected = [r for r in held if r['observable'].endswith('_'+field)]
        if len(selected) != 12:
            raise ValueError('all twelve v17 heldout states required for every field')
        errors = np.array([r['normalized_residual'] for r in selected], float)
        if not np.isfinite(errors).all():
            raise ValueError('finite validation residuals required')
        result.append(dict(field=field, count=len(errors),
            normalized_rms=float(np.sqrt(np.mean(errors**2))),
            normalized_max=float(np.max(abs(errors))),
            count_within_declared_scale=int(np.sum(abs(errors) <= 1.)),
            scales_are_experimental_uncertainties=False))
    if len(held) != 48:
        raise ValueError('unrecognized heldout observation')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve earlier comparison')
    if len({d.name for d in args.models}) != len(args.models):
        raise ValueError('model names must be distinct')
    began = time.perf_counter()
    source, prior, states = source_and_targets()
    observations, provenance = normal_development_observations(source, prior, states)
    residuals, summary, grouped, budget, curves, bindings = [], [], [], [], [], []
    ratios = np.unique(np.r_[np.linspace(1., 1.2, 81), np.linspace(1.2, 4., 85)])
    source_curve = [source.evaluate((r*source.h, 0., 0.)) for r in ratios]
    for path in args.models:
        data, definition, model, cache_type, _ = load_material(path)
        if definition['source_sha256'] != source.reference.sha256:
            raise ValueError('same source potential/conditions required')
        raw = cache_type(observations).matrix(data['best']['decays'])
        matrix, obs = cubic_metric_problem(raw[:, :8], observations)
        extra = raw[:, 8:].copy()
        extra[2:5] = np.linalg.solve(cubic_to_mode_matrix(), extra[2:5])
        prediction = np.column_stack([matrix, extra])@data['best']['coefficients']
        records = [dict(model=path.name, observable=o.name, role=o.role,
            units=o.units, target=o.target, prediction=p, scale=o.scale,
            normalized_residual=(p-o.target)/o.scale)
            for o, p in zip(obs, prediction)]
        residuals.extend(records)
        groups = grouped_validation(records)
        grouped.extend(dict(model=path.name, **r) for r in groups)
        fit = [r['normalized_residual'] for r in records[5:] if r['role'] == 'fit']
        held = [r['normalized_residual'] for r in records if r['role'] == 'heldout']
        perfect = model.evaluate((model.h, 0., 0.))
        reference = source.evaluate((source.h, 0., 0.))
        summary.append(dict(model=path.name, current_fit_count=len(fit),
            current_fit_loss=float(np.dot(fit, fit)),
            heldout_count=len(held), heldout_normalized_rms=float(np.linalg.norm(held)/np.sqrt(len(held))),
            original_saved_fit_loss=data['best']['squared_loss'],
            original_target_generation='v17' if definition.get('normal_development') else 'v16',
            inherited_pair_control=bool(definition.get('fixed_pair_control')),
            perfect_Haa=perfect.hessian[0, 0], source_perfect_Haa=reference.hessian[0, 0],
            perfect_Hxx=perfect.hessian[1, 1], source_perfect_Hxx=reference.hessian[1, 1],
            perfect_Haa_relative_error=perfect.hessian[0, 0]/reference.hessian[0, 0]-1,
            full_material_accepted=False, physical_yield_validated=False))
        for label, state in [('perfect', states['perfect']), ('source_fault', states['fault']),
                             ('source_saddle', states['saddle']),
                             ('opening1.62', (1.62*model.h, 0., 0.))]:
            value, ref = model.evaluate(state), source.evaluate(state)
            for component, jet in value.components.items():
                budget.append(dict(model=path.name, state=label, component=component,
                    a=state[0], s1=state[1], s2=state[2],
                    energy_eV=jet[0], normal_force_eV_L0=jet[1], Haa=jet[4], Hxx=jet[7],
                    total_Haa=value.hessian[0, 0], source_Haa=ref.hessian[0, 0],
                    evaluation_at_same_fixed_coordinate=True))
        for ratio, reference in zip(ratios, source_curve):
            value = model.evaluate((ratio*model.h, 0., 0.))
            curves.append(dict(model=path.name, a_over_h=ratio,
                energy_J_m2=float(UNITS.energy_to_surface(value.energy)),
                source_energy_J_m2=float(UNITS.energy_to_surface(reference.energy)),
                traction_MPa=float(UNITS.force_to_traction_mpa(value.gradient[0])),
                source_traction_MPa=float(UNITS.force_to_traction_mpa(reference.gradient[0])),
                Haa=value.hessian[0, 0], source_Haa=reference.hessian[0, 0],
                rank1_Haa=value.components['angular_1'][4]))
        bindings.append(dict(model=path.name,
            calibration_sha256=hashlib.sha256((path/'calibration.json').read_bytes()).hexdigest(),
            definition_sha256=hashlib.sha256((path/'definition.json').read_bytes()).hexdigest()))
        print('same-dataset normal response replay', path.name, flush=True)
    for name, rows in [('comparison', summary), ('validation_groups', grouped),
                       ('all_residuals', residuals), ('curvature_components', budget),
                       ('normal_curves', curves)]:
        write_csv(args.out/(name+'.csv'), rows)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
    for ax, field, label in zip(axes, ('energy_J_m2', 'traction_MPa', 'Haa'),
                                ('Energy [J/m²]', 'Normal traction [MPa]', 'Haa [eV/L0²]')):
        for j, path in enumerate(args.models):
            rows = [r for r in curves if r['model'] == path.name]
            if j == 0:
                ax.plot(ratios, [r['source_'+field] for r in rows], 'k--', label='Al99 source (0 K)')
            ax.plot(ratios, [r[field] for r in rows], label=path.name)
        ax.set(xlabel='a/h (rigid opening)', ylabel=label)
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=6)
    fig.suptitle('v17: energies, forces AND curvatures; not a physical yield prediction')
    save_figure(fig, args.out/'normal_response.svg')
    save_json(args.out/'scope.json', dict(completed=True, bindings=bindings,
        source_sha256=source.reference.sha256, target_provenance=provenance,
        same_dataset_replay=True, optimizer_run_by_report=False,
        full_material_accepted=False, production_changed=False,
        physical_seconds=False, physical_Hz=False, elapsed_seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
