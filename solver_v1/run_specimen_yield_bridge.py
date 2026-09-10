"""Execute source-to-specimen strain accounting, not a fitted yield-stress law."""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from .fetch_public_strength_data import bounded_fetch
from .fcc111_geometry import fcc111_geometry_from_b
from .finite_source_reference import line_energy_coefficient, pinned_source_branch, solve_pinned_graph
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor
from .public_wire_validation import maximum_fcc_schmid
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .specimen_slip_kinematics import (
    pinned_swept_area, pinned_population_strain_budget, transport_history,
    axial_slip_strain, uniform_stress_work,
)
from .dislocation_line_kinetics import PinnedLineKinetics
from .strength_validation import stress_at_plastic_strain
from .vector_material_calibration import LENGTH_M
from .yield_reference_data import PIGATO_URL, parse_pigato_yield_xml


ROOT = Path(__file__).resolve().parents[1]
MATERIAL = ROOT/'results/fcc111_active_interface/tangent_calibration_v20/shape_refinement/old_family/calibration.json'
LINE_SCOPE = ROOT/'results/strength_fatigue_kinetics_v18/line_dynamics/scope.json'


def run(output, source_xml):
    if output.exists():
        raise FileExistsError('fresh output directory required; preserve completed work')
    began = time.perf_counter()
    if not source_xml.exists():
        raw = bounded_fetch(PIGATO_URL, 4*1024**2)
        parse_pigato_yield_xml(raw)  # validate before saving a downloaded original
        source_xml.parent.mkdir(parents=True, exist_ok=True)
        source_xml.write_bytes(raw)
    raw = source_xml.read_bytes()
    experimental = parse_pigato_yield_xml(raw)
    material_bytes = MATERIAL.read_bytes()
    material = json.loads(material_bytes)
    # Saved v20 predictions explicitly use the independent C11/C12/C44 metric.
    # This is replay of a previous material result, NOT another optimization.
    if not material['completed'] or not material['best']['strictly_positive_LJ']:
        raise ValueError('completed positive-LJ research reference required')
    predictions = np.asarray(material['best']['predictions'])
    np.testing.assert_allclose(predictions[:5], [0, 3.36, 114, 62, 32], atol=2e-10, rtol=0)
    constants = predictions[2:5]
    basis = fcc111_geometry_from_b(1.).plane_basis_in_stacked_cubic_axes()
    tensor = rotate_elastic_tensor(cubic_elastic_tensor(*(constants*1e9)), basis)
    lines = {}
    for count in (32, 64):
        line = line_energy_coefficient(tensor, [LENGTH_M, 0, 0], samples=count)
        lines[count] = replace(line, coefficients=line.coefficients*np.exp(1j*line.frequencies*np.pi/2))
    line = lines[64]
    angles = np.linspace(0, np.pi, 129)
    angular_error = float(np.max(abs(lines[32].evaluate(angles)-line.evaluate(angles))))
    spans = (.2, 1., 5.)
    orientations = {'100': [1, 0, 0], '110': [1, 1, 0], '111': [1, 1, 1]}
    geometry_status = 'hypothetical pins, R=L, core cutoff=2b; no fitted/measured source population'
    eta = .01  # declared spacing study ONLY; not an Al dislocation density
    budget_rows, sweep_rows, refinement_rows, histories = [], [], [], []
    for span_um in spans:
        L = span_um*1e-6
        logarithm = float(np.log(L/(2*LENGTH_M)))
        fold = pinned_swept_area(line, np.pi/2, span_m=L, outer_log_ratio=logarithm)
        critical = fold['critical_shear_outer_only_Pa']
        def angle_for(shear):
            return brentq(lambda t: pinned_source_branch(line, t, span_m=L,
                outer_log_ratio=logarithm, points=17)['applied_shear_Pa']-shear,
                1e-6, np.pi/2, xtol=2e-13)
        for label, direction in orientations.items():
            selected = maximum_fcc_schmid(direction)
            schmid = selected['schmid_factor']
            m = np.array(selected['slip_direction'], float); m /= np.linalg.norm(m)
            n = np.array(selected['slip_plane_normal'], float); n /= np.linalg.norm(n)
            e = np.array(direction, float); e /= np.linalg.norm(e)
            if (e@m)*(e@n) < 0:
                m = -m  # explicit positive resolved shear under axial tension
            budget = pinned_population_strain_budget(area_per_source_m2=fold['area_m2'],
                span_m=L, burgers_m=LENGTH_M, schmid_factor=schmid, axial_criterion=.002)
            budget_rows.append(dict(span_um=span_um, loading_axis=label, schmid_factor=schmid,
                critical_resolved_shear_MPa=critical/1e6,
                conditional_axial_first_source_MPa=critical/schmid/1e6,
                fold_area_over_L_squared=fold['area_over_span_squared'],
                illustrative_overlap=eta, illustrative_max_axial_bow_strain=eta*budget['axial_strain_per_unit_overlap'],
                **budget, geometry_status=geometry_status))
            # One hypothetical source in an explicitly specified cubic averaging volume.
            volume = L**3/eta
            geometry = dict(burgers_vectors_m=[LENGTH_M*m], plane_normals=[n], specimen_volume_m3=volume)
            for sigma in (0., 2., 5., 10., 15., 20., 30., 50.):
                shear = sigma*1e6*schmid
                area = None
                if shear == 0:
                    area = 0.
                elif shear < critical:
                    area = pinned_swept_area(line, angle_for(shear), span_m=L,
                        outer_log_ratio=logarithm)['area_m2']
                sweep_rows.append(dict(span_um=span_um, loading_axis=label,
                    axial_stress_MPa=sigma, resolved_shear_MPa=shear/1e6,
                    critical_resolved_shear_MPa=critical/1e6,
                    hypothetical_specimen_volume_m3=volume, initial_pinned_line_density_m2=eta/L**2,
                    swept_area_m2=area, axial_bow_strain=None if area is None else schmid*LENGTH_M*area/volume,
                    status='subcritical recoverable bow' if area is not None else 'no subcritical branch; emission not implemented',
                    experimental_yield_MPa=None))
            fractions = np.tile([0, .2, .4, .6, .8, .4, 0, -.4, -.8, -.4, 0], 3)
            areas = np.array([0. if f == 0 else np.sign(f)*pinned_swept_area(line,
                angle_for(abs(f)*critical), span_m=L, outer_log_ratio=logarithm)['area_m2'] for f in fractions])
            history = transport_history(areas[:, None], **geometry)
            axial = axial_slip_strain(history['signed_slip_strain_tensor'], e)
            stresses = fractions*critical/schmid
            midpoint = .5*(stresses[1:]+stresses[:-1])[:, None, None]*np.outer(e, e)
            work = uniform_stress_work(midpoint, np.diff(areas)[:, None], **geometry)
            # Only the increasing first segment goes to a loading-branch criterion.
            criterion = stress_at_plastic_strain(stresses[:5]/1e6, axial[:5], .002,
                                                 maximum_bracket_width=.0001)
            assert not criterion['resolved']  # bow does not reach even this offset demand
            for index, f in enumerate(fractions):
                histories.append(dict(span_um=span_um, loading_axis=label, static_step=index,
                    axial_stress_MPa=float(stresses[index]/1e6), signed_area_m2=areas[index],
                    axial_slip_strain=axial[index], forward_area_m2=history['forward_area_m2'][index, 0],
                    backward_area_m2=history['backward_area_m2'][index, 0],
                    gross_area_m2=history['gross_area_m2'][index, 0],
                    area_balance_residual_m2=history['area_balance_residual_m2'][index, 0],
                    work_identity_max_residual_J=float(np.max(abs(work['residual_J']))),
                    physical_hold_performed=False, residual_plasticity_certified=False))
        for angle in (.4, 1., 1.4):
            exact = pinned_swept_area(line, angle, span_m=L, outer_log_ratio=logarithm)
            tighter = pinned_swept_area(line, angle, span_m=L, outer_log_ratio=logarithm, relative_tolerance=1e-11)
            for segments in (32, 64, 128, 256):
                discrete = solve_pinned_graph(line, span_m=L, outer_log_ratio=logarithm,
                    shear_Pa=exact['applied_shear_Pa'], segments=segments)
                if not discrete['converged']:
                    raise RuntimeError('independent line-equilibrium check failed')
                refinement_rows.append(dict(span_um=span_um, endpoint_angle=angle, segments=segments,
                    shear_MPa=exact['applied_shear_Pa']/1e6, analytic_area_m2=exact['area_m2'],
                    discrete_area_m2=discrete['swept_area_m2'],
                    area_error_over_L_squared=abs(discrete['swept_area_m2']-exact['area_m2'])/L**2,
                    quadrature_refinement_change_m2=abs(exact['area_m2']-tighter['area_m2']),
                    dimensionless_force_residual=discrete['dimensionless_force_residual']))
    # Independent numerical time hold of the EXISTING linear pinned-line model.
    # Display dimensionless t/tau; the reused drag does NOT grant an a/s clock.
    scope_bytes = LINE_SCOPE.read_bytes()
    scope = json.loads(scope_bytes)
    L = 1e-6; logratio = np.log(L/(2*LENGTH_M))
    kinetic = PinnedLineKinetics(scope['conditional_drag_Pa_s'],
        float(line.stiffness(0)*logratio), LENGTH_M, L)
    geometry = dict(burgers_vectors_m=[[LENGTH_M, 0, 0]], plane_normals=[[0, 1, 0]], specimen_volume_m3=L**3/eta)
    direction = np.array([1., 1., 0.])/np.sqrt(2)
    hold_rows = []
    for steps in (64, 128, 256):
        scaled = np.linspace(0, 12*np.pi, 6*steps+1)
        hold = np.arange(1, int(np.ceil(24/(2*np.pi/steps)))+1)*(2*np.pi/steps)
        times = np.r_[scaled, scaled[-1]+hold]
        shear = np.r_[.01e6*np.sin(scaled), np.zeros(len(hold))]
        for segments in (32, 64, 128):
            result = kinetic.integrate_backward_euler(times*kinetic.slowest_seconds, shear, segments=segments)
            areas = result['mean_bow_m']*L
            history = transport_history(areas[:, None], **geometry)
            axial = axial_slip_strain(history['signed_slip_strain_tensor'], direction)
            hold_rows.append(dict(steps_per_period=steps, segments=segments,
                hold_time_over_line_tau=float(hold[-1]), max_axial_bow_strain=float(max(abs(axial))),
                final_axial_bow_strain=float(axial[-1]), final_to_peak=float(abs(axial[-1])/max(abs(axial))),
                max_slope=float(max(result['maximum_slope'])), final_gross_area_m2=float(history['gross_area_m2'][-1, 0]),
                permanent_slip_generated=result['permanent_slip_generated'],
                physical_time_scope='conditional v18 line-only; t/tau shown; not production a/s',
                experimental_yield_MPa=None))
    for name, rows in [('reported_yield_references', experimental), ('strain_budget', budget_rows),
                       ('stress_sweep', sweep_rows), ('static_cycle_transport', histories),
                       ('spatial_refinement', refinement_rows), ('unload_hold', hold_rows)]:
        write_csv(output/f'{name}.csv', rows)
    assessment = dict(completed=True, elapsed_seconds=time.perf_counter()-began,
        material_path=MATERIAL.relative_to(ROOT).as_posix(), material_sha256=hashlib.sha256(material_bytes).hexdigest(),
        material_accepted=material['best']['material_accepted'], cubic_GPa=constants,
        material_temperature_K=0, experimental_source_url=PIGATO_URL,
        source_xml_sha256=hashlib.sha256(raw).hexdigest(), source_line_scope_sha256=hashlib.sha256(scope_bytes).hexdigest(),
        angular32_to64_coefficient_error_J_m=angular_error, source_spans_um=spans,
        geometry_status=geometry_status, hypothetical_overlap=eta, target_axial_criterion=.002,
        minimum_required_overlap=min(r['required_overlap_parameter'] for r in budget_rows),
        maximum_required_overlap=max(r['required_overlap_parameter'] for r in budget_rows),
        actual_yield_prediction_MPa=None, core_validated=False, source_population_calibrated=False,
        emission_and_interactions_implemented=False, strength_reduction_factor_fitted=False,
        production_changed=False, production_M_a_phys=None, production_M_s_phys=None, production_t0_seconds=None,
        conclusion='MPa first-source scale is not proof yield; independent recoverable bows fail the strain/persistence budget')
    save_json(output/'assessment.json', assessment)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    for label in orientations:
        rows = [r for r in budget_rows if r['loading_axis'] == label]
        axes[0].loglog(spans, [r['conditional_axial_first_source_MPa'] for r in rows], 'o-', label='['+label+']')
        axes[1].loglog(spans, [r['required_overlap_parameter'] for r in rows], 'o-', label='['+label+']')
    axes[0].set(xlabel='Hypothetical pin span [um]', ylabel='Outer-only first-source axial stress [MPa]', title='Not specimen yield')
    axes[1].axhline(1, color='gray', linestyle='--', label='spacing comparable to span')
    axes[1].set(xlabel='Hypothetical pin span [um]', ylabel='Required source overlap for axial strain 0.002', title='Recoverable-bow strain demand')
    for ax in axes:
        ax.set_xticks(spans, labels=['0.2', '1', '5'])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=.25); ax.legend(fontsize=7)
    fig.savefig(output/'source_vs_proof_strain.png', dpi=160)
    plt.close(fig)
    print(json.dumps(assessment, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--source-xml', type=Path, required=True, help='cached original, fetched if absent; keep outside committed results')
    args = parser.parse_args()
    run(args.out, args.source_xml)
