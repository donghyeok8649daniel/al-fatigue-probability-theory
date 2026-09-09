"""Execute, not replay, the matched full-vector registry/static-stress audit.

Saved candidate unchanged. No optimization of material parameters. Source EAM
is a benchmark only. Uniform-interface saddle energies are NOT finite-defect
activation energies, and the static first fold is NOT experimental yield.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from .reference_eam_targets import MishinRigidFCCReference, SOURCE_URL
from .run_low_stress_cyclic_diagnostic import build_surface, write_csv
from .vector_interface_reference import FullRegistryInterface, MishinVectorInterfaceReference
from .vector_registry_audit import stationary_state, relax_at_x, saddle_downhill_endpoints


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/fcc111_active_interface/vector_registry_v9'


def save_json(path, value):
    def plain(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=plain, allow_nan=False)+'\n', encoding='utf-8')


def point_row(name, label, q, out, units):
    eigen = np.linalg.eigvalsh(out.hessian)
    row = dict(model=name, state=label, a_L0=q[0], ux_L0=q[1], uy_L0=q[2],
               energy_eV_cell=out.energy, energy_J_m2=float(units.energy_to_surface(out.energy)),
               gradient_norm_eV_L0=float(np.linalg.norm(out.gradient)),
               eigenvalue_min_eV_L0sq=eigen[0], eigenvalue_middle_eV_L0sq=eigen[1],
               eigenvalue_max_eV_L0sq=eigen[2])
    for i, axis in enumerate(('a', 'x', 'y')):
        row['force_'+axis+'_eV_L0'] = out.gradient[i]
        row['traction_'+axis+'_MPa'] = float(units.force_to_traction_mpa(out.gradient[i]))
        for j, axis2 in enumerate(('a', 'x', 'y')):
            row['H_'+axis+axis2+'_eV_L0sq'] = out.hessian[i, j]
    return row


def static_points(model, name, units, tau):
    rows = []; roots = {}
    for label, u in [('perfect', np.zeros(2)), ('direct_half', [.5, 0.]),
                     ('partial_saddle_guess', .6*tau), ('intrinsic_fault', tau), ('on_top', 2*tau)]:
        q = np.r_[model.h, u]
        rows.append(point_row(name, label+'_fixed_a', q, model.evaluate(q), units))
        if label == 'direct_half':
            continue  # not a full stationary state; DO NOT disguise omitted force
        expected = {'perfect': 0, 'partial_saddle_guess': 1, 'intrinsic_fault': 0, 'on_top': 2}[label]
        start = np.r_[(1.2 if label == 'on_top' else 1.03)*model.h, u]
        result = stationary_state(model, start, expected_index=expected)
        if not result['valid']:
            raise RuntimeError(f'{name} {label}: stationary residual/index validation failed')
        roots[label] = result
        rows.append(point_row(name, label+'_stationary', result['q'], result['evaluation'], units))
    saddle = roots['partial_saddle_guess']
    for i, endpoint in enumerate(saddle_downhill_endpoints(model, saddle)):
        rows.append(point_row(name, f'saddle_downhill_{i}', endpoint['q'], endpoint['evaluation'], units))
    return rows, roots


def relaxed_registry_curve(model, name, units, n, *, normal_mpa=0.):
    """Fixed x, relax a/y; actual local branch and independent first-fold solve."""
    rows = []; states = []; z = None
    fn = float(units.traction_mpa_to_force(normal_mpa))
    for i, x in enumerate(np.linspace(0., .5, n)):
        state = relax_at_x(model, x, z, normal_force=fn)
        q = state['q']; out = state['evaluation']; z = q[[0, 2]]
        states.append(state)
        row = point_row(name, 'relaxed_x_branch', q, out, units)
        row.update(samples=n, normal_load_MPa=normal_mpa,
                   schur_x_eV_L0sq=state['schur_curvature'],
                   transverse_force_residual_eV_L0=state['transverse_force_residual'],
                   eliminated_eigenvalue_min_eV_L0sq=state['eliminated_eigenvalues'][0])
        rows.append(row)
    crossings = [i for i in range(1, n) if states[i-1]['schur_curvature'] > 0 >= states[i]['schur_curvature']]
    if not crossings:
        raise RuntimeError('first static fold not bracketed; no inferred spinodal')
    i = crossings[0]; left, right = states[i-1], states[i]
    def at(x):
        t = (x-left['q'][1])/(right['q'][1]-left['q'][1])
        guess = ((1-t)*left['q']+t*right['q'])[[0, 2]]
        return relax_at_x(model, x, guess, normal_force=fn)
    root_x = brentq(lambda x: at(x)['schur_curvature'], left['q'][1], right['q'][1],
                    xtol=2e-11, rtol=2e-12)
    state = at(root_x); out = state['evaluation']
    fold = point_row(name, 'first_uniform_interface_fold', state['q'], out, units)
    fold.update(samples=n, normal_load_MPa=normal_mpa,
        schur_x_eV_L0sq=state['schur_curvature'],
        transverse_force_residual_eV_L0=state['transverse_force_residual'],
        eliminated_eigenvalue_min_eV_L0sq=state['eliminated_eigenvalues'][0],
        interpretation='athermal rigid-half interface fold; NOT experimental yield or finite-defect threshold')
    return rows, fold


def stress_probes(model, name, units):
    """Prescribed independent normal and two shear components, MPa; static."""
    protocols = {
        'pure_x_shear': [(0, 0, 0), (0, 25, 0), (0, 50, 0), (0, 150, 0), (0, 0, 0), (0, -50, 0), (0, 0, 0)],
        'mixed_traction': [(0, 0, 0), (25, 25, 15), (50, 50, 25), (150, 100, -50), (0, 0, 0)],
        'compression_shear': [(0, 0, 0), (-50, 50, 15), (-150, 100, 25), (0, 0, 0)],
    }
    rows = []
    for protocol, loads in protocols.items():
        q = np.array([model.h, 0., 0.]); baseline = None
        for step, load in enumerate(loads):
            force = units.traction_mpa_to_force(load)
            result = stationary_state(model, q, force=force, expected_index=0)
            if not result['valid']:
                raise RuntimeError(f'low-stress branch failed for {name}, {protocol}, {load}')
            q = result['q']
            if baseline is None:
                baseline = q.copy()
            row = point_row(name, protocol, q, result['evaluation'], units)
            row.update(step=step, normal_load_MPa=load[0], shear_x_MPa=load[1], shear_y_MPa=load[2],
                       equilibrium_residual_eV_L0=result['force_residual'],
                       du_norm_L0=float(np.linalg.norm(q[1:]-baseline[1:])),
                       da_from_initial_L0=q[0]-baseline[0],
                       static_unload=bool(step and not any(load)), physical_time_available=False)
            rows.append(row)
    return rows


def derivative_rows(model, name):
    rows = []
    for a, x, y in [(1.017, .173, .079), (1.21, -.217, .341), (1.51, .307, -.071)]:
        q = np.array([a*model.h, x, y]); out = model.evaluate(q)
        for step in (2e-5, 1e-5):
            p = [model.evaluate(q+step*e) for e in np.eye(3)]
            m = [model.evaluate(q-step*e) for e in np.eye(3)]
            grad = np.array([(u.energy-v.energy)/(2*step) for u, v in zip(p, m)])
            H = np.column_stack([(u.gradient-v.gradient)/(2*step) for u, v in zip(p, m)])
            rows.append(dict(model=name, a_L0=q[0], ux_L0=x, uy_L0=y, step_L0=step,
                gradient_max_abs_error_eV_L0=float(np.max(abs(grad-out.gradient))),
                hessian_max_abs_error_eV_L0sq=float(np.max(abs(H-out.hessian))),
                gradient_relative_norm_error=float(np.linalg.norm(grad-out.gradient)/np.linalg.norm(out.gradient)),
                hessian_relative_norm_error=float(np.linalg.norm(H-out.hessian)/np.linalg.norm(out.hessian))))
    return rows


def direct_rows(model):
    rows = []
    for q in ([1.09*model.h, .237, .087], [model.h, *model.geometry.tau]):
        expected = model.evaluate(q)
        face = copy.copy(model.face)
        face._active_delta = lambda k, s: face._baseline_delta(k)+np.array(q[1:])+s*face.direction
        for radius, layers in ((24, 24), (48, 48), (72, 72)):
            out = face.direct_reference(q[0], 0., radial_index=radius, layers=layers)
            e = out.energy
            for amplitude, plane in model.moments:
                moment = copy.copy(plane.invariant); moment.interface = face
                e += amplitude*moment.direct_value(q[0], 0., radius=radius, layers=layers)
            rows.append(dict(a_L0=q[0], ux_L0=q[1], uy_L0=q[2], direct_radial_index=radius,
                direct_layers=layers, direct_energy_eV_cell=e, reciprocal_energy_eV_cell=expected.energy,
                abs_error_eV_cell=abs(e-expected.energy), relative_error=abs(e-expected.energy)/abs(expected.energy)))
    return rows


def plots(folder, grid, curves, point_rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for name, marker in [('analytic_candidate', '-'), ('Mishin_source_only', '--')]:
        sel = [r for r in curves if r['model'] == name and r['samples'] == max(t['samples'] for t in curves)]
        axes[0].plot([r['ux_L0'] for r in sel], [r['uy_L0'] for r in sel], marker, label=name)
        axes[1].plot([r['ux_L0'] for r in sel], [r['energy_J_m2'] for r in sel], marker, label=name)
        axes[2].plot([r['ux_L0'] for r in sel], [r['traction_x_MPa']/1000 for r in sel], marker, label=name)
    axes[0].plot([0, .5], [0, 0], ':', label='old constrained direct line')
    axes[0].set_ylabel('relaxed registry y / L0')
    axes[1].set_ylabel('interface energy [J/m2]')
    axes[2].set_ylabel('required uniform shear traction [GPa]')
    for ax in axes:
        ax.set_xlabel('prescribed registry x / L0'); ax.grid(alpha=.2)
    axes[0].legend(fontsize=7); fig.suptitle('Full vector local branch; static research, NOT Al yield or fatigue')
    fig.tight_layout(); fig.savefig(folder/'relaxed_registry_and_traction.png', dpi=150); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, name in zip(axes, ('analytic_candidate', 'Mishin_source_only')):
        sel = [r for r in grid if r['model'] == name]
        colors = ax.tricontourf([r['ux_L0'] for r in sel], [r['uy_L0'] for r in sel],
                               [r['energy_J_m2'] for r in sel], levels=25)
        fig.colorbar(colors, ax=ax, label='fixed-spacing GSF [J/m2]')
        ax.plot([0, .5, 1], [0, np.sqrt(3)/6, 0], 'w--', lw=1)
        ax.set(xlabel='registry x / L0', ylabel='registry y / L0', title=name)
        ax.set_aspect('equal')
    fig.tight_layout(); fig.savefig(folder/'vector_gamma_surface.png', dpi=150); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--grid', type=int, default=21)
    parser.add_argument('--curve-samples', type=int, default=41)
    args = parser.parse_args(); folder = args.output
    if (folder/'scientific_status.json').exists():
        raise FileExistsError('completed result preserved; use a new --output for an independent rerun')
    if args.grid < 5 or args.curve_samples < 9:
        raise ValueError('at least 5x5 vector samples and 9 branch samples required')
    started = time.perf_counter()
    surface, units, meta = build_surface(tolerance=2e-11)
    analytic = FullRegistryInterface(surface)
    source = MishinVectorInterfaceReference(MishinRigidFCCReference(), units.length_scale_m/1e-10)
    save_json(folder/'run_manifest.json', dict(completed=False, parameter_metadata=meta,
        source_url=SOURCE_URL, source_sha256=source.reference.sha256,
        source_status='matched 0-K rigid Mishin-EAM benchmark, NOT experimental truth',
        coordinate_order=['a', 'ux', 'uy'], source_code_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__), Path(__file__).with_name('vector_interface_reference.py'),
             Path(__file__).with_name('vector_registry_audit.py')]}))
    all_points = []; all_curves = []; folds = []; all_stress = []; grid = []; derivative = []; barriers = []
    for name, model in [('analytic_candidate', analytic), ('Mishin_source_only', source)]:
        print(name, 'stationary states / downhill connectivity', flush=True)
        points, roots = static_points(model, name, units, analytic.geometry.tau)
        all_points.extend(points)
        p, f, s = (roots[k]['evaluation'].energy for k in ('perfect', 'intrinsic_fault', 'partial_saddle_guess'))
        direct = model.evaluate([model.h, .5, 0.]).energy
        barriers.append(dict(model=name, forward_saddle_J_m2=float(units.energy_to_surface(s-p)),
            reverse_from_fault_J_m2=float(units.energy_to_surface(s-f)), intrinsic_fault_relaxed_J_m2=float(units.energy_to_surface(f-p)),
            constrained_direct_midpoint_J_m2=float(units.energy_to_surface(direct)),
            direct_over_connected_saddle=direct/(s-p), interpretation='per-interface-area static barriers, not finite defect activation'))
        write_csv(folder/'stationary_points.csv', all_points); write_csv(folder/'barrier_comparison.csv', barriers)
        for n in (args.curve_samples, 2*args.curve_samples-1):
            print(name, 'relaxed x branch', n, flush=True)
            rows, fold = relaxed_registry_curve(model, name, units, n)
            all_curves.extend(rows); folds.append(fold)
            write_csv(folder/'relaxed_registry_branches.csv', all_curves); write_csv(folder/'uniform_interface_folds.csv', folds)
        print(name, 'low-MPa mixed tractions and static unload', flush=True)
        all_stress.extend(stress_probes(model, name, units)); write_csv(folder/'stress_and_unload.csv', all_stress)
        derivative.extend(derivative_rows(model, name)); write_csv(folder/'derivative_validation.csv', derivative)
        for i, xi in enumerate(np.linspace(0., 1., args.grid)):
            for eta in np.linspace(0., 1., args.grid):
                u = xi*analytic.geometry.a1+eta*analytic.geometry.a2
                q = np.r_[model.h, u]; val = model.evaluate(q)
                grid.append(dict(model=name, fraction_a1=xi, fraction_a2=eta, a_L0=model.h,
                    ux_L0=u[0], uy_L0=u[1], energy_J_m2=float(units.energy_to_surface(val.energy)),
                    normal_force_eV_L0=val.gradient[0], registry_x_force_eV_L0=val.gradient[1],
                    registry_y_force_eV_L0=val.gradient[2],
                    min_registry_hessian_eV_L0sq=float(np.linalg.eigvalsh(val.hessian[1:, 1:])[0])))
            if i % 5 == 0:
                print(name, 'registry grid row', i+1, '/', args.grid, flush=True)
        write_csv(folder/'vector_energy_grid.csv', grid)
    print('independent real-space tail checks', flush=True)
    direct = direct_rows(analytic); write_csv(folder/'direct_lattice_validation.csv', direct)
    tight = FullRegistryInterface(build_surface(tolerance=2e-13)[0]); refinement = []
    for state in [r for r in all_points if r['model'] == 'analytic_candidate' and 'stationary' in r['state']]:
        q = [state['a_L0'], state['ux_L0'], state['uy_L0']]
        coarse, fine = analytic.evaluate(q), tight.evaluate(q)
        refinement.append(dict(state=state['state'], tolerance_coarse=2e-11, tolerance_fine=2e-13,
            energy_error_eV_cell=abs(coarse.energy-fine.energy),
            gradient_max_error_eV_L0=float(max(abs(coarse.gradient-fine.gradient))),
            hessian_max_error_eV_L0sq=float(np.max(abs(coarse.hessian-fine.hessian)))))
    write_csv(folder/'reciprocal_neighborhood_refinement.csv', refinement)
    plots(folder, grid, all_curves, all_points)
    summary = dict(completed=True, wall_seconds=time.perf_counter()-started, no_refit=True,
        parameter_metadata=meta, source_sha256=source.reference.sha256,
        barriers=barriers, first_uniform_folds=folds,
        stationary_states=len(all_points), stress_states=len(all_stress), grid_points=len(grid),
        max_unload_displacement_L0=max(r['du_norm_L0'] for r in all_stress if r['static_unload']),
        max_stress_force_residual_eV_L0=max(r['equilibrium_residual_eV_L0'] for r in all_stress),
        full_registry_validated_numerically=True, material_calibration_accepted=False,
        production_registered=False, finite_defect_activation_validated=False,
        physical_seconds_available=False, physical_Hz_available=False,
        static_unload_is_not_dynamic_hold=True,
        overall_status='full vector static reference checked; Al material/core/kinetics and fatigue gates remain closed')
    save_json(folder/'scientific_status.json', summary)
    print(json.dumps(summary, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item(), indent=2), flush=True)


if __name__ == '__main__':
    main()
