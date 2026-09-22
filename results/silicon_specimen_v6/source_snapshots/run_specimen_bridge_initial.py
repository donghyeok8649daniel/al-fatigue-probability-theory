"""Run new Si cleavage and anisotropic K-boundary statics; no material fit.

The source-bound SW/Tersoff engines and each model's v5 relaxed elastic
constants remain distinct. All dimensions and boundary work are recorded.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import time
from types import SimpleNamespace
import zipfile

import numpy as np

from solver_v1.silicon_atomistic_reference import LammpsSilicon
from solver_v1.silicon_crack_research import RelaxedCoordinates, relax_atoms
from solver_v1.silicon_specimen_research import (
    AnisotropicModeI, atomic_crack_boundary, specimen_stress,
)
from .run_atomistic_controls import checksum, save_json, strip_box
from .run_static_probe import diamond111_basis, EV_A2_TO_J_M2, EV_A3_TO_GPA


V5 = Path('results/silicon_atomistic_v5')
SOURCE_FILES = ['solver_v1/silicon_specimen_research.py',
                'results/silicon_wafer_feasibility/run_specimen_bridge.py',
                'solver_v1/silicon_atomistic_reference.py',
                'solver_v1/silicon_crack_research.py']


def source_geometry_audit(archive):
    """Read metadata only; no interpolation of differently strained DFT cells."""
    expected = '1d3efe976c53bcd2889c0556f2d600c26bf4603522b5c92fd8633dff7b7ab5c2'
    if checksum(archive) != expected:
        raise ValueError('unexpected Cambridge archive checksum')
    with zipfile.ZipFile(archive) as z:
        raw = z.read('gp_iter6_sparse9k.xml.xyz')
    lines = iter(raw.decode().splitlines())
    groups = {}
    frame = 0
    for line in lines:
        count = int(line)
        header = next(lines)
        if re.search(r'config_type=["\']?decohesion', header):
            cell = np.fromstring(re.search(r'Lattice="([^"]+)"', header)[1], sep=' ').reshape(3, 3)
            area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
            groups.setdefault(count, []).append(dict(frame=frame, cell_rows_A=cell.tolist(), area_A2=area,
                normal_period_A=float(abs(np.linalg.det(cell))/area), xc_declared='xc_functional=' in header))
        for _ in range(count):
            next(lines)
        frame += 1
    return dict(archive_sha256=expected, total_frames=frame,
        member_sha256=hashlib.sha256(raw).hexdigest(),
        groups=[dict(atoms=count, rows=rows,
                     relative_area_span=(max(r['area_A2'] for r in rows)-min(r['area_A2'] for r in rows))/np.mean([r['area_A2'] for r in rows]))
                for count, rows in groups.items()],
        used_as_cohesive_curve=False,
        reason='in-plane cell changes along all three series; XC absent. Training configurations are not a fixed-area rigid opening curve.')


def cleavage_geometry(lattice, periods):
    basis, cell = diamond111_basis()
    cell = lattice*cell
    r = (basis[None]+np.arange(periods)[:, None, None]*np.array([0., 0., 1.]))@cell
    r = r.reshape(-1, 3)
    cell[2] *= periods
    normal = np.ones(3)/np.sqrt(3)
    height = r@normal
    planes = np.unique(np.round(height, 8))
    gap = np.diff(planes)
    cuts = np.flatnonzero(gap > .5*(min(gap)+max(gap)))
    cut = cuts[np.argmin(abs(.5*(planes[cuts]+planes[cuts+1])-np.mean(height)))]
    upper = height > .5*(planes[cut]+planes[cut+1])
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    return r, cell, normal, upper, area


def cleavage(engine, lattice):
    """One periodic cut: delta U/area is work for TWO newly created surfaces.

    Traction differentiates both positions and cell. The cell term cannot be
    replaced by the net force on a periodic half-slab.
    """
    rows, depth_checks = [], []
    for periods in [3, 6, 9]:
        r, cell, n, upper, area = cleavage_geometry(lattice, periods)
        reference = engine.evaluate(r, cell)

        def evaluate(opening):
            rr, cc = r.copy(), cell.copy()
            rr[upper] += opening*n
            cc[2] += opening*n
            value = engine.evaluate(rr, cc)
            length = float(cc[2]@n)
            gradient_n = value.gradient@n
            reaction = (abs(np.linalg.det(cc))*(n@value.stress@n)/length
                        + np.sum(gradient_n*(upper.astype(float)-(rr@n)/length)))
            return dict(opening_A=float(opening),
                        work_J_m2=float(np.sum(value.site_energy-reference.site_energy)/area*EV_A2_TO_J_M2),
                        traction_GPa=float(reaction/area*EV_A3_TO_GPA))

        selected = [evaluate(delta) for delta in [0., .5, 1., 2., 4., 6.]]
        depth_checks.append(dict(periods=periods, atoms=len(r), area_A2=area, rows=selected))
        if periods == 6:
            rows = [evaluate(d) for d in np.linspace(0., 6., 121)]
            derivative_checks = []
            for delta in [.25, .6, 1.1, 1.7, 2.7]:
                exact = evaluate(delta)['traction_GPa']
                for step in [1e-4, 5e-5]:
                    numerical = (evaluate(delta+step)['work_J_m2']-evaluate(delta-step)['work_J_m2'])/(2*step)*10.
                    derivative_checks.append(dict(opening_A=delta, step_A=step,
                                                  analytic_GPa=exact, finite_difference_GPa=numerical,
                                                  abs_error_GPa=abs(exact-numerical)))
    mismatch = max(abs(a['work_J_m2']-b['work_J_m2']) for a, b in zip(depth_checks[0]['rows'], depth_checks[-1]['rows']))
    if mismatch > 1e-8 or max(r['abs_error_GPa'] for r in derivative_checks) > 2e-4:
        raise AssertionError('cleavage derivative/depth check failed')
    peak = max(rows, key=lambda r:r['traction_GPa'])
    return dict(convention='rigid 0 K (111) shuffle, one periodic cut; two created surfaces; no reconstruction',
        rows=rows, depth_checks=depth_checks, depth_work_max_error_J_m2=mismatch,
        derivative_checks=derivative_checks,
        separated_work_J_m2=rows[-1]['work_J_m2'], plateau_delta_J_m2=rows[-1]['work_J_m2']-rows[-21]['work_J_m2'],
        sampled_peak_traction_GPa=peak['traction_GPa'], peak_opening_A=peak['opening_A'],
        is_measured_toughness=False)


def observables(boundary, positions, cutoff):
    bonds = boundary.bonds
    middle = boundary.reference[bonds].mean(axis=1)
    delta = positions[bonds[:, 1]]-positions[bonds[:, 0]]
    delta[:, 2] -= boundary.period_A*np.round(delta[:, 2]/boundary.period_A)
    valid = ~np.any(boundary.fixed[bonds], axis=1)
    x = np.unique(np.round(middle[valid, 0], 7))
    columns = []
    for xx in x:
        select = valid & np.isclose(middle[:, 0], xx, atol=2e-7, rtol=0)
        lengths = np.linalg.norm(delta[select], axis=1)
        columns.append(dict(x_A=float(xx), gaps_A=delta[select, 1].tolist(), lengths_A=lengths.tolist(),
                            original_pairs_beyond_cutoff=int(np.sum(lengths > cutoff))))
    opened = [c['x_A'] for c in columns if c['original_pairs_beyond_cutoff'] == len(c['lengths_A'])]
    return dict(columns=columns, furthest_fully_open_original_column_A=max(opened) if opened else None,
        cutoff_A=cutoff,
        meaning='reversible geometry of original crossing pairs, not damage, an irreversible crack tip, or a rupture rate')


def load_case(engine, elastic, boundary, K, initial, output, *, tolerance=2e-6):
    prescribed = boundary.displaced(elastic, K)
    origin, cell = strip_box(prescribed, boundary.period_A)
    # Include ample vacuum for relaxed atoms while preserving front periodicity.
    origin[:2] -= 12.
    cell[0, 0] += 24.
    cell[1, 1] += 24.
    model = SimpleNamespace(evaluate=lambda r:engine.evaluate(r-origin, cell, periodic=(False, False, True)))
    coordinates = RelaxedCoordinates(prescribed, boundary.fixed)
    if initial is None:
        initial = prescribed.copy()
    else:
        initial = initial.copy()
        initial[boundary.fixed] = prescribed[boundary.fixed]
    positions, info = relax_atoms(model, coordinates, initial=initial, tolerance=tolerance, maxiter=3500)
    value = model.evaluate(positions)
    unit = boundary.displaced(elastic, 1.)-boundary.reference
    # U includes all atoms, including fixed-atom site energies. Envelope
    # derivative: only the prescribed coordinates move after free equilibration.
    work = float(np.sum(value.gradient[boundary.fixed]*unit[boundary.fixed]))
    row = dict(K_MPa_sqrt_m=float(K), nominal_far_field_G_J_m2=float(elastic.energy_release_J_m2(K)),
               radius_A=boundary.radius_A, grip_width_A=boundary.grip_width_A,
               atoms=len(positions), free_atoms=int(np.sum(~boundary.fixed)),
               front_period_A=boundary.period_A, relaxation=info,
               net_force_eV_A=np.sum(-value.gradient, axis=0).tolist(),
               grip_force_eV_A=np.sum(-value.gradient[boundary.fixed], axis=0).tolist(),
               fixed_position_error_A=float(np.max(abs(positions[boundary.fixed]-prescribed[boundary.fixed]))),
               dU_dK_eV_per_MPa_sqrt_m=work,
               material_calibrated=False, stability_certified=False,
               continuum_boundary='leading anisotropic mode I; finite radius; no flexible-tip correction')
    np.savez_compressed(output.with_suffix('.npz'), positions=positions, prescribed=prescribed,
                         reference=boundary.reference, fixed=boundary.fixed, bonds=boundary.bonds,
                         site_energy=value.site_energy, gradient=value.gradient, cell=cell, origin=origin)
    return positions, row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/silicon_specimen_v6'))
    parser.add_argument('--radii', nargs='+', type=float, default=[28., 40., 56.])
    parser.add_argument('--factors', nargs='+', type=float, default=[.8, 1., 1.2])
    parser.add_argument('--front-repeats', type=int, default=4)
    parser.add_argument('--skip-source-audit', action='store_true')
    args = parser.parse_args()
    target = args.output/'specimen_bridge.json'
    if target.exists() or (args.output/'checkpoint.json').exists():
        raise FileExistsError('refusing to overwrite an existing run')
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    baseline = json.loads((V5/'atomistic_controls.json').read_text())
    result = dict(scope='new static Si reference calculations with specimen loading bridge; no material/kinetic calibration',
        models={}, cases=[], input_v5_sha256=checksum(V5/'atomistic_controls.json'),
        code_sha256={p:checksum(p) for p in SOURCE_FILES}, front_repeats=args.front_repeats)
    if not args.skip_source_audit:
        result['dft_training_geometry_audit'] = source_geometry_audit(Path('.cache/si-atomistic-v5/Si_PRX_GAP.zip'))
    potentials = dict(original_sw_1985=Path('results/silicon_wafer_feasibility/source_Si.sw'),
                      tersoff_1989=V5/'sources/SiC.tersoff')
    for label, old in baseline['models'].items():
        bulk = old['bulk']
        stiffness = [r for r in bulk['rows'] if r.get('relaxed') and 'C11_GPa' in r][-1]
        elastic = AnisotropicModeI(*[stiffness[k] for k in ['C11_GPa', 'C12_GPa', 'C44_GPa']])
        md = old['metadata']
        with LammpsSilicon(md['style'], potentials[label], sha256=md['potential_sha256']) as engine:
            separated = cleavage(engine, bulk['lattice_A'])
            kg = float(elastic.griffith_K(separated['separated_work_J_m2']))
            info = dict(engine=engine.metadata(), bulk_input=stiffness, lattice_A=bulk['lattice_A'],
                        H_GPa_inv=elastic.H_GPa_inv, effective_plane_strain_modulus_GPa=1/elastic.H_GPa_inv,
                        roots=[[p.real, p.imag] for p in elastic.roots], cleavage=separated,
                        griffith_K_MPa_sqrt_m=kg,
                        griffith_interpretation='energetic equality using rigid separation; not K+ or experimental K_Ic',
                        contour_checks=[dict(radius_A=r, points=n, J_J_m2=elastic.contour_J(kg, radius_A=r, points=n))
                                        for r in [8., 32., 128.] for n in [512, 1024]],
                        specimen_examples=[dict(half_crack_um=a, geometry='infinite central through crack', Y=1.,
                                                griffith_nominal_stress_MPa=float(specimen_stress(kg, a*1e-6, geometry_factor=1.)))
                                           for a in [.1, 1., 10., 100.]])
            result['models'][label] = info
            print(json.dumps(dict(event='cleavage', model=label, G_J_m2=separated['separated_work_J_m2'],
                                  KG=kg, H=elastic.H_GPa_inv)), flush=True)
            cutoff = 3.77118 if md['style'] == 'sw' else 3.2
            for radius in args.radii:
                boundary = atomic_crack_boundary(bulk['lattice_A'], radius_A=radius,
                                                 front_repeats=args.front_repeats, grip_width_A=8.)
                previous, previous_k = None, 0.
                for factor in args.factors:
                    k = factor*kg
                    if previous is not None:
                        previous = previous+(boundary.displaced(elastic, k)-boundary.displaced(elastic, previous_k))
                    name = f'{label}_R{radius:g}_K{factor:g}'
                    positions, row = load_case(engine, elastic, boundary, k, previous, args.output/name)
                    row.update(model=label, factor_of_rigid_griffith=factor, artifact=name+'.npz',
                               front=observables(boundary, positions, cutoff))
                    result['cases'].append(row)
                    save_json(args.output/(name+'.json'), row)
                    save_json(args.output/'checkpoint.json', result)
                    print(json.dumps(dict(event='relaxed_boundary', model=label, radius_A=radius, factor=factor,
                                          force_converged=row['relaxation']['force_converged'],
                                          force_max=row['relaxation']['free_gradient_max_eV_A'],
                                          iterations=row['relaxation']['iterations'],
                                          open_column=row['front']['furthest_fully_open_original_column_A'])), flush=True)
                    # Preserve failures, but do not propagate a failed force solution.
                    previous = positions if row['relaxation']['force_converged'] else None
                    previous_k = k
    result.update(complete=True, elapsed_seconds=time.perf_counter()-started,
                  material_calibrated=False, kinetic_calibrated=False, physical_clock_available=False,
                  equilibrium_cases=sum(c['relaxation']['force_converged'] for c in result['cases']),
                  failed_force_cases=sum(not c['relaxation']['force_converged'] for c in result['cases']))
    save_json(target, result)
    print(json.dumps(dict(event='complete', elapsed_seconds=result['elapsed_seconds'],
                          cases=len(result['cases']), failed=result['failed_force_cases'])), flush=True)


if __name__ == '__main__':
    main()
