"""Recompute retained v6 states, check independent engines and draw figures."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np
from solver_v1.silicon_atomistic_reference import LammpsSilicon
from solver_v1.silicon_crack_research import SpatialSW
from solver_v1.silicon_specimen_research import specimen_stress
from .run_atomistic_controls import checksum, save_json
from .run_static_probe import source_parameters


def isolated_bond_safe_ase(potential):
    """ASE 3.26 raises at zeta=0, n<1 before multiplying by grad(zeta)=0.

    For a lone bond there are no active third-neighbor terms. The TOTAL bond
    order force is zero, although db/dzeta considered alone is singular. Pure
    Si here has n=.78734>1/2, so the smooth radial cutoff has this zero-force
    limit too. This local subclass changes that exact-zero case only. It does
    not modify the installed ASE source. Dimer force differences are tested.
    """
    from ase.calculators.tersoff import Tersoff
    class IsolatedBondTersoff(Tersoff):
        zero_zeta_calls = 0
        def _calc_bij_d(self, zeta, beta, n):
            if zeta == 0.:
                if not n > .5:
                    raise ValueError('zero-environment limit not audited for this exponent')
                self.zero_zeta_calls += 1
                return 0.
            return super()._calc_bij_d(zeta, beta, n)
    return IsolatedBondTersoff.from_lammps(potential)


def plot(root, baseline, refined):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for (model, info), color, label in zip(baseline['models'].items(), ['#166d98', '#bc5932'], ['Original SW', 'Tersoff 1989']):
        rr = info['cleavage']['rows']
        axes[0, 0].plot([r['opening_A'] for r in rr], [r['work_J_m2'] for r in rr], color=color, label=label, lw=2)
        sizes = np.geomspace(.1, 100., 200)
        axes[0, 1].loglog(sizes, specimen_stress(info['griffith_K_MPa_sqrt_m'], sizes*1e-6, geometry_factor=1.), color=color, lw=2)
        selected = [r for r in refined['cases'] if r['model'] == model and r['front_repeats'] == 4]
        for radius, marker, style in [(28., 'o', ':'), (40., 's', '--'), (56., '^', '-')]:
            rr = [r for r in selected if r['radius_A'] == radius]
            gaps = [np.mean(min((c for c in r.get('front_after_polish', r['front'])['columns'] if c['x_A'] > 0),
                               key=lambda c:c['x_A'])['gaps_A']) for r in rr]
            axes[1, 0].plot([r['factor_of_rigid_griffith'] for r in rr], gaps, color=color, ls=style,
                            marker=marker, label=f'{label}, R={radius:g} Å')
        for front, marker in [(4, 'o'), (8, 's')]:
            rr = [r for r in refined['cases'] if r['model'] == model and r['radius_A'] == 28
                  and r['front_repeats'] == front and 'stationary_check' in r]
            axes[1, 1].plot([r['factor_of_rigid_griffith'] for r in rr],
                            [r['stationary_check']['spectra'][-1]['eigenvalues_eV_A2'][0] for r in rr],
                            marker=marker, color=color, ls='-' if front == 4 else '--', label=f'{label}, front={front}')
    axes[0, 0].set(xlabel='Rigid interface opening (Å)', ylabel='Work for two surfaces (J/m²)',
                   title='A. Same atomistic energies, no strength scaling', xlim=(0, 4.))
    axes[0, 0].legend(frameon=False)
    axes[0, 1].set(xlabel='Half-crack length c (µm)', ylabel='Nominal Griffith stress (MPa)', title='B. Infinite central crack: Y=1')
    axes[0, 1].text(.04, .05, 'Energy equality only; not measured failure stress', transform=axes[0, 1].transAxes, fontsize=9)
    axes[1, 0].set(xlabel='Applied K / rigid-separation Griffith K', ylabel='First ahead bond normal gap (Å)',
                   title='C. Independent atoms / optical boundary')
    axes[1, 0].legend(frameon=False, fontsize=8, ncol=2)
    axes[1, 1].axhline(0., color='#566', lw=.8)
    axes[1, 1].set(xlabel='Applied K / rigid-separation Griffith K', ylabel='Lowest free-atom curvature (eV/Å²)',
                   title='D. All free Cartesian directions tested')
    axes[1, 1].legend(frameon=False, fontsize=8)
    for ax in axes.ravel():
        ax.grid(alpha=.18)
    fig.suptitle('Si specimen loading bridge — static references, not material calibration', fontsize=14)
    fig.savefig(root/'specimen_bridge.png', dpi=180)
    fig.savefig(root/'specimen_bridge.pdf')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('results/silicon_specimen_v6'))
    args = parser.parse_args()
    root = args.root
    if (root/'validation.json').exists():
        raise FileExistsError('refusing to overwrite completed validation')
    started = time.perf_counter()
    base_dir, refined_dir = root/'verified_affine_boundary', root/'verified_internal_boundary'
    initial = json.loads((base_dir/'specimen_bridge.json').read_text())
    refined = json.loads((refined_dir/'refinement.json').read_text())
    for data in [initial, refined]:
        assert data.get('complete'), 'calculation incomplete'
        for path, expected in data['code_sha256'].items():
            assert checksum(path) == expected, 'executed source changed: '+path
    files = [(base_dir/r['artifact'], r) for r in initial['cases']]
    files += [(refined_dir/r['artifact'], r) for r in refined['cases']]
    assert len(set(p for p, _ in files)) == len(files), 'duplicate state paths'
    records, independent, dimers = [], [], []
    parameters, _ = source_parameters()
    for label, info in initial['models'].items():
        md = info['engine']
        potential = (Path('results/silicon_wafer_feasibility/source_Si.sw') if md['style'] == 'sw'
                     else Path('results/silicon_atomistic_v5/sources/SiC.tersoff'))
        with LammpsSilicon(md['style'], potential, sha256=md['potential_sha256']) as engine:
            if md['style'] == 'tersoff':
                from ase import Atoms
                for distance in [2.35, 2.8, 3.1]:
                    xyz = np.array([[10., 10., 10.], [10.+distance, 10., 10.]])
                    atom = Atoms('Si2', positions=xyz, cell=np.eye(3)*25, pbc=False)
                    atom.calc = isolated_bond_safe_ase(potential)
                    value = engine.evaluate(xyz, np.eye(3)*25, periodic=(False, False, False))
                    energy_error = abs(atom.get_potential_energy()-value.energy)
                    force_error = float(np.max(abs(atom.get_forces()+value.gradient)))
                    h = 1e-5
                    xp, xm = xyz.copy(), xyz.copy()
                    xp[1, 0] += h
                    xm[1, 0] -= h
                    atom.set_positions(xp)
                    ep = atom.get_potential_energy()
                    atom.set_positions(xm)
                    em = atom.get_potential_energy()
                    derivative_error = abs((ep-em)/(2*h)-value.gradient[1, 0])
                    assert max(energy_error, force_error) < 1e-9 and derivative_error < 1e-6
                    dimers.append(dict(distance_A=distance, energy_error_eV=energy_error,
                                       force_error_eV_A=force_error, derivative_error_eV_A=derivative_error))
            for path, row in files:
                if row['model'] != label:
                    continue
                checked_path = refined_dir/row['stationary_check']['output'] if 'stationary_check' in row else path
                if path.parent == base_dir and not row['relaxation']['force_converged']:
                    corrected = next(r for r in refined['original_failed_force_corrections'] if r['input'] == path.name)
                    checked_path = refined_dir/corrected['output']
                with np.load(checked_path) as data:
                    saved = {k:data[k].copy() for k in data.files}
                value = engine.evaluate(saved['positions']-saved['origin'], saved['cell'], periodic=(False, False, True))
                de = abs(float(np.sum(value.site_energy-saved['site_energy'])))
                df = float(np.max(abs(value.gradient-saved['gradient'])))
                fixed = float(np.max(abs(saved['positions'][saved['fixed']]-saved['prescribed'][saved['fixed']])))
                residual = float(np.max(abs(value.gradient[~saved['fixed']])))
                assert de < 1e-7 and df < 1e-8 and fixed < 1e-12 and residual < 2e-6, 'saved state/force mismatch'
                records.append(dict(file=checked_path.relative_to(root).as_posix(), atoms=len(saved['positions']),
                                    energy_error_eV=de, force_error_eV_A=df, fixed_error_A=fixed,
                                    actual_free_residual_eV_A=residual))
                if path.parent == refined_dir and row['radius_A'] == 28 and row['front_repeats'] == 4 and row['factor_of_rigid_griffith'] >= 1.:
                    if md['style'] == 'sw':
                        other = SpatialSW(parameters, front_period=row['front_period_A']).evaluate(saved['positions'], representation='direct')
                        de = abs(other.energy-value.energy)
                        df = float(np.max(abs(other.gradient-value.gradient)))
                        name = 'independent explicit SW triple sum'
                    else:
                        from ase import Atoms
                        atoms = Atoms('Si'*len(saved['positions']), positions=saved['positions']-saved['origin'], cell=saved['cell'], pbc=(False, False, True))
                        atoms.calc = isolated_bond_safe_ase(potential)
                        de = abs(float(atoms.get_potential_energy())-value.energy)
                        df = float(np.max(abs(atoms.get_forces()+value.gradient)))
                        name = 'ASE Tersoff with audited zero-environment chain derivative'
                    assert de < 2e-7 and df < 1e-8, 'independent implementation mismatch'
                    independent.append(dict(model=label, file=checked_path.name, implementation=name, energy_error_eV=float(de), force_error_eV_A=df))
    stable = [r['stationary_check'] for r in refined['cases'] if 'spectra' in r.get('stationary_check', {})]
    assert len(stable) == 8 and all(s['full_q_stability_certified'] for s in stable), 'stability controls incomplete or failed'
    assert all(r['force_converged'] for r in refined['original_failed_force_corrections']), 'force failures unresolved'
    for work in refined['boundary_work']:
        assert work['relative_error'] < 2e-5 and all(r['force_converged'] for r in work['force_checks']), 'boundary work failed'
    plot(root, initial, refined)
    save_json(root/'validation.json', dict(passed=True, actual_recomputed_states=records,
        independent_engine_comparisons=independent, isolated_bond_checks=dimers, checked_distinct_state_paths=len(files),
        full_free_atom_spectra=len(stable), work_checks=refined['boundary_work'],
        final_pipeline_cases=len(initial['cases'])+len(refined['cases']), initial_failed_force_count=initial['failed_force_cases'],
        source_hashes_verified=True, elapsed_seconds=time.perf_counter()-started,
        material_calibrated=False, kinetic_calibrated=False, physical_clock_available=False))
    print(json.dumps(dict(event='validation_complete', states=len(records), independent=len(independent), elapsed_seconds=time.perf_counter()-started)), flush=True)


if __name__ == '__main__':
    main()
