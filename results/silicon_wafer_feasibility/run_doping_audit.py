"""Reproduce v7 doping/elasticity diagnostics, with no external atomistic run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from solver_v1.silicon_doping_research import (
    directional_young_GPa, dopant_site_statistics, excess_electrons_for_material_volume,
    load_elastic_samples, silicon_site_density_cm3,
)
from solver_v1.silicon_specimen_research import AnisotropicModeI, specimen_K


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'results/silicon_doping_v7/sources'


def dump_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(name):
    with (SOURCE/name).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def run(output):
    output = Path(output)
    products = ['elastic_crack_fields.csv', 'chemical_site_counts.csv', 'audit.json', 'doping_audit.png']
    if any((output/name).exists() for name in products):
        raise FileExistsError('preserve existing audit; choose a new output directory')
    output.mkdir(parents=True, exist_ok=True)
    samples = load_elastic_samples()
    elastic, occupancy, j_errors = [], [], []
    nominal_K = float(specimen_K(100, 10e-6, geometry_factor=1.))
    reference = samples['B0.6']
    uncertainty = np.array([.3, .1, .2])
    for sample in samples.values():
        for temperature in (-40., 25., 85.):
            c = sample.constants_GPa(temperature)
            field = sample.crack_field(temperature)
            ref = reference.crack_field(temperature)
            j_exact = float(field.energy_release_J_m2(nominal_K))
            j4096 = field.contour_J(nominal_K, radius_A=40., points=4096)
            j8192 = field.contour_J(nominal_K, radius_A=80., points=8192)
            j_errors.extend([abs(j4096/j_exact-1), abs(j8192/j_exact-1)])
            # Fixed Wsep hypothetical ratio. No Wsep or K_Ic is assigned.
            same_work_ratio = float(np.sqrt(ref.H_GPa_inv/field.H_GPa_inv))
            corners = [AnisotropicModeI(*(c+np.array(sign)*uncertainty)).H_GPa_inv
                       for sign in itertools.product([-1, 1], repeat=3)] if temperature == 25 else []
            elastic.append(dict(sample_id=sample.sample_id, species=sample.species,
                                carrier_average_cm3=sample.carrier_average_cm3, temperature_C=temperature,
                                c11_GPa=c[0], c12_GPa=c[1], c44_GPa=c[2],
                                E100_GPa=directional_young_GPa(*c, [1, 0, 0]),
                                E110_GPa=directional_young_GPa(*c, [1, 1, 0]),
                                E111_GPa=directional_young_GPa(*c, [1, 1, 1]),
                                H_GPa_inv=field.H_GPa_inv, J_at_100MPa_halfcrack10um_J_m2=j_exact,
                                J_contour_4096=j4096, J_contour_8192=j8192,
                                constant_Wsep_K_ratio_to_B06=same_work_ratio,
                                H_c0_corner_min=min(corners) if corners else '',
                                H_c0_corner_max=max(corners) if corners else ''))
    # Existing atom count is factual geometry, not a new doped calculation.
    v6_paths = [ROOT/f'results/silicon_specimen_v6/verified_affine_boundary/original_sw_1985_R{r}_K1.json'
                for r in (28, 56)]
    counts = [json.loads(path.read_text())['atoms'] for path in v6_paths]
    lattice = 5.431  # nominal experimental lattice, Noda supplement Table S1
    for count in counts:
        for density in np.logspace(14, 21, 29):
            s = dopant_site_statistics(density, lattice_A=lattice, site_count=count)
            occupancy.append(dict(site_count=count, chemical_density_cm3=density, **s))
    published = read_csv('noda_2023_ideal_strength.csv')
    wafer = read_csv('chen_2025_wafer.csv')
    strength_changes = {kind: 100*(float(next(r for r in published if r['carrier_type'] == kind and float(r['excess_carrier_cm3']) == 5e21)['ideal_strength_LDA_GPa'])/20.98-1)
                        for kind in ('electron', 'hole')}
    room = [r for r in elastic if r['temperature_C'] == 25]
    summary = dict(
        status='research reference; doping-sensitive elasticity connected; doped fracture/kinetics unavailable',
        calculations=dict(elastic_fields=len(elastic), independent_J_integrals=len(j_errors),
                          chemical_occupancy_scenarios=len(occupancy), c0_sensitivity_corners=7*8,
                          new_DFT=0, new_MD=0, new_atomistic_energy_evaluations=0),
        geometry=dict(plane='111', propagation='11-2', front='1-10', kinematics='plane strain',
                      specimen='infinite plate central through crack', half_crack_m=10e-6,
                      remote_stress_MPa=100, Y=1, K_MPa_sqrt_m=nominal_K),
        validation=dict(max_relative_contour_J_error=max(j_errors),
                        doped_fracture_calibrated=False, physical_clock_calibrated=False),
        room_temperature=dict(H_min_GPa_inv=min(r['H_GPa_inv'] for r in room),
                              H_max_GPa_inv=max(r['H_GPa_inv'] for r in room),
                              same_Wsep_ratio_min=min(r['constant_Wsep_K_ratio_to_B06'] for r in room),
                              same_Wsep_ratio_max=max(r['constant_Wsep_K_ratio_to_B06'] for r in room),
                              reference='B0.6, doped; not intrinsic silicon',
                              assumption='Wsep held equal only for sensitivity; no doped strength assigned'),
        atom_count_scale=dict(lattice_A=lattice, silicon_sites_cm3=silicon_site_density_cm3(lattice),
                              v6_site_counts=counts,
                              at_1e15_7600=dopant_site_statistics(1e15, lattice_A=lattice, site_count=7600),
                              at_1e19_7600=dopant_site_statistics(1e19, lattice_A=lattice, site_count=7600),
                              eight_atom_cell_excess_electrons_at_1e19=excess_electrons_for_material_volume(1e19, material_volume_A3=lattice**3)),
        published_ideal_strength_change_at_5e21_percent=strength_changes,
        wafer_characteristic_strength_range_GPa=[min(float(r['characteristic_strength_GPa']) for r in wafer), max(float(r['characteristic_strength_GPa']) for r in wafer)],
        limitations=['No chemical concentration to carrier conversion', 'No fitted activation fraction',
                     'No DFT ideal-strength extrapolation or use as crack cutoff',
                     'No doped SW/Tersoff atomistic material', 'Published wafer fit is not a PDE input',
                     'Electronic fast-equilibrium closure and mobility require physical validation'],
        source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(SOURCE.glob('*')) if p.is_file()} |
                      {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in v6_paths})
    if max(j_errors) > 1e-9:
        raise RuntimeError('independent J-integral check failed')
    dump_csv(output/'elastic_crack_fields.csv', elastic)
    dump_csv(output/'chemical_site_counts.csv', occupancy)
    (output/'audit.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    plot(output, room, occupancy, published, wafer)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def plot(output, room, occupancy, published, wafer):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    colors = {'B': '#ab3d6e', 'P': '#2867ad', 'As': '#15806d', 'Sb': '#79613b'}
    ax = axes[0, 0]
    for species in ('B', 'As', 'P'):
        rows = sorted([r for r in room if r['species'] == species], key=lambda r: r['carrier_average_cm3'])
        ax.plot([r['carrier_average_cm3']/1e19 for r in rows], [100*(r['constant_Wsep_K_ratio_to_B06']-1) for r in rows], 'o--', color=colors[species], label=species)
    ax.axhline(0, color='.6', lw=.7)
    ax.set(xlabel='Published carrier estimate ($10^{19}$ cm$^{-3}$)', ylabel='K ratio change relative to B0.6 (%)',
           title='A  Measured elasticity only, 25 C\nSame separation work assumed; not strength')
    ax.legend(title='Dopant')
    ax = axes[0, 1]
    for kind, color in [('electron', '#2867ad'), ('hole', '#ab3d6e')]:
        rows = [r for r in published if r['carrier_type'] == kind]
        ax.plot([float(r['excess_carrier_cm3'])/1e21 for r in rows], [float(r['ideal_strength_LDA_GPa']) for r in rows], 'o-', label=kind, color=color)
    ax.set(xlabel='Imposed excess carriers ($10^{21}$ cm$^{-3}$)', ylabel='Published ideal strength (GPa)',
           title='B  Published charge-only DFT [111]\nHomogeneous tension; no chemical dopant atoms')
    ax.legend()
    ax = axes[1, 0]
    for count in sorted({r['site_count'] for r in occupancy}):
        rows = [r for r in occupancy if r['site_count'] == count]
        ax.loglog([r['chemical_density_cm3'] for r in rows], [r['expected_count'] for r in rows], label=f'{count:,} sites')
    ax.axhline(1, color='.5', linestyle=':')
    ax.set(xlabel='Chemical dopant atoms (cm$^{-3}$)', ylabel='Expected dopant count in atomistic domain',
           title='C  Finite atomic domain\nUniform substitution assumption; no rounding')
    ax.legend()
    ax = axes[1, 1]
    for i, row in enumerate(wafer):
        ax.errorbar(i, float(row['characteristic_strength_GPa']), yerr=float(row['standard_error_GPa']), fmt='o', color=colors[row['species']], capsize=4)
    ax.set_xticks(range(len(wafer)), [r['sample_id'] for r in wafer], rotation=20)
    ax.set(ylabel='Published characteristic strength (GPa)', ylim=(1.8, 2.35),
           title='D  Published wafer test, n=50 per condition\nError bars: reported SE; not raw observations')
    fig.suptitle('Silicon doping: separate carrier physics, elasticity and specimen evidence', fontsize=15)
    fig.supxlabel('A: Jaakkola et al. (2014), Tables I/III.  B: Noda et al. (2023), Table S5.  D: Chen et al. (2025), Table III.\nConnecting lines guide the eye. Published data are not new DFT or wafer measurements.', fontsize=9)
    fig.savefig(output/'doping_audit.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output)
