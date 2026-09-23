"""Aggregate completed force audits with count-weighted errors and radius sensitivity.

No model call, refit, error-based selection, oxidation-state assignment or kinetics.
The condensed/extended groups and isolated small systems remain separately visible.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re
from pathlib import Path
import numpy as np


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def main(args):
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.results/'surface_environments/environment_force_errors.csv'
    rows = read_csv(source)
    baseline = json.loads((args.results/'oxidized_surfaces_max320/summary.json').read_text(encoding='utf-8'))
    if not baseline['complete']:
        raise ValueError('complete force audit required')
    aggregated = {}
    group_sums = {}
    for row in rows:
        radius = int(row['radius_set']); group = row['group']; label = row['coordination_label']
        count = int(row['atoms']); mse = float(row['force_component_RMSE_eV_A'])**2
        mae = float(row['force_component_MAE_eV_A'])
        if label in ('H', 'O'):
            category = label
        else:
            match = re.fullmatch(r'Si_O(\d+)_Si(\d+)_H(\d+)', label)
            if match is None:
                raise ValueError(label)
            oxygen = int(match.group(1))
            category = 'Si: O=0' if oxygen == 0 else ('Si: O=1-3' if oxygen < 4 else 'Si: O>=4')
        for key, target in [((radius, group, category), aggregated), ((radius, group), group_sums)]:
            values = target.setdefault(key, dict(atoms=0, sq=0., absolute=0.))
            values['atoms'] += count
            values['sq'] += 3*count*mse
            values['absolute'] += 3*count*mae
    max_group_difference = 0.
    for radius in range(3):
        for group in baseline['groups']:
            values = group_sums[(radius, group['group'])]
            if values['atoms'] != group['atoms']:
                raise ValueError('lost atoms during aggregation')
            difference = abs(np.sqrt(values['sq']/(3*values['atoms'])) - group['force_component_RMSE_eV_A'])
            max_group_difference = max(max_group_difference, float(difference))
    if max_group_difference > 1e-12:
        raise ValueError('group force errors not preserved')
    output = []
    for (radius, group, category), values in sorted(aggregated.items()):
        total = group_sums[(radius, group)]
        output.append(dict(radius_set=radius, group=group, category=category, atoms=values['atoms'],
            atom_fraction=values['atoms']/total['atoms'],
            squared_error_fraction=values['sq']/total['sq'] if total['sq'] else 0.,
            force_component_RMSE_eV_A=float(np.sqrt(values['sq']/(3*values['atoms']))),
            force_component_MAE_eV_A=values['absolute']/(3*values['atoms'])))
    with (args.output/'chemical_classes.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0])); writer.writeheader(); writer.writerows(output)
    net = read_csv(args.results/'surface_environments/net_force_projection_diagnostic.csv')
    drift = []
    for group in baseline['groups']:
        subset = [r for r in net if r['observed_group'] == group['group']]
        raw_sq = sum(3*int(r['atoms'])*float(r['raw_force_component_RMSE_eV_A'])**2 for r in subset)
        uniform_sq = sum(3*int(r['atoms'])*float(r['uniform_error_RMS_eV_A'])**2 for r in subset)
        drift.append(dict(group=group['group'], atoms=group['atoms'],
            raw_force_component_RMSE_eV_A=group['force_component_RMSE_eV_A'],
            uniform_error_RMS_eV_A=float(np.sqrt(uniform_sq/(3*group['atoms']))),
            uniform_squared_error_fraction=uniform_sq/raw_sq if raw_sq else 0.))
    summary = dict(source_csv_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        maximum_group_RMSE_replay_difference=max_group_difference,
        oxidized_surface_chemical_classes=[r for r in output if r['group']=='Si_O_H_vacuum1'],
        group_uniform_error_projection=drift,
        category_meaning='geometric oxygen-neighbor counts; not oxidation state, electronic charge, or crack origin',
        original_force_errors_modified=False, new_potential_calls=0, new_DFT=0, new_MD=0,
        material_approved=False, initiation_probability=None)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), layout='constrained')
    big = [g for g in baseline['groups'] if not g['group'].startswith('small_')]
    labels = ['Si\nno vacuum label', 'Si/O\nno vacuum label', 'Si/H\nvacuum label', 'Si/O/H\nvacuum label']
    x = np.arange(len(big))
    bars = axes[0].bar(x, [g['force_component_RMSE_eV_A'] for g in big], color=['#617a8f','#5d9597','#ae8665','#b86149'])
    axes[0].bar_label(bars, fmt='%.3f', padding=3)
    axes[0].set_xticks(x, labels); axes[0].set_ylim(0,.40)
    axes[0].set_ylabel('Force component RMSE (eV / Angstrom)')
    axes[0].set_title(f'All {sum(g["atoms"] for g in big):,} atoms in the four extended groups')
    small = [g for g in baseline['groups'] if g['group'].startswith('small_')]
    axes[1].bar(np.arange(5), [g['force_component_RMSE_eV_A'] for g in small], color='#8495a7')
    axes[1].set_xticks(np.arange(5), ['O','Si','H','Si/H','Si/O']); axes[1].set_ylabel('Force component RMSE (eV / Angstrom)')
    axes[1].set_title(f'Small systems retained separately: {sum(g["atoms"] for g in small):,} atoms')
    fig.suptitle('MACE / public CP2K force audit: resource selection N <= 320, no error filtering')
    fig.savefig(args.output/'surface_force_groups.png', dpi=180); plt.close(fig)
    categories = ['Si: O=0','Si: O=1-3','Si: O>=4','O','H']
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), layout='constrained')
    lookup = {(r['radius_set'], r['category']):r for r in output if r['group']=='Si_O_H_vacuum1'}
    for radius, label in enumerate(['Si-O 1.8 A','Si-O 2.0 A','Si-O 2.2 A']):
        subset = [lookup.get((radius,c)) for c in categories]
        axes[0].plot(np.arange(5), [r['force_component_RMSE_eV_A'] if r else np.nan for r in subset], 'o-', label=label)
    axes[0].set_xticks(np.arange(5), categories, rotation=15); axes[0].set_ylabel('Force component RMSE (eV / Angstrom)')
    axes[0].legend(fontsize=9); axes[0].set_title('Radius sensitivity of geometric classes')
    central = [lookup[(1,c)] for c in categories]
    x = np.arange(5)
    axes[1].bar(x-.18, [100*r['atom_fraction'] for r in central], .36, label='Atom fraction', color='#547e9b')
    axes[1].bar(x+.18, [100*r['squared_error_fraction'] for r in central], .36, label='Squared-error fraction', color='#b77753')
    axes[1].set_xticks(x, categories, rotation=15); axes[1].set_ylabel('Fraction of Si/O/H group (%)'); axes[1].legend(fontsize=9)
    axes[1].set_title('Population and error contribution, central radii')
    fig.suptitle('Oxidized Si/O/H vacuum group: 423 frames, 106,148 atoms\nNeighbor counts do not identify charge, oxidation state or crack initiation.')
    fig.savefig(args.output/'surface_chemical_environments.png', dpi=180); plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True); main(parser.parse_args())
