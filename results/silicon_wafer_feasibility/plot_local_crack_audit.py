"""Plot computed atomistic coordinates and energies, not illustrative risk fields."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    root = Path('results/silicon_local_crack_v2')
    summary = read(root/'summary.json')
    base = root/'front4'
    geometry = np.load(base/'geometry.npz')
    states = np.load(base/'localized/stationary_positions.npz')
    bonds = geometry['crossing_bonds']
    middle = geometry['reference'][bonds].mean(axis=1)
    initial_gaps = states['initial'][bonds[:, 1], 1]-states['initial'][bonds[:, 0], 1]
    opened_gaps = states['opened_minimum'][bonds[:, 1], 1]-states['opened_minimum'][bonds[:, 0], 1]
    blue, orange = '#176caa', '#d16b26'
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                         'savefig.facecolor':'white'})
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6), layout='constrained')
    ax = axes[0, 0]
    scatter = ax.scatter(middle[:, 0], middle[:, 2], c=opened_gaps-initial_gaps,
                         cmap='viridis', s=95, edgecolors='#333333', linewidths=.35)
    fig.colorbar(scatter, ax=ax, label='Change in normal bond gap (Å)', shrink=.85)
    ax.set(xlabel='Propagation coordinate x (Å)', ylabel='Front coordinate z (Å)',
           title='A. One newly opened bond at an existing crack front')
    ax.set_aspect('equal', adjustable='box')
    ax = axes[0, 1]
    for label, title, color in [('localized', 'One bond, independent atoms', blue),
                               ('coherent', 'Four bonds moving together', orange)]:
        path = read(base/label/'summary.json')
        rows = read(base/label/'profile.json')
        endpoint = path['states']['opened_minimum']
        points = [(row['q_A'], row['relative_energy_eV']) for row in rows if row['q_A'] <= endpoint['q_A']]
        points += [(path['states'][key]['q_A'], path['states'][key]['relative_energy_eV'])
                   for key in ('saddle', 'opened_minimum')]
        points.sort()
        ax.plot(*np.array(points).T, color=color, label=title, linewidth=2)
        saddle = path['states']['saddle']
        ax.scatter(saddle['q_A'], saddle['relative_energy_eV'], marker='*', color=color, s=110)
    ax.set(xlabel='Selected bond normal separation a (Å)', ylabel='Total energy above initial (eV)',
           title='B. Same potential and fixed grips; different paths')
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=.18)
    ax = axes[1, 0]
    lengths = [case['geometry']['front_period_A'] for case in summary['cases']]
    for label, title, color in [('localized', 'Localized first step', blue),
                               ('coherent', 'Whole-front control', orange)]:
        barriers = [case['paths'][label]['forward_barrier_eV'] for case in summary['cases']]
        ax.plot(lengths, barriers, 'o-', color=color, label=title, linewidth=2)
        for x, y in zip(lengths, barriers):
            ax.annotate(f'{y:.3f}', (x, y), xytext=(0, 9), textcoords='offset points', ha='center', color=color)
    ax.set(xlabel='Periodic front length (Å)', ylabel='Stationary energy barrier (eV)',
           title='C. Front-length comparison: local and collective paths', ylim=(0, 7.5))
    ax.legend(frameon=False, fontsize=9, loc='center left')
    ax.grid(alpha=.18)
    ax = axes[1, 1]
    sequence = read(root/'sequence/summary.json')
    previous_energy = 0.
    for step in sequence['steps']:
        index = step['step']
        folder = base/'localized' if index == 1 else root/f'sequence/step{index}'
        path = read(folder/'summary.json')
        rows = read(folder/'profile.json')
        qstart = rows[0]['q_A']
        qend = path['states']['opened_minimum']['q_A']
        points = [(row['q_A'], row['relative_energy_eV']) for row in rows if row['q_A'] <= qend]
        points += [(path['states'][key]['q_A'], path['states'][key]['relative_energy_eV'])
                   for key in ('saddle', 'opened_minimum')]
        points.sort()
        q, e = np.array(points).T
        ax.plot(index-1+(q-qstart)/(qend-qstart), e+previous_energy, color=blue, linewidth=2)
        ax.scatter(index, step['endpoint_energy_from_initial_eV'], color=blue, zorder=4)
        previous_energy = step['endpoint_energy_from_initial_eV']
    coherent = summary['cases'][0]['paths']['coherent']['forward_barrier_eV']
    maximum = sequence['maximum_static_path_energy_above_initial_eV']
    ax.axhline(coherent, color=orange, linestyle='--', label=f'Together: {coherent:.3f} eV')
    ax.axhline(maximum, color=blue, linestyle=':', label=f'Sequential path: {maximum:.3f} eV')
    ax.set(xlabel='Successive local opening coordinates (one segment per bond)',
           ylabel='Total energy above the same initial state (eV)',
           title='D. Verified local steps reach the same final state', xticks=range(5), ylim=(-.08, 2.6))
    ax.legend(frameon=False, fontsize=9, loc='upper left')
    ax.grid(alpha=.18)
    fig.suptitle('Si crack-front localization — original SW, 0 K, finite fixed grips\n'
                 'Computed static energies; no wafer-strength, fatigue-life or physical-time calibration', fontsize=14)
    fig.savefig(root/'local_crack_audit.png', dpi=170)
    fig.savefig(root/'local_crack_audit.svg')
    plt.close(fig)


if __name__ == '__main__':
    main()
