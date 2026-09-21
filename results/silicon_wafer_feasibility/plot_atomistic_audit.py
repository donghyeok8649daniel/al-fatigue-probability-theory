"""Plot actual v5 static comparisons; no interpolated physical probability."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path('results/silicon_atomistic_v5')
    audit = json.loads((root/'dft_material_audit.json').read_text())
    controls = json.loads((root/'atomistic_controls.json').read_text())
    polish = json.loads((root/'reconstruction_polish.json').read_text())
    models = ['original_sw_1985', 'tersoff_1989']
    colors = ['#1768ac', '#d16b20']
    labels = ['Original SW', 'Tersoff 1989']
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), constrained_layout=True)
    categories = [('dia', 'PW91', 'Diamond\n489'), ('crack_111_1-10', 'PW91', 'Crack (111)\n10'),
                  ('crack_110_1-10', 'PW91', 'Crack (110)\n7'),
                  ('surface_111', 'PW91', 'Surface (111)\n47'),
                  ('surface_111_pandey', 'PBE', 'Pandey*\n50')]
    x = np.arange(len(categories))
    for j, (model, color, label) in enumerate(zip(models, colors, labels)):
        values = [next(r['component_rmse_eV_A'] for r in audit['group_results']
                       if (r['model'], r['config_type'], r['declared_xc']) == (model, kind, xc))
                  for kind, xc, _ in categories]
        axes[0].bar(x+(j-.5)*.34, values, width=.34, label=label, color=color)
    axes[0].set_xticks(x, [c[2] for c in categories], fontsize=9)
    axes[0].set_ylabel('Force component RMSE (eV / Å)')
    axes[0].set_title('A  Same structures vs archived DFT', loc='left', fontsize=11)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].text(.02, .94, '*PBE; others PW91\nCounts are source frames', transform=axes[0].transAxes,
                 va='top', fontsize=8, color='#555555')
    axes[0].set_ylim(0, 1.5)
    for i, (model, color, label) in enumerate(zip(models, colors, labels)):
        rows = controls['models'][model]['bulk']['rows']
        a = [next(r['C44_GPa'] for r in reversed(rows) if r.get('relaxed') == relaxed and 'C44_GPa' in r)
             for relaxed in [False, True]]
        bars = axes[1].bar(np.arange(2)+(i-.5)*.34, a, width=.34, color=color, label=label)
        axes[1].bar_label(bars, fmt='%.2f', fontsize=9, padding=3)
    axes[1].set_xticks([0, 1], ['Affine basis', 'Internally relaxed basis'])
    axes[1].set_ylabel('$C_{44}$ (GPa)')
    axes[1].set_ylim(0, 140)
    axes[1].set_title('B  Same potential, different relaxation', loc='left', fontsize=11)
    axes[1].text(.02, .96, '0 K, each model\'s own lattice equilibrium', transform=axes[1].transAxes,
                 va='top', fontsize=8, color='#555555')
    sw = [r['relaxed_delta_vs_relaxed_state00_eV'] for r in
          controls['models'][models[0]]['reconstruction']['rows']]
    tf = [r['energy_vs_state00_eV'] for r in polish['rows']]
    for j, (values, color, label) in enumerate(zip([sw, tf], colors, labels)):
        # Discrete minima only: connecting lines would suggest an energy path
        # or barrier that was not calculated for Tersoff.
        axes[2].plot(np.arange(4)+(j-.5)*.07, values, 'o', color=color, label=label)
    axes[2].axhline(0, color='#aaaaaa', lw=.8)
    axes[2].set_xticks(np.arange(4), ['state00', 'state04', 'state01', 'state02'])
    axes[2].set_ylabel('Energy relative to own state00 (eV)')
    axes[2].set_title('C  Fixed-gap reconstructed minima', loc='left', fontsize=11)
    axes[2].text(.03, .05, 'Identical grips / gap; conditional stability only', transform=axes[2].transAxes,
                 fontsize=8, color='#555555')
    axes[2].set_ylim(-3.6, .5)
    for axis in axes:
        axis.spines[['top', 'right']].set_visible(False)
        axis.grid(axis='y', alpha=.18)
        axis.set_axisbelow(True)
    figure.suptitle('Si atomistic foundation: numerical agreement does not establish material accuracy', fontsize=14)
    figure.savefig(root/'atomistic_audit.png', dpi=160)
    figure.savefig(root/'atomistic_audit.pdf')
    plt.close(figure)


if __name__ == '__main__':
    main()
