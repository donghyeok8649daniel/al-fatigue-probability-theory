"""Plot completed static control data. No new physical prediction is made."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
curves = np.genfromtxt(root / 'rigid_cleavage.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
grid = np.genfromtxt(root / 'shuffle_energy_grid.csv', delimiter=',', names=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, (left, right) = plt.subplots(1, 2, figsize=(11.4, 4.5), layout='constrained')
fig.suptitle('Si(111) feasibility: rigid original-SW control at 0 K', fontsize=14)
for cut, color in [('shuffle', '#176aac'), ('glide', '#c75d20')]:
    rows = curves[curves['cut'] == cut]
    left.plot(rows['opening_angstrom'], rows['energy_J_m2'], color=color, label=cut, linewidth=2.2)
left.set(xlabel='Added opening (angstrom)', ylabel='Separation work (J / m$^2$)',
         title='Two possible cuts through the diamond lattice', xlim=(0, 4), ylim=(0, 9))
left.legend(frameon=False)
left.grid(alpha=.18)
opening = np.unique(grid['opening_angstrom'])
slip = np.unique(grid['slip_angstrom'])
energy = grid['energy_eV_surface_cell'].reshape(len(opening), len(slip))
image = right.pcolormesh(slip, opening, energy, shading='nearest', cmap='magma')
right.set(xlabel='In-plane slip (angstrom)', ylabel='Added opening (angstrom)',
          title='Shuffle opening / slip energy surface')
fig.colorbar(image, ax=right, label='Energy (eV / surface cell)')
fig.text(.5, -.025, 'Rigid, unreconstructed slabs. These are not measured wafer strengths or failure probabilities.',
         ha='center', fontsize=9)
fig.savefig(root / 'static_probe.png', dpi=170, bbox_inches='tight')
fig.savefig(root / 'static_probe.svg', bbox_inches='tight')
plt.close(fig)
