"""Feed the Si static table to the unchanged SG kernel; numerical reuse only.

Arbitrary mobilities and an imposed reflecting box: no Si kinetics, absorbing
crack boundary, physical time, specimen probability, or calibration is claimed.
"""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import sys

import numpy as np

root = Path(__file__).resolve().parent
sys.path.insert(0, str(root.parents[1]))
from solver_v1.probability_pde_2d import Grid2D, _sg_generator_2d, _implicit_step_2d

table_path = root / 'shuffle_energy_grid.csv'
table = np.genfromtxt(table_path, delimiter=',', names=True)
opening, slip = np.unique(table['opening_angstrom']), np.unique(table['slip_angstrom'])
energy = table['energy_eV_surface_cell'].reshape(len(opening), len(slip))
grid = Grid2D(opening, slip, float(opening[1] - opening[0]), float(slip[1] - slip[0]))
model = SimpleNamespace(p=SimpleNamespace(kT=.02, mobility_a=1., mobility_s=.05))
generator = _sg_generator_2d(energy, model, grid)
column_error = float(np.max(abs(np.asarray(generator.sum(axis=0)))))
offdiag = generator.tocoo()
minimum_rate = float(offdiag.data[offdiag.row != offdiag.col].min())
gibbs = np.exp(-(energy - energy.min()) / model.p.kT)
gibbs /= gibbs.sum() * grid.cell_volume
stationary_error = float(np.max(abs(generator @ gibbs.ravel())))
scale = max(1., float(abs(generator).sum(axis=0).max()) * float(gibbs.max()))
relative_stationary_error = stationary_error / scale
assert column_error < 1e-10
assert minimum_rate >= 0.
assert relative_stationary_error < 1e-12

density = np.zeros_like(energy)
density[len(opening) // 2, len(slip) // 2 + 2] = 1 / grid.cell_volume
mass_error = 0.
minimum_density = float('inf')
correction = 0.
for _ in range(20):
    density, repair, minimum = _implicit_step_2d(density, energy, model, grid, .005, 1e-12)
    mass_error = max(mass_error, abs(float(density.sum() * grid.cell_volume) - 1.))
    minimum_density = min(minimum_density, minimum)
    correction += repair
assert mass_error < 1e-12
assert minimum_density >= 0.
assert correction == 0.
source = root.parents[1] / 'solver_v1/probability_pde_2d.py'
result = {'scope': 'numerical SG reuse on imposed reflecting Si-energy box; no physical kinetics or crack probability',
    'table_sha256': hashlib.sha256(table_path.read_bytes()).hexdigest(),
    'existing_kernel_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'grid': list(energy.shape), 'arbitrary_mobility_a': 1., 'arbitrary_mobility_s': .05,
    'numerical_kT_eV': .02, 'dt_arbitrary_time': .005, 'steps': 20,
    'max_generator_column_sum': column_error, 'minimum_offdiagonal_rate': minimum_rate,
    'gibbs_stationarity_max_absolute_residual': stationary_error,
    'gibbs_stationarity_scaled_residual': relative_stationary_error,
    'max_mass_error': mass_error, 'minimum_density': minimum_density,
    'negative_mass_correction': correction, 'verification': 'PASS'}
(root / 'sg_reuse_summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
