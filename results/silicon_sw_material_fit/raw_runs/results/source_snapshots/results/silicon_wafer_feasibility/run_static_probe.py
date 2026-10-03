"""Independent Si feasibility calculation, not a production or calibrated model.

Pure Si, original Stillinger-Weber, rigid unreconstructed (111) cleavage at 0 K.
The comparison is a geometry/unit/control test, NOT validation of SW fracture.
No MD, probability PDE, mobility fit, wafer solver or existing module is changed.
Run from the repository root: python results/silicon_wafer_feasibility/run_static_probe.py
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
from pathlib import Path
import time
from urllib.request import urlopen

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree


PARAMETER_URL = 'https://raw.githubusercontent.com/lammps/lammps/develop/potentials/Si.sw'
PARAMETER_SHA256 = 'c6d3a7d26db28ee8e3feaf2f8059607de9a4accb45d9719fcabe59488ecf1cad'
EV_A2_TO_J_M2 = 16.02176634
EV_A3_TO_GPA = 160.2176634
ROOT = Path(__file__).resolve().parent


def source_parameters():
    cached = ROOT / 'source_Si.sw'
    if not cached.exists():
        with urlopen(PARAMETER_URL, timeout=30) as response:
            data = response.read()
        cached.write_bytes(data)
    data = cached.read_bytes()
    if hashlib.sha256(data).hexdigest() != PARAMETER_SHA256:
        raise ValueError('source file differs from the audited original-SW parameters')
    tokens = ' '.join(line.split('#')[0] for line in data.decode().splitlines()).split()
    if tokens[:3] != ['Si', 'Si', 'Si'] or len(tokens) != 14:
        raise ValueError('unexpected source parameter format')
    keys = 'epsilon sigma a lambda gamma costheta0 A B p q tol'.split()
    params = dict(zip(keys, map(float, tokens[3:])))
    return params, hashlib.sha256(data).hexdigest()


def pair(r, p, derivative=False):
    r = np.asarray(r, float)
    x = r / p['sigma']
    out = np.zeros_like(x)
    active = x < p['a']
    y = x[active]
    v = p['B'] * y**(-p['p']) - y**(-p['q'])
    exponential = np.exp(1 / (y - p['a']))
    if derivative:
        dv = (-p['p'] * p['B'] * y**(-p['p'] - 1)
              + p['q'] * y**(-p['q'] - 1))
        out[active] = p['epsilon'] * p['A'] / p['sigma'] * exponential * (
            dv - v / (y - p['a'])**2)
    else:
        out[active] = p['epsilon'] * p['A'] * v * exponential
    return out


def diamond111_basis():
    """Six atoms in a primitive surface cell and three-bilayer [111] period."""
    cell = np.array([[.5, -.5, 0], [.5, 0, -.5], [1, 1, 1]])
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    diamond = np.concatenate([fcc, fcc + .25])
    translations = np.array(list(itertools.product(range(-2, 4), repeat=3)))
    positions = (translations[:, None, :] + diamond).reshape(-1, 3)
    fractional = positions @ np.linalg.inv(cell)
    keep = np.all((fractional >= -1e-10) & (fractional < 1 - 1e-10), axis=1)
    result = fractional[keep]
    assert result.shape == (6, 3)
    return result, cell


def analytic_shuffle_hessian(bond, p):
    """Independent bond/angle expansion for rigid shuffle cleavage at zero shift."""
    x, gap = bond / p['sigma'], bond / p['sigma'] - p['a']
    v = p['B'] * x**(-p['p']) - x**(-p['q'])
    dv = -p['p'] * p['B'] * x**(-p['p'] - 1) + p['q'] * x**(-p['q'] - 1)
    ddv = p['p'] * (p['p'] + 1) * p['B'] * x**(-p['p'] - 2) - p['q'] * (p['q'] + 1) * x**(-p['q'] - 2)
    haa = p['epsilon'] * p['A'] / p['sigma']**2 * np.exp(1 / gap) * (
        ddv - 2 * dv / gap**2 + v * (2 / gap**3 + 1 / gap**4))
    # Each endpoint has three tetrahedral bonds; the transverse dyadic sum
    # is 4/3 per endpoint. The second derivative contributes another factor 2.
    hss = float(pair(bond, p, derivative=True)) / bond + 16 / 3 * p['epsilon'] * p['lambda'] * np.exp(
        2 * p['gamma'] / gap) / bond**2
    return np.diag([haa, hss])


class RigidSlab:
    def __init__(self, p, lattice, cut_kind, repeats=2, periods=3, image_shell=2):
        basis, cell = diamond111_basis()
        cell *= lattice
        copies = np.array(list(itertools.product(range(repeats), range(repeats), range(periods))))
        self.positions = ((copies[:, None, :] + basis) @ cell).reshape(-1, 3)
        self.cell = cell * np.array([repeats, repeats, periods])[:, None]
        self.invcell = np.linalg.inv(self.cell)
        self.normal = np.ones(3) / np.sqrt(3)
        self.slip = cell[0] / np.linalg.norm(cell[0])
        self.period = np.linalg.norm(cell[0])
        self.cells = repeats**2
        self.area = np.linalg.norm(np.cross(cell[0], cell[1]))
        heights = np.unique(np.round(self.positions @ self.normal, 10))
        spacing = np.diff(heights)
        threshold = .5 * (spacing.min() + spacing.max())
        candidates = np.flatnonzero(spacing > threshold if cut_kind == 'shuffle' else spacing < threshold)
        cut = candidates[np.argmin(abs(candidates + .5 - .5 * (len(heights) - 1)))]
        self.gap = float(spacing[cut])
        self.upper = self.positions @ self.normal > .5 * (heights[cut] + heights[cut + 1])
        images = np.array(list(itertools.product(range(-image_shell, image_shell + 1), repeat=2)))
        self.images = images @ self.cell[:2]
        self.p = p
        self.reference = self.total_energy(self.positions)

    def total_energy(self, positions):
        fractional = positions @ self.invcell
        fractional[:, :2] %= 1.
        positions = fractional @ self.cell
        replicas = (self.images[:, None, :] + positions).reshape(-1, 3)
        cutoff = self.p['a'] * self.p['sigma']
        tree = cKDTree(replicas)
        neighbors = tree.query_ball_point(positions, cutoff)
        pair_energy = 0.
        angular_energy = 0.
        for center, ids in zip(positions, neighbors):
            vectors = replicas[ids] - center
            lengths = np.linalg.norm(vectors, axis=1)
            active = (lengths > 1e-8) & (lengths < cutoff)
            vectors, lengths = vectors[active], lengths[active]
            pair_energy += .5 * np.sum(pair(lengths, self.p))
            left, right = np.triu_indices(len(lengths), k=1)
            cosine = np.einsum('ij,ij->i', vectors[left], vectors[right]) / (lengths[left] * lengths[right])
            envelope = np.exp(self.p['gamma'] * self.p['sigma'] * (
                1 / (lengths[left] - cutoff) + 1 / (lengths[right] - cutoff)))
            angular_energy += self.p['epsilon'] * self.p['lambda'] * np.sum(
                envelope * (cosine - self.p['costheta0'])**2)
        return float(pair_energy + angular_energy)

    def energy(self, opening, slip=0.):
        positions = self.positions.copy()
        positions[self.upper] += opening * self.normal + slip * self.slip
        return (self.total_energy(positions) - self.reference) / self.cells


def hessian(energy, step):
    zero = energy(0., 0.)
    haa = (energy(step, 0.) + energy(-step, 0.) - 2 * zero) / step**2
    hss = (energy(0., step) + energy(0., -step) - 2 * zero) / step**2
    has = (energy(step, step) - energy(step, -step) - energy(-step, step)
           + energy(-step, -step)) / (4 * step**2)
    return np.array([[haa, has], [has, hss]])


def elasticity_probe():
    # Hopcroft, Nix, Kenny (2010), DOI 10.1109/JMEMS.2009.2039697.
    # Room-temperature experimental reference; separate from the 0 K SW calculation.
    c11, c12, c44 = 165.7, 63.9, 79.6  # GPa
    c = np.zeros((6, 6))
    c[:3, :3] = c12
    c[np.arange(3), np.arange(3)] = c11
    c[np.arange(3, 6), np.arange(3, 6)] = c44
    compliance = np.linalg.inv(c)
    moduli = {}
    for label, direction in [('100', [1, 0, 0]), ('110', [1, 1, 0]), ('111', [1, 1, 1])]:
        e = np.array(direction, float)
        e /= np.linalg.norm(e)
        mixed = e[0]**2 * e[1]**2 + e[1]**2 * e[2]**2 + e[2]**2 * e[0]**2
        inverse_e = compliance[0, 0] - 2 * (
            compliance[0, 0] - compliance[0, 1] - .5 * compliance[3, 3]) * mixed
        moduli[label] = float(1 / inverse_e)
    normal = np.ones(3) / np.sqrt(3)
    m1, m2 = np.array([1., -1., 0.]) / np.sqrt(2), np.array([1., 1., -2.]) / np.sqrt(6)
    tractions = {}
    for label, e in [('normal_111', normal), ('axis_001', np.array([0., 0., 1.]))]:
        traction = np.outer(e, e) @ normal  # unit applied axial stress
        tractions[label] = dict(zip(['normal', 'shear1', 'shear2'], map(float,
            [normal @ traction, m1 @ traction, m2 @ traction])))
    assert abs(tractions['normal_111']['shear1']) < 1e-14
    assert abs(tractions['normal_111']['shear2']) < 1e-14
    return {'reference_C_GPa': {'C11': c11, 'C12': c12, 'C44': c44},
            'young_modulus_GPa': moduli, 'resolved_traction_per_unit_axial_stress': tractions}


def main():
    start = time.perf_counter()
    p, checksum = source_parameters()
    optimum = minimize_scalar(lambda r: float(pair(r, p)), bounds=(2.2, 2.5),
                              method='bounded', options={'xatol': 1e-14})
    bond = float(optimum.x)
    lattice = 4 * bond / np.sqrt(3)
    summary = {'scope': '0 K rigid original-SW control; not relaxed Si fracture or a calibration',
        'source_url': PARAMETER_URL, 'source_sha256': checksum, 'parameters': p,
        'lattice_angstrom': lattice, 'bond_angstrom': bond,
        'bulk_energy_eV_atom': float(2 * pair(bond, p)), 'elasticity_separate_RT_reference': elasticity_probe()}
    rows = []
    for kind, expected_bonds in [('shuffle', 1), ('glide', 3)]:
        slab = RigidSlab(p, lattice, kind)
        expected_sep = -expected_bonds * float(pair(bond, p))
        sep = slab.energy(6.)
        assert abs(sep - expected_sep) < 1e-10
        hs = hessian(slab.energy, 1e-3)
        hf = hessian(slab.energy, 5e-4)
        h_error = float(np.max(abs(hf - hs)))
        assert np.all(np.linalg.eigvalsh(hf) > 0)
        assert h_error < 1e-4
        periodic_error = max(abs(slab.energy(.3, slip + slab.period) - slab.energy(.3, slip))
                             for slip in [0., .23 * slab.period, .51 * slab.period])
        assert periodic_error < 1e-10
        repeat_errors = []
        for repeats, periods, images in [(1, 2, 2), (2, 4, 2), (3, 3, 2), (2, 3, 3)]:
            control = RigidSlab(p, lattice, kind, repeats, periods, images)
            error = max(abs(control.energy(a, s) - slab.energy(a, s))
                        for a, s in [(0., 0.), (.3, .21 * slab.period), (1., 0.), (6., 0.)])
            repeat_errors.append({'in_plane_repeats': repeats, 'normal_periods': periods,
                                  'image_shell': images, 'max_energy_error_eV_cell': error})
            assert error < 1e-10
        detail = {'initial_gap_angstrom': slab.gap, 'surface_cell_area_angstrom2': slab.area,
            'slip_period_angstrom': slab.period, 'separation_eV_cell': sep,
            'separation_J_m2_two_created_surfaces': sep / slab.area * EV_A2_TO_J_M2,
            'hessian_eV_angstrom2': hf.tolist(), 'hessian_step_difference': h_error,
            'periodicity_max_error_eV_cell': periodic_error, 'finite_slab_and_image_checks': repeat_errors}
        if kind == 'shuffle':
            analytic_h = analytic_shuffle_hessian(bond, p)
            analytic_h_error = float(np.max(abs(analytic_h - hf)))
            assert analytic_h_error < 1e-5
            openings = np.linspace(0., 4., 41)
            analytic_error = max(abs(slab.energy(a) - float(pair(bond + a, p) - pair(bond, p)))
                                 for a in openings)
            assert analytic_error < 1e-10
            peak = minimize_scalar(lambda r: -float(pair(r, p, derivative=True)),
                bounds=(bond + 1e-4, p['a'] * p['sigma'] - 1e-5), method='bounded',
                options={'xatol': 1e-12})
            detail.update({'analytic_shuffle_curve_max_error_eV_cell': analytic_error,
                'analytic_hessian_eV_angstrom2': analytic_h.tolist(),
                'analytic_hessian_max_error_eV_angstrom2': analytic_h_error,
                'rigid_ideal_peak_traction_GPa': -float(peak.fun) / slab.area * EV_A3_TO_GPA,
                'opening_at_peak_angstrom': float(peak.x - bond)})
        for a in np.linspace(0., 4., 81):
            w = slab.energy(a)
            rows.append([kind, float(a), w, w / slab.area * EV_A2_TO_J_M2])
        summary[kind] = detail
        print(kind, json.dumps(detail), flush=True)
    shuffle = RigidSlab(p, lattice, 'shuffle')
    with (ROOT / 'shuffle_energy_grid.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['opening_angstrom', 'slip_angstrom', 'energy_eV_surface_cell'])
        for a in np.linspace(0., 2.5, 11):
            for s in np.linspace(-.5 * shuffle.period, .5 * shuffle.period, 17):
                writer.writerow([a, s, shuffle.energy(a, s)])
    with (ROOT / 'rigid_cleavage.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['cut', 'opening_angstrom', 'energy_eV_surface_cell', 'energy_J_m2'])
        writer.writerows(rows)
    summary['elapsed_seconds'] = time.perf_counter() - start
    summary['verification'] = 'PASS: analytical shuffle curve and Hessian, broken-bond limit, periodicity, positive local Hessian, h refinement, slab/area/image controls, traction orientation'
    (ROOT / 'static_probe_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ['lattice_angstrom', 'bulk_energy_eV_atom', 'elapsed_seconds', 'verification']}, indent=2))


if __name__ == '__main__':
    main()
