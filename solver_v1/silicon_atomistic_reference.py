"""Optional, source-bound atomistic Si references; never a production model.

Coordinates/cell rows are in Angstrom, total energies in eV, gradients in
eV/Angstrom, and tensile stresses in eV/Angstrom**3. LAMMPS evaluates every
periodic image, including cells shorter than twice the potential cutoff.
The original coordinate frame and atom ordering are preserved on return.
Neither this adapter nor Newtonian atomic masses supplies a probability clock.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
import hashlib
from pathlib import Path

import numpy as np


EV_A3_TO_BAR = 1602176.634


def periodic_basis_audit(cell, periodic):
    """Explicitly reindex a left-handed fully periodic lattice, without moving R.

    The integer basis operation diag(1,1,-1) preserves the periodic lattice.
    Return its signs for the provenance record. This is NOT a Cartesian mirror.
    Mixed/free boundary geometry is rejected instead of reinterpreted.
    """
    cell = np.asarray(cell, float)
    flags = np.asarray(periodic)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
        raise ValueError('finite (3,3) cell required')
    if flags.shape != (3,) or flags.dtype != np.dtype(bool):
        raise ValueError('three boolean periodic flags required')
    signs = np.ones(3, int)
    if np.linalg.det(cell) < 0:
        if not np.all(flags):
            raise ValueError('left-handed nonperiodic cell requires explicit geometry review')
        signs[2] = -1
    changed = signs[:, None]*cell
    restricted_cell(changed)  # Also reject singular/ill-conditioned geometry.
    return changed, signs


def restricted_cell(cell):
    """Return B, Q with B = cell @ Q lower triangular and det(Q) = +1.

    Cell vectors are ROWS. Positions transform as r_new = r @ Q and forces
    return as f = f_new @ Q.T. Reflections and singular cells are rejected.
    """
    c = np.asarray(cell, float)
    if (c.shape != (3, 3) or not np.all(np.isfinite(c))
            or np.linalg.det(c) <= 0 or np.linalg.cond(c) > 1e10):
        raise ValueError('finite, nonsingular right-handed (3,3) cell required')
    q, r = np.linalg.qr(c.T)
    signs = np.where(np.diag(r) < 0, -1., 1.)
    q = q * signs
    b = c @ q
    b[np.triu_indices(3, 1)] = 0.
    return b, q


@dataclass
class ReferenceEvaluation:
    energy: float
    gradient: np.ndarray
    site_energy: np.ndarray
    stress: np.ndarray


class LammpsSilicon:
    """Serial SW or ordinary Tersoff reference using an explicit potential file.

    Import of the optional engine occurs only on construction. A required
    checksum prevents a different Si parameterization from being substituted.
    One instance must not be shared between concurrent threads/processes.
    Species are pure Si only: dopants and electronic charge are not supported.
    """
    def __init__(self, style, potential, *, sha256):
        if style not in ('sw', 'tersoff'):
            raise ValueError('only explicitly audited sw/tersoff styles supported')
        potential = Path(potential).resolve()
        actual_hash = hashlib.sha256(potential.read_bytes()).hexdigest()
        if actual_hash != sha256:
            raise ValueError('potential file checksum differs from declared source')
        from lammps import lammps
        self.lmp = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
        self.style, self.potential, self.sha256 = style, potential, actual_hash
        self._key = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if self.lmp is not None:
            self.lmp.close()
            self.lmp = None

    def metadata(self):
        library = Path(self.lmp.lib._name)
        return dict(style=self.style, potential_filename=self.potential.name,
                    potential_sha256=self.sha256, lammps_version=int(self.lmp.version()),
                    library_filename=library.name,
                    library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
                    units='metal: eV, Angstrom; no dynamics', pure_silicon=True)

    def evaluate(self, positions, cell, *, periodic=(True, True, True)):
        positions, cell = np.asarray(positions, float), np.asarray(cell, float)
        periodic = np.asarray(periodic)
        if (positions.ndim != 2 or positions.shape[1] != 3 or not len(positions)
                or not np.all(np.isfinite(positions))):
            raise ValueError('nonempty finite (atoms,3) positions required')
        if periodic.shape != (3,) or periodic.dtype != np.dtype(bool):
            raise ValueError('three boolean periodic flags required')
        box, rotation = restricted_cell(cell)
        fractional = np.linalg.solve(cell.T, positions.T).T
        fractional[:, periodic] %= 1.
        if np.any((fractional[:, ~periodic] <= 0) | (fractional[:, ~periodic] >= 1)):
            raise ValueError('atoms must lie strictly inside nonperiodic cell faces')
        current = np.ascontiguousarray(fractional @ box, dtype=np.float64)
        n = len(positions)
        key = (n, cell.tobytes(), periodic.tobytes())
        lmp = self.lmp
        if lmp is None:
            raise RuntimeError('reference engine has been closed')
        if key != self._key:
            lmp.command('clear')
            lmp.command('units metal')
            lmp.command('atom_style atomic')
            lmp.command('atom_modify map array sort 0 0.0')
            lmp.command('boundary ' + ' '.join('p' if b else 'f' for b in periodic))
            lmp.command('box tilt large')
            lx, ly, lz = np.diag(box)
            xy, xz, yz = box[1, 0], box[2, 0], box[2, 1]
            lmp.command(f'region box prism 0 {lx:.17g} 0 {ly:.17g} 0 {lz:.17g} '
                        f'{xy:.17g} {xz:.17g} {yz:.17g} units box')
            lmp.command('create_box 1 box')
            made = lmp.create_atoms(n, list(range(1, n+1)), [1]*n,
                                   current.ravel().tolist())
            if made != n:
                raise RuntimeError('LAMMPS did not create every supplied atom')
            lmp.command('mass 1 28.0855')  # Required by engine; run 0 only.
            lmp.command(f'pair_style {self.style}')
            lmp.command(f'pair_coeff * * "{self.potential.as_posix()}" Si')
            lmp.command('neighbor 0.6 bin')
            lmp.command('neigh_modify every 1 delay 0 check yes')
            lmp.command('compute atomic_energy all pe/atom')
            lmp.command('thermo_style custom pe pxx pyy pzz pxy pxz pyz')
            lmp.command('thermo_modify norm no')
            self._key = key
        else:
            pointer = current.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            lmp.scatter_atoms('x', 1, 3, pointer)
        lmp.command('run 0 post no')
        if int(lmp.get_natoms()) != n:
            raise RuntimeError('atom count changed during run-0 evaluation')
        energy = float(lmp.get_thermo('pe'))
        f = np.ctypeslib.as_array(lmp.gather_atoms('f', 1, 3), shape=(3*n,)).reshape(n, 3).copy()
        ids = lmp.numpy.extract_atom('id')[:n].copy()-1
        site = np.empty(n)
        site[ids] = lmp.numpy.extract_compute('atomic_energy', 1, 1)[:n]
        xx, yy, zz, xy, xz, yz = [float(lmp.get_thermo(k)) for k in
                                  ('pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz')]
        stress = -np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])/EV_A3_TO_BAR
        result = ReferenceEvaluation(energy, -f @ rotation.T, site,
                                     rotation @ stress @ rotation.T)
        if not (np.isfinite(energy) and np.all(np.isfinite(result.gradient))
                and np.all(np.isfinite(result.site_energy))
                and np.all(np.isfinite(result.stress))):
            raise FloatingPointError('nonfinite atomistic reference result')
        if abs(float(np.sum(site))-energy) > 1e-8*max(1., abs(energy)):
            raise RuntimeError('per-site energy does not sum to total energy')
        return result


def force_metrics(predicted, reference):
    """Cartesian component RMSE, not the RMS vector norm (sqrt(3) different)."""
    predicted, reference = np.asarray(predicted, float), np.asarray(reference, float)
    if (predicted.shape != reference.shape or predicted.ndim != 2
            or predicted.shape[1] != 3 or not len(predicted)
            or not np.all(np.isfinite(predicted)) or not np.all(np.isfinite(reference))):
        raise ValueError('matching finite (atoms,3) forces required')
    diff = predicted-reference
    return dict(components=int(diff.size), squared_error=float(np.sum(diff**2)),
                component_rmse_eV_A=float(np.sqrt(np.mean(diff**2))),
                max_vector_error_eV_A=float(np.max(np.linalg.norm(diff, axis=1))),
                reference_component_rms_eV_A=float(np.sqrt(np.mean(reference**2))))
