"""Regression for decimal load labels: atomic states must never overwrite."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from results.silicon_wafer_feasibility.run_specimen_bridge import load_case, source_cutoff
from .silicon_specimen_research import AnisotropicModeI, atomic_crack_boundary


class StorageTests(unittest.TestCase):
    def test_decimal_cases_are_distinct_and_existing_state_is_protected(self):
        class StationaryEngine:
            def evaluate(self, positions, cell, **kwargs):
                return SimpleNamespace(energy=0., gradient=np.zeros_like(positions),
                                       site_energy=np.zeros(len(positions)))
        engine = StationaryEngine()
        boundary = atomic_crack_boundary(5.43, radius_A=20., front_repeats=1)
        elastic = AnisotropicModeI(151.4, 76.4, 56.4)
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory)/f'state_K{k:g}.npz' for k in [.8, 1., 1.2]]
            for path, k in zip(paths, [.8, 1., 1.2]):
                load_case(engine, elastic, boundary, k, None, path)
            self.assertEqual(len(list(Path(directory).glob('*.npz'))), 3)
            with np.load(paths[0]) as a, np.load(paths[1]) as b, np.load(paths[2]) as c:
                self.assertGreater(np.max(abs(c['positions']-b['positions'])), .01)
                self.assertGreater(np.max(abs(a['positions']-b['positions'])), .01)
            original = paths[1].read_bytes()
            with self.assertRaises(FileExistsError):
                load_case(engine, elastic, boundary, 1.2, None, paths[1])
            self.assertEqual(paths[1].read_bytes(), original)
            with self.assertRaises(ValueError):
                load_case(engine, elastic, boundary, 1.2, None, Path(directory)/'state_K1.2')

    def test_actual_source_cutoffs(self):
        root = Path(__file__).resolve().parents[1]
        self.assertAlmostEqual(source_cutoff('sw', root/'results/silicon_wafer_feasibility/source_Si.sw'), 3.77118)
        self.assertAlmostEqual(source_cutoff('tersoff', root/'results/silicon_atomistic_v5/sources/SiC.tersoff'), 3.0)


if __name__ == '__main__':
    unittest.main()
