from pathlib import Path

import numpy as np

from simulations.fem_tension_ui import load_fem_history
from simulations.fvm1d_solver import run


def test_fvm_cylinder_tracks_poisson_diameter_with_zero_transverse_stress(tmp_path: Path):
    run(
        elements=4,
        length_m=0.05,
        area_m2=np.pi * 0.006**2 / 4.0,
        young_pa=70.0e9,
        diameter_m=0.006,
        poisson_ratio=0.33,
        tensile_axis=(0.0, 2.0, 0.0),
        stress_mean_mpa=70.0,
        stress_amplitude_mpa=0.0,
        frequency_hz=1.0,
        cycles=1,
        steps_per_cycle=4,
        outdir=tmp_path,
    )
    nodes, elements = load_fem_history(tmp_path)

    np.testing.assert_allclose(elements["strain"], 1.0e-3)
    np.testing.assert_allclose(elements["transverse_strain"], -0.33e-3)
    np.testing.assert_allclose(elements["diameter_m"], 0.006 * (1.0 - 0.33e-3))
    np.testing.assert_allclose(elements["transverse_stress_pa"], 0.0)
    np.testing.assert_allclose(nodes["position_x_m"], 0.0)
    np.testing.assert_allclose(nodes["position_z_m"], 0.0)
    np.testing.assert_allclose(nodes["position_y_m"], nodes["x_m"])
    np.testing.assert_allclose(nodes["displacement_y_m"], nodes["displacement_m"])

    metadata = (tmp_path / "metadata.csv").read_text(encoding="utf-8")
    assert "tensile_axis_y,1" in metadata
    assert "transverse_applied_stress_pa,0" in metadata
