import numpy as np

from solver_v1.model import ModelParams
from solver_v1.probability_tt_6d import Grid6DParams, compress_initial_gibbs_6d


def test_n3_gibbs_tensor_tt_prototype_preserves_mass_and_accuracy():
    result = compress_initial_gibbs_6d(
        model_params=ModelParams(
            n_cells=3,
            kT=0.03,
            mobility_a=1.0,
            mobility_s=0.05,
            chi_axial_projection=0.20,
        ),
        grid_params=Grid6DParams(n_a=4, n_s=5, s_wells=1, a_upper=1.45),
        relative_tolerance=1.0e-9,
    )
    assert abs(float(result["exact_mass"]) - 1.0) < 1.0e-12
    assert abs(float(result["reconstructed_mass"]) - 1.0) < 1.0e-6
    assert float(result["relative_frobenius_error"]) < 2.0e-8
    assert float(result["negative_mass"]) < 1.0e-7
    assert max(result["tt_ranks"]) > 1


def test_n3_tt_prototype_reports_real_compression_metrics():
    result = compress_initial_gibbs_6d(
        model_params=ModelParams(n_cells=3, kT=0.04),
        grid_params=Grid6DParams(n_a=4, n_s=4, s_wells=1, a_upper=1.42),
        relative_tolerance=1.0e-6,
    )
    assert int(result["dense_storage"]) == 4**6
    assert int(result["tt_storage"]) > 0
    assert np.isfinite(float(result["compression_ratio"]))
