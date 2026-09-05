import numpy as np

from solver_v1.tensor_train import relative_frobenius_error, tt_svd


def test_tt_svd_exact_reconstruction_without_truncation():
    x = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
    tt = tt_svd(x, relative_tolerance=0.0)
    xr = tt.reconstruct()
    assert xr.shape == x.shape
    assert relative_frobenius_error(x, xr) < 1.0e-12


def test_tt_rank_greater_than_one_represents_correlation():
    i = np.linspace(-1.0, 1.0, 7)
    j = np.linspace(-1.0, 1.0, 7)
    ii, jj = np.meshgrid(i, j, indexing="ij")
    correlated = np.exp(-4.0 * (ii - jj) ** 2)
    tt = tt_svd(correlated, relative_tolerance=1.0e-10)
    assert max(tt.ranks) > 1
    assert relative_frobenius_error(correlated, tt.reconstruct()) < 1.0e-8


def test_tt_compresses_simple_low_rank_tensor():
    a = np.linspace(1.0, 2.0, 8)
    b = np.linspace(0.5, 1.5, 9)
    c = np.linspace(2.0, 3.0, 10)
    tensor = a[:, None, None] * b[None, :, None] * c[None, None, :]
    tt = tt_svd(tensor, relative_tolerance=1.0e-12)
    assert tt.compression_ratio > 1.0
    assert relative_frobenius_error(tensor, tt.reconstruct()) < 1.0e-10
