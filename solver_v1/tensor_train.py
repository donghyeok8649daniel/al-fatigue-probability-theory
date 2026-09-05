"""Minimal tensor-train utilities for compressed probability representations.

These routines are numerical-compression tools only.  They do not impose a
product closure on the physics: ranks are allowed to exceed one and carry
cross-coordinate correlation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TensorTrain:
    cores: tuple[np.ndarray, ...]
    shape: tuple[int, ...]

    @property
    def ranks(self) -> tuple[int, ...]:
        if not self.cores:
            return (1,)
        return (1,) + tuple(int(core.shape[2]) for core in self.cores[:-1]) + (1,)

    @property
    def storage(self) -> int:
        return int(sum(core.size for core in self.cores))

    @property
    def dense_storage(self) -> int:
        return int(np.prod(self.shape, dtype=np.int64))

    @property
    def compression_ratio(self) -> float:
        return float(self.dense_storage / self.storage)

    def reconstruct(self) -> np.ndarray:
        if not self.cores:
            return np.empty(self.shape)
        result = self.cores[0]
        for core in self.cores[1:]:
            result = np.tensordot(result, core, axes=([-1], [0]))
        # Shape after contractions: (1,n1,n2,...,nd,1)
        result = np.squeeze(result, axis=(0, -1))
        return np.asarray(result).reshape(self.shape)


def _select_rank(
    singular_values: np.ndarray,
    *,
    local_error_budget: float,
    max_rank: int | None,
) -> int:
    full_rank = singular_values.size
    cap = full_rank if max_rank is None else min(full_rank, int(max_rank))
    if cap < 1:
        raise ValueError("max_rank must be positive when supplied")
    if local_error_budget <= 0.0:
        return cap

    sq = singular_values**2
    tail = np.concatenate((np.cumsum(sq[::-1])[::-1], np.array([0.0])))
    for rank in range(1, cap + 1):
        if tail[rank] <= local_error_budget**2:
            return rank
    return cap


def tt_svd(
    tensor: np.ndarray,
    *,
    relative_tolerance: float = 1.0e-8,
    max_rank: int | None = None,
) -> TensorTrain:
    """Compress a dense tensor with TT-SVD.

    ``relative_tolerance`` bounds the accumulated Frobenius truncation target
    in the standard TT-SVD sense.  Setting it to zero performs no tolerance-
    based truncation (subject to ``max_rank`` if supplied).
    """

    array = np.asarray(tensor, dtype=float)
    if array.ndim < 1:
        raise ValueError("tensor must have at least one dimension")
    if not np.all(np.isfinite(array)):
        raise ValueError("tensor contains non-finite values")
    if relative_tolerance < 0.0:
        raise ValueError("relative_tolerance must be nonnegative")

    shape = tuple(int(n) for n in array.shape)
    dimension = array.ndim
    norm = float(np.linalg.norm(array.ravel()))
    if norm == 0.0:
        cores = []
        left_rank = 1
        for axis, n in enumerate(shape):
            right_rank = 1
            cores.append(np.zeros((left_rank, n, right_rank), dtype=float))
            left_rank = right_rank
        return TensorTrain(tuple(cores), shape)

    local_budget = 0.0
    if relative_tolerance > 0.0 and dimension > 1:
        local_budget = relative_tolerance * norm / np.sqrt(dimension - 1.0)

    unfolding = array.copy()
    cores: list[np.ndarray] = []
    left_rank = 1

    for axis in range(dimension - 1):
        n_axis = shape[axis]
        unfolding = unfolding.reshape(left_rank * n_axis, -1)
        u, singular_values, vh = np.linalg.svd(unfolding, full_matrices=False)
        rank = _select_rank(
            singular_values,
            local_error_budget=local_budget,
            max_rank=max_rank,
        )
        u = u[:, :rank]
        singular_values = singular_values[:rank]
        vh = vh[:rank, :]
        cores.append(u.reshape(left_rank, n_axis, rank))
        unfolding = singular_values[:, None] * vh
        left_rank = rank

    cores.append(unfolding.reshape(left_rank, shape[-1], 1))
    return TensorTrain(tuple(cores), shape)


def relative_frobenius_error(reference: np.ndarray, approximation: np.ndarray) -> float:
    reference = np.asarray(reference, dtype=float)
    approximation = np.asarray(approximation, dtype=float)
    denominator = float(np.linalg.norm(reference.ravel()))
    numerator = float(np.linalg.norm((reference - approximation).ravel()))
    if denominator == 0.0:
        return 0.0 if numerator == 0.0 else np.inf
    return numerator / denominator
