"""Exact finite-free-domain analytic Hessian assembly, not a new operator.

H=sum_i D_i^T e_i,zz D_i + sum_(i,R) B_iR^T L_iR B_iR.
Only the free degrees of freedom are assembled; all affected site energies
remain. This independently accelerates dense lowest-mode checks on small and
moderate research cores. The matrix-free implementation remains the default.
"""
import numpy as np


def assemble_core_hessian(free_ids, energy_ids, destination, gradients, local, site_hessian):
    free_ids, energy_ids, destination = map(np.asarray, (free_ids, energy_ids, destination))
    count = len(free_ids)
    total_sites = max(int(free_ids.max()), int(energy_ids.max()), int(destination.max()))+1
    mapping = np.full(total_sites, -1, int); mapping[free_ids] = np.arange(count)
    source, target = mapping[energy_ids], mapping[destination]
    sites, neighbors, channels, _ = gradients.shape
    # Recover local channel Hessians through their existing analytic action.
    channel_hessians = np.stack([site_hessian(np.broadcast_to(np.eye(channels)[j], (sites, channels)))
                                for j in range(channels)], axis=2)
    matrix = np.zeros((3*count, 3*count))
    diagonal = np.zeros((count, 3, 3))
    valid_source = source >= 0
    np.add.at(diagonal, source[valid_source], local[valid_source].sum(axis=1))
    valid_target = target >= 0
    np.add.at(diagonal, target[valid_target], local[valid_target])
    for number in range(count):
        matrix[3*number:3*number+3, 3*number:3*number+3] += diagonal[number]
    both = valid_target & valid_source[:, None]
    src = np.broadcast_to(source[:, None], target.shape)[both]
    dst = target[both]
    row = 3*src[:, None]+np.arange(3)
    col = 3*dst[:, None]+np.arange(3)
    np.add.at(matrix, (row[:, :, None], col[:, None, :]), -local[both])
    np.add.at(matrix, (col[:, None, :], row[:, :, None]), -local[both])
    for site in range(sites):
        indices = target[site][target[site] >= 0]
        if source[site] >= 0:
            indices = np.r_[indices, source[site]]
        indices = np.unique(indices)
        if not len(indices):
            continue
        inverse = {value: j for j, value in enumerate(indices)}
        derivative = np.zeros((channels, 3*len(indices)))
        for neighbor in np.flatnonzero(target[site] >= 0):
            start = 3*inverse[target[site, neighbor]]
            derivative[:, start:start+3] += gradients[site, neighbor]
        if source[site] >= 0:
            start = 3*inverse[source[site]]
            derivative[:, start:start+3] -= gradients[site].sum(axis=0)
        dofs = (3*indices[:, None]+np.arange(3)).ravel()
        matrix[np.ix_(dofs, dofs)] += derivative.T@channel_hessians[site]@derivative
    # Do not symmetrize an implementation error away. The caller can inspect
    # the independent antisymmetric residual before any eigensolver.
    return matrix
