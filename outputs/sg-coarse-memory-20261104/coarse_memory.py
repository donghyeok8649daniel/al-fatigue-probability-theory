"""Finite-state aggregation audit in the column-vector convention."""
from __future__ import annotations

import numpy as np


def aggregation(labels: np.ndarray, lifting_weights: np.ndarray):
    """Return summation P and a positive lifting R with P R = I."""
    labels = np.asarray(labels, dtype=int)
    weights = np.asarray(lifting_weights, dtype=float)
    groups = np.unique(labels)
    if not np.array_equal(groups, np.arange(groups.size)):
        raise ValueError("Labels must be consecutive and nonnegative")
    if weights.shape != labels.shape or np.any(weights < 0):
        raise ValueError("Nonnegative lifting weights must match labels")
    p = np.zeros((groups.size, labels.size))
    r = np.zeros((labels.size, groups.size))
    for group in groups:
        cells = labels == group
        total = float(weights[cells].sum())
        if total <= 0:
            raise ValueError("Every aggregate needs positive lifting mass")
        p[group, cells] = 1.0
        r[cells, group] = weights[cells] / total
    return p, r


def exact_memory(steps, p, r, initial, *, retain_initial=True):
    """Eliminate hidden states into an ordered, nonstationary history sum."""
    liftings = [r] * (len(steps) + 1) if isinstance(r, np.ndarray) and r.ndim == 2 else r
    if len(liftings) != len(steps) + 1:
        raise ValueError("Moving liftings require one lifting at each time node")
    q = [np.eye(initial.size) - lifting @ p for lifting in liftings]
    b = [p @ t @ liftings[n] for n, t in enumerate(steps)]
    h = [p @ t @ q[n] for n, t in enumerate(steps)]
    c = [q[n + 1] @ t @ liftings[n] for n, t in enumerate(steps)]
    d = [q[n + 1] @ t @ q[n] for n, t in enumerate(steps)]
    hidden_initial = q[0] @ initial
    history = [p @ initial]
    initial_terms = []
    memory_terms = []
    for n in range(len(steps)):
        noise = hidden_initial.copy()
        for k in range(n):
            noise = d[k] @ noise
        initial_term = h[n] @ noise
        memory_term = np.zeros(p.shape[0])
        for j in range(n):
            contribution = c[j] @ history[j]
            for k in range(j + 1, n):
                contribution = d[k] @ contribution
            memory_term += h[n] @ contribution
        updated = b[n] @ history[n] + memory_term
        if retain_initial:
            updated += initial_term
        history.append(updated)
        initial_terms.append(initial_term)
        memory_terms.append(memory_term)
    return np.asarray(history), np.asarray(initial_terms), np.asarray(memory_terms)


def markov_audit(steps, fine_path, p, r, valid_states):
    """Compute finite-horizon TV certificates and actual aggregate errors."""
    coarse = p @ fine_path[0]
    posterior = 0.0
    uniform = 0.0
    records = []
    for n, t in enumerate(steps):
        lifting = r if isinstance(r, np.ndarray) and r.ndim == 2 else r[n]
        b = p @ t @ lifting
        defect = p @ t - b @ p
        posterior += 0.5 * float(np.abs(defect @ fine_path[n]).sum())
        # The cemetery column has zero defect. Valid surviving mass may be <1.
        alive = valid_states[n].copy()
        alive[-1] = False
        epsilon = 0.5 * float(np.max(np.abs(defect[:, alive]).sum(axis=0)))
        uniform += epsilon * float(fine_path[n, :-1].sum())
        coarse = b @ coarse
        truth = p @ fine_path[n + 1]
        records.append({
            "step": n + 1,
            "cdf_truth": float(fine_path[n + 1, -1]),
            "cdf_markov": float(coarse[-1]),
            "cdf_error": float(abs(truth[-1] - coarse[-1])),
            "aggregate_tv_error": 0.5 * float(np.abs(truth - coarse).sum()),
            "posterior_tv_bound": min(1.0, posterior),
            "uniform_tv_bound": min(1.0, uniform),
            "one_step_tv_defect": epsilon,
        })
    return records


def cdf_bracket(steps, p, labels, initial_macro, valid_states):
    """Bellman envelope over unresolved cells, including a singleton cemetery."""
    lower = np.zeros(p.shape[0])
    upper = np.zeros(p.shape[0])
    lower[-1] = upper[-1] = 1.0
    for n in range(len(steps) - 1, -1, -1):
        aggregate_columns = p @ steps[n]
        lower_cells = lower @ aggregate_columns
        upper_cells = upper @ aggregate_columns
        next_lower = lower.copy()
        next_upper = upper.copy()
        for group in range(p.shape[0]):
            cells = (labels == group) & valid_states[n]
            if not np.any(cells):
                # Its true input mass is zero. The post-mask transition also
                # gives it zero incoming mass at this node.
                next_lower[group], next_upper[group] = 0.0, 1.0
                continue
            next_lower[group] = np.min(lower_cells[cells])
            next_upper[group] = np.max(upper_cells[cells])
        lower, upper = next_lower, next_upper
    return float(lower @ initial_macro), float(upper @ initial_macro)


def closure_witness(step, p, labels, valid):
    """Find two admissible microstates in one block with distinct next laws."""
    columns = p @ step
    best = {"diameter_tv": 0.0}
    for group in range(p.shape[0] - 1):
        cells = np.flatnonzero((labels == group) & valid)
        for position, left in enumerate(cells):
            for right in cells[position + 1:]:
                gap = 0.5 * float(np.abs(columns[:, left] - columns[:, right]).sum())
                if gap > best["diameter_tv"]:
                    best = {
                        "group": int(group), "left": int(left), "right": int(right),
                        "diameter_tv": gap,
                        "any_memoryless_column_minimax_error_lower_bound": gap / 2,
                        "left_next_macro": columns[:, left].tolist(),
                        "right_next_macro": columns[:, right].tolist(),
                    }
    return best
