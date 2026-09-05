# Deterministic normal Theory Core solver

The canonical desktop path uses `deterministic_normal.py`, which implements the active finite generalized-LJ normal chain and deterministic initial-measure push-forward.

## What it solves

The active state is the finite set of normal spacings and spacing rates

$$
(\lambda_1,\ldots,\lambda_M,c_1,\ldots,c_M).
$$

It integrates the conservative spacing equations from the active theory under $q(\tau)=\sigma_n(t)/E$. A weighted discrete initial measure is pushed forward without thermal noise, mobility, diffusion, resampling, or a named probability distribution. The default is the single ideal state $\mu_0=\delta_{(\lambda=1,c=0)}$.

The local probability is the exact $1/M$ spatial counting measure, extended by the declared weights when a nontrivial discrete $\mu_0$ is supplied. First passage occurs when a spacing reaches

$$
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)},
$$

where $\phi''(\lambda_c)=0$.

## Important scope

This is a **dimensionless conservative reference solver**, not a calibrated aluminum fatigue-life predictor.

- The generalized-LJ exponents are $m=12.19$ and $n=6$.
- A physical specimen initial measure $\mu_0$ and its correlation scale remain open.
- Registry/slip $s$, residual plasticity, and irreversible $G3$ are not active without a consistently derived inertia/evolution law and irreversible microscopic mechanism.
- Specimen-scale characteristic length/area aggregation is deliberately deferred.
- The atomic-to-laboratory time bridge remains open.

## Run

From the repository root:

```bash
python -m pytest solver_v1/test_solver.py
```

Outputs are written to `solver_v1/output/`.

## Theory-to-code map

- `deterministic_normal.py`: canonical finite-chain mechanics, exact delta-supported phase-space measure, discrete-$\mu_0$ push-forward, and first passage.
- `model.py`, `solver.py`, and `run_demo.py`: historical two-row stochastic v4 screening implementation. They are retained for provenance but are not imported by the canonical desktop solve.

## Current screening parameters

- $M$ is the user-selected number of chain spacings.
- $m=12.19$, $n=6$.
- $\mu_0=\delta_{(\lambda=1,c=0)}$ by default.
- No stochastic trajectory count exists.

The normalized interaction remains a model assumption; quantitative pure-Al fatigue calibration is not claimed.
