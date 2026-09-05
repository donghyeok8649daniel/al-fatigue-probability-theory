# Theory Core v1 solver

This package contains numerical implementations of the paper's correlated configurational fatigue-initiation framework.

## Canonical production direction

The production fatigue-probability solver is defined by the probability density itself,

\[
P_N(\mathbf q,t),\qquad
\mathbf q=(a_1,\ldots,a_N,s_1,\ldots,s_N),
\]

and by direct deterministic solution of the many-body Smoluchowski equation,

\[
\partial_t P_N
=\nabla_{\mathbf q}\cdot\left[
\mathbf M\left(P_N\nabla_{\mathbf q}\mathcal G_N
+k_BT\nabla_{\mathbf q}P_N\right)\right].
\]

Crack initiation is probability mass absorbed through the mechanically defined opening dividing surface.  The production probability is **not** defined as a Monte Carlo fraction.

The initial finite-temperature distribution is obtained from the correlated interaction energy as a conditional Gibbs measure in the declared intact initial basin, not from an imposed Gaussian spacing law or a product closure.

## Numerical implementations in this package

### `probability_pde_2d.py` -- N=1 deterministic gold standard

Directly evolves `P(a,s,t)` with conservative Scharfetter--Gummel finite-volume fluxes.  No RNG or trajectory counting is used.

Primary checks are Gibbs normalization, probability conservation, positivity/CFL behaviour, compression sign handling, absorbing first passage, survival monotonicity, and grid/time-step convergence.

### `probability_pde_4d.py` -- N=2 dense correlated reference

Directly evolves

\[
P_2(a_1,a_2,s_1,s_2,t)
\]

on a deliberately small four-dimensional tensor grid.  The energy is the interacting two-cell energy from `model.py`; no product closure is imposed.  The solver records cross-cell covariances and an explicit L1 discrepancy between the full joint density and the product of its one-cell marginals.

This solver is a convergence/reference tool, not a scalable production implementation.

### `tensor_train.py` -- numerical compression utility

Implements TT-SVD and reconstruction diagnostics.  Tensor rank is allowed to exceed one.  A rank-one product state is never imposed as a physical assumption; higher TT ranks carry cross-coordinate correlation.

### `probability_tt_6d.py` -- N=3 TT initial-state prototype

Constructs the full correlated six-dimensional Gibbs initial density on a small verification grid,

\[
P_3(a_1,a_2,a_3,s_1,s_2,s_3,0),
\]

then compresses it with TT-SVD and reports ranks, storage, compression ratio, reconstruction error, mass error, and negative reconstructed mass.

This is the first Layer-C prototype.  It does **not yet** time-integrate the six-dimensional Smoluchowski equation.

### `solver.py` -- Euler--Maruyama reference implementation

Integrates stochastic trajectories of the correlated state.  It is retained only as a **reference/cross-validation mechanism implementation**, not as the canonical production probability estimator.  Finite-ensemble first-passage fractions must not be presented as the final continuum probability law.

### `model.py` -- shared interaction and opening mechanics

Contains the two-row LJ geometry, correlated interaction energy, macroscopic strain bridge, periodic configurational wells, and the local normal-opening saddle/barrier lookup used by all current numerical solvers.

## Development sequence

The numerical hierarchy is now explicit:

1. `N=1`: 2D direct probability PDE gold standard;
2. `N=2`: 4D dense correlated probability reference;
3. `N=3`: 6D tensor-train / sparse-grid compressed production development.

The next production step is a tensor-train time integrator for the six-dimensional Smoluchowski operator, with mass, positivity, equilibrium and first-passage behaviour checked against the N=1/N=2 reference solvers before any UI integration.

See `PROBABILITY_PDE_ROADMAP.md` for the detailed validation plan.

## Why the N=3 PDE needs compression

For `N=3`, the density depends on six coordinates:

\[
(a_1,a_2,a_3,s_1,s_2,s_3).
\]

A full grid with `m` points per coordinate stores `m^6` values.  At `m=41`, one scalar field already contains more than 4.75 billion doubles, so a dense six-dimensional finite-volume grid is not practical.

The production `N=3` solver will therefore use a validated compressed representation, with adaptive sparse-grid and tensor-train approaches compared against lower-dimensional gold standards.  Compression is numerical only: it must not impose the physical product closure

\[
P_N=\prod_i P_i.
\]

## Scientific scope

All current LJ parameters, mobilities, thermal scale, and axial projection coefficients remain dimensionless mechanism-screening quantities.  None of these numerical implementations is yet a calibrated pure-Al fatigue-life predictor.

Quantitative aluminum prediction still requires:

- an Al-specific EAM/MEAM or validated energy landscape;
- mobility/time-scale calibration;
- a physically derived axial configurational bridge;
- characteristic correlation length/area for specimen-scale aggregation;
- experimental validation.

## Historical stochastic screening demo

The historical mechanism demo remains available:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m solver_v1.run_demo
```

It is useful only as a reference mechanism check while the direct probability-PDE solver is developed.
