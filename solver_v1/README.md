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

The initial finite-temperature distribution is to be obtained from the correlated interaction energy (conditional Gibbs measure in the declared intact initial basin), not from an imposed Gaussian spacing law or a product closure.

## Numerical implementations in this package

### `probability_pde_2d.py` -- deterministic probability reference

This is the new `N=1`, two-coordinate `(a,s)` gold-standard solver.  It evolves the density directly with a conservative Scharfetter--Gummel finite-volume discretization of the Smoluchowski flux.  It contains no RNG or trajectory counting.

Its purpose is to validate:

- Gibbs initial normalization;
- conservative probability transport;
- positivity/CFL behaviour;
- compression versus tensile-opening sign handling;
- absorbing first passage;
- survival monotonicity;
- grid/time-step convergence.

After the `N=1` reference is validated, an `N=2` dense correlated reference will be used to validate the eventual compressed `N=3` six-dimensional solver.

See `PROBABILITY_PDE_ROADMAP.md` for the sparse-grid / tensor-train development plan.

### `solver.py` -- Euler--Maruyama reference implementation

`solver.py` integrates stochastic trajectories of the same correlated state.  It is retained as a **reference/cross-validation mechanism implementation**, not as the canonical production probability estimator.  Its finite ensemble first-passage fractions must not be presented as the final continuum probability law.

### `model.py` -- shared interaction and opening mechanics

`model.py` contains the two-row LJ geometry, correlated interaction energy, macroscopic strain bridge, periodic configurational wells, and the local normal-opening saddle/barrier lookup used by the numerical solvers.

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

All current LJ parameters, mobilities, thermal scale, and axial projection coefficients remain dimensionless mechanism-screening quantities.  Neither the stochastic reference solver nor the probability-PDE development branch is a calibrated pure-Al fatigue-life predictor.

Quantitative aluminum prediction still requires:

- an Al-specific EAM/MEAM or validated energy landscape;
- mobility/time-scale calibration;
- a physically derived axial configurational bridge;
- characteristic correlation length/area for specimen-scale aggregation;
- experimental validation.

## Existing stochastic screening demo

The historical mechanism demo remains available:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m solver_v1.run_demo
```

It is useful only as a reference mechanism check while the direct probability-PDE solver is developed.
