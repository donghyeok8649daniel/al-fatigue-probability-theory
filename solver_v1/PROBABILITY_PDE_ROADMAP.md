# Deterministic probability-PDE solver roadmap

## Purpose

The production fatigue-probability solver is not defined by Monte Carlo trajectory counts.  Its state is the probability density itself,

\[
P_N(\mathbf q,t),\qquad
\mathbf q=(a_1,\ldots,a_N,s_1,\ldots,s_N),
\]

and the governing equation is the correlated many-body Smoluchowski equation

\[
\partial_t P_N
=\nabla_{\mathbf q}\!\cdot\!\left[
\mathbf M\left(P_N\nabla_{\mathbf q}\mathcal G_N
+k_B T\nabla_{\mathbf q}P_N\right)\right].
\]

Crack initiation is probability mass that leaves the intact basin through the mechanically defined opening dividing surface.  No random-number sampling is required to define or solve this PDE.

The existing Euler--Maruyama ensemble remains a reference/cross-validation implementation only.

## Why a dense six-dimensional grid is not acceptable

For the minimum correlated proof case, `N=3`, the state has six coordinates:

\[
(a_1,a_2,a_3,s_1,s_2,s_3).
\]

A tensor-product grid with `m` points per coordinate stores `m^6` density values.  For example:

- `m=21`: 85,766,121 cells;
- `m=31`: 887,503,681 cells;
- `m=41`: 4,750,104,241 cells.

One double-precision scalar field at `m=41` already requires about 38 GB before storing energy, fluxes, masks, temporary arrays, or linear-system workspaces.  Therefore the production `N=3` solver must compress the probability representation.

## What a sparse grid means

A sparse grid does **not** change the physics.  It changes where the PDE is represented numerically.

A full tensor grid refines every combination of all coordinates.  A hierarchical sparse grid refines mainly those multi-dimensional basis functions that contribute appreciably.  For sufficiently smooth functions, this can reduce the number of active degrees of freedom by orders of magnitude.

Advantages:

- deterministic, no Monte Carlo resolution floor;
- adaptive refinement can target the Gibbs core, configurational-well transitions, and the crack dividing surface;
- natural for quadrature and interpolation.

Risks for this project:

- naive sparse interpolation can produce negative density;
- strict mass conservation is harder than on a finite-volume grid;
- a moving absorbing crack surface is awkward on a generic sparse collocation grid.

For that reason a sparse grid should not be introduced before the probability-flux and first-passage implementation is validated in lower dimension.

## What a tensor / tensor-train representation means

On a structured grid the density is a `2N`-way array,

\[
P[i_1,\ldots,i_{2N}].
\]

A tensor-train (TT) approximation writes this high-dimensional array as a product of low-rank cores,

\[
P[i_1,\ldots,i_d]
\approx
G_1[i_1]G_2[i_2]\cdots G_d[i_d],\qquad d=2N.
\]

Instead of storing `m^d` numbers, storage is roughly `O(d m r^2)` when the TT rank is `r`.  Correlation is **not** removed: the TT ranks carry correlation between coordinate groups.  A rank-one tensor would correspond to a product closure, but the solver must allow ranks to grow and must never impose rank one as a physical assumption.

Advantages:

- keeps a structured coordinate grid, useful for conservative derivative operators;
- can represent a genuinely correlated `P_N` without storing the full tensor;
- scales much better to the six-dimensional `N=3` state when ranks remain moderate.

Risks:

- low-rank truncation is a numerical approximation that needs convergence tests;
- positivity is not automatic;
- the interacting drift field and moving absorbing boundary can increase ranks.

## Development order

The production implementation will be built in three validation layers.

### Layer A -- 2D finite-volume gold standard (`N=1`)

State: `(a,s)`.

Use a conservative Scharfetter--Gummel finite-volume discretization of the Smoluchowski flux.  This gives a deterministic probability calculation with no RNG and treats drift and thermal diffusion in a thermodynamically consistent flux form.

Required checks:

1. Gibbs initial density normalizes to one;
2. probability mass is conserved when no absorbing opening boundary is reached;
3. survival is monotone non-increasing when absorption is active;
4. compression does not spuriously create tensile-opening first passage;
5. athermal opening loss appears when the mechanically defined opening basin disappears;
6. grid/time-step convergence is documented.

### Layer B -- 4D correlated validation (`N=2`)

Use a small dense grid only for convergence/reference cases.  This is the first test that cell-cell correlation is represented explicitly.

### Layer C -- 6D production representation (`N=3`)

Compare two compressed representations against Layers A/B:

- adaptive hierarchical sparse-grid finite volume / Galerkin;
- tensor-train structured-grid representation.

The selected method must pass mass, positivity, equilibrium, first-passage, and reduced-dimension convergence checks before it is connected to the UI as the production probability solver.

## Initial probability

The default finite-temperature initial state is a conditional correlated Gibbs density in the intact, principal configurational basin,

\[
P_N(\mathbf q,0)
=Z_0^{-1}
\exp\!\left[-\beta\left(\mathcal G_N(\mathbf q;\sigma_{\rm pre})-\mathcal G_{\min}\right)\right]
\mathbf 1_{\Omega_{\rm intact}}
\prod_i\mathbf 1_{[-b/2,b/2)}(s_i).
\]

`G_min` is only a numerical energy offset.  It does not change the probability law.

The canonical baseline is `sigma_pre = 0` unless a physical preload protocol is explicitly declared.

## First passage

The intact survival mass is

\[
S(t)=\int_{\Omega_b(t)}P_N(\mathbf q,t)\,d\mathbf q,
\qquad
P_{\rm init}(t)=1-S(t).
\]

For a fixed dividing surface,

\[
-\dot S(t)=\int_{\Gamma_c}\mathbf J\cdot\mathbf n\,dS.
\]

Numerically, mass that crosses the opening dividing surface is absorbed and accumulated as initiated probability.  Failed mass is never renormalized back into the intact density.

## Implementation status on `probability-pde-solver-v1`

### Implemented

- `probability_pde_2d.py`: N=1 direct deterministic Smoluchowski PDE reference using conservative Scharfetter--Gummel fluxes.
- `test_probability_pde_2d.py`: Gibbs normalization, zero-load mass conservation, compression/no-tensile-opening check, spinodal absorption, and monotone survival checks.
- `probability_pde_4d.py`: N=2 dense correlated Smoluchowski PDE reference on `(a1,a2,s1,s2)` without product closure.
- `test_probability_pde_4d.py`: N=2 normalization, explicit non-product correlation diagnostic, mass conservation, compression check, and monotone survival checks.
- `tensor_train.py`: rank-adaptive TT-SVD compression utility; rank one is not imposed.
- `probability_tt_6d.py`: N=3 six-dimensional correlated Gibbs initial-density construction and TT compression prototype on small verification grids.
- `test_tensor_train.py` and `test_probability_tt_6d.py`: reconstruction, nontrivial rank, mass, error, and positivity diagnostics.

### Not yet implemented

- TT time integration of the six-dimensional Smoluchowski operator;
- positivity-preserving TT truncation / correction;
- TT representation of the moving absorbing opening boundary and first-passage flux;
- sparse-grid competitor;
- systematic grid/time/rank convergence tables;
- calibrated aluminum energy/mobility data.

The N=2 dense solver and N=3 TT prototype are reference-development code.  They must not be connected to the UI as production physics until the lower-dimensional tests pass locally and the compressed time integrator reproduces those references.

## Scientific status

This branch is numerical-method development.  The current LJ parameters and mobilities remain dimensionless mechanism parameters.  The compressed high-dimensional solver must not be described as a calibrated pure-Al lifetime predictor until Al-specific energy/mobility calibration and experimental validation are completed.
