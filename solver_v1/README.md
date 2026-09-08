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

Crack initiation is probability mass absorbed through the mechanically defined opening dividing surface. The production probability is **not** defined as a Monte Carlo fraction.

The initial finite-temperature distribution is obtained from the correlated interaction energy as a conditional Gibbs measure in the declared intact initial basin, not from an imposed Gaussian spacing law or a product closure.

## Canonical infinite-lattice LJ kernel

The lower row is now represented analytically by the Poisson-summed infinite-lattice kernel rather than by a fixed finite image cutoff. For one upper cell,

\[
W(a,s)=4\epsilon\left[\sigma^{12}S_6(a,s)-\sigma^6S_3(a,s)\right],
\]

with

\[
S_p(a,s)=\sum_{n=-\infty}^{\infty}[a^2+((n+1/2)b-s)^2]^{-p}.
\]

Poisson summation gives a reciprocal-lattice cosine series whose coefficients contain modified Bessel functions `K_{p-1/2}`. The energy and its exact first derivatives with respect to `a` and `s` are implemented in `lattice_bessel.py`. These derivatives generate the conservative drift term in the probability PDE.

`model.py` uses this Bessel kernel as the canonical lower-row interaction. The historical `lower_images` parameter remains only for API compatibility and no longer controls production lower-row physics. A large finite direct sum is retained only as a verification reference.

See `BESSEL_LATTICE_DERIVATION.md` for the derivation and phase convention.

## Numerical implementations in this package

### `lattice_bessel.py` -- canonical infinite lower-row interaction

Evaluates the exact Poisson/Bessel representation of the staggered infinite lower lattice with tolerance-controlled reciprocal-space truncation. It returns `W`, `dW/da`, and `dW/ds`. The Bessel representation is an energy kernel, not a separate empirical dynamics law.

### `probability_pde_2d.py` -- N=1 deterministic gold standard

Directly evolves `P(a,s,t)` with conservative Scharfetter--Gummel finite-volume fluxes. No RNG or trajectory counting is used.

Primary checks are Gibbs normalization, probability conservation, positivity/CFL behaviour, compression sign handling, absorbing first passage, survival monotonicity, and grid/time-step convergence.

An optional backward-Euler mode uses the same conservative
Scharfetter--Gummel generator for stiff timestep and fast/slow diagnostics; the
validated default remains explicit. Linear relaxation, harmonic response,
stable $a^*(s,f)$, and finite-temperature conditional fast-$a$ references are
implemented in `dynamics_diagnostics.py` and described in
[`DYNAMICS_FAST_A_DIAGNOSTICS.md`](DYNAMICS_FAST_A_DIAGNOSTICS.md).

The finite-temperature elimination derivation, moving-boundary correction,
truncated-Gibbs/QSD distinction, and experimental slow-$s$ solver are in
[`FAST_A_REDUCTION.md`](FAST_A_REDUCTION.md) and `reduced_fast_a.py`. The full
2D solver remains the reference, and the reduced module does not provide a
validated crack probability.

Local absorbed-mass bookkeeping, configurational well fluxes, convergence
floors, and the strictly post-processing statistical correlation-area layer
are documented in
[`PROBABILITY_AND_PLASTICITY_RESOLUTION.md`](PROBABILITY_AND_PLASTICITY_RESOLUTION.md).

The analytic LJ--EAM construction and its pure-Al target/identifiability audit
are documented separately in
[`ANALYTIC_LJ_EAM_HYBRID.md`](ANALYTIC_LJ_EAM_HYBRID.md) and
[`ALUMINUM_ANALYTIC_CALIBRATION.md`](ALUMINUM_ANALYTIC_CALIBRATION.md). The
result is explicitly best-feasible and practically non-identifiable; it is not
labeled a calibrated aluminum potential.

The exact mobility/time nondimensionalization, uncalibrated kinetic-data
status, calibration routes, and physical-time UI gate are documented in
[`PHYSICAL_TIME_AND_MOBILITY.md`](PHYSICAL_TIME_AND_MOBILITY.md).

### `probability_pde_4d.py` -- N=2 dense correlated reference

Directly evolves

\[
P_2(a_1,a_2,s_1,s_2,t)
\]

on a deliberately small four-dimensional tensor grid. The energy is the interacting two-cell energy from `model.py`; no product closure is imposed. The solver records cross-cell covariances and an explicit L1 discrepancy between the full joint density and the product of its one-cell marginals.

This solver is a convergence/reference tool, not a scalable production implementation.

### `tensor_train.py` -- numerical compression utility

Implements TT-SVD and reconstruction diagnostics. Tensor rank is allowed to exceed one. A rank-one product state is never imposed as a physical assumption; higher TT ranks carry cross-coordinate correlation.

### `probability_tt_6d.py` -- N=3 TT initial-state prototype

Constructs the full correlated six-dimensional Gibbs initial density on a small verification grid,

\[
P_3(a_1,a_2,a_3,s_1,s_2,s_3,0),
\]

then compresses it with TT-SVD and reports ranks, storage, compression ratio, reconstruction error, mass error, and negative reconstructed mass.

This is the first Layer-C prototype. It does **not yet** time-integrate the six-dimensional Smoluchowski equation.

### `solver.py` -- Euler--Maruyama reference implementation

Integrates stochastic trajectories of the same correlated state. It is retained only as a historical/reference cross-check and is **not** the canonical production probability estimator. Finite-ensemble first-passage fractions must not be presented as the final continuum probability law.

### `model.py` -- shared interaction and opening mechanics

Contains the two-row geometry, analytic infinite-lattice lower-row interaction, explicit upper-cell correlations, macroscopic strain bridge, periodic configurational wells, and the local normal-opening saddle/barrier lookup used by all current numerical solvers.

## Development sequence

The numerical hierarchy is explicit:

1. `N=1`: 2D direct probability PDE gold standard;
2. `N=2`: 4D dense correlated probability reference;
3. `N=3`: 6D tensor-train / sparse-grid compressed production development.

The next production step is a tensor-train time integrator for the six-dimensional Smoluchowski operator, with mass, positivity, equilibrium and first-passage behaviour checked against the N=1/N=2 reference solvers before final UI integration.

See `PROBABILITY_PDE_ROADMAP.md` for the detailed validation plan and `RESULT_FIELDS.md` for user-facing output meanings.

See `CONFIGURATIONAL_PLASTICITY.md` for the exact registry-well population and
flux balances, the separation of interwell plastic flow from selective crack
absorption, the candidate crystallographic interpretation of $\chi$, and the
current (conservative) numerical classification of the LJ slip mechanism.

See `PLASTICITY_AND_STRAIN_AUDIT.md` for the production survivor-normalized
strain audit, deterministic grid/time/interface refinement, unload/zero-stress
hold study, barrier comparison, exact desktop energy-model path, and the
separate crack/plasticity resolution classifications.

See `ANALYTIC_LJ_EAM_HYBRID.md` for the optional analytic LJ--EAM hybrid. The
verified Poisson/Bessel LJ pair term remains the base interaction; an
exponential environment-density lattice sum and differentiable embedding term
are added in a separate `AnalyticLJEAM` model. Its default parameter set is a
mathematical sensitivity example, not calibrated aluminum data.

See `FCC111_FULL_STACK_DERIVATION.md` for the separate full FCC(111)
infinite-plane-stack static reference. It derives the triangular reciprocal
lattice, exact ABC phases, two-dimensional Poisson/Bessel plane kernels,
zeta(4)/zeta(10) layer means, and the rule that the embedding function is
applied once after summing the complete environment density. This research
model does not replace the reduced production energy surfaces.

See `ALUMINUM_FULL_FCC_CALIBRATION.md` and
`FCC111_ACTIVE_INTERFACE_DERIVATION.md` for the deterministic bulk refit and
the localized half-crystal opening/slip construction. The independent-mode
audit rejects the initial six-observation rank-five claim and adds the missing
homogeneous shear information. A new linear-embedding candidate fits the bulk
targets with strong pair/embedding cancellation, but its interface is not
physically validated. Another locally stable candidate has negative separation
work. Neither is connected to the production probability PDE.

Run the actual optimization with `python -m solver_v1.run_full_fcc_calibration_audit`.
Current tables are in the `audited_v2/` result directories; the earlier tables
are explicitly superseded historical evidence. `SOLVER_VALIDATION_GATES.md`
records the required physical checks before any specimen UI or PDE promotion.

## Why the N=3 PDE needs compression

For `N=3`, the density depends on six coordinates:

\[
(a_1,a_2,a_3,s_1,s_2,s_3).
\]

A full grid with `m` points per coordinate stores `m^6` values. At `m=41`, one scalar field already contains more than 4.75 billion doubles, so a dense six-dimensional finite-volume grid is not practical.

The production `N=3` solver will therefore use a validated compressed representation, with adaptive sparse-grid and tensor-train approaches compared against lower-dimensional gold standards. Compression is numerical only: it must not impose the physical product closure

\[
P_N=\prod_i P_i.
\]

## Scientific scope

All current LJ parameters, mobilities, thermal scale, and axial projection coefficients remain dimensionless mechanism-screening quantities. None of these numerical implementations is yet a calibrated pure-Al fatigue-life predictor.

Quantitative aluminum prediction still requires:

- an Al-specific EAM/MEAM or validated energy landscape;
- mobility/time-scale calibration;
- a physically derived axial configurational bridge;
- externally calibrated statistical correlation area for specimen-scale aggregation;
- experimental validation.

## Historical stochastic screening demo

The historical mechanism demo remains available:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m solver_v1.run_demo
```

It is useful only as a reference mechanism check while the direct probability-PDE solver is developed.
