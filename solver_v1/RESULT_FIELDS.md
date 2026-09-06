# Solver result-field meanings

This file defines the user-facing meaning of the probability-PDE outputs. The UI should use these meanings directly and should not reinterpret the fields as Monte Carlo statistics.

## Load calibration used before characteristic-scale work

The current probability solver deliberately does **not** introduce a characteristic length, area, or volume. The user-entered signed normal stress is first reduced by the user/material Young modulus,

\[
q(t)=\frac{\sigma(t)}{E},
\]

but `q` is not the LJ generalized force itself. Define the extension gradient and pristine local Hessian by

\[
\mathbf c=\begin{bmatrix}1\\ \chi\end{bmatrix},\qquad
H_0=\begin{bmatrix}W_{aa}&W_{as}\\W_{as}&W_{ss}\end{bmatrix}_{(a_0,0)}.
\]

The relaxed total-axial-strain calibration is

\[
\boxed{f^*(t)=\kappa_{\rm axial}q(t)},\qquad
\boxed{\kappa_{\rm axial}=\frac{a_0}{\mathbf c^T H_0^{-1}\mathbf c}}.
\]

This is a tangent calibration of the force coordinate only. It guarantees that the **relaxed total axial strain** reproduces \(\varepsilon\simeq\sigma/E\) at infinitesimal signed load while the actual PDE continues to use the full nonlinear Bessel-LJ energy. The frozen-registry scale \(a_0W_{aa}(a_0,0)\) remains available as a diagnostic and is not the canonical load mapping.

The implementation is `TwoRowLJ.force_from_sigma_over_E(...)` and `cyclic_load_from_sigma_over_E(...)`.

## Core probability fields

### `time`
Dimensionless model-time coordinate used by the probability PDE. It is not yet calibrated to laboratory seconds for pure Al.

### `force`
Current dimensionless generalized normal load \(f^*\) passed into the effective energy landscape. It should be produced from \(\sigma/E\) through the tangent conversion above, not by setting `force = sigma/E` directly.

### `survival`
Current intact probability mass,

\[
S(t)=\int_{\Omega_b(t)}P_N(\mathbf q,t)\,d\mathbf q.
\]

For a correct absorbing first-passage calculation, physical survival is non-increasing. The UI should retain enough numerical precision to reveal tiny changes near one.

### `initiation_probability`
Cumulative crack-initiation probability,

\[
P_{\rm init}(t)=1-S(t),
\]

or, in the preferred production accounting, the cumulative probability mass absorbed through the crack-opening boundary. It is not a fraction of randomly sampled trajectories.

### `first_passage_flux`
Instantaneous rate at which intact probability mass crosses the crack-opening dividing surface. For a fixed boundary,

\[
-\dot S(t)=\int_{\Gamma_c}\mathbf J\cdot\mathbf n\,dS.
\]

This can vary strongly within each load cycle even though survival itself should not increase.

## Mechanical/configurational fields

The unwrapped configurational coordinate is decomposed exactly as

\[
s=bn+\xi,\qquad n\in\mathbb Z,\quad -b/2\le\xi<b/2.
\]

For N=1 the survivor-conditioned strain decomposition is

\[
\boxed{\varepsilon=\varepsilon_a+\varepsilon_\xi+\varepsilon_p}
\]

with

\[
\varepsilon_a=\left\langle\frac{a-a_0}{a_0}\right\rangle_S,
\qquad
\varepsilon_\xi=\left\langle\frac{\chi\xi}{a_0}\right\rangle_S,
\qquad
\varepsilon_p=\left\langle\frac{\chi b n}{a_0}\right\rangle_S.
\]

### `strain`
Total survivor-conditioned axial strain \(\varepsilon\). This is the sum of the three fields below.

### `normal_strain`
Fast normal-opening contribution \(\varepsilon_a\). The entered Young modulus calibrates the relaxed total strain, not this component in isolation.

### `intrawell_strain`
Reversible/anelastic registry contribution \(\varepsilon_\xi\) from motion inside the current configurational well.

### `plastic_strain`
Signed permanent well-index contribution \(\varepsilon_p\). This is computed from the PDE probability carried by nonzero well indices. It must **not** be hard-coded to zero. It may nevertheless remain zero or extremely small for a particular load if the solved PDE carries no appreciable probability across a well boundary.

### `mean_well_index`
Survivor-conditioned signed mean \(\langle n\rangle_S\). Its sign distinguishes positive and negative accumulated configurational registry changes.

### `plastic_well_activity`
Survivor-conditioned \(\langle|n|\rangle_S\). This is a nonnegative diagnostic of inter-well rearrangement. It is not the same quantity as signed plastic strain and is not an empirical damage variable.

### `cov_a12`
For the N=2 dense reference solver, covariance between the two normal-spacing coordinates `a1` and `a2` among survivors. Nonzero values indicate explicit cross-cell correlation.

### `cov_s12`
For the N=2 dense reference solver, covariance between configurational coordinates `s1` and `s2` among survivors.

### `product_closure_l1_error`
For the N=2 dense reference solver,

\[
\int |P_2-P_1^{(1)}P_1^{(2)}|\,d\mathbf q.
\]

A nonzero value measures how strongly the actual correlated joint density differs from an independent product closure.

## Numerical-convergence diagnostics

### `s_truncation_boundary_mass`
Probability mass close to the artificial numerical truncation boundary in the `s` direction. It is a convergence warning, not a physical observable. If it becomes appreciable, the state-space domain must be enlarged.

### `a_lower_boundary_mass`
Probability mass near the lower numerical `a` boundary. This should remain negligible for a well-resolved calculation.

### `a_upper_boundary_mass`
Probability mass near the upper numerical `a` boundary. This is separate from the mechanically defined crack-opening dividing surface and should remain negligible unless the computational domain is too small.

## Probability bookkeeping diagnostics

The N=1 PDE and desktop adapter expose:

- `intact_probability_mass`: direct numerical integral of the current intact density;
- `cumulative_absorbed_mass`: accumulated first-passage mass through the crack boundary;
- `mass_balance_residual = intact_probability_mass + cumulative_absorbed_mass - 1`;
- `negative_mass_correction`: largest roundoff-level negative mass removed so far;
- `minimum_density`: smallest uncorrected density encountered so positivity loss is visible.

These fields are necessary when displaying rare probabilities near machine/numerical tolerance, because they distinguish physical first-passage probability from numerical mass drift.

## UI plotting behaviour

Plots should preserve scientific-axis offset notation rather than forcing survival to the full 0--1 range. The user should be able to zoom and pan interactively:

- mouse wheel: zoom about the cursor;
- horizontal/vertical constrained zoom where supported;
- drag: pan;
- double click or Home: autoscale/reset.

Each result field should remember its own last view range when the user switches between fields.
