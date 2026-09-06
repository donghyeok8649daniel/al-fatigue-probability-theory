# Solver result-field meanings

This file defines the user-facing meaning of the probability-PDE outputs. The UI should use these meanings directly and should not reinterpret the fields as Monte Carlo statistics.

## Core probability fields

### `time`
Physical/model time coordinate used by the probability PDE. In the present dimensionless mechanism solver this is not yet calibrated to seconds for pure Al.

### `force`
Current dimensionless normal loading value passed into the effective free-energy landscape. In a calibrated model this corresponds to the stress/work term used in `G_N`.

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

### `strain`
Survivor-conditioned macroscopic axial strain obtained by averaging the microscopic strain bridge over the intact probability density. For the current model,

\[
\varepsilon(\mathbf q)=\frac{1}{Na_0}\sum_i[(a_i-a_0)+\chi s_i].
\]

### `plastic_well_activity`
For the N=1 reference solver, the survivor-conditioned activity of the configurational well index associated with

\[
s=bn+\xi.
\]

This is a diagnostic of inter-well configurational rearrangement; it is not a scalar empirical fatigue-damage variable.

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

## Recommended production-only diagnostics

The UI integration should additionally expose the following once the probability bookkeeping is upgraded:

- `raw_intact_mass`: direct numerical integral of the current intact density;
- `cumulative_absorbed_mass`: accumulated first-passage mass through the crack boundary;
- `mass_balance_residual = raw_intact_mass + cumulative_absorbed_mass - 1`.

These fields are necessary when displaying rare probabilities near machine/numerical tolerance, because they distinguish physical first-passage probability from numerical mass drift.

## UI plotting behaviour

Plots should preserve scientific-axis offset notation rather than forcing survival to the full 0--1 range. The user should be able to zoom and pan interactively:

- mouse wheel: zoom about the cursor;
- horizontal/vertical constrained zoom where supported;
- drag: pan;
- double click or Home: autoscale/reset.

Each result field should remember its own last view range when the user switches between fields.
