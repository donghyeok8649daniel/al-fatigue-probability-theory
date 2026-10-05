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

The implementation uses the selected energy surface's
`force_from_sigma_over_E(...)` together with
`cyclic_load_from_sigma_over_E(...)`. Each surface recomputes its own $a_0$,
Hessian, and relaxed axial kappa. Result metadata records the exact model ID,
Python class, parameter source, calibration status, and kappa; a `TwoRowLJ`
result must not be labeled as the Al-target hybrid.

## Core probability fields

### `time`
Dimensionless model-time coordinate used by the probability PDE. It is not yet calibrated to laboratory seconds for pure Al.

### `force`
Current dimensionless generalized normal load \(f^*\) passed into the effective energy landscape. It should be produced from \(\sigma/E\) through the tangent conversion above, not by setting `force = sigma/E` directly.

### `survival`
Physical local survival defined from accumulated opening absorption,

$$
S_{\mathrm{local}}(t)=1-A_{\mathrm{abs}}(t).
$$

It is non-increasing by construction. The independently integrated intact
density is `intact_probability_mass`, not this physical bookkeeping field.

### `initiation_probability`
Cumulative crack-initiation probability,

$$
P_{\mathrm{init,local}}(t)=A_{\mathrm{abs}}(t)
=1-S_{\mathrm{local}}(t).
$$

Only mass removed by the opening absorbing mask enters $A_{\mathrm{abs}}$.
Numerical positivity repair and conservative-solve roundoff do not. It is not
a fraction of randomly sampled trajectories.

### `first_passage_flux`
Record-interval average rate at which probability is removed at the opening
dividing surface. For records $t_{k-1}<t_k$,

$$
\Phi_k=\frac{A_{\mathrm{abs}}(t_k)-A_{\mathrm{abs}}(t_{k-1})}
{t_k-t_{k-1}}.
$$

The initial mask removal is reported separately and is not assigned an
artificial finite rate. `integrated_first_passage_flux` accumulates later
discrete absorption, while `flux_consistency_residual` exposes their mismatch.

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

### `well_populations`
Absolute survivor masses $P_n(t)$ in each configurational partition
$s\in[(n-1/2)b,(n+1/2)b)$. At least $n=-1,0,+1$ are retained on the default
three-well domain.

### `interwell_net_flux`
Signed finite-volume configurational flux through each represented well
boundary. Positive sign points toward increasing $s$.

### `interwell_gross_flux`
Sum of the two nonnegative one-way SG rates at a well boundary. It detects
bidirectional thermal crossing even when signed net transfer is nearly zero.

### `interwell_forward_flux` and `interwell_backward_flux`
The two nonnegative one-way SG interface rates. With positive direction toward
increasing $s$,

$$
J_{+}=\frac{J_{\mathrm{gross}}+J_{\mathrm{net}}}{2},\qquad
J_{-}=\frac{J_{\mathrm{gross}}-J_{\mathrm{net}}}{2}.
$$

They expose whether small net flow comes from genuinely small activity or from
large forward/backward cancellation.

### `well_population_balance_residual`
Difference between directly integrated $P_n(t)$ and the population reconstructed
from cumulative left/right boundary fluxes and opening absorption. This is a
numerical diagnostic, not plastic strain.

### `unnormalized_registry_moment`
$M_n=\sum_n nP_n$, formed from unnormalized intact well masses. This quantity
can change by both interwell transport and selective opening absorption.

### `accumulated_net_registry_transfer`
Time integral of the signed interwell interface fluxes, summed over represented
interfaces. Positive sign means accumulated transport toward increasing $s$.

### `cumulative_interwell_forward_transfer`,
`cumulative_interwell_backward_transfer`, and
`cumulative_interwell_gross_transfer`
Interface-resolved time integrals of the corresponding one-way/gross SG rates.
Their interface sums are exposed to the UI as
`cumulative_forward_registry_activity`,
`cumulative_backward_registry_activity`, and
`cumulative_gross_registry_activity`. These are probability-traffic measures;
they are never multiplied by specimen area and do not replace signed plastic
strain.

### `absorbed_registry_moment`
$\sum_n nA_n(t)$, the well-index moment of probability removed by the opening
boundary. It is a crack-loss term, not plastic flow.

### `registry_moment_balance_residual`
Numerical residual of

$$
M_n(t)-M_n(0)-\sum_n\int_0^t\mathcal J_{n+1/2}\,dt'
+\sum_n nA_n(t)=0.
$$

### `net_interwell_registry_rate`
$\sum_n\mathcal J_{n+1/2}$ for the unnormalized intact density. It measures
directional interwell registry transport.

### `gross_interwell_activity_rate`
Sum of all SG forward-plus-backward one-way crossing rates. It is not
$|\mathcal J_{\rm net}|$ and can be nonzero for symmetric thermal hopping.

### `net_plastic_flow_rate`
The survivor-conditioned axial rate contribution

$$
\frac{\chi b}{a_0S}\sum_n\mathcal J_{n+1/2}.
$$

This is the model plastic-flow term. A numerical nonzero is still subject to
grid, timestep, interface-alignment, and domain convergence tests.

### `selective_opening_plastic_rate`
The change in survivor-conditioned $\epsilon_p$ caused solely by selective
opening absorption,

$$
\frac{\chi b}{a_0S}
\left(-\sum_n n\dot A_n+\langle n\rangle_S\sum_n\dot A_n\right).
$$

It is exposed specifically so it cannot be misclassified as plastic flow.

### `gross_configurational_slip_activity`
Axially scaled survivor-conditioned gross SG interwell activity. It records
forward-plus-backward configurational traffic and is distinct from signed net
plastic flow.

The complete derivation is in `CONFIGURATIONAL_PLASTICITY.md`. The current
survivor-normalization, energy-model-path, refinement, unload/hold, and barrier
audit is in `PLASTICITY_AND_STRAIN_AUDIT.md`.

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
- `absorbed_mass_increment`: opening mass removed since the previous record;
- `initial_absorbed_mass`: mass removed by the opening mask at the initial instant;
- `integrated_first_passage_flux`: accumulated post-initial discrete opening flux;
- `flux_consistency_residual = cumulative_absorbed_mass - initial_absorbed_mass - integrated_first_passage_flux`;
- `mass_balance_residual = intact_probability_mass + cumulative_absorbed_mass - 1`;
- `negative_mass_correction`: largest roundoff-level negative mass removed so far;
- `cumulative_negative_mass_correction`: accumulated numerical positivity repair;
- `minimum_density`: smallest uncorrected density encountered so positivity loss is visible.

These fields are necessary when displaying rare probabilities near machine/numerical tolerance, because they distinguish physical first-passage probability from numerical mass drift.

The local rare-event floor is the maximum observed scale from mass residual,
accumulated numerical repair, flux/absorption discrepancy, and grid/time/
integrator variation. A single PDE run supplies only a lower bound and does not
certify a signal for specimen amplification.

## Specimen aggregation fields

$A_c$ is a statistical correlation area and
$A_{\mathrm{stressed}}$ is the effective specimen surface area under the
modeled loading. They are distinct from the atomic geometry and from the
effective free energy $\mathcal F_{\mathrm{eff}}$. Neither enters the
constitutive force mapping.

For the independent-equivalent-region approximation,

$$
N_{\mathrm{eff}}=\frac{A_{\mathrm{stressed}}}{A_c},\qquad
\log S_{\mathrm{spec}}=N_{\mathrm{eff}}\log(1-P_{\mathrm{init,local}}).
$$

`specimen_initiation_probability` and `specimen_survival_probability` are
exposed as physical fields only when the local signal is convergence-certified.
Otherwise they remain unavailable and `specimen_probability_extrapolation` is
diagnostic only.

## Time-basis metadata

- `model_time`: canonical nondimensional solver time;
- `physical_time_seconds`: available only with validated kinetic mobility;
- `plot_time`: selected display basis without altering `model_time`;
- `time_basis`: `model` or `physical`;
- `time_unit` and `frequency_unit`: explicitly `model time` and
  `cycles / model time` unless a validated physical calibration is loaded;
- `t0_seconds`, `frequency_hz`, `physical_period_seconds`, and
  `physical_duration_seconds`: null in uncalibrated model-time mode;
- `physical_M_a`, `physical_M_s`, and `kinetic_calibration_source`: provenance
  for enabled physical-time mode.

`per_cycle_diagnostics` contains each cycle's absorbed opening mass, sampled
minimum opening barrier, peak first-passage flux, and end survival. It helps
identify first-cycle Gibbs-tail depletion but does not alone certify an event.

## UI plotting behaviour

`kinetic_calibration` stores the full supplied, model-bound calibration when
physical mode is used (including temperature, mobilities and parameter
fingerprint). It is null in model mode. Loading a new calibration never changes
the old result's `time_basis` or arrays.

`state_coordinates`, `probability_state_dimension`, `loading_mode`,
`independent_shear_input`, `orientation_input_active`, `vector_registry_pde`,
and `spatial_specimen_solver` describe the actual production path. Current UI
state is `[a,s]`; the four latter feature flags are false. Separate tensor
traction projection/research mixed-load experiments do not change these flags.
The former orientation widget was unused and is now explicitly disabled.

Plots should preserve scientific-axis offset notation rather than forcing survival to the full 0--1 range. The user should be able to zoom and pan interactively:

- mouse wheel: zoom about the cursor;
- horizontal/vertical constrained zoom where supported;
- drag: pan;
- double click or Home: autoscale/reset.

Each result field should remember its own last view range when the user switches between fields.
