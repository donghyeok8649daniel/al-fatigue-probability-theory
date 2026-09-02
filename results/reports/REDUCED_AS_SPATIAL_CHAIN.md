# Deterministic reduced $(a_i,s_i)$ spatial-chain verification

## Purpose

This run tests Milestone 21 without FCC geometry, stochastic forcing, a named PDF,
Boltzmann equilibrium, Fokker--Planck dynamics, or damping.  The only spatial
extension is the explicitly declared reduced-model assumption

$$
V_M=\sum_i U_0(a_i,s_i),\qquad a_i=x_{i+1}-x_i,
$$

with the same active multiplicity-free row/layer energy $U_0(a,s)$ used in
`theory-core`.

The numerical evaluator uses the direct $(k,p)$ sum with `kmax=12`, `pmax=30`.
These are numerical truncations, not physical interaction cutoffs.

## Normalized test setup

- cells: $M=48$;
- $m=12$, $n=6$, $b=\sigma_{LJ}=\epsilon_{LJ}=1$;
- stable reference found numerically at $a_0=0.9919601754$ and $s_0/b=0.5$;
- the $s_0=0.5$ phase is only a registry-origin convention; shifting
  $\tilde s=s-b/2$ places the same well at zero;
- velocity-Verlet step: $dt=0.004$;
- cyclic frequency: $\omega=0.35$;
- normal boundary drive after one-cycle ramp:

$$
Q_a(t)=0.2+0.2\sin(\omega t).
$$

Two cases were evaluated.

1. `normal_only_5cycle`: $q_s(t)=0$.
2. `registry_driven_5cycle`: $q_s(t)=0.8\sin(\omega t)$ after the same ramp.
   This is a **normalized mechanism diagnostic**, not a calibrated Al loading
   or Schmid-law claim.
3. `registry_driven_10cycle`: the second case extended to ten cycles as a
   persistence check.

## Result 1: boundary loading alone generates a nontrivial spacing density

The normal-only run starts exactly from

$$
a_i(0)=a_0,\qquad s_i(0)=s_0,
$$

with all velocities zero.  Nevertheless the right-boundary force produces a
spatially nonuniform wave field.  Over five cycles,

$$
\max_t \operatorname{Var}_i(a_i)=9.4562\times10^{-5},
$$

with

$$
0.972760\le a_i\le1.023305.
$$

Therefore the empirical measure in the normal coordinate broadens from its
initial delta **without random forcing**.

The corresponding global spacing-velocity covariance diagnostic reaches

$$
\max_t\widehat\Theta^{\rm glob}_{aa}
=1.1955\times10^{-5}.
$$

The hat and `glob` label are essential: this is the finite-cell unconditional
velocity covariance, not the conditional field
$\Theta_{aa}(a,s,t)=\operatorname{Var}(\dot a\mid a,s)$.

## Result 2: exact registry symmetry blocks spontaneous $s$ spreading

In the same normal-only run,

$$
\max_t\operatorname{Var}_i(s_i)
=1.41\times10^{-19},
$$

and the observed registry range remains

$$
0.49999999897\le s_i\le0.50000000007.
$$

Similarly,

$$
\max_t|\widehat\Theta^{\rm glob}_{as}|
=4.29\times10^{-12},
$$

$$
\max_t\widehat\Theta^{\rm glob}_{ss}
=2.35\times10^{-18}.
$$

These values are numerical zero at the truncation/integration level.  Thus the
simulation confirms the Milestone-21 symmetry statement: if the system starts
exactly at the symmetric registry well and no physical registry drive or other
symmetry-breaking mechanism is present, normal spatial heterogeneity alone does
**not** generate a meaningful $s$ distribution.

This is a falsifiable restriction of the current reduced model, not a defect to
be hidden with random noise.

## Result 3: once registry symmetry is physically displaced, the existing
## $U_0(a,s)$ coupling transfers diversity into $s$

The declared registry-driven diagnostic gives, over five cycles,

$$
\max_t\operatorname{Var}_i(s_i)
=3.8972\times10^{-5},
$$

at cycle $4.7126$, while

$$
\max_t\operatorname{Var}_i(a_i)
=8.7061\times10^{-5}.
$$

The registry range becomes

$$
0.448875\le s_i\le0.534739.
$$

The global velocity-covariance diagnostics become genuinely coupled:

$$
\max_t\widehat\Theta^{\rm glob}_{aa}=9.5510\times10^{-6},
$$

$$
\max_t|\widehat\Theta^{\rm glob}_{as}|=7.0778\times10^{-6},
$$

$$
\max_t\widehat\Theta^{\rm glob}_{ss}=6.1267\times10^{-6}.
$$

Thus a nonzero cross-velocity covariance can emerge deterministically from
boundary-generated spatial diversity plus the already existing mixed
$U_0(a,s)$ dependence.  No correlated stochastic forcing was inserted.

At the snapshot of maximum registry variance, the finite-cell state covariance
matrix is approximately

$$
\operatorname{Cov}(a,s)=
\begin{bmatrix}
8.2527\times10^{-5} & -5.6581\times10^{-5}\\
-5.6581\times10^{-5} & 3.8972\times10^{-5}
\end{bmatrix}.
$$

Its eigenvalues are

$$
1.2164\times10^{-7},\qquad1.2138\times10^{-4},
$$

with ratio about $1.00\times10^{-3}$.  Therefore the cloud is **not a broad
full-rank 2D blob yet**; it is a thin, strongly correlated state manifold.  The
state correlation at this snapshot is about $-0.9977$.

This matters for the exact density-shape equation: the divided
$\boldsymbol\Theta^{-1}$ form must not automatically be used merely because the
global covariance is nonzero.  A local conditional covariance field can remain
ill-conditioned or singular when the deterministic state cloud is nearly
lower-dimensional.

## Result 4: mean intrinsic energy is generated directly from the state cloud

For every snapshot the finite-cell G2 counterpart is

$$
\bar U_M(t)=\frac1M\sum_i
\left[U_0(a_i,s_i)-U_0(a_0,s_0)\right].
$$

The maximum values observed in the five-cycle runs are

$$
\max_t\bar U_M=0.0211378
$$

for normal-only loading and

$$
\max_t\bar U_M=0.0348324
$$

for the registry-driven diagnostic.

The ten-cycle registry-driven check reaches

$$
\max_t\bar U_M=0.0450945,
$$

while the peak registry variance rises to

$$
5.3830\times10^{-5}.
$$

This confirms the intended chain

$$
\text{boundary mechanics}\to\{a_i,s_i\}\to P_M(a,s,t)\to\bar U_M(t)
$$

without assigning an analytic PDF family.

## Result 5: conservative work/energy balance is numerically satisfied

Because the model contains no damping or irreversible mechanism, the relevant
trajectory-level validation is

$$
\Delta E_{\rm mech}=W_{\rm ext}.
$$

The final relative residuals are

- normal-only, five cycles: $5.90\times10^{-8}$;
- registry-driven, five cycles: $4.75\times10^{-8}$;
- registry-driven, ten cycles: $4.00\times10^{-7}$.

Therefore the observed state spreading is not numerical artificial
dissipation.  It is conservative spatial redistribution.  Consequently this
simulation still does **not** provide G3; a physical irreversible mechanism is
required before positive cumulative hysteretic dissipation is claimed.

## What this simulation establishes

The current reduced model now has a demonstrated deterministic source for the
normal part of the empirical state density:

$$
\boxed{
\text{cyclic boundary force}
\to\text{wave/spatial nonuniformity}
\to P_M(a,t).
}
$$

It also demonstrates a conditional route for the joint state:

$$
\boxed{
\text{registry displacement}
+U_{as}
+\text{normal spatial diversity}
\to P_M(a,s,t),\ \widehat\Theta^{\rm glob}_{as}\ne0.
}
$$

But it simultaneously shows a limitation: with exact registry symmetry and no
physical $s$-drive, the active model stays effectively one-dimensional in
state space.  Therefore a two-coordinate fatigue theory cannot obtain a broad
$s$ distribution merely by declaring the coordinate.  A physically justified
registry-driving or symmetry-breaking mechanism must be supplied.

## Probability-scale caution

The generated

$$
P_M(a,s,t)=\frac1M\sum_i\delta(a-a_i)\delta(s-s_i)
$$

is a **spatial one-point empirical measure** over cells in one deterministic
realization.  It supports G1/G2 and local tail diagnostics.  It is not by itself
specimen-to-specimen crack probability, and cell states are strongly
correlated.  No independent-cell product survival law is used here.

## Files

- simulation: `simulations/run_reduced_as_spatial_chain.py`
- summary: `results/data/reduced_as_spatial_chain/summary.csv`
- raw peak cloud: `results/data/reduced_as_spatial_chain/peak_state_registry_driven.csv`
