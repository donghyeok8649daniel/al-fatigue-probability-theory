# Master 1D $P$–$u$–$\Theta$ formulation

## Status

This document freezes the active mathematical mainline for the present paper:

$$
\boxed{
\text{1D nonlinear LJ chain}
\longrightarrow
F(\lambda,c,\tau)
\longrightarrow
\{P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)\}
\longrightarrow
\{\bar a,\bar U,\text{first passage}\}.
}
$$

It is intentionally **normal-only**. Registry/slip $s$, FCC geometry,
Smoluchowski closure, Boltzmann equilibrium, Gaussian/Weibull probability
families, empirical damage laws, and ad-hoc damping are not part of this
mainline.

Every equation below is labelled by status:

- **MODEL** — defining assumption of the reduced physical model;
- **DEFINITION** — mathematical definition;
- **EXACT** — exact identity under the stated finite chain / projected measure;
- **CONDITIONAL** — exact only when an explicitly stated additional condition is met;
- **OPEN** — physics or closure not yet derived.

The objective is not to make the reduced probability description look closed
when it is not. The finite microscopic chain is closed. Its one-point
probability projection is exact but hierarchical.

---

# I. Physical and nondimensional variables

Let physical time be $t$ and define the atomic mechanical time scale

$$
\boxed{
t_0=\sqrt{\frac{m_a a_0}{EA_0}}.
}
$$

Here $m_a$ is the represented atomic/repeat mass, $a_0$ the equilibrium normal
spacing, $E$ the reference Young modulus, and $A_0$ the effective 1D reference
area used in the present stress-to-force calibration.

Define nondimensional time

$$
\boxed{\tau=\frac{t}{t_0}.}
$$

Unless explicitly marked otherwise, dots in the microscopic equations below
mean derivatives with respect to $\tau$.

Let

$$
\boxed{
U_{\rm ref}=EA_0a_0,
\qquad
F_{\rm ref}=EA_0.
}
$$

For an applied normal stress $\sigma_n(t)$, the present chain mapping is

$$
\boxed{
q(\tau)=\frac{F_{\rm ext}}{EA_0}=\frac{\sigma_n(t)}{E}.
}
$$

For sinusoidal loading,

$$
\boxed{
\sigma_n(t)=\sigma_m+\sigma_a\sin(2\pi f t),
}
$$

and the nondimensional angular frequency is

$$
\boxed{
\omega^*=2\pi f t_0.
}
$$

**Status:** DEFINITION under the current $A_0,E,a_0$ calibration bridge. The
mapping is not a statement that the 1D LJ chain is an exact atomistic model of
aluminium.

---

# II. Microscopic 1D generalized-LJ chain

Consider $M+1$ atomic nodes $x_0,\ldots,x_M$ in nondimensional position units,
with the left node fixed,

$$
\boxed{x_0(\tau)=0.}
$$

Define the $M$ normalized spacings

$$
\boxed{
\lambda_i=x_i-x_{i-1}>0,
\qquad i=1,\ldots,M.
}
$$

The physical spacing is

$$
\boxed{a_i=a_0\lambda_i.}
$$

## II.1 Normalized generalized-LJ energy

The active nearest-neighbour normalized interaction is

$$
\boxed{
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
\qquad m>n>1.
}
$$

Its first two derivatives are

$$
\boxed{
\phi'(\lambda)
=
\frac{\lambda^{-n-1}-\lambda^{-m-1}}{m-n},
}
$$

and

$$
\boxed{
\phi''(\lambda)
=
\frac{(m+1)\lambda^{-m-2}-(n+1)\lambda^{-n-2}}{m-n}.
}
$$

By construction,

$$
\boxed{
\phi'(1)=0,
\qquad
\phi''(1)=1.
}
$$

The finite-chain configurational energy is

$$
\boxed{
V^*(\boldsymbol\lambda)=\sum_{i=1}^{M}\phi(\lambda_i).
}
$$

The physical energy corresponding to $V^*$ is $U_{\rm ref}V^*$.

**Status:** MODEL. The generalized LJ form is the adopted interaction law.
Everything derived from it below is exact only within this stated model.

## II.2 Closed node equations

For identical unit masses in nondimensional variables, the interior atomic
nodes satisfy

$$
\boxed{
\ddot x_j
=
\phi'(\lambda_{j+1})-\phi'(\lambda_j),
\qquad j=1,\ldots,M-1.
}
$$

The loaded right end satisfies

$$
\boxed{
\ddot x_M
=-\phi'(\lambda_M)+q(\tau).
}
$$

These equations are the closed microscopic dynamics used by the active finite
chain.

**Status:** EXACT under the MODEL above and the stated boundary conditions.

## II.3 Spacing equations

For bulk spacings $i=2,\ldots,M-1$,

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
}
$$

The boundary spacings obey different equations:

$$
\boxed{
\ddot\lambda_1
=
\phi'(\lambda_2)-\phi'(\lambda_1),
}
$$

and

$$
\boxed{
\ddot\lambda_M
=q(\tau)+\phi'(\lambda_{M-1})-2\phi'(\lambda_M).
}
$$

Therefore the bulk formula must not be silently applied at the boundaries.

**Status:** EXACT.

---

# III. Exact energy balance and the spacing-coordinate mass metric

The dimensionless total mechanical energy is

$$
\boxed{
E_{\rm mech}^*
=
T^*+V^*,
\qquad
T^*=\frac12\sum_{j=1}^{M}\dot x_j^2.
}
$$

The external power balance is

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}
=q(\tau)\dot x_M.
}
$$

Thus the present baseline is conservative except for prescribed external work.

Because

$$
x_j=\sum_{k=1}^{j}\lambda_k,
$$

let $\mathbf L$ be the lower-triangular matrix

$$
L_{jk}=\begin{cases}
1,&k\le j,\\
0,&k>j.
\end{cases}
$$

Then

$$
\boxed{
\boldsymbol x=\mathbf L\boldsymbol\lambda,
\qquad
\dot{\boldsymbol x}=\mathbf L\boldsymbol c,
\qquad
c_i=\dot\lambda_i,
}
$$

and therefore

$$
\boxed{
T^*=\frac12\boldsymbol c^T\mathbf G_\lambda\boldsymbol c,
\qquad
\mathbf G_\lambda=\mathbf L^T\mathbf L.
}
$$

Explicitly,

$$
\boxed{
(G_\lambda)_{k\ell}
=M-\max(k,\ell)+1.
}
$$

This is important: the spacing rates are not independent unit-mass generalized
coordinates. Exact total kinetic energy contains cross-correlations between
$c_i$ and $c_j$.

**Status:** EXACT kinematics.

---

# IV. Probability meaning without assuming thermal randomness

For a deterministic finite chain at a fixed time, choose a spacing index $I$
uniformly from the represented set. Then $\lambda_I$ and $c_I$ are random
variables with respect to this **spatial counting measure**. No independence
between different spacings is implied.

Define the finite empirical phase-space measure

$$
\boxed{
F_M(\lambda,c,\tau)
=
\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(\tau)]
\delta[c-c_i(\tau)].
}
$$

It is normalized:

$$
\boxed{
\iint F_M\,dc\,d\lambda=1.
}
$$

The empirical spacing density is

$$
\boxed{
P_M(\lambda,\tau)
=
\int F_M(\lambda,c,\tau)\,dc
=
\frac1M\sum_i\delta[\lambda-\lambda_i(\tau)].
}
$$

Thus the probability density is generated by mechanics and spatial sampling;
it is not selected from a named distribution family.

A smooth $F$ or $P$ below means a continuum/coarse representation of these
empirical measures. Gaussian kernels used in numerical verification are only
estimators; they are not physical Gaussian-distribution assumptions.

**Status:** DEFINITION.

---

# V. Exact empirical phase-space transport

Define the finite acceleration flux

$$
\boxed{
\mathcal G_M(\lambda,c,\tau)
=
\frac1M\sum_i\ddot\lambda_i(\tau)
\delta[\lambda-\lambda_i(\tau)]
\delta[c-c_i(\tau)].
}
$$

Distributional differentiation gives

$$
\boxed{
\partial_\tau F_M
+\partial_\lambda(cF_M)
+\partial_c\mathcal G_M
=0.
}
$$

**Status:** EXACT identity for the empirical measure.

If a smooth representation exists, define

$$
\boxed{
A(\lambda,c,\tau)
=
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c],
}
$$

so that

$$
\mathcal G=AF.
$$

Then

$$
\boxed{
\partial_\tau F
+\partial_\lambda(cF)
+\partial_c(AF)
=0.
}
$$

This equation is exact as a projected identity when $A$ is the true conditional
acceleration. It is **not an autonomous closure** because $A$ can depend on
hidden neighbouring spacings and velocities through their conditional
statistics.

No Markov, Langevin, white-noise, Boltzmann, or Fokker--Planck assumption is
required for this identity.

---

# VI. Complete raw velocity-moment hierarchy

Define raw velocity-moment densities

$$
\boxed{
R_r(\lambda,\tau)
=\int_{-\infty}^{\infty}c^rF(\lambda,c,\tau)\,dc,
\qquad r=0,1,2,\ldots
}
$$

and acceleration-moment sources

$$
\boxed{
B_r(\lambda,\tau)
=\int c^{r-1}A(\lambda,c,\tau)F(\lambda,c,\tau)\,dc,
\qquad r\ge1.
}
$$

Assuming the required moments exist and $c^rAF\to0$ at velocity-space
boundaries, multiplication of the phase-space equation by $c^r$ gives

$$
\boxed{
\partial_\tau R_r
+\partial_\lambda R_{r+1}
=rB_r,
\qquad r\ge0,
}
$$

with the right-hand side zero for $r=0$.

This is the exact one-point hierarchy.

**Status:** EXACT under the stated moment/boundary regularity.

The first raw moments are

$$
R_0=P,
$$

$$
R_1=Pu,
$$

$$
R_2=P(u^2+\Theta),
$$

and

$$
R_3=P(u^3+3u\Theta+C_3),
$$

where

$$
\boxed{
u(\lambda,\tau)=\mathbb E[c\mid\lambda]}
$$

is the conditional mean spacing rate,

$$
\boxed{
\Theta(\lambda,\tau)
=\operatorname{Var}(c\mid\lambda)
=\mathbb E[(c-u)^2\mid\lambda],
}
$$

and

$$
\boxed{
C_3(\lambda,\tau)
=\mathbb E[(c-u)^3\mid\lambda].
}
$$

---

# VII. Exact $P$ and $u$ equations

The $r=0$ equation is

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0.
}
$$

Define the spacing-space current

$$
\boxed{J=Pu.}
$$

Then

$$
\boxed{
\partial_\tau P=-\partial_\lambda J.
}
$$

The $r=1$ equation is

$$
\boxed{
\partial_\tau(Pu)
+\partial_\lambda\left[P(u^2+\Theta)\right]
=P\mathcal A,
}
$$

where

$$
\boxed{
\mathcal A(\lambda,\tau)
=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda].
}
$$

Using continuity,

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta),
}
$$

with

$$
\boxed{
D_\tau u=\partial_\tau u+u\partial_\lambda u.
}
$$

**Status:** EXACT.

---

# VIII. Exact $\Theta$ evolution — corrected general form

The $r=2$ raw moment equation is

$$
\partial_\tau[P(u^2+\Theta)]
+\partial_\lambda[P(u^3+3u\Theta+C_3)]
=2P\mathbb E[c\ddot\lambda\mid\lambda].
$$

Define

$$
\boxed{
\Psi(\lambda,\tau)
=\operatorname{Cov}(c,\ddot\lambda\mid\lambda)
=\mathbb E[(c-u)\ddot\lambda\mid\lambda].
}
$$

Combining the raw second-moment equation with the continuity and mean-velocity
equations gives

$$
\boxed{
D_\tau\Theta
+2\Theta\,\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi.
}
$$

This is the correct general one-point second-central-moment equation for the
projected LJ chain.

The shorter equation

$$
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)=0
$$

is only valid when

$$
\boxed{\Psi=0.}
$$

That is automatic if acceleration is deterministic at fixed $\lambda$ and
$\tau$, but it is **not automatic in the actual spatial LJ chain**, because
$\ddot\lambda_i$ depends on neighbouring spacings.

**Status:** EXACT general equation; zero-source version is CONDITIONAL.

Therefore a closed evolution for $(P,u,\Theta)$ has not been derived. Exact
prediction introduces at least $C_3$, $\Psi$, and neighbour-conditioned
statistics.

---

# IX. Full nonlinear LJ acceleration enters through neighbour joint states

For bulk spacings,

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}).
$$

Define the ordered neighbour joint densities

$$
\boxed{
P_2^+(\lambda,\lambda',\tau)
}
$$

for central/right-neighbour spacings and

$$
\boxed{
P_2^-(\lambda,\lambda',\tau)
}
$$

for central/left-neighbour spacings, with central marginal $P$.

Then

$$
\boxed{
\mathcal A_{\rm bulk}(\lambda,\tau)
=m_+(\lambda,\tau)+m_-(\lambda,\tau)-2\phi'(\lambda),
}
$$

where

$$
\boxed{
m_+(\lambda,\tau)
=\frac1P\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)\,d\lambda',
}
$$

and

$$
\boxed{
m_-(\lambda,\tau)
=\frac1P\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)\,d\lambda'.
}
$$

No neighbour-independence approximation is made.

To expose the source in the $\Theta$ equation, define central-velocity/neighbour
joint densities

$$
F_2^+(\lambda,c,\lambda',\tau),
\qquad
F_2^-(\lambda,c,\lambda',\tau).
$$

At fixed central $\lambda$, the term proportional to the central force
$-2\phi'(\lambda)$ drops out of the covariance because
$\mathbb E[c-u\mid\lambda]=0$. Therefore

$$
\boxed{
\Psi_{\rm bulk}
=\frac1P\iint
(c-u)\phi'(\lambda')
\left[F_2^++F_2^-\right]
\,dc\,d\lambda'.
}
$$

This is the precise mechanical origin of the acceleration-covariance source in
the $\Theta$ equation.

Boundary spacings have their own acceleration laws and must be included with
separate boundary statistics if they are retained in the one-point average.

**Status:** EXACT for bulk spacings.

---

# X. Exact instantaneous shape relation for $P$

From

$$
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta),
$$

expand

$$
\frac1P\partial_\lambda(P\Theta)
=\partial_\lambda\Theta
+\Theta\partial_\lambda\ln P.
$$

Hence

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta.
}
$$

Where $P>0$ and $\Theta>0$ are smooth,

$$
\boxed{
\partial_\lambda\ln P
=\frac{\mathcal A-D_\tau u}{\Theta}
-\partial_\lambda\ln\Theta.
}
$$

Integrating at fixed $\tau$ gives

$$
\boxed{
P(\lambda,\tau)
=\frac{C(\tau)}{\Theta(\lambda,\tau)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}
{\Theta(\eta,\tau)}\,d\eta
\right].
}
$$

For a normalized intact density on the full positive-spacing domain,
$C(\tau)$ is fixed by

$$
\boxed{
\int_0^\infty P(\lambda,\tau)\,d\lambda=1.
}
$$

**Status:** EXACT instantaneous shape representation under smoothness,
$P>0$, $\Theta>0$, and existence of the required moments.

It is not an independent evolution law because $D_\tau u$ is itself obtained
from the same first-moment balance. It is primarily a reconstruction,
consistency, and shape constraint.

## X.1 Degenerate regime

If

$$
\Theta(\lambda,\tau)=0,
$$

the divided shape formula is invalid. The undivided first-moment relation and
continuity equation remain valid.

For the ideal homogeneous initial state

$$
\boxed{
\lambda_i(0)=1,
\qquad
c_i(0)=0,
}
$$

one has

$$
\boxed{
F_M(\lambda,c,0)=\delta(\lambda-1)\delta(c),
\qquad
P_M(\lambda,0)=\delta(\lambda-1),
\qquad
\Theta=0.
}
$$

Thus the smooth shape formula is not used to invent initial broadening.
Broadening is generated by the actual finite-chain dynamics and boundary-wave
propagation.

---

# XI. Same-force history dependence

Let $\tau_L$ and $\tau_U$ be loading and unloading times such that

$$
q(\tau_L)=q(\tau_U)=q^*,
$$

with

$$
\dot q(\tau_L)>0,
\qquad
\dot q(\tau_U)<0.
$$

Define the second-order reduced descriptor

$$
\boxed{
\mathcal R_2(\tau)=\{P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)\}.
}
$$

If

$$
\boxed{
\mathcal R_2(\tau_L)\ne\mathcal R_2(\tau_U),
}
$$

then no memoryless map

$$
\mathcal R_2=\mathcal S[q(\tau)]
$$

can describe that trajectory.

The current nonlinear-chain numerical test shows this same-force
non-retracing. Therefore $(P,u,\Theta)$ is a **history-bearing reduced
descriptor**.

However, because its exact evolution introduces $C_3$, $\Psi$, and neighbour
joint statistics, $\mathcal R_2$ is not claimed to be a closed Markov state.

---

# XII. G1 — mean spacing

Using the normalized spacing density,

$$
\boxed{
\bar\lambda(\tau)
=\int_0^\infty\lambda P(\lambda,\tau)\,d\lambda.
}
$$

The physical mean spacing is

$$
\boxed{
\bar a(t)=a_0\bar\lambda(t/t_0).
}
$$

From continuity,

$$
\boxed{
\frac{d\bar\lambda}{d\tau}
=-[\lambda J]_{0}^{\infty}
+\int_0^\infty J\,d\lambda.
}
$$

When spacing-space boundary fluxes vanish,

$$
\boxed{
\frac{d\bar\lambda}{d\tau}
=\int Pu\,d\lambda
=\mathbb E[c].
}
$$

**Status:** DEFINITION plus EXACT moment identity.

---

# XIII. G2 — mean intrinsic configurational energy

For the active nearest-neighbour chain, use the same intrinsic energy that
generates the microscopic force. Define the one-bond reference-subtracted
energy

$$
\boxed{
\Delta\phi(\lambda)=\phi(\lambda)-\phi(1).
}
$$

Then the mean intrinsic configurational energy per represented spacing is

$$
\boxed{
\bar U(\tau)
=U_{\rm ref}
\int_0^\infty
\Delta\phi(\lambda)P(\lambda,\tau)\,d\lambda.
}
$$

For the finite empirical density over all $M$ spacings,

$$
\boxed{
V_{\rm phys}-M U_{\rm ref}\phi(1)
=M\bar U.
}
$$

Thus G2 is exactly consistent with the same nearest-neighbour energy used in
the dynamics.

Its rate is

$$
\boxed{
\frac{1}{U_{\rm ref}}
\frac{d\bar U}{d\tau}
=-[\Delta\phi J]_{0}^{\infty}
+\int_0^\infty\phi'(\lambda)J\,d\lambda.
}
$$

When boundary fluxes vanish,

$$
\boxed{
\frac{d\bar U}{d\tau}
=U_{\rm ref}\int\phi'(\lambda)P u\,d\lambda.
}
$$

This is a configurational-energy moment. It is not the full mechanical energy,
which also contains the nonlocal kinetic term of Section III.

**Status:** DEFINITION and EXACT under the active nearest-neighbour chain.

A long-range/zeta energy must not be inserted into G2 and simultaneously called
mechanically exact unless the equations of motion are derived consistently from
that same long-range energy.

---

# XIV. G3 — irreversible hysteresis energy: mathematical status

The fixed observable definition remains

$$
\boxed{
E_{\rm hyst}(t)
=\int_0^t\dot D_{\rm irr}(t')\,dt',
\qquad
\dot D_{\rm irr}\ge0.
}
$$

The present conservative LJ chain contains no irreversible force, so its exact
baseline is

$$
\boxed{
\dot D_{\rm irr}=0,
\qquad
E_{\rm hyst}=0.
}
$$

This does **not** mean the reduced state retraces at fixed force; inertia and
wave propagation can produce non-retracing while total mechanical energy is
conserved apart from external work.

To show what a mechanically consistent irreversible extension would require,
let additional nondimensional irreversible node forces $r_j^{\rm irr}$ be
introduced explicitly. Then

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}
=q\dot x_M
+\sum_{j=1}^{M}r_j^{\rm irr}\dot x_j.
}
$$

If the mechanism satisfies

$$
\sum_jr_j^{\rm irr}\dot x_j\le0,
$$

define

$$
\boxed{
\dot D_{\rm irr}^*
=-\sum_jr_j^{\rm irr}\dot x_j\ge0.
}
$$

Then

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}
=q\dot x_M-\dot D_{\rm irr}^*.
}
$$

Over one cycle,

$$
\boxed{
W_{\rm ext}^{\rm cyc}
=\Delta E_{\rm mech}^{\rm cyc}
+D_{\rm irr}^{\rm cyc}.
}
$$

Therefore a loop area can be identified with irreversible dissipation only when
the stored mechanical-energy change is separately accounted for. In the
current conservative simulations,

$$
D_{\rm irr}^{\rm cyc}=0
$$

and cycle work equals the mechanical-energy change during transients.

The older current/mobility expression for $\dot D_{\rm irr}$ is retained only
as a **CONDITIONAL Smoluchowski-type realization** requiring an overdamped
mobility law and its associated assumptions. It is not part of the exact
current mainline.

**Status:** G3 observable DEFINED; physical irreversible mechanism OPEN.

---

# XV. G4 — normalization and first-passage survival

## XV.1 Nonabsorbing density

Before imposing an initiation sink,

$$
\boxed{
\int_0^\infty P(\lambda,\tau)\,d\lambda=1.
}
$$

An instantaneous tail diagnostic is

$$
\boxed{
Q_c(\tau)
=\int_{\lambda_c}^{\infty}P(\lambda,\tau)\,d\lambda.
}
$$

This is **not** cumulative crack-initiation probability because probability in
a nonabsorbing computation can return below $\lambda_c$.

## XV.2 Mechanical instability threshold

The current operational local initiation threshold is the first loss of
positive tangent stiffness,

$$
\boxed{
\phi''(\lambda_c)=0.
}
$$

Hence

$$
\boxed{
\lambda_c
=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
}
$$

For $m=12.19$, $n=6$,

$$
\lambda_c\approx1.1077715386.
$$

This is an operational local stability-loss criterion within the ideal chain;
it is not a claim that a macroscopic free crack surface is already fully
formed at that instant.

## XV.3 Exact finite empirical first passage

For each represented spacing define the first-hitting time

$$
\boxed{
\tau_i^c
=\inf\{\tau\ge0:\lambda_i(\tau)\ge\lambda_c\}.
}
$$

Define the alive indicator

$$
\chi_i(\tau)=\mathbf 1_{\{\tau<\tau_i^c\}}.
$$

The finite local survivor fraction is

$$
\boxed{
S_M(\tau)
=\frac1M\sum_i\chi_i(\tau).
}
$$

The cumulative local first-passage fraction is

$$
\boxed{
F_{{\rm ci},M}^{\rm local}(\tau)
=1-S_M(\tau).
}
$$

In the distributional sense,

$$
\boxed{
-\frac{dS_M}{d\tau}
=\frac1M\sum_i\delta(\tau-\tau_i^c).
}
$$

**Status:** EXACT finite empirical first-passage definition.

## XV.4 Smooth kinetic survivor equation

Let $F_b(\lambda,c,\tau)$ be the survivor phase-space subdensity on

$$
0<\lambda<\lambda_c.
$$

Because the dynamics is second order, the correct kinetic absorbing condition
is a **no-inflow condition** at the right boundary: once a trajectory has left
the intact domain, it is not reintroduced. Thus

$$
\boxed{
F_b(\lambda_c,c,\tau)=0
\quad\text{for incoming }c<0.
}
$$

Outgoing states with $c>0$ carry escape flux. Assuming no loss through the
lower boundary,

$$
\boxed{
S(\tau)
=\int_0^{\lambda_c}\int_{-\infty}^{\infty}
F_b(\lambda,c,\tau)\,dc\,d\lambda,
}
$$

and

$$
\boxed{
j_{\rm esc}(\tau)
=\int_0^\infty
cF_b(\lambda_c^-,c,\tau)\,dc\ge0.
}
$$

Therefore

$$
\boxed{
\frac{dS}{d\tau}=-j_{\rm esc}(\tau).
}
$$

The local cumulative initiation and hazard are

$$
\boxed{
F_{\rm ci}^{\rm local}=1-S,
}
$$

and, when $S>0$,

$$
\boxed{
h(\tau)=\frac{j_{\rm esc}(\tau)}{S(\tau)}
=-\frac{d}{d\tau}\ln S(\tau).
}
$$

**Status:** EXACT kinetic first-passage balance under the stated absorbing and
lower-boundary conditions.

## XV.5 Survivor-conditioned moments

The survivor subdensity is not normalized to one. Therefore moments over it
must be labelled correctly. Define the normalized conditional survivor density

$$
\boxed{
\widehat P_b(\lambda,\tau)
=\frac{P_b(\lambda,\tau)}{S(\tau)}.
}
$$

Then

$$
\boxed{
\bar\lambda_{\rm surv}
=\frac1S\int_0^{\lambda_c}\lambda P_b\,d\lambda,
}
$$

and

$$
\boxed{
\bar U_{\rm surv}
=\frac{U_{\rm ref}}{S}
\int_0^{\lambda_c}\Delta\phi(\lambda)P_b\,d\lambda.
}
$$

Multiplication by the time-only factor $1/S$ does not change the spatial
logarithmic slope, so the local shape identity has the same form for
$\widehat P_b$ in the interior when the conditional fields are defined over
survivors.

---

# XVI. Local probability is not automatically specimen crack probability

For a single deterministic chain, define the specimen first-initiation time

$$
\boxed{
\tau_{\rm spec}^c
=\min_i\tau_i^c.
}
$$

This is a deterministic event time for that microscopic realization.

The one-point survivor fraction $S_M$ is a fraction of represented local
spacings that have not yet crossed. In general,

$$
\boxed{
1-S_M(\tau)
\ne
\Pr(\tau_{\rm spec}^c\le\tau).
}
$$

To define a specimen-level probability, an ensemble of microscopic/specimen
realizations $\omega$ is required:

$$
\boxed{
S_{\rm spec}(\tau)
=\Pr_\omega\!\left[
\min_i\tau_i^c(\omega)>\tau
\right].
}
$$

No independence product such as $S_{\rm spec}=S_{\rm local}^{N}$ is adopted.
The strong spatial correlations already observed in the chain prevent such a
step without a separately validated statistical-length model.

**Status:** specimen-probability bridge OPEN.

---

# XVII. What is mathematically closed and what is not

## XVII.1 Closed exactly

The finite microscopic state

$$
\boxed{
\mathbf Z(\tau)
=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
}
$$

with the LJ node equations and boundary forcing is a closed deterministic
system.

From any solution $\mathbf Z(\tau)$, the following are exactly generated:

$$
F_M,
\quad
P_M,
\quad
u,
\quad
\Theta,
\quad
C_3,
\quad
\Psi,
\quad
P_2^\pm,
\quad
F_2^\pm,
\quad
\bar a,
\quad
\bar U,
\quad
\tau_i^c.
$$

## XVII.2 Exact but not autonomously closed after projection

The one-point equations

$$
\partial_\tau P+\partial_\lambda(Pu)=0,
$$

$$
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta),
$$

and

$$
D_\tau\Theta
+2\Theta u_{,\lambda}
+\frac1P\partial_\lambda(PC_3)
=2\Psi
$$

are exact but not a closed three-field PDE system because

$$
\boxed{
\{P,u,\Theta\}
\longrightarrow
\{C_3,\Psi,P_2^\pm,F_2^\pm,\ldots\}.
}
$$

No arbitrary closure is inserted in the active theory.

## XVII.3 Physically unresolved

The following are not yet supplied by the current conservative mainline:

1. a mechanically justified irreversible mechanism giving $\dot D_{\rm irr}>0$;
2. a laboratory-fatigue-time-scale mechanism that preserves or converts the
   microscopic history dependence over many Hz-scale cycles;
3. a validated specimen-level probability bridge from correlated local first
   passage;
4. experimental validation/calibration of the predicted distribution and
   initiation statistics.

These are physical open problems, not missing algebraic terms hidden inside the
present equations.

---

# XVIII. Assumption ledger

The active master formulation assumes only the following at its base level:

1. one-dimensional normal motion;
2. identical finite chain masses after nondimensionalization;
3. nearest-neighbour generalized-LJ configurational energy
   $V^*=\sum_i\phi(\lambda_i)$;
4. fixed left boundary and prescribed right-end normal force;
5. the current stress mapping $q=\sigma/E$ through the calibrated $A_0$ scale;
6. an empirical/spatial probability measure over represented spacings;
7. smooth $F,P$ only when invoking continuum moment equations;
8. sufficient velocity-moment decay to remove $c$-boundary terms;
9. $P>0$ and $\Theta>0$ only where the divided density-shape formula is used;
10. local tangent-stiffness loss $\phi''(\lambda_c)=0$ as the operational
    initiation threshold.

The active master formulation does **not** assume:

- Boltzmann/Gibbs equilibrium;
- Gaussian, Weibull, lognormal, or any named $P$ family;
- independence of neighbouring spacings;
- Markovian stochastic dynamics;
- Langevin noise or white noise;
- Fokker--Planck/Smoluchowski closure;
- empirical fatigue damage accumulation;
- viscous damping;
- an independent FEM element probability product;
- FCC geometry or registry/slip coordinates.

---

# XIX. Final compact equation set

The active paper-level mathematical core can be written as the following chain.

### Microscopic mechanics

$$
\boxed{
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
}
$$
for bulk spacings, with the actual boundary laws retained separately.

### Empirical phase-space state

$$
\boxed{
F_M(\lambda,c,\tau)
=\frac1M\sum_i\delta(\lambda-\lambda_i)\delta(c-c_i)
}
$$

### Exact projected transport

$$
\boxed{
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
}
$$

### Reduced fields

$$
\boxed{
P=\int Fdc,
\qquad
u=\mathbb E[c\mid\lambda],
\qquad
\Theta=\operatorname{Var}(c\mid\lambda)
}
$$

### Continuity

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0
}
$$

### Mean-flow balance

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta)
}
$$

### Exact $P$ shape

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta
}
$$

### Exact $\Theta$ balance

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi
}
$$

### Full LJ conditional acceleration

$$
\boxed{
\mathcal A
=m_++m_- -2\phi'(\lambda)
}
$$

with $m_\pm$ generated from neighbour joint densities.

### G1

$$
\boxed{
\bar a=a_0\int\lambda P\,d\lambda
}
$$

### G2

$$
\boxed{
\bar U
=U_{\rm ref}\int[\phi(\lambda)-\phi(1)]P\,d\lambda
}
$$

### G3

$$
\boxed{
E_{\rm hyst}=\int\dot D_{\rm irr}\,dt,
\qquad
\dot D_{\rm irr}\ge0
}
$$

with $\dot D_{\rm irr}=0$ in the present conservative baseline and the physical
irreversible mechanism still OPEN.

### G4 / first passage

$$
\boxed{
S(\tau)=\iint_{\lambda<\lambda_c}F_b\,dc\,d\lambda,
\qquad
\dot S=-j_{\rm esc},
\qquad
F_{\rm ci}^{\rm local}=1-S
}
$$

with

$$
\boxed{
\phi''(\lambda_c)=0.
}
$$

This is the mathematically complete active 1D formulation. What remains is not
an omitted algebraic step; it is the physical derivation of irreversibility,
fatigue-scale memory, specimen-scale probability, and experimental validation.
