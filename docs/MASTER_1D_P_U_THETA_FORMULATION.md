# Master 1D $P$–$u$–$\Theta$ formulation

## 0. Active scope

The active mathematical mainline is

$$
\boxed{
\text{1D nonlinear LJ chain}
\to F_M(\lambda,c,\tau)
\to \{P,u,\Theta\}
\to \{\bar a,\bar U,\text{first passage}\}.
}
$$

This document uses five labels:

- **MODEL**: adopted physical model;
- **DEFINITION**: mathematical definition;
- **EXACT**: exact identity under the stated model;
- **CONDITIONAL**: exact only under an additional stated condition;
- **OPEN**: physical ingredient not yet derived.

No Boltzmann/Gibbs equilibrium, Gaussian/Weibull PDF family, Fokker--Planck or
Smoluchowski closure, white noise, empirical damage law, FCC geometry, or slip
coordinate is required by the active 1D formulation.

---

## 1. Physical scaling

Let $a_0$ be the reference spacing, $m_a$ the represented atomic/repeat mass,
$E$ the reference Young modulus, and $A_0$ the current effective 1D reference
area. Define

$$
\boxed{
t_0=\sqrt{\frac{m_a a_0}{EA_0}},
\qquad
\tau=\frac{t}{t_0}.
}
$$

The force and energy scales are

$$
\boxed{
F_{\rm ref}=EA_0,
\qquad
U_{\rm ref}=EA_0a_0.
}
$$

The current stress-to-chain-force mapping is

$$
\boxed{
q(\tau)=\frac{F_{\rm ext}}{EA_0}=\frac{\sigma_n(t)}{E}.
}
$$

For

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(2\pi f t),
$$

the nondimensional angular frequency is

$$
\boxed{\omega^*=2\pi f t_0.}
$$

**Status:** DEFINITION under the current calibration bridge.

---

## 2. Microscopic generalized-LJ chain

Take $M+1$ nodes $x_0,\ldots,x_M$ in units of $a_0$, with

$$
\boxed{x_0(\tau)=0.}
$$

Define normalized spacings

$$
\boxed{
\lambda_i=x_i-x_{i-1}>0,
\qquad i=1,\ldots,M,
}
$$

and physical spacings $a_i=a_0\lambda_i$.

The normalized generalized-LJ energy is

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

Hence

$$
\boxed{
\phi'(\lambda)
=\frac{\lambda^{-n-1}-\lambda^{-m-1}}{m-n},
}
$$

$$
\boxed{
\phi''(\lambda)
=\frac{(m+1)\lambda^{-m-2}-(n+1)\lambda^{-n-2}}{m-n},
}
$$

with

$$
\boxed{\phi'(1)=0,\qquad \phi''(1)=1.}
$$

The dimensionless configurational energy is

$$
\boxed{
V^*(\boldsymbol\lambda)=\sum_{i=1}^{M}\phi(\lambda_i).
}
$$

**Status:** MODEL.

### 2.1 Exact node equations

For identical unit masses after nondimensionalization,

$$
\boxed{
\ddot x_j
=\phi'(\lambda_{j+1})-\phi'(\lambda_j),
\qquad j=1,\ldots,M-1,
}
$$

and the loaded end obeys

$$
\boxed{
\ddot x_M=-\phi'(\lambda_M)+q(\tau).
}
$$

Dots in Sections 2--15 denote derivatives with respect to $\tau$ unless
physical time is written explicitly.

**Status:** EXACT under the MODEL and boundary conditions.

### 2.2 Exact spacing equations

For $i=2,\ldots,M-1$,

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
}
$$

The boundary spacings are different:

$$
\boxed{
\ddot\lambda_1
=\phi'(\lambda_2)-\phi'(\lambda_1),
}
$$

$$
\boxed{
\ddot\lambda_M
=q(\tau)+\phi'(\lambda_{M-1})-2\phi'(\lambda_M).
}
$$

**Status:** EXACT. The bulk equation must not be used at the boundaries without
these corrections.

---

## 3. Exact mechanical-energy structure

The dimensionless mechanical energy is

$$
\boxed{
E_{\rm mech}^*=T^*+V^*,
\qquad
T^*=\frac12\sum_{j=1}^{M}\dot x_j^2.
}
$$

The exact external-power balance is

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}=q(\tau)\dot x_M.
}
$$

Thus the baseline chain is conservative apart from prescribed boundary work.

### 3.1 Spacing-coordinate mass metric

Since

$$
x_j=\sum_{k=1}^{j}\lambda_k,
$$

define $L_{jk}=1$ for $k\le j$ and $0$ otherwise. With
$c_i=\dot\lambda_i$,

$$
\boxed{
\boldsymbol x=\mathbf L\boldsymbol\lambda,
\qquad
\dot{\boldsymbol x}=\mathbf L\boldsymbol c.
}
$$

Therefore

$$
\boxed{
T^*=\frac12\boldsymbol c^T\mathbf G_\lambda\boldsymbol c,
\qquad
\mathbf G_\lambda=\mathbf L^T\mathbf L,
}
$$

where, using one-based indices,

$$
\boxed{
(G_\lambda)_{k\ell}=M-\max(k,\ell)+1.
}
$$

Thus spacing rates are not independent unit-mass velocities. Exact total
kinetic energy requires cross-spacing rate correlations.

**Status:** EXACT kinematics.

---

## 4. Probability meaning

For a deterministic chain at fixed time, sample a spacing index uniformly from
the represented spacings. The resulting probability is a **spatial empirical
counting measure**; it does not assume thermal randomness or independence.

Define

$$
\boxed{
F_M(\lambda,c,\tau)
=\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(\tau)]
\delta[c-c_i(\tau)].
}
$$

Then

$$
\boxed{
\iint F_M\,dc\,d\lambda=1
}
$$

and

$$
\boxed{
P_M(\lambda,\tau)
=\int F_M\,dc
=\frac1M\sum_i\delta[\lambda-\lambda_i(\tau)].
}
$$

A smooth $F$ or $P$ is a continuum/coarse representation of the empirical
measure. Numerical KDE smoothing is an estimator only, not a physical Gaussian
PDF assumption.

**Status:** DEFINITION.

---

## 5. Exact empirical phase-space transport

Define the empirical acceleration flux

$$
\boxed{
\mathcal G_M(\lambda,c,\tau)
=\frac1M\sum_i\ddot\lambda_i
\delta(\lambda-\lambda_i)\delta(c-c_i).
}
$$

Distributional differentiation gives

$$
\boxed{
\partial_\tau F_M
+\partial_\lambda(cF_M)
+\partial_c\mathcal G_M=0.
}
$$

**Status:** EXACT empirical identity.

For a smooth representation define

$$
\boxed{
A(\lambda,c,\tau)
=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c].
}
$$

Then $\mathcal G=AF$ and

$$
\boxed{
\partial_\tau F
+\partial_\lambda(cF)
+\partial_c(AF)=0.
}
$$

This is an exact projected identity when $A$ is the true conditional
acceleration. It is not an autonomous one-point closure because $A$ retains
hidden-neighbour information.

---

## 6. Complete raw moment hierarchy

Define

$$
\boxed{
R_r(\lambda,\tau)=\int c^rF(\lambda,c,\tau)\,dc,
\qquad r=0,1,2,\ldots
}
$$

and, for $r\ge1$,

$$
\boxed{
B_r(\lambda,\tau)
=\int c^{r-1}A(\lambda,c,\tau)F(\lambda,c,\tau)\,dc.
}
$$

If the required moments exist and the velocity-boundary terms vanish,

$$
\boxed{
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r,
\qquad r\ge0,
}
$$

with zero right-hand side for $r=0$.

**Status:** EXACT.

Define

$$
\boxed{
u(\lambda,\tau)=\mathbb E[c\mid\lambda]}
$$

as the conditional mean spacing rate,

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
C_3(\lambda,\tau)=\mathbb E[(c-u)^3\mid\lambda].
}
$$

Then

$$
\boxed{R_0=P,}
$$

$$
\boxed{R_1=Pu,}
$$

$$
\boxed{R_2=P(u^2+\Theta),}
$$

$$
\boxed{R_3=P(u^3+3u\Theta+C_3).}
$$

---

## 7. Exact $P$ and $u$ equations

The zeroth moment gives

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0.
}
$$

Define

$$
\boxed{J=Pu.}
$$

The first moment gives

$$
\boxed{
\partial_\tau(Pu)
+\partial_\lambda[P(u^2+\Theta)]
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

## 8. Exact $\Theta$ equation

The second raw moment balance is

$$
\boxed{
\partial_\tau[P(u^2+\Theta)]
+\partial_\lambda[P(u^3+3u\Theta+C_3)]
=2P\mathbb E[c\ddot\lambda\mid\lambda].
}
$$

Define the conditional spacing-rate/acceleration covariance

$$
\boxed{
\Psi(\lambda,\tau)
=\operatorname{Cov}(c,\ddot\lambda\mid\lambda)
=\mathbb E[(c-u)\ddot\lambda\mid\lambda].
}
$$

Combining the zeroth, first, and second moment equations gives

$$
\boxed{
D_\tau\Theta
+2\Theta\,\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi.
}
$$

**Status:** EXACT general second-central-moment equation.

The shorter equation with zero right-hand side is **CONDITIONAL** on

$$
\boxed{\Psi=0.}
$$

That condition is not automatic in the LJ spacing chain because
$\ddot\lambda_i$ depends on neighbouring spacings. Setting either $C_3=0$ or
$\Psi=0$ without derivation is a closure assumption.

---

## 9. Exact nonlinear-LJ neighbour terms

For a bulk spacing,

$$
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}).
$$

Let $P_2^+(\lambda,\lambda',\tau)$ and
$P_2^-(\lambda,\lambda',\tau)$ be ordered central/right and central/left
neighbour joint densities with central marginal $P$. Define

$$
\boxed{
m_+
=\frac1P\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)\,d\lambda',
}
$$

$$
\boxed{
m_-
=\frac1P\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)\,d\lambda'.
}
$$

Then

$$
\boxed{
\mathcal A_{\rm bulk}=m_++m_- -2\phi'(\lambda).
}
$$

No neighbour-independence assumption is used.

For the $\Theta$ source, let
$F_2^+(\lambda,c,\lambda',\tau)$ and
$F_2^-(\lambda,c,\lambda',\tau)$ include the central spacing rate. Since
$\mathbb E[c-u\mid\lambda]=0$, the central force term drops from the covariance,
and

$$
\boxed{
\Psi_{\rm bulk}
=\frac1P\iint
(c-u)\phi'(\lambda')
[F_2^++F_2^-]\,dc\,d\lambda'.
}
$$

Boundary spacings require their own acceleration statistics if retained in the
one-point average.

**Status:** EXACT for the bulk chain.

---

## 10. Exact instantaneous shape identity for $P$

From the mean-flow balance,

$$
\frac1P\partial_\lambda(P\Theta)
=\partial_\lambda\Theta
+\Theta\partial_\lambda\ln P.
$$

Thus

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

At fixed $\tau$,

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

For a normalized nonabsorbing density, $C(\tau)$ is fixed by

$$
\boxed{
\int_0^\infty P(\lambda,\tau)\,d\lambda=1.
}
$$

**Status:** EXACT instantaneous shape representation under smoothness,
positivity, and moment-existence conditions. It is a reconstruction/constraint,
not an independent evolution law.

### 10.1 Degenerate case

At $\Theta=0$, the divided shape formula is invalid. The undivided moment
balance and continuity equation remain valid.

For the ideal homogeneous initial state,

$$
\boxed{
\lambda_i(0)=1,
\qquad
c_i(0)=0,
}
$$

so

$$
\boxed{
F_M(\lambda,c,0)=\delta(\lambda-1)\delta(c),
\qquad
P_M(\lambda,0)=\delta(\lambda-1).
}
$$

No artificial initial PDF width is introduced.

---

## 11. What $\Theta$ means mechanically

Exactly,

$$
\boxed{
\mathbb E[c^2\mid\lambda]=u^2+\Theta.
}
$$

Therefore $\Theta$ is the conditional spacing-rate dispersion that is lost by
reducing the phase-space state to $P$ alone. It can distinguish loading and
unloading populations at the same spacing.

However, because of Section 3.1,

$$
\frac12(u^2+\Theta)
$$

is only a local spacing-rate quadratic moment. It is not by itself the exact
chain kinetic-energy density. Exact total kinetic energy requires the metric
$\mathbf G_\lambda$ and cross-spacing correlations.

---

## 12. Same-force history dependence

Let $\tau_L$ and $\tau_U$ satisfy

$$
q(\tau_L)=q(\tau_U)=q^*,
$$

with $\dot q(\tau_L)>0$ and $\dot q(\tau_U)<0$. Define

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

then there is no memoryless map $\mathcal R_2=\mathcal S[q(\tau)]$ for that
trajectory. The current finite-chain numerical test gives this non-retracing.

Thus $(P,u,\Theta)$ is a **history-bearing reduced descriptor**. It is not
claimed to be an autonomous closed Markov state, because $C_3$, $\Psi$ and
neighbour joint states enter its exact evolution.

---

## 13. G1: mean spacing

Define

$$
\boxed{
\bar\lambda(\tau)=\int_0^\infty\lambda P(\lambda,\tau)\,d\lambda,
\qquad
\bar a(t)=a_0\bar\lambda(t/t_0).
}
$$

From continuity,

$$
\boxed{
\frac{d\bar\lambda}{d\tau}
=-[\lambda J]_0^\infty+\int_0^\infty J\,d\lambda.
}
$$

If the spacing-space boundary flux vanishes,

$$
\boxed{
\frac{d\bar\lambda}{d\tau}
=\int_0^\infty Pu\,d\lambda
=\mathbb E[c].
}
$$

**Status:** G1 DEFINITION plus EXACT moment identity.

---

## 14. G2: mean intrinsic configurational energy

Use the same nearest-neighbour energy that generates the microscopic forces:

$$
\boxed{
\Delta\phi(\lambda)=\phi(\lambda)-\phi(1).
}
$$

The mean intrinsic configurational energy per represented spacing is

$$
\boxed{
\bar U(\tau)
=U_{\rm ref}\int_0^\infty
\Delta\phi(\lambda)P(\lambda,\tau)\,d\lambda.
}
$$

If $P_M$ includes all $M$ spacings,

$$
\boxed{
V_{\rm phys}-MU_{\rm ref}\phi(1)=M\bar U.
}
$$

The rate is

$$
\boxed{
\frac1{U_{\rm ref}}\frac{d\bar U}{d\tau}
=-[\Delta\phi J]_0^\infty
+\int_0^\infty\phi'(\lambda)J\,d\lambda.
}
$$

For vanishing spacing-space boundary flux,

$$
\boxed{
\frac{d\bar U}{d\tau}
=U_{\rm ref}\int_0^\infty\phi'(\lambda)Pu\,d\lambda.
}
$$

This is configurational energy, not total mechanical energy.

**Status:** G2 DEFINITION and EXACT under the active nearest-neighbour chain.

A long-range/zeta energy must not be substituted into G2 while retaining
nearest-neighbour equations of motion and then called mechanically exact.

---

## 15. G3: irreversible hysteresis energy

The fixed observable is

$$
\boxed{
E_{\rm hyst}(t)=\int_0^t\dot D_{\rm irr}(t')\,dt',
\qquad
\dot D_{\rm irr}\ge0.
}
$$

The current conservative baseline contains no irreversible force, hence

$$
\boxed{
\dot D_{\rm irr}=0,
\qquad
E_{\rm hyst}=0
}
$$

for that baseline.

If a future physical irreversible node force $r_j^{\rm irr}$ is derived, then

$$
\frac{dE_{\rm mech}^*}{d\tau}
=q\dot x_M+\sum_jr_j^{\rm irr}\dot x_j.
$$

If

$$
\sum_jr_j^{\rm irr}\dot x_j\le0,
$$

define

$$
\boxed{
\dot D_{\rm irr}^*=-\sum_jr_j^{\rm irr}\dot x_j\ge0.
}
$$

Then

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}
=q\dot x_M-\dot D_{\rm irr}^*,
}
$$

and over one cycle

$$
\boxed{
W_{\rm ext}^{\rm cyc}
=\Delta E_{\rm mech}^{\rm cyc}+D_{\rm irr}^{\rm cyc}.
}
$$

Therefore same-force non-retracing or a transient loop does not by itself prove
irreversible dissipation.

**Status:** G3 observable DEFINED; physical $r_j^{\rm irr}$ is OPEN.

---

## 16. G4 and local first passage

The operational local instability threshold is

$$
\boxed{\phi''(\lambda_c)=0,}
$$

so

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

### 16.1 Nonabsorbing tail

$$
\boxed{
Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)\,d\lambda
}
$$

is an instantaneous instability-tail diagnostic, not cumulative first passage.

### 16.2 Exact finite empirical first passage

Define

$$
\boxed{
\tau_i^c=\inf\{\tau\ge0:\lambda_i(\tau)\ge\lambda_c\}.
}
$$

With

$$
\chi_i(\tau)=\mathbf1_{\{\tau<\tau_i^c\}},
$$

the local survivor fraction is

$$
\boxed{
S_M(\tau)=\frac1M\sum_i\chi_i(\tau),
}
$$

and

$$
\boxed{
F_{{\rm ci},M}^{\rm local}=1-S_M.
}
$$

Distributionally,

$$
\boxed{
-\frac{dS_M}{d\tau}
=\frac1M\sum_i\delta(\tau-\tau_i^c).
}
$$

**Status:** EXACT finite empirical definition.

### 16.3 Smooth kinetic absorbing boundary

Let $F_b(\lambda,c,\tau)$ be the survivor phase-space subdensity for
$0<\lambda<\lambda_c$. The interior equation is the same projected kinetic
transport. At the right boundary impose no inflow from the failed side:

$$
\boxed{
F_b(\lambda_c,c,\tau)=0
\quad\text{for incoming }c<0.
}
$$

Outgoing $c>0$ states give

$$
\boxed{
j_{\rm esc}(\tau)
=\int_0^\infty cF_b(\lambda_c^-,c,\tau)\,dc\ge0.
}
$$

Assuming no lower-boundary loss,

$$
\boxed{
S(\tau)
=\int_0^{\lambda_c}\int_{-\infty}^{\infty}F_b\,dc\,d\lambda,
}
$$

$$
\boxed{
\frac{dS}{d\tau}=-j_{\rm esc},
\qquad
F_{\rm ci}^{\rm local}=1-S.
}
$$

For $S>0$,

$$
\boxed{
h_\tau=\frac{j_{\rm esc}}S=-\frac{d}{d\tau}\ln S,
\qquad
h_t=\frac{h_\tau}{t_0}.
}
$$

**Status:** EXACT kinetic first-passage balance under the stated boundary
conditions.

### 16.4 Survivor-conditioned observables

The survivor spacing marginal

$$
P_b(\lambda,\tau)=\int F_b\,dc
$$

has mass $S$, not one. Define

$$
\boxed{
\widehat P_b=\frac{P_b}{S}.
}
$$

Then

$$
\boxed{
\bar\lambda_{\rm surv}
=\frac1S\int_0^{\lambda_c}\lambda P_b\,d\lambda,
}
$$

$$
\boxed{
\bar U_{\rm surv}
=\frac{U_{\rm ref}}S
\int_0^{\lambda_c}\Delta\phi(\lambda)P_b\,d\lambda.
}
$$

---

## 17. Local versus specimen probability

For one deterministic chain realization,

$$
\boxed{
\tau_{\rm spec}^c=\min_i\tau_i^c.
}
$$

The quantity $1-S_M$ is a local spatial first-passage fraction. It is not, in
general, specimen-to-specimen crack probability.

A specimen ensemble $\omega$ would require

$$
\boxed{
S_{\rm spec}(\tau)
=\Pr_\omega\left[\min_i\tau_i^c(\omega)>\tau\right].
}
$$

No independent-cell product is assumed.

**Status:** specimen-probability bridge OPEN.

Also, once $\tau_{\rm spec}^c$ is reached, the intact pre-crack chain is no
longer a physical post-initiation propagation model. Continued local crossings
past the first specimen event are mathematical diagnostics unless a post-crack
model is added.

---

## 18. Exact closure status

The full microscopic state

$$
\boxed{
\mathbf Z(\tau)
=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
}
$$

is a closed deterministic state under the chain equations.

From it one can generate exactly

$$
F_M,\ P_M,\ u,\ \Theta,\ C_3,\ \Psi,\ P_2^\pm,\ F_2^\pm,
\ \bar a,\ \bar U,\ \tau_i^c.
$$

The projected three-field system is exact but not autonomously closed:

$$
\boxed{
\{P,u,\Theta\}
\to
\{C_3,\Psi,P_2^\pm,F_2^\pm,\ldots\}.
}
$$

This is not a missing algebraic derivation. It is the exact hierarchy produced
by the reduced projection.

---

## 19. Assumption ledger

The active base model assumes:

1. one-dimensional normal motion;
2. identical masses after nondimensionalization;
3. nearest-neighbour generalized-LJ energy $V^*=\sum_i\phi(\lambda_i)$;
4. fixed left boundary and prescribed right-end normal force;
5. current calibration map $q=\sigma/E$;
6. spatial empirical probability over represented spacings;
7. smooth continuum fields only where moment equations are used;
8. sufficient velocity-space decay for moment integration by parts;
9. $P>0$ and $\Theta>0$ only where the divided shape formula is used;
10. $\phi''(\lambda_c)=0$ as the operational local initiation threshold.

It does not assume equilibrium statistics, named PDF families, neighbour
independence, stochastic Markov dynamics, white noise, Fokker--Planck,
Smoluchowski, empirical fatigue damage, viscous damping, independent FEM
elements, FCC geometry, or registry slip.

---

## 20. Final paper-level equation set

Microscopic bulk mechanics:

$$
\boxed{
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}).
}
$$

Empirical state:

$$
\boxed{
F_M=\frac1M\sum_i\delta(\lambda-\lambda_i)\delta(c-c_i).
}
$$

Projected transport:

$$
\boxed{
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0.
}
$$

Reduced fields:

$$
\boxed{
P=\int F\,dc,
\qquad
u=\mathbb E[c\mid\lambda],
\qquad
\Theta=\operatorname{Var}(c\mid\lambda).
}
$$

Continuity:

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0.
}
$$

Mean-flow balance:

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta).
}
$$

Density shape:

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta.
}
$$

Second-central-moment balance:

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi.
}
$$

Bulk LJ acceleration:

$$
\boxed{
\mathcal A=m_++m_- -2\phi'(\lambda).
}
$$

G1:

$$
\boxed{
\bar a=a_0\int\lambda P\,d\lambda.
}
$$

G2:

$$
\boxed{
\bar U=U_{\rm ref}\int[\phi(\lambda)-\phi(1)]P\,d\lambda.
}
$$

G3:

$$
\boxed{
E_{\rm hyst}=\int\dot D_{\rm irr}\,dt,
\qquad
\dot D_{\rm irr}\ge0,
}
$$

with $\dot D_{\rm irr}=0$ for the current conservative baseline and the
physical irreversible mechanism OPEN.

G4 / local first passage:

$$
\boxed{
\phi''(\lambda_c)=0,
\quad
\dot S=-j_{\rm esc},
\quad
F_{\rm ci}^{\rm local}=1-S.
}
$$

This is the complete active 1D mathematical formulation. The remaining
problems are physical rather than hidden algebra: irreversible mechanism,
laboratory-fatigue time-scale bridge, specimen-level probability bridge, and
experimental validation.
