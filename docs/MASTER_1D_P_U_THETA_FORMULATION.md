# Master derivation — active 1D $P$–$u$–$\Theta$ theory

This document is the canonical differential derivation of the active normal-only model.
이 문서는 현재 1D normal-only 모델의 **기준 미분형 유도문서**이다.

For navigation and definitions, use:

- `../README_EQUATION_INDEX.md`
- `EQUATION_SUMMARY_1D_P_U_THETA.md`
- `VARIABLE_INDEX_1D_P_U_THETA.md`
- `AUXILIARY_SYMBOL_INDEX_1D.md`
- `MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`

The conditional mean spacing-rate symbol is **always** $u$.

---

# 1. Physical scaling

Let $a_0$ be the equilibrium spacing, $m_a$ the represented microscopic mass,
$E$ the reference Young modulus, and $A_0$ the effective 1D reference area.
Define

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

The active stress-to-chain-force bridge is

$$
\boxed{
q(\tau)=\frac{F_{\rm ext}(t)}{EA_0}=\frac{\sigma_n(t)}{E}.
}
$$

For sinusoidal loading,

$$
\boxed{
\sigma_n(t)=\sigma_m+\sigma_a\sin(2\pi f t),
\qquad
\omega^*=2\pi f t_0.
}
$$

**Status:** DEFINITION under the current calibration bridge.

---

# 2. Microscopic 1D generalized-LJ chain

Take $M+1$ nodes $x_0,\ldots,x_M$ with

$$
\boxed{x_0(\tau)=0.}
$$

Define the $M$ normalized spacings

$$
\boxed{
\lambda_i=x_i-x_{i-1}>0,
\qquad i=1,\ldots,M,
}
$$

and physical spacings

$$
\boxed{a_i=a_0\lambda_i.}
$$

Dots below denote $d/d\tau$.

## 2.1 Active interaction energy

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

Its derivatives are

$$
\boxed{
\phi'(\lambda)
=\frac{\lambda^{-n-1}-\lambda^{-m-1}}{m-n},
}
$$

$$
\boxed{
\phi''(\lambda)
=\frac{(m+1)\lambda^{-m-2}-(n+1)\lambda^{-n-2}}{m-n}.
}
$$

Thus

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

**Status:** $\phi$ is MODEL; the derivatives are EXACT consequences of that model.

## 2.2 Closed node equations

For identical unit masses after nondimensionalization,

$$
\boxed{
\ddot x_j
=\phi'(\lambda_{j+1})-\phi'(\lambda_j),
\qquad j=1,\ldots,M-1,
}
$$

and at the loaded end

$$
\boxed{
\ddot x_M=-\phi'(\lambda_M)+q(\tau).
}
$$

Therefore the bulk spacing equation is

$$
\boxed{
\ddot\lambda_i
=\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}),
\qquad i=2,\ldots,M-1.
}
$$

The boundary spacings satisfy

$$
\boxed{
\ddot\lambda_1=\phi'(\lambda_2)-\phi'(\lambda_1),
}
$$

$$
\boxed{
\ddot\lambda_M
=q(\tau)+\phi'(\lambda_{M-1})-2\phi'(\lambda_M).
}
$$

**Status:** EXACT under the active chain model and boundary law.

---

# 3. Exact mechanical energy

Define

$$
\boxed{
T^*=\frac12\sum_{j=1}^{M}\dot x_j^2,
\qquad
E_{\rm mech}^*=T^*+V^*.
}
$$

Differentiating and using the equations of motion gives

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}=q(\tau)\dot x_M.
}
$$

Hence the active baseline is conservative except for prescribed external work.

## 3.1 Spacing-coordinate mass metric

Since

$$
x_j=\sum_{k=1}^{j}\lambda_k,
$$

write

$$
\boxed{
\boldsymbol x=\mathbf L\boldsymbol\lambda,
\qquad
\dot{\boldsymbol x}=\mathbf L\boldsymbol c,
\qquad
\boldsymbol c=\dot{\boldsymbol\lambda}.
}
$$

Then

$$
\boxed{
T^*=\frac12\boldsymbol c^T\mathbf G_\lambda\boldsymbol c,
\qquad
\mathbf G_\lambda=\mathbf L^T\mathbf L,
}
$$

with

$$
\boxed{
(G_\lambda)_{k\ell}=M-\max(k,\ell)+1.
}
$$

Thus local spacing rates are not independent unit-mass velocities.

---

# 4. Mechanics-generated empirical probability state

For one deterministic chain, define the empirical phase-space measure

$$
\boxed{
F_M(\lambda,c,\tau)
=\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(\tau)]
\delta[c-c_i(\tau)],
}
$$

where

$$
\boxed{c_i=\dot\lambda_i.}
$$

Its spacing marginal is

$$
\boxed{
P_M(\lambda,\tau)
=\int F_M(\lambda,c,\tau)dc
=\frac1M\sum_i\delta[\lambda-\lambda_i(\tau)].
}
$$

No named probability family is assumed.

---

# 5. Exact empirical phase-space transport

Define the empirical acceleration flux

$$
\boxed{
\mathcal G_M(\lambda,c,\tau)
=\frac1M\sum_i\ddot\lambda_i
\delta(\lambda-\lambda_i)
\delta(c-c_i).
}
$$

Differentiating $F_M$ distributionally gives

$$
\boxed{
\partial_\tau F_M
+\partial_\lambda(cF_M)
+\partial_c\mathcal G_M=0.
}
$$

For a smooth projected representation define

$$
\boxed{
A(\lambda,c,\tau)
=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c].
}
$$

Then

$$
\boxed{
\mathcal G=AF
}
$$

and

$$
\boxed{
\partial_\tau F
+\partial_\lambda(cF)
+\partial_c(AF)=0.
}
$$

This is an exact projected identity when $A$ is the true mechanics-generated conditional acceleration. It is not an autonomous one-point constitutive closure.

---

# 6. Exact raw moment hierarchy

Define

$$
\boxed{
R_r(\lambda,\tau)=\int_{-\infty}^{\infty}c^rF(\lambda,c,\tau)dc,
\qquad r=0,1,2,\ldots
}
$$

and for $r\ge1$

$$
\boxed{
B_r(\lambda,\tau)
=\int c^{r-1}A(\lambda,c,\tau)F(\lambda,c,\tau)dc.
}
$$

Multiply the projected phase-space equation by $c^r$ and integrate over $c$.
Assuming the required moments exist and the velocity-boundary terms vanish,

$$
\boxed{
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r.
}
$$

**Status:** EXACT.

---

# 7. Reduced fields $P,u,\Theta$

Define

$$
\boxed{
P(\lambda,\tau)=R_0=\int Fdc,
}
$$

$$
\boxed{
 u(\lambda,\tau)=\mathbb E[c\mid\lambda]
=\frac{R_1}{P},
}
$$

$$
\boxed{
\Theta(\lambda,\tau)
=\operatorname{Var}(c\mid\lambda)
=\mathbb E[(c-u)^2\mid\lambda],
}
$$

$$
\boxed{
C_3(\lambda,\tau)=\mathbb E[(c-u)^3\mid\lambda],
}
$$

$$
\boxed{
\mathcal A(\lambda,\tau)
=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda],
}
$$

and

$$
\boxed{
\Psi(\lambda,\tau)
=\operatorname{Cov}(c,\ddot\lambda\mid\lambda)
=\mathbb E[(c-u)\ddot\lambda\mid\lambda].
}
$$

The first raw moments are

$$
\boxed{
R_0=P,
\qquad
R_1=Pu,
\qquad
R_2=P(u^2+\Theta),
}
$$

$$
\boxed{
R_3=P(u^3+3u\Theta+C_3).
}
$$

---

# 8. Exact continuity and mean-flow equations

The $r=0$ equation gives

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0.
}
$$

Define the spacing-space probability current

$$
\boxed{J=Pu.}
$$

The $r=1$ equation gives

$$
\boxed{
\partial_\tau(Pu)
+\partial_\lambda[P(u^2+\Theta)]
=P\mathcal A.
}
$$

Using continuity,

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta),
}
$$

where

$$
\boxed{
D_\tau=\partial_\tau+u\partial_\lambda.
}
$$

**Status:** EXACT.

---

# 9. Correct exact $\Theta$ equation

The $r=2$ raw equation is

$$
\boxed{
\partial_\tau[P(u^2+\Theta)]
+\partial_\lambda[P(u^3+3u\Theta+C_3)]
=2P\mathbb E[c\ddot\lambda\mid\lambda].
}
$$

Because

$$
\mathbb E[c\ddot\lambda\mid\lambda]
=u\mathcal A+\Psi,
$$

combining the $r=0,1,2$ balances gives

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi.
}
$$

The shorter form

$$
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)=0
$$

is only valid under the extra condition

$$
\boxed{\Psi=0.}
$$

For the actual spatial LJ chain, $\Psi=0$ is not automatic because $\ddot\lambda_i$ depends on neighbouring spacings.

---

# 10. Exact neighbour-conditioned acceleration

For a bulk spacing,

$$
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}).
$$

Let $P_2^+(\lambda,\lambda',\tau)$ and $P_2^-(\lambda,\lambda',\tau)$ be ordered central/right and central/left neighbour joint densities. Define

$$
\boxed{
m_+
=\frac1P\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)d\lambda',
}
$$

$$
\boxed{
m_-
=\frac1P\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)d\lambda'.
}
$$

Then

$$
\boxed{
\mathcal A_{\rm bulk}=m_++m_- -2\phi'(\lambda).
}
$$

For $\Psi$, let $F_2^+(\lambda,c,\lambda',\tau)$ and $F_2^-(\lambda,c,\lambda',\tau)$ include the central spacing rate. Since

$$
\mathbb E[c-u\mid\lambda]=0,
$$

the central-force term vanishes inside the covariance and

$$
\boxed{
\Psi_{\rm bulk}
=\frac1P\iint
(c-u)\phi'(\lambda')
[F_2^++F_2^-]dc\,d\lambda'.
}
$$

No neighbour-independence approximation is used.

---

# 11. Exact instantaneous density-shape relation

From the mean-flow equation,

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

On a smooth interval where $P>0$ and $\Theta>0$,

$$
\boxed{
\partial_\lambda\ln P
=\frac{\mathcal A-D_\tau u}{\Theta}
-\partial_\lambda\ln\Theta.
}
$$

Integrating in $\lambda$ gives

$$
\boxed{
P(\lambda,\tau)
=\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}
{\Theta(\eta,\tau)}d\eta
\right].
}
$$

The normalization factor $\mathcal N_P(\tau)$ is determined by

$$
\boxed{
\int_0^\infty P(\lambda,\tau)d\lambda=1.
}
$$

At $\Theta=0$ the divided form is invalid; the undivided moment equations remain valid.

---

# 12. Exact time-integral representations

The continuity equation integrates to

$$
\boxed{
P(\lambda,\tau)
=P_0(\lambda)
-\partial_\lambda
\int_{\tau_0}^{\tau}P(\lambda,s)u(\lambda,s)ds.
}
$$

The first moment integrates to

$$
\boxed{
Pu
=P_0u_0
-\partial_\lambda\int_{\tau_0}^{\tau}P(u^2+\Theta)ds
+\int_{\tau_0}^{\tau}P\mathcal A ds.
}
$$

The second raw moment integrates to

$$
\boxed{
\begin{aligned}
P(u^2+\Theta)(\lambda,\tau)
={}&P_0(u_0^2+\Theta_0)(\lambda)\\
&-\partial_\lambda\int_{\tau_0}^{\tau}
P(u^3+3u\Theta+C_3)ds\\
&+2\int_{\tau_0}^{\tau}P(u\mathcal A+\Psi)ds.
\end{aligned}
}
$$

The full push-forward and characteristic integral forms are given in `MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`.

---

# 13. Meaning of $\Theta$

Exactly,

$$
\boxed{
\mathbb E[c^2\mid\lambda]=u^2+\Theta.
}
$$

Thus $\Theta$ is the conditional spacing-rate dispersion lost when reducing $F$ to $P$ alone.
However,

$$
\frac12(u^2+\Theta)
$$

is not the full chain kinetic-energy density because the exact kinetic energy uses $\mathbf G_\lambda$ and cross-spacing rate correlations.

---

# 14. Same-force history dependence

At two times $\tau_L$ and $\tau_U$ satisfying

$$
q(\tau_L)=q(\tau_U)=q^*,
$$

with

$$
\dot q(\tau_L)>0,
\qquad
\dot q(\tau_U)<0,
$$

define

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

then no memoryless scalar-load map

$$
\mathcal R_2=\mathcal S[q(\tau)]
$$

exists for that trajectory. This establishes dynamic history dependence, not irreversible dissipation.

---

# 15. G1 — mean spacing

Define

$$
\boxed{
\bar\lambda(\tau)=\int_0^\infty\lambda P(\lambda,\tau)d\lambda,
}
$$

$$
\boxed{
\bar a(t)=a_0\bar\lambda(t/t_0).
}
$$

Using continuity,

$$
\boxed{
\frac{d\bar\lambda}{d\tau}
=-[\lambda J]_0^\infty+\int_0^\infty Jd\lambda.
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

---

# 16. G2 — mean intrinsic configurational energy

Define the equilibrium-subtracted interaction energy

$$
\boxed{
\Delta\phi(\lambda)=\phi(\lambda)-\phi(1).
}
$$

Then

$$
\boxed{
\bar U(\tau)
=U_{\rm ref}\int_0^\infty\Delta\phi(\lambda)P(\lambda,\tau)d\lambda.
}
$$

For the empirical density over all $M$ spacings,

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
+\int_0^\infty\phi'(\lambda)Jd\lambda.
}
$$

If boundary flux vanishes,

$$
\boxed{
\frac{d\bar U}{d\tau}
=U_{\rm ref}\int_0^\infty\phi'(\lambda)Pu\,d\lambda.
}
$$

This is configurational energy, not total mechanical energy.

---

# 17. G3 — irreversible hysteresis energy

The fixed observable is

$$
\boxed{
E_{\rm hyst}(t)=\int_0^t\dot D_{\rm irr}(t')dt',
\qquad
\dot D_{\rm irr}\ge0.
}
$$

For the current conservative baseline,

$$
\boxed{
\dot D_{\rm irr}=0,
\qquad
E_{\rm hyst}=0.
}
$$

If a future physical irreversible node force $r_j^{\rm irr}$ is derived, define

$$
\boxed{
\dot D_{\rm irr}^*=-\sum_jr_j^{\rm irr}\dot x_j\ge0
}
$$

whenever the irreversible force performs nonpositive mechanical power. Then

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
=\Delta E_{\rm mech}^{\rm cyc}+D_{\rm irr}^{\rm cyc}.
}
$$

Thus same-force non-retracing alone does not prove $\dot D_{\rm irr}>0$.

---

# 18. G4 — mechanical first passage

The operational local instability threshold is defined by

$$
\boxed{
\phi''(\lambda_c)=0.
}
$$

For the active generalized-LJ form,

$$
\boxed{
\lambda_c
=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
}
$$

For each spacing,

$$
\boxed{
\tau_i^c
=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}.
}
$$

The instantaneous nonabsorbing tail

$$
\boxed{
Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)d\lambda
}
$$

is not cumulative first passage.

For finite empirical trajectories define

$$
\boxed{
\chi_i(\tau)=\mathbf1_{\{\tau<\tau_i^c\}},
\qquad
S_M(\tau)=\frac1M\sum_i\chi_i(\tau).
}
$$

Then

$$
\boxed{
F_{{\rm ci},M}^{\rm local}=1-S_M.
}
$$

For a smooth survivor phase-space subdensity $F_b$ on $0<\lambda<\lambda_c$, impose no inflow from the failed side:

$$
\boxed{
F_b(\lambda_c,c,\tau)=0
\qquad(c<0).
}
$$

The outward first-passage flux is

$$
\boxed{
j_{\rm esc}(\tau)
=\int_0^\infty cF_b(\lambda_c^-,c,\tau)dc.
}
$$

Define

$$
\boxed{
S(\tau)=\int_0^{\lambda_c}\int_{-\infty}^{\infty}F_b\,dc\,d\lambda.
}
$$

Then

$$
\boxed{
\dot S=-j_{\rm esc},
\qquad
F_{\rm ci}^{\rm local}=1-S,
}
$$

and for $S>0$,

$$
\boxed{
h_\tau=\frac{j_{\rm esc}}S=-\frac{d}{d\tau}\ln S,
\qquad
h_t=\frac{h_\tau}{t_0}.
}
$$

The normalized survivor density is

$$
\boxed{
\widehat P_b=\frac{P_b}{S},
\qquad
P_b(\lambda,\tau)=\int F_bdc.
}
$$

---

# 19. Local versus specimen probability

For one realization,

$$
\boxed{
\tau_{\rm spec}^c=\min_i\tau_i^c.
}
$$

The finite local first-passage fraction $1-S_M$ is not automatically a specimen-to-specimen probability.
With an ensemble measure $\mu_0$, the exact specimen survival formula is

$$
\boxed{
S_{\rm spec}(\tau)
=\int
\mathbf1\left[
\max_i\sup_{s\in[\tau_0,\tau]}
\Lambda_i(s;\Gamma_0)<\lambda_c
\right]
\mu_0(d\Gamma_0).
}
$$

Thus the mathematical survival formula exists; the physical construction and calibration of $\mu_0$ and its correlation scale remain open.

---

# 20. Exact closure status

The full microscopic state

$$
\boxed{
\Gamma(\tau)
=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
}
$$

is closed under the finite-chain ODEs.

The projected fields satisfy

$$
\boxed{
\{P,u,\Theta\}
\longrightarrow
\{C_3,\Psi,P_2^\pm,F_2^\pm,\ldots\},
}
$$

so the three-field PDE is exact but hierarchical rather than autonomous.
However, the closed microscopic flow yields exact push-forward integrals for all these quantities; see `MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`.

Therefore

$$
\boxed{
\text{lack of autonomous three-field closure}
\neq
\text{lack of exact mathematical solution representation}.
}
$$

---

# 21. Active assumption ledger

The active model assumes:

1. one-dimensional normal motion;
2. identical masses after nondimensionalization;
3. nearest-neighbour generalized-LJ energy $V^*=\sum_i\phi(\lambda_i)$;
4. fixed left boundary and prescribed right-end normal force;
5. current calibration map $q=\sigma_n/E$;
6. spatial empirical probability over represented spacings;
7. smooth continuum fields only where moment equations are used;
8. sufficient velocity-space decay for moment integration by parts;
9. $P>0$ and $\Theta>0$ only where the divided density-shape formula is used;
10. $\phi''(\lambda_c)=0$ as the operational local initiation threshold.

The active theory does **not** assume Boltzmann/Gibbs equilibrium, Gaussian/Weibull spacing or life distributions, neighbour independence, stochastic Markov dynamics, white noise, Fokker--Planck, Smoluchowski, empirical fatigue damage, arbitrary viscous damping, independent FEM-element probabilities, FCC geometry, or registry slip.

---

# 22. Final governing set

Microscopic mechanics:

$$
\boxed{
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}).
}
$$

Empirical phase-space state:

$$
\boxed{
F_M
=\frac1M\sum_i
\delta(\lambda-\lambda_i)\delta(c-c_i).
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
P=\int Fdc,
\qquad
u_{\rm reserved}\;\text{is not used},
\qquad
u_{\rm reserved}\neq u,
}
$$

and the active mean rate is

$$
\boxed{
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

Variance balance:

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi.
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
\bar U
=U_{\rm ref}\int[\phi(\lambda)-\phi(1)]P\,d\lambda.
}
$$

G3:

$$
\boxed{
E_{\rm hyst}=\int\dot D_{\rm irr}dt,
\qquad
\dot D_{\rm irr}\ge0,
}
$$

with $\dot D_{\rm irr}=0$ in the present conservative baseline.

G4:

$$
\boxed{
\phi''(\lambda_c)=0,
\qquad
\dot S=-j_{\rm esc},
\qquad
F_{\rm ci}^{\rm local}=1-S.
}
$$

This is the active mathematical theory. The remaining blockers are physical: an irreversible G3 mechanism, a physical specimen measure/correlation scale, microscopic-to-laboratory fatigue time-scale bridging, and experimental validation.
