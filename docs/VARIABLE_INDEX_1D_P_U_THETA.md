# Variable index — active 1D $P$–$u$–$\Theta$ formulation

This index accompanies `MASTER_1D_P_U_THETA_FORMULATION.md` and uses its
classification: MODEL, DEFINITION, EXACT, CONDITIONAL, OPEN.

## A. Physical scales and loading

| Symbol | Definition | Meaning | Unit / scale | Status |
|---|---|---|---|---|
| $t$ | physical time | laboratory time | s | DEFINITION |
| $\tau$ | $t/t_0$ | nondimensional time | 1 | DEFINITION |
| $t_0$ | $\sqrt{m_a a_0/(EA_0)}$ | atomic mechanical time scale | s | DEFINITION under calibration |
| $m_a$ | represented atomic/repeat mass | inertia scale | kg | input/calibration |
| $a_0$ | equilibrium spacing | normal length scale | m | input/calibration |
| $E$ | reference Young modulus | stress scale | Pa | empirical input |
| $A_0$ | effective 1D reference area | stress-force bridge | m$^2$ | calibration quantity |
| $U_{\rm ref}$ | $EA_0a_0$ | energy scale | J | DEFINITION |
| $F_{\rm ref}$ | $EA_0$ | force scale | N | DEFINITION |
| $\sigma_n(t)$ | applied normal stress | external loading | Pa | input history |
| $\sigma_m$ | mean stress | cycle mean | Pa | DEFINITION |
| $\sigma_a$ | stress amplitude | cycle amplitude | Pa | DEFINITION |
| $f$ | loading frequency | cycles per second | Hz | input |
| $\omega^*$ | $2\pi f t_0$ | nondimensional angular frequency | 1 | DEFINITION |
| $q(\tau)$ | $F_{\rm ext}/(EA_0)=\sigma_n/E$ | nondimensional end force | 1 | current mapping |

## B. Microscopic chain

| Symbol | Definition | Meaning | Unit | Status |
|---|---|---|---|---|
| $M$ | number of represented spacings | finite chain size | 1 | DEFINITION |
| $x_j$ | node position in units of $a_0$ | atom/node coordinate | 1 | DEFINITION |
| $\lambda_i$ | $x_i-x_{i-1}$ | normalized spacing | 1 | DEFINITION |
| $a_i$ | $a_0\lambda_i$ | physical spacing | m | DEFINITION |
| $c_i$ | $\dot\lambda_i=d\lambda_i/d\tau$ | spacing rate | 1 per $\tau$ | DEFINITION |
| $\ddot\lambda_i$ | $d^2\lambda_i/d\tau^2$ | spacing acceleration | 1 per $\tau^2$ | DEFINITION |
| $m,n$ | LJ exponents, $m>n>1$ | repulsive/attractive exponents | 1 | MODEL parameters |
| $\phi(\lambda)$ | normalized generalized-LJ energy | bond energy | 1 | MODEL |
| $\phi'$ | $d\phi/d\lambda$ | normalized bond force | 1 | EXACT derivative |
| $\phi''$ | $d^2\phi/d\lambda^2$ | normalized tangent stiffness | 1 | EXACT derivative |
| $V^*$ | $\sum_i\phi(\lambda_i)$ | chain configurational energy | 1 | EXACT under MODEL |
| $T^*$ | $\frac12\sum_j\dot x_j^2$ | chain kinetic energy | 1 | EXACT |
| $E_{\rm mech}^*$ | $T^*+V^*$ | total mechanical energy | 1 | EXACT |
| $\mathbf L$ | lower-triangular cumulative-sum matrix | $\boldsymbol x=\mathbf L\boldsymbol\lambda$ | 1 | DEFINITION |
| $\mathbf G_\lambda$ | $\mathbf L^T\mathbf L$ | spacing-coordinate mass metric | 1 | EXACT |

## C. Empirical and smooth phase-space state

| Symbol | Definition | Meaning | Status |
|---|---|---|---|
| $I$ | uniformly sampled spacing index | spatial counting random index | DEFINITION |
| $F_M(\lambda,c,\tau)$ | $M^{-1}\sum_i\delta(\lambda-\lambda_i)\delta(c-c_i)$ | finite empirical phase-space measure | DEFINITION |
| $P_M(\lambda,\tau)$ | $\int F_Mdc$ | finite empirical spacing measure | DEFINITION |
| $\mathcal G_M$ | $M^{-1}\sum_i\ddot\lambda_i\delta\delta$ | empirical acceleration flux | DEFINITION |
| $F(\lambda,c,\tau)$ | smooth representation of $F_M$ | projected phase-space density | DEFINITION |
| $P(\lambda,\tau)$ | $\int Fdc$ | spacing density | DEFINITION |
| $A(\lambda,c,\tau)$ | $E[\ddot\lambda\mid\lambda,c]$ | conditional phase-space acceleration | DEFINITION |
| $J(\lambda,\tau)$ | $Pu$ | spacing-space current | DEFINITION |

## D. Conditional moment fields

| Symbol | Definition | Meaning | Status |
|---|---|---|---|
| $u(\lambda,\tau)$ | $E[c\mid\lambda]$ | conditional mean spacing rate | DEFINITION |
| $\Theta(\lambda,\tau)$ | $\operatorname{Var}(c\mid\lambda)$ | conditional spacing-rate variance | DEFINITION |
| $C_3(\lambda,\tau)$ | $E[(c-u)^3\mid\lambda]$ | third conditional central rate moment | DEFINITION |
| $\mathcal A(\lambda,\tau)$ | $E[\ddot\lambda\mid\lambda]$ | one-point conditional acceleration | DEFINITION |
| $\Psi(\lambda,\tau)$ | $\operatorname{Cov}(c,\ddot\lambda\mid\lambda)$ | acceleration-covariance source in $\Theta$ balance | DEFINITION |
| $D_\tau$ | $\partial_\tau+u\partial_\lambda$ | material derivative in spacing space | DEFINITION |
| $R_r$ | $\int c^rFdc$ | raw velocity-moment density | DEFINITION |
| $B_r$ | $\int c^{r-1}AFdc$ | acceleration-moment source | DEFINITION |

Exact hierarchy:

$$
\boxed{
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r.
}
$$

Exact first three field balances:

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0,
}
$$

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta),
}
$$

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi.
}
$$

## E. Neighbour statistics

| Symbol | Definition | Meaning | Status |
|---|---|---|---|
| $P_2^+(\lambda,\lambda')$ | ordered central/right-neighbour spacing density | neighbour correlation | DEFINITION |
| $P_2^-(\lambda,\lambda')$ | ordered central/left-neighbour spacing density | neighbour correlation | DEFINITION |
| $F_2^+(\lambda,c,\lambda')$ | central spacing/rate + right-neighbour density | source for $\Psi$ | DEFINITION |
| $F_2^-(\lambda,c,\lambda')$ | central spacing/rate + left-neighbour density | source for $\Psi$ | DEFINITION |
| $m_+$ | $P^{-1}\int\phi'(\lambda')P_2^+d\lambda'$ | conditional right-neighbour force mean | EXACT |
| $m_-$ | $P^{-1}\int\phi'(\lambda')P_2^-d\lambda'$ | conditional left-neighbour force mean | EXACT |

Bulk conditional acceleration:

$$
\boxed{
\mathcal A=m_++m_- -2\phi'(\lambda).
}
$$

## F. Density-shape relation

| Symbol | Definition | Meaning |
|---|---|---|
| $\lambda_*$ | arbitrary reference point in a smooth positive-support interval | integration reference |
| $C(\tau)$ | time-dependent integration/normalization constant | density normalization |

Exact relation for $P>0$, $\Theta>0$:

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta,
}
$$

$$
\boxed{
P(\lambda,\tau)
=\frac{C(\tau)}{\Theta(\lambda,\tau)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\mathcal A-D_\tau u}{\Theta}\,d\eta
\right].
}
$$

## G. G1–G4 observables

| Symbol | Definition | Meaning | Status |
|---|---|---|---|
| $\bar\lambda$ | $\int\lambda P d\lambda$ | mean normalized spacing | G1 |
| $\bar a$ | $a_0\bar\lambda$ | mean physical spacing | G1 |
| $\Delta\phi$ | $\phi(\lambda)-\phi(1)$ | reference-subtracted intrinsic bond energy | DEFINITION |
| $\bar U$ | $U_{\rm ref}\int\Delta\phi P d\lambda$ | mean intrinsic configurational energy per spacing | G2 |
| $\dot D_{\rm irr}$ | nonnegative irreversible power once a physical mechanism is supplied | irreversible dissipation rate | G3 OPEN physically |
| $E_{\rm hyst}$ | $\int\dot D_{\rm irr}dt$ | accumulated irreversible energy | G3 |
| $S$ | intact/survivor probability mass | local survival | G4 |
| $F_{\rm ci}^{\rm local}$ | $1-S$ | cumulative local first-passage fraction | G4 |
| $h$ | $-d(\ln S)/dt$ | local initiation hazard | G4 |

The present conservative chain has

$$
\boxed{\dot D_{\rm irr}=0.}
$$

## H. First-passage variables

| Symbol | Definition | Meaning | Status |
|---|---|---|---|
| $\lambda_c$ | root of $\phi''(\lambda_c)=0$ | operational local stability threshold | MODEL-based mechanical criterion |
| $\tau_i^c$ | $\inf\{\tau:\lambda_i\ge\lambda_c\}$ | local first-hitting time | DEFINITION |
| $\chi_i$ | $1_{\{\tau<\tau_i^c\}}$ | local alive indicator | DEFINITION |
| $F_b$ | survivor phase-space subdensity for $\lambda<\lambda_c$ | intact phase-space population | DEFINITION |
| $P_b$ | $\int F_bdc$ | survivor spacing subdensity | DEFINITION |
| $j_{\rm esc}$ | $\int_0^\infty cF_b(\lambda_c^-,c)dc$ | outward first-passage flux | EXACT under absorbing boundary |
| $\widehat P_b$ | $P_b/S$ | normalized conditional survivor density | DEFINITION |
| $\tau_{\rm spec}^c$ | $\min_i\tau_i^c$ | first local initiation in one deterministic specimen realization | DEFINITION |
| $S_{\rm spec}$ | $\Pr_\omega[\tau_{\rm spec}^c>t]$ | specimen survival over a realization ensemble | OPEN bridge |

## I. Critical cautions

1. $\Theta$ is an exact conditional spacing-rate variance, not a fitted damage
   variable.
2. $\frac12(u^2+\Theta)$ is a spacing-rate quadratic moment, not by itself the
   exact chain kinetic-energy density; total kinetic energy uses
   $\mathbf G_\lambda$ and cross-correlations.
3. The exact $\Theta$ equation contains $2\Psi$. Omitting it requires the extra
   condition $\Psi=0$.
4. $(P,u,\Theta)$ is history-bearing but not autonomously closed.
5. $1-S_{\rm local}$ is not automatically specimen crack probability.
6. Same-force non-retracing is not equivalent to irreversible G3 dissipation.
7. The numerical proof-of-principle $\omega^*=0.02$ must not be read as a
   laboratory fatigue frequency without the $t_0$ mapping.
