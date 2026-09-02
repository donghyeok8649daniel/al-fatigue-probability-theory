# Milestone 20 — Exact reduced moment hierarchy and the real closure boundary

## 0. Scope

The mainline remains the reduced fatigue state

$$
\mathbf q=(a,s),
$$

with the current intrinsic lattice energy

$$
U_0(a,s)
=\sum_{k\ge1}\sum_{p\in\mathbb Z}
v_{m,n}\!\left(\sqrt{k^2a^2+(pb+s)^2}\right).
$$

No FCC lattice model, Boltzmann distribution, Gaussian/Weibull family,
Fokker--Planck law, Markov bath, or ad hoc spatial spring coupling is introduced
here.

Milestone 18 derived the exact smooth shape identity

$$
\boldsymbol\Theta\nabla\ln P
=\boldsymbol{\mathcal A}-D_t\mathbf u-\nabla\cdot\boldsymbol\Theta.
$$

Milestone 19 showed how a declared reduced mechanics can supply
$\boldsymbol{\mathcal A}$.  The present milestone identifies exactly what still
prevents a predictive closure.

---

## 1. Reduced phase-space empirical measure

For represented states $\alpha=1,\dots,M$, define

$$
\mathbf q_\alpha=(a_\alpha,s_\alpha),
\qquad
\mathbf v_\alpha=\dot{\mathbf q}_\alpha.
$$

The exact empirical phase-space measure is

$$
\boxed{
F_M(\mathbf q,\mathbf v,t)
=\frac1M\sum_{\alpha=1}^M
\delta(\mathbf q-\mathbf q_\alpha)
\delta(\mathbf v-\mathbf v_\alpha).
}
$$

This is only the same reduced state augmented by its velocities.  It is not a
many-body statistical-mechanics replacement for the project.

Suppose the declared reduced mechanics gives

$$
\dot{\mathbf q}=\mathbf v,
\qquad
\dot{\mathbf v}=\mathbf B(\mathbf q,\mathbf v,t).
$$

Direct differentiation of the empirical measure gives

$$
\boxed{
\partial_tF_M
+\nabla_{\mathbf q}\cdot(\mathbf vF_M)
+\nabla_{\mathbf v}\cdot(\mathbf B F_M)=0.
}
$$

Classification: **EXACT EMPIRICAL TRANSPORT IDENTITY** for the declared reduced
mechanics.  No stochastic process or equilibrium distribution has been
postulated.

The state density used by G1--G4 is the velocity marginal

$$
\boxed{
P(\mathbf q,t)=\int F(\mathbf q,\mathbf v,t)\,d\mathbf v.
}
$$

---

## 2. First three conditional velocity moments

Define

$$
\mathbf u(\mathbf q,t)
=\mathbb E[\mathbf v\mid\mathbf q],
$$

$$
\boxed{
\Theta_{ij}
=\mathbb E[(v_i-u_i)(v_j-u_j)\mid\mathbf q]
}
$$

and the third central moment

$$
\boxed{
C_{ijk}
=\mathbb E[(v_i-u_i)(v_j-u_j)(v_k-u_k)\mid\mathbf q].
}
$$

The first two exact smooth moment equations are

$$
\boxed{
\partial_tP+\partial_j(Pu_j)=0
}
$$

and

$$
\boxed{
D_tu_i
=\mathcal A_i
-\frac1P\partial_j(P\Theta_{ij}),
}
$$

where

$$
\mathcal A_i=\mathbb E[B_i\mid\mathbf q].
$$

Expanding the second equation yields precisely

$$
\boxed{
\boldsymbol\Theta\nabla\ln P
=\boldsymbol{\mathcal A}-D_t\mathbf u
-\nabla\cdot\boldsymbol\Theta.
}
$$

Therefore the Theta density-shape equation and the first velocity-moment
balance are **the same information written in different forms**.

---

## 3. Important correction: the Theta shape formula is not an independent PDE

The representation

$$
\nabla\ln P
=\boldsymbol\Theta^{-1}
[\boldsymbol{\mathcal A}-D_t\mathbf u-\nabla\cdot\boldsymbol\Theta]
$$

is exact wherever $P>0$ and $\boldsymbol\Theta$ is invertible.  It is very useful
for:

1. reconstructing $P$ from independently measured mechanical moment fields;
2. checking compatibility by the two-dimensional curl condition;
3. testing numerical consistency of a proposed closure;
4. connecting the mechanically generated state distribution to G1 and G2.

But if $D_t\mathbf u$ is first calculated from the same momentum equation and
then substituted into the shape formula, the result is an identity.  It does
not generate a new predictive equation for $P$.

Thus the earlier one-dimensional reconstruction success validates the exact
moment identity and its numerical extraction from mechanics.  It does **not**
by itself prove that $P$ has been closed without higher moment information.

---

## 4. Exact evolution equation for Theta

For the general reduced acceleration $\mathbf B(\mathbf q,\mathbf v,t)$, define
its conditional mean $\boldsymbol{\mathcal A}$ and the symmetric
velocity--acceleration covariance source

$$
\Xi_{ij}
=\mathbb E[
(v_i-u_i)(B_j-\mathcal A_j)
+(v_j-u_j)(B_i-\mathcal A_i)
\mid\mathbf q].
$$

Then the exact conditional-covariance equation is

$$
\boxed{
D_t\Theta_{ij}
+\Theta_{kj}\partial_k u_i
+\Theta_{ik}\partial_k u_j
+\frac1P\partial_k(PC_{ijk})
=\Xi_{ij}.
}
$$

This exposes the next member of the hierarchy directly.

### Constant-metric conservative subcase

If the declared reduced mechanics has a constant generalized mass metric and

$$
\mathbf B(\mathbf q,t)
=\mathbf G^{-1}[\mathbf Q(t)-\nabla U_0(\mathbf q)],
$$

then $\mathbf B$ has no conditional velocity fluctuation at fixed $\mathbf q$,
so

$$
\boxed{\boldsymbol\Xi=0.}
$$

Therefore

$$
\boxed{
D_t\boldsymbol\Theta
+(\nabla\mathbf u)\boldsymbol\Theta
+\boldsymbol\Theta(\nabla\mathbf u)^T
+\frac1P\nabla\cdot(P\mathbf C)=0.
}
$$

The third central moment $\mathbf C$ remains.  Hence even this simplest
conservative reduced mechanics does **not** close exactly at $P,\mathbf u,
\boldsymbol\Theta$.

Setting $\mathbf C=0$ would be a closure assumption corresponding to a special
conditional velocity-shape restriction.  It is not adopted here.

---

## 5. One-dimensional reduction

In one coordinate the exact covariance equation becomes

$$
\boxed{
D_t\Theta
+2\Theta\,\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=\Xi.
}
$$

For acceleration independent of velocity, $\Xi=0$.

This confirms that the old one-dimensional $\Theta$ shape law also belongs to a
moment hierarchy: predicting $\Theta$ exactly requires the third conditional
velocity moment unless the full reduced phase-space state is evolved.

---

## 6. Exact non-closure statement for the current mainline

The current $U_0(a,s)$ determines an intrinsic energy landscape and its force
components

$$
-\partial_aU_0,
\qquad
-\partial_sU_0.
$$

It does **not**, by itself, determine:

1. the physical ensemble/initial spread in $(a,s,\dot a,\dot s)$;
2. spatial coupling among several representative states;
3. the third central velocity moment $C_{ijk}$;
4. an irreversible mechanism for G3.

Therefore the exact hierarchy cannot be closed from $U_0$ alone without
additional physical information.

A tempting extension such as

$$
\ddot{\mathbf q}_i
\stackrel{?}{=}
\nabla U_0(\mathbf q_{i+1})
-2\nabla U_0(\mathbf q_i)
+\nabla U_0(\mathbf q_{i-1})
$$

is **not derived** from the current model and must not be introduced merely by
analogy with the historical one-dimensional spacing chain.  The current
$U_0(a,s)$ is a local reduced state energy, not a demonstrated interaction
energy between neighboring probability patches.

---

## 7. The exact no-closure route

There is one way to avoid a moment closure completely: evolve the finite
reduced phase-space ensemble itself,

$$
\{\mathbf q_\alpha(t),\mathbf v_\alpha(t)\}_{\alpha=1}^M,
$$

under the declared reduced mechanics and form

$$
P_M(a,s,t)
=\frac1M\sum_\alpha
\delta[a-a_\alpha(t)]\delta[s-s_\alpha(t)].
$$

Equivalently one can evolve the reduced phase-space transport equation for
$F$.  The trajectory ensemble is usually the simpler numerical realization.

This route uses no Gaussian, Weibull, Boltzmann, Fokker--Planck, or finite-
moment closure.  However, it still requires the **physical origin of the
initial/state diversity** to be specified.  If every represented state has
identical initial $(\mathbf q,\mathbf v)$ and obeys identical deterministic
mechanics, then

$$
P_M(a,s,t)=\delta[a-a(t)]\delta[s-s(t)]
$$

and $\boldsymbol\Theta=0$ for all time.

Thus broadening cannot appear from nowhere.

---

## 8. Energy meaning of Theta

For a constant generalized mass metric $\mathbf G$, the conditional kinetic
energy is

$$
\boxed{
\mathbb E[T\mid\mathbf q]
=\frac12\left[
\mathbf u^T\mathbf G\mathbf u
+\operatorname{tr}(\mathbf G\boldsymbol\Theta)
\right].
}
$$

Hence $\boldsymbol\Theta$ is not merely a probability-width parameter.  It
represents real conditional kinetic-energy dispersion in the reduced state.

G2 remains

$$
\boxed{
\bar U(t)
=\iint\Delta U_0(a,s)P(a,s,t)\,da\,ds.
}
$$

For a conservative reduced ensemble with no boundary flux, the ensemble total
mechanical energy

$$
\bar E_{\rm mech}
=\iint P\left[
\frac12\mathbf u^T\mathbf G\mathbf u
+\frac12\operatorname{tr}(\mathbf G\boldsymbol\Theta)
+U_0
\right]da\,ds
$$

obeys the external-work balance

$$
\boxed{
\frac{d\bar E_{\rm mech}}{dt}
=\iint P\,\mathbf Q\cdot\mathbf u\,da\,ds
}
$$

when the stated constant-metric conservative assumptions hold.

This provides a useful future numerical verification independently of G2.

---

## 9. Consequence for the four governing observables

Nothing here changes the official observables:

$$
G1:\quad \bar a=\iint aP\,da\,ds,
$$

$$
G2:\quad \bar U=\iint\Delta U_0P\,da\,ds,
$$

$$
G3:\quad E_{\rm hyst}=\int\dot D_{\rm irr}\,dt,
$$

$$
G4:\quad \iint P\,da\,ds=1
$$

(or survival mass under an absorbing crack boundary).

The correction is only about **how $P$ is predicted**.  G1, G2, and G4 can be
computed once a valid mechanically generated $P$ exists.  G3 still requires a
physically justified irreversible mechanism and is not generated by purely
conservative moment transport.

---

## 10. Mainline decision after this milestone

The scientifically safe hierarchy is now

$$
\boxed{
U_0(a,s)+\text{declared reduced mechanics}
\rightarrow
F(a,s,v_a,v_s,t)
\rightarrow
P(a,s,t)
\rightarrow
\{\bar a,\bar U,\langle z\rangle,S\}.
}
$$

The $\Theta$ equation remains central as a diagnostic and reduced description,
but it must not be advertised as a standalone closed probability solver.

For computation, the preferred exact/no-closure realization is a finite
ensemble of reduced trajectories.  The unresolved physical task is to identify
the source of distinct represented states/initial conditions and, separately,
the irreversible mechanism needed for G3 and permanent plasticity.
