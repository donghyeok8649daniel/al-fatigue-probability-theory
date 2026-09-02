# Milestone 24 — $P$–$u$–$\Theta$ as a history-bearing reduced state

## 0. Scope

This milestone keeps the active 1D normal generalized-LJ chain and the exact
spacing-space moment structure.  It introduces no Boltzmann law, Fokker--Planck
closure, damping, empirical damage variable, or prescribed PDF family.

The question is:

> At the same externally applied normal force, can the mechanically generated
> state differ between the loading and unloading passages of a cycle?

The answer must be separated into two statements:

1. **state non-retracing / history dependence** of $(P,u,\Theta)$;
2. **irreversible dissipation** required by G3.

They are not the same statement.

---

## 1. Exact state definitions

For local spacing $a_i(t)$ and spacing velocity $c_i=\dot a_i$, define the
finite empirical phase-space measure

$$
F_M(a,c,t)=\frac1M\sum_{i=1}^M
\delta[a-a_i(t)]\delta[c-c_i(t)].
$$

Its spacing marginal is

$$
\boxed{
P_M(a,t)=\int F_M(a,c,t)\,dc
=\frac1M\sum_i\delta[a-a_i(t)].
}
$$

For the smooth one-point representation,

$$
\boxed{u(a,t)=\mathbb E[c\mid a]}
$$

and

$$
\boxed{
\Theta(a,t)=\operatorname{Var}(c\mid a)
=\mathbb E[(c-u)^2\mid a].
}
$$

Thus $\Theta$ is not a fitted memory coefficient.  It is an exact conditional
second central velocity moment of the mechanically generated state.

---

## 2. Exact moment relation retained

The exact first two phase-space moment equations give

$$
\partial_tP+\partial_a(Pu)=0
$$

and

$$
D_tu
=\mathcal A-\frac1P\partial_a(P\Theta),
$$

where

$$
D_tu=\partial_tu+u\partial_a u,
\qquad
\mathcal A(a,t)=\mathbb E[\ddot a_i\mid a_i=a].
$$

Rearranging,

$$
\boxed{
\Theta\,\partial_a\ln P
=\mathcal A-D_tu-\partial_a\Theta.
}
$$

Where $P>0$ and $\Theta>0$ are smooth,

$$
\boxed{
\partial_a\ln P
=\frac{\mathcal A-D_tu}{\Theta}
-\partial_a\ln\Theta.
}
$$

This is an exact shape identity, not an independent probability closure.

---

## 3. Meaning of same-force non-retracing

Let $t_L$ and $t_U$ be two times in one cyclic forcing history such that

$$
Q(t_L)=Q(t_U)=Q^*,
$$

with

$$
\dot Q(t_L)>0,
\qquad
\dot Q(t_U)<0.
$$

If

$$
\boxed{P(a,t_L)\ne P(a,t_U)}
$$

or

$$
\boxed{u(a,t_L)\ne u(a,t_U)}
$$

or

$$
\boxed{\Theta(a,t_L)\ne\Theta(a,t_U),}
$$

then the reduced state is not a single-valued function of the instantaneous
external force alone.

Equivalently,

$$
(P,u,\Theta)\ne\mathcal S[Q(t)]
$$

for any memoryless scalar map $\mathcal S$.

This is a precise reduced-state notion of **history dependence**.

---

## 4. Why $\Theta$ can carry history

At a fixed spacing value $a$, two ensembles can have the same local spacing
population but different incoming/outgoing velocity populations.  Then their
conditional mean velocities and conditional velocity spreads differ:

$$
u_L(a)\ne u_U(a),
\qquad
\Theta_L(a)\ne\Theta_U(a).
$$

Therefore $P(a)$ alone need not distinguish all load/unload states.  The
phase-space information retained by $(u,\Theta)$ can distinguish states whose
instantaneous spacing statistics are similar.

The conditional kinetic contribution is

$$
\boxed{
\mathbb E\!\left[\frac12c^2\middle|a\right]
=\frac12\left[u(a,t)^2+\Theta(a,t)\right]
}
$$

in the current unit-mass normalized 1D chain.  Hence $\Theta$ also has a direct
mechanical energy interpretation.

---

## 5. This does not yet prove irreversible fatigue hysteresis

The active normal chain is conservative.  With external forcing,

$$
\boxed{
\frac{dE_{\rm mech}}{dt}=Q(t)\,\dot x_{\rm end}(t).
}
$$

A non-retracing force--state loop can therefore arise from reversible inertia,
wave propagation, and phase-space redistribution even when no irreversible
mechanism is present.

Consequently

$$
(P_L,u_L,\Theta_L)\ne(P_U,u_U,\Theta_U)
$$

is sufficient to establish **history-dependent reduced state**, but it is not
sufficient to establish

$$
\dot D_{\rm irr}>0.
$$

G3 remains a separate physical requirement.

The correct hierarchy is

$$
\boxed{
\text{same-force non-retracing}
\not\Rightarrow
\text{irreversible dissipation}.
}
$$

---

## 6. Exact evolution of $\Theta$ and the remaining hierarchy

The second central moment obeys

$$
\boxed{
D_t\Theta
+2\Theta\,\partial_a u
+\frac1P\partial_a(PC_3)=0,
}
$$

with

$$
C_3(a,t)=\mathbb E[(c-u)^3\mid a].
$$

Thus $\Theta$ is mathematically defined and dynamically meaningful, but its
standalone predictive evolution still belongs to the exact moment hierarchy.
Setting $C_3=0$ would be a new closure assumption and is not adopted here.

For direct deterministic simulation this is not an obstacle: the microscopic
chain generates $F_M$, hence $P,u,\Theta,C_3,\ldots$ directly.

---

## 7. Consequence for the active mainline

The current 1D mainline is therefore

$$
\boxed{
\text{LJ chain mechanics}
\to F_M(a,c,t)
\to\{P(a,t),u(a,t),\Theta(a,t)\}
\to\{\bar a,\bar U,\text{first passage}\}.
}
$$

The immediate numerical falsification test is

$$
\boxed{
(P,u,\Theta)_{\rm load}(Q^*)
\stackrel{?}{=}
(P,u,\Theta)_{\rm unload}(Q^*)
}
$$

at several identical force levels in the same cycle.

If unequal, the chain-generated distribution state contains load-path history.
The result must still be labeled **dynamic history dependence** until a
mechanically justified irreversible mechanism supplies G3.

A further scale check is mandatory before relating this inertial history to
laboratory fatigue frequencies, because the current atomic normal dynamics and
laboratory cycling have already been shown to have a very large time-scale
separation.
