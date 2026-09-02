# Milestone 24 — $P$–$u$–$\Theta$ as a history-bearing reduced state

## 0. Scope

This milestone keeps the active 1D normal generalized-LJ chain and the exact
spacing-space moment structure. It introduces no Boltzmann law, Fokker--Planck
closure, damping, empirical damage variable, or prescribed PDF family.

The purpose is to distinguish three statements that must not be conflated:

1. exact mathematical definition of $P$, $u$, and $\Theta$;
2. load-path dependence / same-force non-retracing of those reduced fields;
3. irreversible dissipation required by G3.

The first two are established in the current reduced theory and numerical
protocol. The third is not supplied by the present conservative chain.

---

## 1. Exact empirical and smooth state definitions

For represented local spacing $a_i(t)$ and spacing rate

$$
c_i(t)=\dot a_i(t),
$$

define the finite empirical phase-space measure

$$
\boxed{
F_M(a,c,t)=\frac1M\sum_{i=1}^{M}
\delta[a-a_i(t)]\delta[c-c_i(t)].
}
$$

Its spacing marginal is

$$
\boxed{
P_M(a,t)=\int F_M(a,c,t)\,dc
=\frac1M\sum_i\delta[a-a_i(t)].
}
$$

For a smooth one-point representation, wherever $P>0$, define

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

Thus $\Theta$ is not a fitted memory coefficient. It is the exact conditional
second central moment of the mechanically generated spacing-rate population.

---

## 2. Exact phase-space and first two moment equations

Let

$$
A(a,c,t)=\mathbb E[\ddot a_i\mid a_i=a,\dot a_i=c]
$$

be the conditional acceleration field. The smooth projected kinetic equation is

$$
\boxed{
\partial_tF+\partial_a(cF)+\partial_c(AF)=0,
}
$$

provided the conditional field is defined from the underlying mechanics. This
is a projected identity, not an autonomous closure: $A$ can retain information
about hidden neighbouring degrees of freedom.

Integrating over $c$ gives

$$
\boxed{
\partial_tP+\partial_a(Pu)=0.
}
$$

Define the one-point conditional acceleration

$$
\boxed{
\mathcal A(a,t)=\mathbb E[\ddot a_i\mid a_i=a].
}
$$

Multiplying the kinetic equation by $c$ and integrating gives

$$
\boxed{
\partial_t(Pu)
+\partial_a\!\left[P(u^2+\Theta)\right]
=P\mathcal A.
}
$$

Therefore

$$
\boxed{
D_tu
=\mathcal A-\frac1P\partial_a(P\Theta),
\qquad
D_tu=\partial_tu+u\partial_a u.
}
$$

Rearranging yields the exact density-shape identity

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

This is an exact shape identity, not an independent predictive closure.

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

If any of

$$
P(a,t_L)\ne P(a,t_U),
$$

$$
u(a,t_L)\ne u(a,t_U),
$$

or

$$
\Theta(a,t_L)\ne\Theta(a,t_U)
$$

holds, then the reduced descriptor is not a single-valued function of the
instantaneous external force alone. Equivalently,

$$
\boxed{
(P,u,\Theta)\ne\mathcal S[Q(t)]
}
$$

for a memoryless scalar map $\mathcal S$.

This is the precise reduced-state meaning of **history dependence** used here.
It does not imply that $(P,u,\Theta)$ is a closed Markov state.

---

## 4. What $\Theta$ measures — and what it does not

At fixed spacing $a$,

$$
\boxed{
\mathbb E[c^2\mid a]=u(a,t)^2+\Theta(a,t).
}
$$

Hence $\Theta$ measures the unresolved spread of local spacing rates at fixed
spacing. It distinguishes, for example, incoming and outgoing populations that
can have similar spacing statistics.

However, $c_i=\dot a_i$ is a bond-spacing rate, not an independent unit-mass
particle velocity. For a chain of unit-mass moving atoms with the left end
fixed, write the spacing vector as $\boldsymbol a$ and the atomic-position
vector as

$$
\boldsymbol x=\mathbf L\boldsymbol a,
$$

where $L_{jk}=1$ for $k\le j$ and $0$ otherwise. Then the exact kinetic energy is

$$
\boxed{
T=\frac12\dot{\boldsymbol a}^{T}
\mathbf G_a
\dot{\boldsymbol a},
\qquad
\mathbf G_a=\mathbf L^{T}\mathbf L.
}
$$

Thus the actual chain kinetic energy contains cross terms between different
spacing rates. One-point $u$ and $\Theta$ alone do not determine total kinetic
energy. The quantity

$$
\frac12[u^2+\Theta]
$$

is therefore a **local spacing-rate quadratic diagnostic**, not by itself the
exact chain kinetic-energy density.

---

## 5. Exact $\Theta$ evolution: acceleration-covariance source

Define the third conditional central spacing-rate moment

$$
\boxed{
C_3(a,t)=\mathbb E[(c-u)^3\mid a]
}
$$

and the conditional spacing-rate/acceleration covariance

$$
\boxed{
\Psi(a,t)
=\operatorname{Cov}(c,\ddot a_i\mid a_i=a)
=\mathbb E[(c-u)\ddot a_i\mid a_i=a].
}
$$

The exact second-central-moment balance is

$$
\boxed{
D_t\Theta
+2\Theta\,\partial_a u
+\frac1P\partial_a(PC_3)
=2\Psi.
}
$$

This is the general equation for the current projected chain.

The commonly shorter form

$$
D_t\Theta
+2\Theta\,\partial_a u
+\frac1P\partial_a(PC_3)=0
$$

is valid only when

$$
\Psi=\operatorname{Cov}(c,\ddot a\mid a)=0,
$$

for example if the acceleration is deterministic at fixed $a,t$. That
condition is **not automatic** for the actual LJ spacing chain because the
central spacing acceleration depends on neighbouring spacings.

Therefore the predictive hierarchy contains both a third velocity moment and
neighbour-conditioned acceleration information. Setting either $C_3=0$ or
$\Psi=0$ without a mechanical derivation would be a closure assumption.

---

## 6. Conservative history dependence is not G3 dissipation

The active normal chain is conservative apart from external work. With
prescribed end force $Q(t)$,

$$
\boxed{
\frac{dE_{\rm mech}}{dt}=Q(t)\dot x_{\rm end}(t).
}
$$

Consequently, a non-retracing force--state loop can arise from inertia, wave
propagation, and phase-space redistribution without irreversible dissipation.
For a cycle,

$$
\boxed{
W_{\rm ext}^{\rm cyc}
=\Delta E_{\rm mech}^{\rm cyc}
+D_{\rm irr}^{\rm cyc}.
}
$$

In the present conservative baseline,

$$
\boxed{D_{\rm irr}^{\rm cyc}=0,}
$$

so a nonzero cycle work during a transient is stored as a change in mechanical
energy. Therefore

$$
\boxed{
\text{same-force non-retracing}
\not\Rightarrow
\dot D_{\rm irr}>0.
}
$$

G3 remains a separate physical requirement.

---

## 7. Current status of the 1D mainline

The exact finite chain is closed at the microscopic state level. Its one-point
projection is

$$
\boxed{
\text{LJ chain mechanics}
\to F_M(a,c,t)
\to\{P(a,t),u(a,t),\Theta(a,t)\}
\to\{\bar a,\bar U,\text{first passage}\}.
}
$$

The tuple $(P,u,\Theta)$ is a mechanically generated, history-bearing reduced
descriptor. It is not yet an autonomous closed state because exact evolution
introduces $C_3$, $\Psi$, and neighbouring joint statistics.

The same-force numerical test establishes dynamic history dependence of this
reduced descriptor. A further physical time-scale and irreversibility mechanism
is still required before that non-retracing can be identified with laboratory
fatigue hysteresis and long-cycle accumulation.
