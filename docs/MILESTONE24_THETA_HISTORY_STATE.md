# Milestone 24 — $P$–$u$–$\Theta$ as a history-bearing reduced descriptor

## Scope

The active model is the conservative 1D normal generalized-LJ chain. No
Boltzmann law, prescribed PDF family, Fokker--Planck/Smoluchowski closure,
damping, or empirical damage variable is introduced here.

The authoritative full derivation is `MASTER_1D_P_U_THETA_FORMULATION.md`.
This milestone records the history-state interpretation and the corrected
second-central-moment equation.

## Exact definitions

For spacing $a_i(t)$ and spacing rate $c_i=\dot a_i$, define

$$
F_M(a,c,t)=\frac1M\sum_i
\delta[a-a_i(t)]\delta[c-c_i(t)],
$$

$$
P_M(a,t)=\int F_M(a,c,t)\,dc.
$$

For a smooth one-point representation,

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

The exact first two moment equations are

$$
\boxed{
\partial_tP+\partial_a(Pu)=0,
}
$$

$$
\boxed{
\partial_t(Pu)+\partial_a[P(u^2+\Theta)]=P\mathcal A,
}
$$

where

$$
\mathcal A(a,t)=\mathbb E[\ddot a_i\mid a_i=a].
$$

Hence

$$
\boxed{
D_tu=\mathcal A-\frac1P\partial_a(P\Theta),
}
$$

and

$$
\boxed{
\Theta\partial_a\ln P
=\mathcal A-D_tu-\partial_a\Theta.
}
$$

## Correct general $\Theta$ balance

Define

$$
C_3(a,t)=\mathbb E[(c-u)^3\mid a]
$$

and

$$
\boxed{
\Psi(a,t)=\operatorname{Cov}(c,\ddot a\mid a)
=\mathbb E[(c-u)\ddot a\mid a].
}
$$

Then the general exact second-central-moment equation is

$$
\boxed{
D_t\Theta
+2\Theta\partial_a u
+\frac1P\partial_a(PC_3)
=2\Psi.
}
$$

The shorter zero-right-hand-side equation is valid only under the additional
condition $\Psi=0$. That condition is not automatic for the spatial LJ chain,
because a spacing acceleration depends on neighbouring spacings.

Thus $(P,u,\Theta)$ is exact as a reduced descriptor but is not an autonomous
three-field closure; its evolution exposes $C_3$, $\Psi$, and neighbour joint
statistics.

## Same-force history dependence

Let $t_L$ and $t_U$ satisfy

$$
Q(t_L)=Q(t_U)=Q^*,
\qquad
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

holds, then the reduced descriptor is not a memoryless function of the
instantaneous force. The current deterministic-chain numerical test confirms
this non-retracing.

This proves **dynamic history dependence** of the reduced descriptor. It does
not prove irreversible fatigue dissipation.

## Mechanical meaning of $\Theta$

Exactly,

$$
\boxed{
\mathbb E[c^2\mid a]=u^2+\Theta.
}
$$

Thus $\Theta$ is conditional spacing-rate dispersion. It is not, by itself,
the exact chain kinetic-energy density. With spacing-rate vector
$\boldsymbol c$,

$$
\boxed{
T=\frac12\boldsymbol c^T\mathbf G_a\boldsymbol c,
\qquad
\mathbf G_a=\mathbf L^T\mathbf L,
}
$$

so total kinetic energy also depends on cross-spacing rate correlations.

## G3 distinction

For the present conservative chain,

$$
\boxed{
\frac{dE_{\rm mech}}{dt}=Q(t)\dot x_{\rm end},
\qquad
D_{\rm irr}=0.
}
$$

Therefore

$$
\boxed{
\text{same-force non-retracing}
\not\Rightarrow
\dot D_{\rm irr}>0.
}
$$

A physical irreversible mechanism is still required for G3 and long-cycle
fatigue accumulation.
