# Exact Derivations

This file records only identities or equations that are exact under an explicitly stated microscopic model.

## Thermodynamic-limit spacing density

For local spacing variables $a_i(t)$, define

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta(a-a_i(t)),
\qquad
P(a,t)=\lim_{N\to\infty}P_N(a,t).
$$

## Exact continuity equation in spacing space

Differentiating the empirical density under smooth deterministic dynamics gives

$$
\partial_tP+\partial_a(Pv)=0,
$$

where

$$
v(a,t)=\langle \dot a_i\mid a_i=a\rangle.
$$

This equation is exact as a kinematic identity. It is not, by itself, a closed evolution law because $v(a,t)$ depends on unresolved microscopic information.

## Moment identities

For any sufficiently regular moment,

$$
\frac{d}{dt}\langle a^n\rangle
=n\langle a^{n-1}v\rangle.
$$

In particular,

$$
\dot{\bar a}=\langle v\rangle,
$$

and

$$
\frac{d}{dt}\operatorname{Var}(a)=2\operatorname{Cov}(a,v).
$$

## Phase-space lift

Because atomistic mechanics is second order, introduce a joint spacing-velocity density

$$
F(a,c,t),\qquad c=\dot a,
$$

with

$$
P(a,t)=\int F(a,c,t)\,dc.
$$

The projected Liouville form is

$$
\partial_tF+\partial_a(cF)+\partial_c(AF)=0,
$$

where

$$
A(a,c,t)=\langle \ddot a_i\mid a_i=a,\dot a_i=c\rangle.
$$

## Nearest-neighbor chain relation

For $V=\sum_i v(a_i)$ and $a_i=x_{i+1}-x_i$, Newton's equations imply

$$
m\ddot a_i=v'(a_{i+1})-2v'(a_i)+v'(a_{i-1}).
$$

Therefore a one-point density $P(a,t)$ is generally not exactly closed; neighboring-spacing correlations enter the acceleration.

## Exact pair-distance energy hierarchy

The distance to the $k$-th neighbor is

$$
R_i^{(k)}=a_i+a_{i+1}+\cdots+a_{i+k-1}.
$$

If $P_k(r,t)$ is the exact density of $R_i^{(k)}$, then the configurational pair-potential energy per atom is

$$
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
$$

For a uniform lattice $P_k(r)=\delta(r-ka)$, this reduces to

$$
U(a)=\sum_{k=1}^{\infty}v(ka).
$$
