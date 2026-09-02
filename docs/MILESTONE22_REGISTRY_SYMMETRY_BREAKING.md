# Milestone 22 — Registry symmetry breaking without an imposed shear drive

## 0. Scope

This milestone keeps the active reduced model only:

$$
U_0(a,s)=\sum_{k\ge1}\sum_{p\in\mathbb Z}
v_{m,n}\!\left(\sqrt{k^2a^2+(pb+s)^2}\right).
$$

No FCC geometry, Boltzmann distribution, Fokker--Planck closure, damping, random
noise, or empirical damage law is introduced.

The question is narrower:

> Can cyclic **normal** loading move the registry coordinate away from its
> symmetric minimum without inserting an arbitrary external $q_s(t)$?

The answer has three levels that must be kept distinct:

1. exact symmetry invariance;
2. static loss of registry stability;
3. dynamic/parametric amplification of a nonzero physical seed.

---

## 1. Exact invariant registry manifold

For the present periodic row energy, the registry dependence is a cosine
series.  At a symmetry point $s=s_0$ (for the current phase convention the
stable well used numerically is $s_0/b=1/2$),

$$
\boxed{U_s(a,s_0)=0\qquad\text{for every }a.}
$$

With no generalized registry force,

$$
\mu_s\ddot s=-U_s(a,s),
$$

so the initial condition

$$
s(0)=s_0,\qquad \dot s(0)=0
$$

has the exact solution

$$
\boxed{s(t)=s_0.}
$$

Therefore normal loading cannot *deterministically create* a finite registry
perturbation from mathematically exact zero in the perfectly symmetric reduced
model.

This explains the normal-only spatial-chain result in which
$\operatorname{Var}(s)$ stayed at numerical zero.

---

## 2. Linear stability of a small registry perturbation

Let

$$
s=s_0+\xi,\qquad |\xi|\ll b.
$$

Because $U_s(a,s_0)=0$,

$$
U_s(a,s_0+\xi)
=U_{ss}(a,s_0)\xi+O(\xi^2).
$$

Hence

$$
\boxed{
\mu_s\ddot\xi+K_s(a(t))\xi=0,
\qquad K_s(a)=U_{ss}(a,s_0).
}
$$

Thus the normal coordinate changes the **registry stiffness** even though it
does not exert a direct shear force at the symmetry point.

This equation is exact to first order in the registry perturbation.

---

## 3. Static registry instability

A static registry minimum is locally stable when

$$
K_s(a)=U_{ss}(a,s_0)>0.
$$

It becomes marginal at

$$
\boxed{U_{ss}(a_s^*,s_0)=0}
$$

and locally unstable for $U_{ss}<0$.

The normal force-controlled branch has its own local stability condition

$$
U_{aa}(a,s_0)>0,
$$

with marginal normal opening at

$$
\boxed{U_{aa}(a_a^*,s_0)=0.}
$$

For the normalized diagnostic parameters used by the spatial-chain run

$$
m=12,\quad n=6,\quad b=\sigma_{LJ}=\epsilon_{LJ}=1,
$$

and $s_0/b=1/2$, direct $(k,p)$ differentiation gives the converged roots

$$
\boxed{a_a^*=1.130690887}
$$

and

$$
\boxed{a_s^*=1.264187982.}
$$

The registry-curvature zero therefore occurs about $11.81\%$ beyond the normal
curvature zero:

$$
\frac{a_s^*-a_a^*}{a_a^*}\approx0.1181.
$$

At the normal marginal point,

$$
\boxed{U_{ss}(a_a^*,s_0)\approx2.38763>0.}
$$

so the registry direction is still locally stable when the normal direction
already loses local stability.

### Consequence

For this normalized active energy, **static normal opening does not produce a
registry instability before the normal cohesive/opening instability**.
Therefore a proposed mechanism

$$
\text{normal opening}\to U_{ss}<0\to\text{slip before crack}
$$

is rejected for this parameter set.

This is a model result, not a universal materials statement.

---

## 4. Numerical convergence of the two curvature roots

Direct-sum truncation gives

| $k_{\max}$ | $p_{\max}$ | $a_a^*$ | $a_s^*$ |
|---:|---:|---:|---:|
| 20 | 50 | 1.130691244 | 1.264187982 |
| 40 | 100 | 1.130690910 | 1.264187982 |
| 80 | 200 | 1.130690888 | 1.264187982 |
| 120 | 300 | 1.130690887 | 1.264187982 |
| 200 | 500 | 1.130690887 | 1.264187982 |

The ordering $a_a^*<a_s^*$ is therefore not a truncation artifact at these
resolutions.

---

## 5. Cyclic normal loading still modulates registry stiffness

Even while $U_{ss}>0$, cyclic normal motion gives

$$
K_s(t)=U_{ss}(a(t),s_0),
$$

so a small physical registry perturbation obeys

$$
\boxed{
\mu_s\ddot\xi+U_{ss}(a(t),s_0)\xi=0.
}
$$

This is a linear equation with time-dependent stiffness.  If $a(t)$ is
periodic, it is a Hill equation.

This statement does **not** assume a sinusoidal response or a Mathieu model.
The exact periodic $a(t)$ from the spatial chain can be inserted directly and
the monodromy/Floquet multiplier can be computed.

---

## 6. Mathieu form is only a controlled local approximation

If, additionally,

1. the normal oscillation is small around $\bar a$;
2. its first harmonic dominates;
3. higher-order dependence of $U_{ss}$ on $a$ is negligible over the excursion,

then write

$$
a(t)\approx\bar a+A_a\cos\Omega t.
$$

Expanding the registry stiffness,

$$
U_{ss}(a(t),s_0)
\approx K_0+K_1\cos\Omega t,
$$

where

$$
K_0=U_{ss}(\bar a,s_0),
$$

$$
K_1=U_{ass}(\bar a,s_0)A_a.
$$

Then

$$
\boxed{
\ddot\xi+
\left[
\omega_s^2+h\cos\Omega t
\right]\xi=0,
}
$$

with

$$
\omega_s^2=\frac{K_0}{\mu_s},
\qquad
h=\frac{K_1}{\mu_s}.
$$

This is Mathieu-type parametric excitation.  The common principal resonance
near

$$
\Omega\approx2\omega_s
$$

is **not** adopted as a governing law; it is only a small-modulation diagnostic.
For the actual reduced model the preferred test is the full Hill/Floquet
problem using computed $a_i(t)$.

For the normalized reference state $a_0=0.9919601754$,

$$
U_{ss}(a_0,s_0)\approx25.7179,
$$

and

$$
U_{ass}(a_0,s_0)\approx-413.789.
$$

Thus registry stiffness is highly sensitive to normal opening in the current
reduced potential even though the symmetry point remains stationary.

For the previous diagnostic loading frequency $\Omega=0.35$, the simple
principal-resonance estimate would require approximately

$$
\mu_s\approx
\frac{U_{ss}(a_0,s_0)}{(\Omega/2)^2}
\approx8.40\times10^2
$$

in the same normalized units.  Because $\mu_s$ has not yet been physically
fixed, this number is only a scaling diagnostic, not a material prediction.

---

## 7. Seed versus instability must not be conflated

A Hill/Floquet instability can amplify a nonzero perturbation, but the exact
symmetric solution remains an exact solution:

$$
\xi(0)=\dot\xi(0)=0
\quad\Longrightarrow\quad
\xi(t)=0.
$$

Therefore the current ideal deterministic reduced model does not itself define
the initial seed.

Possible *physical* seed categories include:

- a specified finite-temperature microscopic initial state;
- a pre-existing registry imperfection/defect;
- a boundary or geometric asymmetry;
- residual registry motion from a previous loading history.

These are categories, not assumptions adopted here.  In particular, no
Boltzmann/Gaussian seed distribution is introduced by this milestone.

---

## 8. Implication for the probability model

The normal coordinate already has a deterministic spatial source of spread:

$$
Q_a(t)\to\{a_i(t)\}\to P_M(a,t).
$$

For registry, the ideal symmetric $T=0$ reduced model gives

$$
\boxed{P_M(s,t)=\delta(s-s_0)}
$$

unless a nonzero physical seed is present.

Once such a seed exists, the deterministic equation

$$
\mu_s\ddot\xi+U_{ss}(a_i(t),s_0)\xi=0
$$

can amplify or suppress it cell by cell.  If amplification differs spatially,
a nontrivial joint $P_M(a,s,t)$ and nonzero $\Theta_{as}$ can emerge without an
imposed $q_s(t)$.

Therefore the next quantitative task is **not** to invent a shear drive.  It is
to compute the Floquet stability of the registry perturbation using the actual
normal-only spatial-chain trajectories $a_i(t)$ as coefficients, while keeping
$\mu_s$ explicit until its physical mapping is justified.

---

## 9. Current conclusion

For the present reduced energy:

$$
\boxed{
\text{exact symmetry}\Rightarrow s=s_0\text{ remains invariant},
}
$$

$$
\boxed{
\text{static }s\text{-instability occurs after normal instability},
}
$$

but

$$
\boxed{
\text{cyclic normal motion can parametrically modulate the stability of any
physical registry seed through }U_{ss}(a(t),s_0).
}
$$

This is the only internal no-$q_s$ symmetry-breaking route identified so far
that is compatible with the active reduced model.  Whether it is quantitatively
relevant depends on the still-unresolved registry inertia and on the actual
normal-cycle trajectory.
