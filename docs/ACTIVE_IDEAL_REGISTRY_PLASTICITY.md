# Historical multilayer spacing--registry extension

> **STATUS: NON-MAINLINE / HISTORICAL EXTENSION.**
>
> Despite the legacy filename `ACTIVE_IDEAL_REGISTRY_PLASTICITY.md`, this document is
> **not** the active governing theory. The active paper-level baseline is the 1D
> deterministic normal-chain formulation
>
> \[
> \boxed{
> P(\lambda,t),\qquad u(\lambda,t),\qquad \Theta(\lambda,t)
> }
> \]
>
> or equivalently in nondimensional time
>
> \[
> \boxed{
> P(\lambda,\tau),\qquad u(\lambda,\tau),\qquad \Theta(\lambda,\tau),
> \qquad \tau=t/t_0.
> }
> \]
>
> The `U0(a,s)`, Bessel/polylog, unwrapped-registry `z`, and Smoluchowski material
> below are retained as mathematically useful extension/history. They must not be
> cited as the active normal-only evolution law or as a proven low-frequency
> normal-loading slip mechanism without new physical justification.
>
> Current source of truth: `README_EQUATION_INDEX.md`,
> `docs/MASTER_1D_P_U_THETA_FORMULATION.md`, and
> `docs/ASSUMPTIONS.md`.

## Scope of this historical extension

This extension considered repeated **uniaxial tensile loading of a single crystal**
with one normal spacing coordinate `a` and one collective registry coordinate `s`.
A 2D state space `(a,s)` is not a 2D continuum constitutive law. The present active
normal-only paper does not require `s`.

Current checks further constrain the interpretation: for an ideal pure-normal,
perfect-symmetry baseline, registry remains on the symmetric manifold unless a
physical symmetry-breaking mechanism is supplied. Therefore the conservative
baseline reduction is

\[
\boxed{
P(a,s,t)=P(a,t)\,\delta(s-s_0).
}
\]

## Counting convention

For row repeat `b`, row separation `d`, and common registry `s`,

```text
W(d,s) = sum_{p in Z} v_mn(sqrt(d^2+(p b+s)^2))
```

is one row--row kernel. The local extension assumes equally spaced normal layers at
`a,2a,3a,...`, so its intrinsic energy is

```text
U0(a,s) = sum_{k>=1} W(k a,s)
        = sum_{k>=1} sum_{p in Z}
          v_mn(sqrt(k^2 a^2+(p b+s)^2)).
```

There is **no prefactor `k`**. A weighted `sum k W(k a,s)` counts all pairs between
two half-spaces and is a different interface-energy convention. The same
collective/unwrapped `s` appears for every layer; neither `ks` nor `js` is used.

## Exact multilayer sum

With `delta=s/b`, `eta=a/b`, and `q>2`, define

```text
B_q(delta,eta) = sum_p [(p+delta)^2+eta^2]^(-q/2),
H_q(delta,eta) = sum_k B_q(delta,k eta).
```

The retained Mellin--Poisson derivation gives, with `nu=(q-1)/2`,

```text
H_q = sqrt(pi)/Gamma(q/2) [
  Gamma(nu) eta^(1-q) zeta(q-1)
  + 4 sum_{ell>=1} cos(2 pi ell delta) (pi ell/eta)^nu
      Kcal_nu(2 pi ell eta)
],
Kcal_nu(x) = sum_{k>=1} k^(-nu) K_nu(kx).
```

This is an exact convergent representation, not a harmonic approximation. The
absolute double sum requires `q>2`; the single-row `B_q` needs only `q>1`.
Registry differences remove the zero Fourier mode and retain the weaker
exponentially convergent slip-excess structure.

For the well-depth convention

```text
C_mn = m/(m-n) (m/n)^(n/(m-n)),
v_mn = C_mn epsilon_LJ [(sigma_LJ/r)^m-(sigma_LJ/r)^n],
```

the common intrinsic potential is

```text
U0 = C_mn epsilon_LJ [
  (sigma_LJ/b)^m H_m - (sigma_LJ/b)^n H_n
],                     m>n>2.
```

External work is not part of `U0`.

## One-energy normal/slip split

For reference `(a0,s0)`,

```text
Delta U0(a,s) = U0(a,s)-U0(a0,s0)
              = Delta U_n(a)+V_slip(a,s),
Delta U_n     = U0(a,s0)-U0(a0,s0),
V_slip        = U0(a,s)-U0(a,s0).
```

This identity prevents double counting. The old collinear `U_infinity` and
single-row `W` are useful reduced derivations but are not summed as the same total
energy.

The slip excess uses

```text
Delta H_q = H_q(delta,eta)-H_q(delta0,eta)
```

and contains only cosine differences. Therefore `V_slip(a,s0)=0` and
`V_slip(a,s+b)=V_slip(a,s)` exactly. For `s0=0`, each difference is
`cos(2 pi ell delta)-1=-2 sin^2(pi ell delta)`.

## Independently derived 12--6 polylog closure

The half-integer identities are

```text
K_5/2(x)  = sqrt(pi/(2x)) exp(-x) (1+3/x+3/x^2),
K_11/2(x) = sqrt(pi/(2x)) exp(-x)
             (1+15/x+105/x^2+420/x^3+945/x^4+945/x^5).
```

Multiplication by `k^(-nu)` shows directly that

```text
Kcal_5/2(x) = sqrt(pi/(2x)) [
  Li_3(e^-x)+3/x Li_4(e^-x)+3/x^2 Li_5(e^-x)],

Kcal_11/2(x) = sqrt(pi/(2x)) [
  Li_6(e^-x)+15/x Li_7(e^-x)+105/x^2 Li_8(e^-x)
  +420/x^3 Li_9(e^-x)+945/x^4 Li_10(e^-x)
  +945/x^5 Li_11(e^-x)].
```

The orders start at `q/2` because the square-root prefactor contributes one
additional `k^(-1/2)`. Tests compare these expressions with both the direct double
sum and the unsimplified Bessel--Lambert series.

## Historical conditional Smoluchowski closure — not active

The earlier extension considered the single applied stress

```text
sigma(t)=sigma_m+sigma_a sin(omega t),
Q_a=A0 sigma(t),       Q_s=A0 M sigma(t),
```

and then adopted the conditional reduced law

```text
partial_t P = -partial_a J_a-partial_s J_s,
J_a = -M_a [P(partial_a U0-Q_a)+kBT partial_a P],
J_s = -M_s [P(partial_s U0-Q_s)+kBT partial_s P].
```

These equations require an eliminated isothermal bath, overdamped/fast velocity
relaxation, Markov reduction, mobilities, and an Einstein fluctuation--dissipation
relation. Those assumptions have not been established for the active deterministic
normal-chain theory. Consequently these equations are **historical conditional
closure material**, not the current governing law.

Likewise the historical dissipation expression

```text
dot(D)_irr = integral integral [J_a^2/(M_a P)+J_s^2/(M_s P)] da ds >=0
```

is not the active G3 unless the bath reduction and fluctuation--dissipation structure
are physically re-derived.

## Historical observables and plasticity extension

The extension used

```text
G1  bar(a) = integral integral a P da ds.
G2  bar(U) = integral integral [U0(a,s)-U0(a0,s0)] P da ds.
G3  E_hyst(t) = integral_0^t dot(D)_irr dt.
G4  integral integral P da ds = 1
    (or S(t)<=1 with an absorbing fracture boundary).
```

The labels G1--G4 are retained project-wide, but the active 1D paper evaluates G1,
G2 and first passage from the deterministic normal-chain `P(lambda,t)` state, while
nonzero G3 remains physically OPEN.

Registry was represented by

```text
s=s0+z b+tilde(s),       z in Z,
p_z=integral_{W_z} P da ds,
gamma_p=(b/h_slip)<z>,   epsilon_p=M gamma_p.
```

Finite-rate intrawell lag is recoverable. A residual `Delta<z> != 0` after unloading
and relaxation was the intended plasticity criterion. Current normal-only checks do
not establish the physical symmetry-breaking mechanism required to generate such a
transition.

## Historical first-passage extension

At a given tensile drive, the normal escape boundary was written as

```text
partial_a U0(a^dagger,s)=Q_a(t),
partial_a^2 U0(a^dagger,s)<0.
```

For a moving graph `a=a^dagger(s,t)`, Reynolds transport gives the relative outflux
(periodic/no-flux registry edges)

```text
-dS/dt = integral [J_a-J_s partial_s a^dagger
                   -P partial_t a^dagger] ds.
```

This remains useful extension mathematics, but the active 1D baseline uses the
normal-chain kinetic first-passage construction and its own mechanically consistent
threshold.

## Present interpretation

This document now records a **mathematically defined but non-mainline extension**.
It does not determine the active paper's probability evolution, irreversible G3,
or a proven low-frequency slip mechanism. Any future return to `(a,s)` must state
and validate the added symmetry-breaking/irreversible physics explicitly.
