# Active multilayer spacing--registry fatigue theory

> **Geometry correction (Milestone 20).** The historical reduced driving pair
> `Q_a=A0 sigma`, `Q_s=A0 M sigma` must not be treated as a single literal
> FCC plane-opening/plane-slip virtual-work identity when `M != 0`.  If `a` is
> the normal opening of the same slip plane that carries `s`, exact virtual
> work gives `Q_a=A0 sigma (l.n)^2` and `Q_s=A0 M sigma`.  If `a` remains the
> loading-axis spacing, the coupled FCC energy must instead be re-derived for
> loading-axis deformation plus internal slip.  See
> `MILESTONE20_FCC_GEOMETRY_CONSISTENCY.md`.  The probability observables G1--G4
> are unchanged by this correction.

## Scope

The active fundamental model describes repeated **uniaxial tensile loading of
a single crystal** with one normal spacing coordinate `a` and one declared
crystallographic slip coordinate `s`.  A 2D state space `(a,s)` is not a 2D
continuum constitutive law.  No independent shear-fatigue input, multiaxial
criterion, fitted Peierls sinusoid, or EAM/DFT lookup surface is active.

## Counting convention

For row repeat `b`, row separation `d`, and common registry `s`,

```text
W(d,s) = sum_{p in Z} v_mn(sqrt(d^2+(p b+s)^2))
```

is one row--row kernel.  The local fatigue state assumes equally spaced normal
layers at `a,2a,3a,...`, so its intrinsic energy is

```text
U0(a,s) = sum_{k>=1} W(k a,s)
        = sum_{k>=1} sum_{p in Z}
          v_mn(sqrt(k^2 a^2+(p b+s)^2)).
```

There is **no prefactor `k`**.  A weighted `sum k W(k a,s)` counts all pairs
between two half-spaces and is a different interface-energy convention.  The
same collective/unwrapped `s` appears for every layer; neither `ks` nor `js`
is used.

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

This is an exact convergent representation, not a harmonic approximation.
The absolute double sum requires `q>2`; the single-row `B_q` needs only
`q>1`.  Registry differences remove the zero Fourier mode and retain the
weaker exponentially convergent slip-excess structure.

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

This identity prevents double counting.  The old collinear `U_infinity` and
single-row `W` are useful reduced derivations but are not summed as the active
fundamental energy.

The slip excess uses

```text
Delta H_q = H_q(delta,eta)-H_q(delta0,eta)
```

and contains only cosine differences.  Therefore `V_slip(a,s0)=0` and
`V_slip(a,s+b)=V_slip(a,s)` exactly.  For `s0=0`, each difference is
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
additional `k^(-1/2)`.  Tests compare these expressions with both the direct
double sum and the unsimplified Bessel--Lambert series.

## Uniaxial driving: geometry status

The historical reduced model wrote

```text
sigma(t)=sigma_m+sigma_a sin(omega t),
Q_a=A0 sigma(t),
Q_s=A0 M sigma(t).
```

The `Q_s` expression is the signed Schmid projection for a declared tensile
axis, slip-plane normal, and in-plane slip direction.  Milestone 20 shows that
if `a` and `s` are interpreted literally as normal and tangential relative
coordinates of that same plane patch, virtual work instead gives

```text
Q_a=A0 sigma(t) (l.n)^2,
Q_s=A0 sigma(t) (l.n)(l.d)=A0 M sigma(t).
```

Hence `Q_a=A0 sigma` and nonzero `Q_s=A0 M sigma` are not simultaneously exact
for one literal plane-relative embedding.  The project must choose between:

1. **slip-plane coordinates:** use the projected `Q_a` above and validate/extend
   `U0(a,s)` against FCC stacking; or
2. **loading-axis spacing plus internal slip:** retain the axial meaning of `a`
   but re-derive the coupled FCC energy rather than reading the current row
   geometry literally.

Until that choice is completed, the old pair is retained only as a historical
reduced closure and is not a foundational identity.

## Historical reduced Smoluchowski closure

The earlier implementation used

```text
partial_t P = -partial_a J_a-partial_s J_s,
J_a = -M_a [P(partial_a U0-Q_a)+kBT partial_a P],
J_s = -M_s [P(partial_s U0-Q_s)+kBT partial_s P].
```

This remains a **specific overdamped/Markov/isothermal closure**, not the
fundamental definition of `P`.  The newer mainline probability foundation is
the exact empirical transport/moment hierarchy and the Theta-based shape
identity.  The Smoluchowski form must only be used when its assumptions are
explicitly justified.

## The four governing equations

The official observables remain

```text
G1  bar(a) = integral integral a P da ds.

G2  bar(U) = integral integral [U0(a,s)-U0(a0,s0)] P da ds.

G3  E_hyst(t) = integral_0^t dot(D)_irr dt.

G4  integral integral P da ds = 1
    (or S(t)<=1 with an absorbing fracture boundary).
```

For the historical Smoluchowski closure only,

```text
dot(D)_irr = integral integral [J_a^2/(M_a P)+J_s^2/(M_s P)] da ds >=0.
```

That quadratic expression is not promoted as a closure-independent universal
law.  `bar(U)` can decrease after a jump into an equivalent registry well
because `U0(a,s+b)=U0(a,s)`.  `E_hyst` is cumulative irreversible/hysteretic
dissipation and must be tied to a physically justified irreversible mechanism.

## Plasticity and crack initiation

Registry is never folded before its well population is measured:

```text
s=s0+z b+tilde(s),       z in Z,
p_z=integral_{W_z} P da ds,
gamma_p=(b/h_slip)<z>,   epsilon_p=M gamma_p.
```

Finite-rate intrawell lag is recoverable.  A well crossing is a slip event; the
strong plasticity criterion additionally requires residual well-index change
after unloading and an explicitly justified relaxation/irreversible mechanism.

At a given tensile drive, the normal escape boundary is the outer root

```text
partial_a U0(a^dagger,s)=Q_a(t),
partial_a^2 U0(a^dagger,s)<0.
```

The correct `Q_a` depends on the coordinate interpretation described above.
For a moving graph `a=a^dagger(s,t)`, Reynolds transport gives the relative
outflux (periodic/no-flux registry edges)

```text
-dS/dt = integral [J_a-J_s partial_s a^dagger
                   -P partial_t a^dagger] ds.
```

At an ideal absorbing boundary `P=0`, the boundary-motion term vanishes.
Crack initiation is this first-passage probability `1-S`, not an arbitrary
spacing or accumulated-energy threshold.

## Current numerical status and limitations

The direct `(k,p)` sum, the Fourier--Bessel representation, the
Bessel--Lambert form, and the 12--6 polylog closure are tested at several
registries and normal separations.  Numerical truncations must still be
refined for each new parameter regime.

The theory is a mathematically derived reduced single-slip mechanism.  It
does not yet determine `A0`, the finite energy--mass patch mapping, the physical
origin/evolution of the full conditional covariance `Theta`, the irreversible
mechanism for G3, `h_slip`, dislocation storage/hardening, or a quantitative
active-slip selection rule.  EAM/DFT remains only a future quantitative-
aluminum validation/extension; it does not replace the current generalized-LJ
governing potential without an explicit model change.
