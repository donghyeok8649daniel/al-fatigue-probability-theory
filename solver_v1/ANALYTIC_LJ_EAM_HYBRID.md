# Analytic LJ--EAM hybrid

## 1. Purpose and status

This model is a minimal many-body extension of the existing analytic two-row
Lennard-Jones model. The LJ pair interaction is deliberately retained. An
EAM-style environment density and embedding energy are added:

$$
E=E_{\mathrm{LJ,pair}}+E_{\mathrm{embedding}},
$$

$$
E_{\mathrm{LJ,pair}}
=\frac12\sum_i\sum_{j\ne i}\phi_{\mathrm{LJ}}(r_{ij}),
\qquad
E_{\mathrm{embedding}}=\sum_iF(\varrho_i).
$$

This is an **analytic LJ--EAM hybrid**, not a standard, exact, or calibrated
aluminum EAM potential. EAM pair/embedding decompositions have gauge freedom,
so adding an embedding term to LJ defines a new hybrid model even when the
individual expressions resemble an EAM parameterization.

The probability dynamics are unchanged. The full two-coordinate
Smoluchowski equation, opening first passage, interwell flux, strain
decomposition, and absorbed-mass bookkeeping use the selected energy surface.

## 2. LJ pair base

The pair potential remains

$$
\phi_{\mathrm{LJ}}(r)
=4\epsilon_{\mathrm{LJ}}
\left[
\left(\frac{\sigma_{\mathrm{LJ}}}{r}\right)^{12}
-\left(\frac{\sigma_{\mathrm{LJ}}}{r}\right)^6
\right].
$$

For the uniform infinite two-row geometry, the registry-dependent cross-row
pair energy is exactly the existing Poisson/Bessel surface
$W_{\mathrm{LJ}}(a,s)$. Its analytic first and second derivatives are not
changed or refitted. Same-row LJ pair energies are independent of $a$ and $s$
for a fixed row spacing $b$ and are therefore an additive cell constant in the
present reduced energy.

## 3. Exponential environment-density kernel

Use a distinct environment-density symbol $\varrho$, not the probability
density $P(a,s,t)$. The site contribution is

$$
f_\varrho(r)
=\varrho_0\exp\left[-\beta_\varrho
\left(\frac{r}{r_e}-1\right)\right]
=C_\varrho e^{-\kappa r},
$$

where

$$
C_\varrho=\varrho_0e^{\beta_\varrho},
\qquad
\kappa=\frac{\beta_\varrho}{r_e}.
$$

$\varrho_0$ carries the selected electron/environment-density unit,
$\beta_\varrho$ is dimensionless, $r_e$ is a length, and $\kappa$ has inverse
length units. In the repository's dimensionless implementation, lengths are
expressed in the same unit as $b$, and energies in the LJ energy unit.

## 4. Two-row geometry and same-row density

For an upper reference atom and lower-row site $n$,

$$
x_n(s)=\left(n+\frac12\right)b-s,
\qquad
r_n(a,s)=\sqrt{a^2+x_n(s)^2}.
$$

The cross-row density is the untruncated infinite sum

$$
\varrho_\perp(a,s)
=C_\varrho\sum_{n=-\infty}^{\infty}
e^{-\kappa\sqrt{a^2+[(n+1/2)b-s]^2}}.
$$

The same-row reference atom is excluded. The remaining sites give the exact
geometric series

$$
\varrho_\parallel
=2C_\varrho\sum_{n=1}^{\infty}e^{-\kappa nb}
=\boxed{\frac{2C_\varrho}{e^{\kappa b}-1}}.
$$

## 5. Fourier transform of the exponential radial kernel

Let

$$
g(x)=e^{-\kappa\sqrt{a^2+x^2}},
\qquad
\widehat g(k)=\int_{-\infty}^{\infty}g(x)e^{-ikx}\,dx.
$$

First use the standard transform

$$
\int_{-\infty}^{\infty}
\frac{e^{-\kappa\sqrt{a^2+x^2}}}{\sqrt{a^2+x^2}}
e^{-ikx}\,dx
=2K_0(aq),
\qquad q=\sqrt{\kappa^2+k^2}.
$$

Differentiation with respect to $\kappa$ is permitted by exponential
domination. Since the derivative of the integrand is $-g(x)e^{-ikx}$ and
$K_0'(z)=-K_1(z)$,

$$
-\widehat g(k)
=\frac{\partial}{\partial\kappa}\left[2K_0(aq)\right]
=-\frac{2a\kappa}{q}K_1(aq).
$$

Therefore

$$
\boxed{
\widehat g(k)=\frac{2a\kappa}{\sqrt{\kappa^2+k^2}}
K_1\left(a\sqrt{\kappa^2+k^2}\right)
}.
$$

Independent adaptive quadrature of the original Fourier integral verifies
this identity for zero and nonzero wave number.

## 6. Poisson-summed cross-row density

For the shift $d=b/2-s$, Poisson summation gives

$$
\sum_{n\in\mathbb Z}g(nb+d)
=\frac1b\sum_{m\in\mathbb Z}
\widehat g(k_m)e^{ik_md},
\qquad
k_m=\frac{2\pi m}{b}.
$$

Because $\widehat g$ is real and even, define

$$
q_m=\sqrt{\kappa^2+k_m^2},
\qquad
H_m(a)=\frac{2a\kappa}{q_m}K_1(aq_m),
\qquad
\theta_m=k_m\left(\frac b2-s\right).
$$

The canonical reciprocal representation is

$$
\boxed{
\varrho_\perp(a,s)
=\frac{C_\varrho}{b}H_0(a)
+\frac{2C_\varrho}{b}
\sum_{m=1}^{\infty}H_m(a)\cos\theta_m
}.
$$

Equivalently, the requested coefficients are

$$
C_0(a)=\frac{2C_\varrho a}{b}K_1(a\kappa),
$$

$$
C_m(a)=\frac{4C_\varrho a\kappa}{bq_m}K_1(aq_m),
\quad m\ge1,
$$

so that

$$
\varrho_\perp=C_0+\sum_{m\ge1}C_m
\cos\left[2\pi m\left(\frac12-\frac{s}{b}\right)\right].
$$

No finite real-space neighbor cutoff is used in canonical evaluation.

## 7. Analytic density derivatives

The stable reusable amplitude identities are

$$
H_m'(a)=-2\kappa aK_0(aq_m),
$$

$$
H_m''(a)=2\kappa
\left[aq_mK_1(aq_m)-K_0(aq_m)\right].
$$

They follow from
$d[zK_1(z)]/dz=-zK_0(z)$ and avoid expanding fragile combinations of many
Bessel derivatives.

The same-row term is constant in $a,s$. Hence

$$
\varrho_a
=\frac{C_\varrho}{b}H_0'
+\frac{2C_\varrho}{b}\sum_{m\ge1}H_m'\cos\theta_m,
$$

$$
\varrho_s
=\frac{2C_\varrho}{b}\sum_{m\ge1}k_mH_m\sin\theta_m,
$$

$$
\varrho_{aa}
=\frac{C_\varrho}{b}H_0''
+\frac{2C_\varrho}{b}\sum_{m\ge1}H_m''\cos\theta_m,
$$

$$
\varrho_{as}
=\frac{2C_\varrho}{b}\sum_{m\ge1}k_mH_m'\sin\theta_m,
$$

$$
\varrho_{ss}
=-\frac{2C_\varrho}{b}\sum_{m\ge1}k_m^2H_m\cos\theta_m.
$$

These formulas establish both periodicity in $s$ and the expected symmetry
$\varrho_s(a,0)=0$.

## 8. Equality of upper and lower site densities

The upper reference atom sees horizontal offsets
$(n+1/2)b-s$. A lower reference atom sees the translated upper row with the
opposite registry sign. Because the radial kernel is even, the substitution
$n\mapsto-n-1$ maps one infinite sum exactly to the other. Both rows also have
the same spacing $b$ and therefore the same $\varrho_\parallel$. Thus

$$
\varrho_{\mathrm{upper}}(a,s)
=\varrho_{\mathrm{lower}}(a,s)
=\varrho(a,s)
=\varrho_\parallel+\varrho_\perp(a,s).
$$

This equality is specific to the present uniform infinite two-row geometry;
it must not be assumed for a surface, defect, finite row, or inequivalent
multi-site environment.

## 9. Unit-cell counting and hybrid energy

Choose one upper and one lower atom per translational cell. The cross-row LJ
sum counts each upper--lower pair exactly once per cell and equals the existing
$W_{\mathrm{LJ}}(a,s)$. Pair double counting has already been removed by the
standard $1/2$ in the full pair energy. There are two embedding sites with the
same density, so the variable cell energy is

$$
\boxed{
W_{\mathrm{hybrid}}(a,s)
=W_{\mathrm{LJ}}(a,s)+2F(\varrho(a,s))
}.
$$

The factor 2 is therefore a derived site multiplicity, not an adjustable
coefficient. Fixed same-row pair energy may be included as an additive
constant but has no effect on forces, Hessians, or barriers.

An optional referenced energy is

$$
W_{\mathrm{hybrid,ref}}(a,s)
=W_{\mathrm{hybrid}}(a,s)-W_{\mathrm{hybrid}}(a_{0,h},0).
$$

This changes no derivative or barrier.

## 10. General embedding derivatives

For any twice differentiable scalar $F(\varrho)$,

$$
W_a=W_{\mathrm{LJ},a}+2F'(\varrho)\varrho_a,
$$

$$
W_s=W_{\mathrm{LJ},s}+2F'(\varrho)\varrho_s,
$$

$$
W_{aa}=W_{\mathrm{LJ},aa}
+2\left[F''(\varrho)\varrho_a^2+F'(\varrho)\varrho_{aa}\right],
$$

$$
W_{as}=W_{\mathrm{LJ},as}
+2\left[F''(\varrho)\varrho_a\varrho_s
+F'(\varrho)\varrho_{as}\right],
$$

$$
W_{ss}=W_{\mathrm{LJ},ss}
+2\left[F''(\varrho)\varrho_s^2+F'(\varrho)\varrho_{ss}\right].
$$

If future geometry produces inequivalent site densities, these expressions
generalize to a sum over sites; replacing them by one averaged density would
not be exact.

## 11. Minimal analytic embedding family

The first differentiable candidate is

$$
F(\varrho)=-A\sqrt{\frac{\varrho}{\varrho_{\mathrm{ref}}}},
$$

$$
F'(\varrho)
=-\frac{A}{2\sqrt{\varrho_{\mathrm{ref}}}\sqrt{\varrho}},
$$

$$
F''(\varrho)
=\frac{A}{4\sqrt{\varrho_{\mathrm{ref}}}\varrho^{3/2}}.
$$

$A$ has energy units and $\varrho_{\mathrm{ref}}$ has the same density unit
as $\varrho_0$. `AnalyticEmbedding` permits replacement by another analytic
family without changing the lattice-density implementation.

## 12. Equilibrium and relaxed axial calibration

Embedding generally shifts the reference spacing. It is recomputed from

$$
W_{\mathrm{hybrid},a}(a_{0,h},0)=0
$$

using a bracketed high-precision root. The original LJ reference remains

$$
a_{0,\mathrm{LJ}}=0.7713438268704838.
$$

For the clearly hypothetical sensitivity set

$$
\varrho_0=1,\quad\beta_\varrho=3,\quad r_e=1,
\quad A=0.08,\quad\varrho_{\mathrm{ref}}=1,
$$

the hybrid values are

$$
a_{0,h}=0.7693829741307685,
\qquad
|W_a(a_{0,h},0)|=3.07\times10^{-14},
$$

$$
H_{0,h}\approx
\begin{bmatrix}
126.840164 & -1.94\times10^{-14}\\
-1.94\times10^{-14} & 52.169960
\end{bmatrix}.
$$

Both eigenvalues are positive. With
$\mathbf c=(1,\chi)^T$, each energy surface has its own calibration

$$
\kappa_{\mathrm{axial}}
=\frac{a_0}{\mathbf c^TH_0^{-1}\mathbf c}.
$$

The unchanged LJ value is 86.2929648874, while this hypothetical hybrid gives
88.9391933042. Stationary solutions at $\sigma/E=\mp10^{-5}$ give tangent
ratios 0.99992506 and 1.00007494. No area factor enters either mapping.

## 13. Reciprocal truncation and numerical verification

The reciprocal series is exponentially convergent for $a>0$. Production
evaluation continues until consecutive value/gradient/Hessian mode
contributions and a geometric tail estimate lie below a declared tolerance,
subject only to a safety ceiling. It returns both `modes_used` and
`estimated_tail_absolute`.

Across 24 deterministic validation states spanning $a$, $s$, $b$,
$\beta_\varrho$, $r_e$, and $\varrho_0$:

- maximum reciprocal-versus-900-image real-space absolute error:
  $2.67\times10^{-15}$;
- maximum relative error: $7.53\times10^{-16}$;
- maximum reciprocal modes used at tolerance $10^{-14}$: 14;
- maximum first-density-derivative absolute finite-difference error:
  $6.45\times10^{-10}$;
- maximum second-density-derivative absolute error:
  $9.59\times10^{-7}$;
- maximum hybrid-energy gradient absolute error: $1.05\times10^{-8}$;
- maximum hybrid Hessian relative error: $6.29\times10^{-7}$.

Relative derivative errors close to symmetry zeros are ill-conditioned, so
absolute and scaled-relative errors are both tested.

## 14. Landscape sensitivity, not calibration

For the same hypothetical set, comparison at identical dimensionless force is:

| $f^*$ | LJ configuration barrier | hybrid configuration barrier | LJ opening barrier at config saddle | hybrid opening barrier at config saddle |
|---:|---:|---:|---:|---:|
| 0.5 | 0.59494 | 0.62476 | 0.92170 | 0.99671 |
| 2.0 | 0.26614 | 0.29241 | 0.32320 | 0.36537 |
| 3.5 | 0.03474 | 0.05045 | 0.10720 | 0.12432 |
| 4.0 | 0.000363 | 0.006890 | 0.09492 | 0.09668 |

The LJ configurational spinodal is bracketed near $4.025$--$4.030$. The
hypothetical hybrid moves it to approximately $4.18$--$4.20$. The principal
normal-opening instability moves from 5.252 to 5.456, while the minimum opening
threshold over one registry period moves from 3.753 to 3.926.

Thus the embedding strengthens both registry corrugation and normal cohesion
for this parameter choice. It does not establish slip before opening: the
minimum registry-dependent opening threshold remains below the configurational
spinodal. Event ordering remains a coupled numerical/first-passage question,
not an empirical yield criterion.

## 15. Fast-normal free-energy reduction

With generalized energy

$$
G(a,s;f)=W_{\mathrm{hybrid}}(a,s)
-f[(a-a_{0,h})+\chi s],
$$

the existing fixed-boundary reduction remains

$$
\mathcal F_{\mathrm{eff}}(s,f)
=-k_BT\ln\int e^{-G(a,s;f)/(k_BT)}\,da,
$$

$$
\partial_s\mathcal F_{\mathrm{eff}}
=\langle G_s\rangle_{a|s,f}.
$$

The implementation now obtains local energy and derivatives through the common
energy-surface interface, so the embedding term is included. A hybrid numerical
check gives absolute gradient-identity error $3.86\times10^{-9}$. For moving
bound-basin limits, the previously derived Leibniz boundary correction remains
mandatory and unchanged.

## 16. Limiting cases

The implementation verifies:

- $A\to0$ recovers `TwoRowLJ` exactly;
- $\varrho_0\to0$ recovers `TwoRowLJ` exactly;
- increasing $\beta_\varrho$ localizes the relative neighbor contribution;
- $\varrho_\perp\to0$ exponentially as $a\to\infty$;
- $\varrho(a,s+b)=\varrho(a,s)$ and
  $W_{\mathrm{hybrid}}(a,s+b)=W_{\mathrm{hybrid}}(a,s)$;
- the Hessian is symmetric and the declared pristine state is stable for the
  hypothetical sensitivity set.

## 17. Mathematical status

**Exact analytic under the declared infinite two-row geometry**

- same-row geometric density sum;
- Fourier-transform identity;
- infinite Poisson representation;
- LJ pair lattice identity already present;
- analytic density, embedding, gradient, and Hessian formulas;
- two-site embedding multiplicity.

**Semi-analytic**

- tolerance-truncated reciprocal Bessel evaluation with returned diagnostics;
- high-precision equilibrium roots;
- minimum/saddle and spinodal continuation;
- bound-basin quadrature for $\mathcal F_{\mathrm{eff}}$.

**Numerical**

- finite-volume probability evolution;
- opening first-passage time integration;
- parameter sensitivity and event-ordering studies.

The tolerance-truncated reciprocal series is not described as a finite closed
form.

## 18. Parameters and future constraints

| Group | Parameter | Role | Possible future constraint |
|---|---|---|---|
| LJ pair | $\epsilon_{\mathrm{LJ}}$ | pair energy scale | cohesive/decohesion energy |
| LJ pair | $\sigma_{\mathrm{LJ}}$ | pair length scale | equilibrium spacing |
| density | $\varrho_0$ | density amplitude | atomistic site-density convention |
| density | $\beta_\varrho$ | radial decay | environment/vacancy response |
| density | $r_e$ | reference decay length | lattice geometry |
| embedding | $A$ | many-body energy scale | cohesive and elastic response |
| embedding | $\varrho_{\mathrm{ref}}$ | density reference | chosen EAM gauge/convention |
| geometry | $b$ | row periodicity | declared crystal plane/direction |
| geometry | $\chi$ | axial slip projection | crystallographic orientation and plane spacing |

Future aluminum work must jointly constrain equilibrium spacing, cohesion,
elastic constants, vacancy/environment response, generalized stacking-fault
corrugation, and decohesion. A published EAM/MEAM/GSF comparison is required
before quantitative FCC-Al claims. Parameters must not be tuned merely to make
interwell slip precede opening.

## 19. Current limitations

- Only the uniform N=1 two-row cell has exact equal-site density bookkeeping.
- Surface, vacancy, finite-row, and multi-site environments may have
  inequivalent $\varrho_i$ and require a new counting derivation.
- The square-root family is illustrative and uncalibrated.
- EAM gauge freedom prevents identification with an existing aluminum
  potential without a complete matched parameterization.
- No empirical hardening, yield cutoff, characteristic volume, Monte Carlo
  estimator, or lifetime fit has been introduced.
- The statistical correlation area $A_c$ remains solely specimen-probability
  post-processing and never appears in $f^*=\kappa_{\mathrm{axial}}\sigma/E$.
