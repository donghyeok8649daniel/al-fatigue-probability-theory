# Finite-temperature fast-normal-coordinate reduction

## Canonical status

The full $P(a,s,t)$ Smoluchowski equation remains the N=1 reference. The
reduction in `reduced_fast_a.py` is experimental and is validated here only for
strain and configurational dynamics when opening loss is negligible. It does
not provide a production-valid reduced crack probability.

The load mapping remains

$$
f^*=\kappa_{\mathrm{axial}}\frac{\sigma}{E}.
$$

No geometric multiplier enters this mapping.

## Full equation and scale separation

The implemented full equation is

$$
\partial_t P
=\partial_a\!\left[M_a\left(PG_a+k_BT\,\partial_aP\right)\right]
+\partial_s\!\left[M_s\left(PG_s+k_BT\,\partial_sP\right)\right],
$$

The configurational marginal used below is denoted $p_s(s,t)$.

Let $\varepsilon_m=M_s/M_a$ and use the slow time $T=M_st$. Then

$$
\partial_TP
=\varepsilon_m^{-1}\mathcal L_aP+\mathcal L_sP.
$$

On a fixed normal interval with reflecting boundaries, the nullspace of
$\mathcal L_a$ is one-dimensional. The leading expansion is

$$
P_0(a,s,T)=p_s(s,T)\rho_{\mathrm{eq}}(a\mid s,f),
$$

with

$$
\rho_{\mathrm{eq}}(a\mid s,f)
=\frac{\exp[-G(a,s;f)/k_BT]}{Z_a(s,f)},
\qquad
Z_a(s,f)=\int_{a_L}^{a_U}\exp[-G(a,s;f)/k_BT],da.
$$

Define the effective configurational free energy

$$
\mathcal F_{\mathrm{eff}}(s,f)=-k_BT\ln Z_a(s,f).
$$

The symbol denotes a free energy only.

## Fixed-domain solvability condition

For fixed $a_L,a_U$, differentiation under the integral gives

$$
\partial_s Z_a
=-\frac{1}{k_BT}\int_{a_L}^{a_U}
G_s\exp[-G/k_BT],da,
$$

and therefore

$$
\partial_s\mathcal F_{\mathrm{eff}}
=\left\langle G_s\right\rangle_{\rho_{\mathrm{eq}}}.
$$

Integrating the order-one equation over $a$, using the fast no-flux boundary
condition and this identity, gives

$$
\partial_Tp_s
=\partial_s\!\left[
p_s\,\partial_s\mathcal F_{\mathrm{eff}}
+k_BT\,\partial_sp_s
\right].
$$

Returning to the repository's model time:

$$
\partial_tp_s
=\partial_s\!\left[
M_s\left(
p_s\,\partial_s\mathcal F_{\mathrm{eff}}
+k_BT\,\partial_sp_s
\right)
\right].
$$

The experimental solver applies a conservative Scharfetter--Gummel flux to
this equation.

## Moving bound-basin derivative

For

$$
Z_b(s,f)=\int_{a_L(s,f)}^{a^{\dagger}(s,f)}
\exp[-G(a,s;f)/k_BT],da,
$$

Leibniz differentiation at fixed $f$ gives

$$
\begin{aligned}
\partial_s Z_b={}&
e^{-G(a^{\dagger},s;f)/k_BT}\,\partial_sa^{\dagger}
-e^{-G(a_L,s;f)/k_BT}\,\partial_sa_L \\
&-\frac{1}{k_BT}\int_{a_L}^{a^{\dagger}}
G_s e^{-G/k_BT},da.
\end{aligned}
$$

Thus

$$
\begin{aligned}
\partial_s\mathcal F_b
={}&\left\langle G_s\right\rangle_b \\
&-\frac{k_BT}{Z_b}
\left[
e^{-G(a^{\dagger},s;f)/k_BT}\,\partial_sa^{\dagger}
-e^{-G(a_L,s;f)/k_BT}\,\partial_sa_L
\right].
\end{aligned}
$$

The saddle condition $G_a(a^{\dagger},s;f)=0$ does not make the boundary term
zero. For the present potential,

$$
\partial_sa^{\dagger}
=-\frac{G_{as}(a^{\dagger},s;f)}
{G_{aa}(a^{\dagger},s;f)}.
$$

At $s=0.17,f^*=3.5$, the numerical values are

| Term | Value |
|---|---:|
| $\langle G_s\rangle_b$ | 0.0115262637608 |
| moving-boundary correction | $-2.79096866\times10^{-5}$ |
| predicted $\partial_s\mathcal F_b$ | 0.0114983540742 |
| centered finite difference | 0.0114983542243 |

At higher barriers the correction is exponentially small. For example, at
$s=0.17,f^*=1.5$ its magnitude is about $2.9\times10^{-17}$. It cannot be
dropped as an identity; its smallness is a metastable asymptotic statement.

## Equilibrium, truncated Gibbs, and QSD

Three conditional objects have different meanings.

1. On a fixed reflecting interval, normalized Gibbs is the exact fast null
   state and produces the potential-of-mean-force equation above.
2. Gibbs truncated at the instantaneous opening saddle is a normalized
   metastable approximation. It is not an eigenfunction of an absorbing fast
   generator and need not vanish at the dividing surface.
3. With a reflecting lower boundary and absorbing opening boundary, the fast
   operator has no zero eigenvalue. Its leading right eigenfunction is a QSD,
   and its eigenvalue supplies conditional survival decay.

For the discrete fast Scharfetter--Gummel operator at $s=0.18$:

| $f^*$ | QSD escape rate for $M_a=1$ | QSD/truncated-Gibbs L1 difference |
|---:|---:|---:|
| 2.2 | $8.20\times10^{-9}$ | $3.95\times10^{-9}$ |
| 3.5 | 0.0370430 | 0.0114882 |

The truncated distribution is excellent at the higher barrier and visibly
different when the barrier falls. Therefore the actual absorbing singular
limit is the QSD, not truncated Gibbs.

If the killed fast operator has a well-separated principal eigenvalue and its
escape is slow on the conditional relaxation scale, projection may produce a
local loss term of the form

$$
\partial_tp_s=\mathcal L_{mathrm{slow}}p_s-k_{\mathrm{open}}(s,f)p_s.
$$

The drift then generally requires the left and right eigenfunctions of the
killed operator, not only $\mathcal F_b$. If the spectral separation fails or
the boundary moves on the fast relaxation scale, memory and nonlocal coupling
can remain. No such loss term is implemented because it has not yet been
validated against the full first-passage flux.

## Point equilibrium and Laplace approximation

The stable mechanical point is only a reference:

$$
G_a(a^*,s;f)=0,
\qquad
G_{aa}(a^*,s;f)>0.
$$

Because repository coordinates are nondimensionalized by $b$, the integration
measure $da$ is dimensionless. For an interior minimum far from the boundaries,

$$
Z_a
\sim e^{-G(a^*,s;f)/k_BT}
\sqrt{\frac{2\pi k_BT}{G_{aa}(a^*,s;f)}},
$$

so, up to an additive constant independent of $s$ and $f$,

$$
\mathcal F_{\mathrm{eff}}
\sim G(a^*,s;f)
+\frac{k_BT}{2}
\ln\!\left[
\frac{G_{aa}(a^*,s;f)}{2\pi k_BT}
\right].
$$

At $s=0.08,f^*=0.6$:

| $k_BT$ | $\langle a\rangle-a^*$ | exact-minus-Laplace free energy |
|---:|---:|---:|
| 0.02 | 0.00178744 | $-1.46458\times10^{-4}$ |
| 0.01 | 0.000883540 | $-3.62990\times10^{-5}$ |
| 0.005 | 0.000439303 | $-9.03632\times10^{-6}$ |

Both errors decrease systematically with temperature.

## Conditional strain

For either a justified equilibrium conditional or QSD, define

$$
\overline a(s,f)=\int a\rho(a\mid s,f),da.
$$

With $s=bn+\xi$ and $\xi\in[-b/2,b/2)$,

$$
\epsilon_{\mathrm{normal}}
=\int p_s(s,t)\frac{\overline a(s,f(t))-a_0}{a_0},ds,
$$

$$
\epsilon_{\mathrm{intrawell}}
=\int p_s(s,t)\frac{\chi\xi}{a_0},ds,
$$

$$
\epsilon_{\mathrm{plastic}}
=\int p_s(s,t)\frac{\chi bn}{a_0},ds,
$$

and their sum is the total axial strain. At the current $k_BT=0.02$, the
principal-well reduced equilibrium has a centered small-stress tangent of
1.06953 with respect to $\sigma/E$. The difference from the zero-temperature
calibration is finite-temperature anharmonic susceptibility, not a change to
$\kappa_{\mathrm{axial}}$.

The potential retains the exact tilt relation

$$
\mathcal F_{\mathrm{eff}}(s+b,f)
=\mathcal F_{\mathrm{eff}}(s,f)-f\chi b,
$$

while $\overline a(s+b,f)=\overline a(s,f)$.

## Quadrature verification

At $s=0.11,f^*=0.8$ on the fixed normal interval:

| Gauss order | $\mathcal F_{\mathrm{eff}}$ | $\langle a\rangle$ |
|---:|---:|---:|
| 32 | -1.8782178633 | 0.8154798828 |
| 64 | -1.8781945908 | 0.8155295426 |
| 128 | -1.8781945908 | 0.8155295426 |
| 256 | -1.8781945908 | 0.8155295426 |

Adding a constant $C$ to $G$ shifts $\mathcal F_{\mathrm{eff}}$ by $C$ while
leaving the conditional density, mean, and effective gradient unchanged.

## Resolved Case A singular limit

Case A uses -900 to -100 MPa. The full reference uses 81 by 91 cells and
backward Euler with the same SG generator. The reduced calculation uses a fixed
normal interval because compression creates no tensile opening saddle.

Reduced strain range: -0.00812405 to 0.00185782. Its maximum same-stress
loading/unloading difference is $1.37\times10^{-5}$.

At $M_a=1$, the resolved full-2D strain range is -0.00599902 to
-0.0000511318. Its normal, intrawell, and plastic ranges are respectively
-0.00538584 to 0.000553641, -0.000613772 to -0.000603799, and roundoff zero.

| $M_a$ | final marginal L1 | mean-$s$ error | max normal-strain error | max total-strain error | full-PDE loop difference |
|---:|---:|---:|---:|---:|---:|
| 1 | $4.1154\times10^{-4}$ | $2.0628\times10^{-6}$ | 0.00396935 | 0.00396908 | 0.00460139 |
| 2 | $2.0114\times10^{-4}$ | $1.0163\times10^{-6}$ | 0.00275879 | 0.00275857 | 0.00440200 |
| 5 | $4.1002\times10^{-5}$ | $2.1168\times10^{-7}$ | 0.00132789 | 0.00132770 | 0.00250803 |
| 10 | $7.0173\times10^{-6}$ | $2.5945\times10^{-8}$ | 0.000691525 | 0.000691355 | 0.00135857 |
| 20 | $4.6353\times10^{-6}$ | $2.9681\times10^{-8}$ | 0.000350321 | 0.000350158 | 0.000702826 |

The marginal and strain errors decrease until spatial/quadrature error begins
to dominate.

## Resolved Case B singular limit

Case B uses -1100 to 2900 MPa. Fixed-domain Gibbs is invalid here because the
tensile tilt makes the detached side of the reflecting interval dominate; it
produced strain above one and was rejected. The comparison therefore uses the
instantaneous bound-basin truncated conditional, with the QSD limitation stated
above.

Reduced strain range: -0.00723109 to 0.06709759. Its maximum same-stress
loading/unloading difference is 0.0001530.

At $M_a=1$, the resolved full-2D strain range is 0.00632932 to 0.04029748.
Its normal range is 0.00474062 to 0.03866022, its intrawell range is
0.00157526 to 0.00164377, and its plastic contribution is below
$1.4\times10^{-12}$.

| $M_a$ | final marginal L1 | mean-$s$ error | max normal-strain error | max total-strain error | full-PDE loop difference | survival |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.00273996 | $2.5097\times10^{-5}$ | 0.0335437 | 0.0335527 | 0.0286227 | 1.0000000000 |
| 2 | 0.00163602 | $1.6262\times10^{-5}$ | 0.0260219 | 0.0260281 | 0.0324236 | 0.9999999999 |
| 5 | 0.000552997 | $6.5202\times10^{-6}$ | 0.0155398 | 0.0155440 | 0.0224763 | 0.9999999997 |
| 10 | 0.000256527 | $2.5505\times10^{-6}$ | 0.00924257 | 0.00924502 | 0.0157359 | 0.9999999721 |
| 20 | 0.000167349 | $9.2786\times10^{-7}$ | 0.00502074 | 0.00502240 | 0.00942401 | 0.9999994552 |

The configurational marginal converges monotonically. The remaining total
strain error is dominated by finite normal relaxation and decreases with
$M_a$. The reduced same-stress loop is much smaller than the full loop.

In both cases intrawell errors decrease with the marginal error. Well-index
plastic strain remains at roundoff-to-$10^{-12}$ scale; no irreversible
interwell plasticity is established. Most full-PDE hysteresis at the present
period is mechanism-one normal phase lag. A small reversible configurational
loop remains after conditional elimination.

## Recommendation

The evidence supports status **B**:

- the full two-dimensional PDE remains canonical/reference;
- the experimental reduction is validated for strain and slow configurational
  dynamics in the tested negligible-opening-loss regimes;
- fixed-domain Gibbs must not be used across a tensile opening saddle;
- truncated Gibbs is only a metastable approximation to the QSD;
- reduced crack probability remains not production-valid.

Broader promotion requires direct validation of a QSD-projected loss and drift
against the full moving-boundary first-passage flux.
