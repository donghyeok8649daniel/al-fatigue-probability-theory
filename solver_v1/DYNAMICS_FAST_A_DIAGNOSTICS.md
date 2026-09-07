# Model time and fast-normal-coordinate diagnostics

## Scope and status

This document separates four questions which must not be conflated:

1. the verified static stress/strain calibration;
2. the finite-mobility dynamics of the current full two-dimensional PDE;
3. the fast-normal-coordinate reference at fixed configurational state;
4. a future conversion from model time to physical seconds.

The calculations below are diagnostics. They do not replace the canonical N=1
probability PDE, change the Bessel/Poisson-summed energy, tune a mobility, or
calibrate physical fatigue frequency.

## Static calibration

At the pristine state, let

$$
q_0=(a_0,0),
\qquad
\mathbf c=(1,\chi)^T,
\qquad
H_0=\nabla_q^2 W(q_0).
$$

The verified relaxed tangent mapping remains

$$
f^*=\kappa_{\mathrm{axial}}\frac{\sigma}{E},
\qquad
\kappa_{\mathrm{axial}}
=\frac{a_0}{\mathbf c^T H_0^{-1}\mathbf c}.
$$

| Quantity | Value |
|---|---:|
| $a_0$ | 0.7713438268704838 |
| $W_{aa}$ | 122.78696473236548 |
| $W_{as}$ | $-1.8791267641884867\times10^{-14}$ |
| $W_{ss}$ | 50.347580146252554 |
| $\kappa_{\mathrm{axial}}$ | 86.29296488740997 |

This is a zero-frequency, fully relaxed tangent statement. It does not require
the finite-frequency PDE strain to equal $\sigma/E$ at every instant.

## Current finite-mobility linear diagnostic

Linearizing the current overdamped two-coordinate model gives

$$
\dot{\delta q}=-M H_0\delta q+M\mathbf c\,\delta f(t),
\qquad
M=\mathrm{diag}(M_a,M_s).
$$

The positive relaxation rates are the eigenvalues of the full symmetric matrix

$$
R=M^{1/2}H_0M^{1/2}.
$$

The implementation does not assume pure-$a$ or pure-$s$ eigenvectors. For
$M_a=1$ and $M_s=0.05$:

| Quantity | Value in model-time units |
|---|---:|
| $\lambda_{\mathrm{fast}}$ | 122.78696473236548 |
| $\lambda_{\mathrm{slow}}$ | 2.517379007312628 |
| $\tau_{\mathrm{fast}}$ | 0.00814418698417756 |
| $\tau_{\mathrm{slow}}$ | 0.39723855529705393 |

For

$$
\delta f(t)=\mathrm{Re}\!\left[\widehat f e^{i\omega t}\right],
$$

the current finite-mobility model has

$$
\widehat{\delta q}
=\left(i\omega I+MH_0\right)^{-1}M\mathbf c\,\widehat f,
$$

and

$$
G(\omega)
=\frac{\widehat\epsilon}{\widehat\sigma/E}
=\frac{\kappa_{\mathrm{axial}}}{a_0}
\mathbf c^T\left(i\omega I+MH_0\right)^{-1}M\mathbf c.
$$

Consequently $G(0)=1$. The API is in `dynamics_diagnostics.py`.

The desktop load period is 0.04 model time, so its model frequency is 25
cycles per model-time unit and $\omega=157.07963267948966$. Thus

$$
\omega\tau_{\mathrm{fast}}=1.279285899947692,
\qquad
\omega\tau_{\mathrm{slow}}=62.39808635219237.
$$

This is not the quasi-static regime of the current finite-$M_a$ equations.

## Analytic response compared with the probability PDE

The table uses reduced-stress amplitude $10^{-4}$, discards two transient
cycles, and fits the next two cycles. The PDE uses the same
Scharfetter--Gummel spatial operator with backward Euler for this stiff
diagnostic, 100 steps per cycle, and a refined one-well grid. Frequency is
reported by $De_s=\omega\tau_{\mathrm{slow}}$, not in hertz.

| $De_s$ | $\omega$ | analytic $\lvert G\rvert$ | analytic phase | PDE amplitude ratio | PDE phase |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.2517379007 | 0.999173 | -0.6118 deg | 1.066870 | -1.1237 deg |
| 1 | 2.5173790073 | 0.957259 | -3.7803 deg | 0.984926 | -4.7044 deg |
| 10 | 25.1737900731 | 0.895223 | -12.1267 deg | 0.910460 | -13.9672 deg |
| 62.3981 | 157.0796326795 | 0.562256 | -52.0733 deg | 0.514764 | -54.9444 deg |

The probability PDE includes finite-temperature and grid-resolved effects
absent from the point-Hessian transfer function. Agreement is within about 8.5
percent in amplitude and 3 degrees in phase over this diagnostic range.

At the desktop model frequency, backward-Euler refinement gave:

| Steps/cycle | $\Delta t$ | amplitude ratio | phase | maximum mass residual |
|---:|---:|---:|---:|---:|
| 50 | 0.0008 | 0.507925 | -53.7068 deg | $1.6\times10^{-15}$ |
| 100 | 0.0004 | 0.514764 | -54.9444 deg | $2.9\times10^{-14}$ |
| 200 | 0.0002 | 0.518335 | -55.5634 deg | $6.6\times10^{-14}$ |

The response converges under time-step refinement. Its difference from the
static tangent is therefore model dynamics, finite-temperature response, and
spatial discretization rather than an unresolved time step.

## Fast normal mechanical equilibrium

For fixed $s$ and $f^*$, the point mechanical reference follows

$$
W_a(a^*,s)-f^*=0,
\qquad
W_{aa}(a^*,s)>0.
$$

Given the current configurational marginal $P_s(s,t)$, its low-temperature
strain reference is

$$
\epsilon_{a^*}(t)
=\int
\frac{a^*(s,f^*(t))-a_0+\chi s}{a_0}
P_s(s,t)\,ds.
$$

For the finite-temperature Smoluchowski PDE, the exact adiabatic limit is a
conditional equilibrium rather than a forced delta distribution:

$$
P_{\mathrm{fast}}(a,s,t)
=P_s(s,t)
\frac{\exp[-G(a,s;f^*(t))/k_BT]}
{Z_a(s;f^*(t))}.
$$

The normal partition integral is over the intact normal basin. This projection
preserves $P_s$ exactly. The point branch is its low-temperature Laplace
approximation. Both references are implemented, and loss of a stable/intact
normal branch is reported rather than silently renormalized.

If an adiabatic one-dimensional slow equation is later validated, its exact
finite-temperature free energy is

$$
A_{\mathrm{eff}}(s;f^*)=-k_BT\ln Z_a(s;f^*),
$$

not merely $G(a^*,s;f^*)$ except in the low-temperature approximation.

## Case A: compression-only cycle

The stress range is -900 to -100 MPa.

| Response | Total-strain range |
|---|---:|
| Current desktop preview, full PDE, first cycle | -0.0094821 to -0.0090780 |
| Stable $a^*$ at fixed $s=0$ | -0.0110099 to -0.00130855 |
| Fully relaxed stable $(a,s)$ branch | -0.0119233 to -0.00143370 |
| Stable $a^*$ averaged over initial $P_s$ | -0.00973449 to -0.00003358 |
| Refined finite-temperature conditional fast-$a$ | -0.00811349 to 0.00186643 |
| Refined full PDE, $M_a=1$, first cycle | -0.00626829 to 0.00016230 |

The positive upper endpoint in finite-temperature refined references is an
anharmonic thermal mean shift relative to zero-temperature $a_0$; it does not
reverse the compressive stress.

The desktop preview has 21 normal cells over a broad domain. Its normal-cell
width is about 0.046, while the stable point branch moves only about 0.0075 in
Case A. It therefore under-resolves the constitutive displacement.

On a refined grid, increasing only $M_a$ gives this mathematical singular-limit
test. $M_s$, energy, load, and temperature remain fixed.

| $M_a$ | full-PDE strain range | conditional fast-$a$ range | maximum deviation |
|---:|---:|---:|---:|
| 1 | -0.0060885 to 0.0000301 | -0.0081135 to 0.0018664 | 0.0040209 |
| 2 | -0.0072441 to 0.0009367 | -0.0081134 to 0.0018668 | 0.0027884 |
| 5 | -0.0079161 to 0.0016129 | -0.0081138 to 0.0018673 | 0.0013312 |
| 10 | -0.0080484 to 0.0017827 | -0.0081140 to 0.0018676 | 0.0006872 |

This monotone convergence directly shows that finite normal mobility causes
much of the remaining lag. It is not a calibration or selection of $M_a$.

## Case B: nonlinear mechanism stress test

The stress range is -1100 to 2900 MPa.

| Response | Total-strain range |
|---|---:|
| Current desktop preview, full PDE, first cycle | 0.00614993 to 0.02718052 |
| Stable $a^*$ at fixed $s=0$ | -0.01324835 to 0.05615692 |
| Fully relaxed stable $(a,s)$ branch | -0.01431324 to 0.07389977 |
| Stable $a^*$ averaged over initial Case-B $P_s$ | -0.00887200 to 0.06077012 |
| Refined finite-temperature conditional fast-$a$ | -0.00723242 to 0.06706601 |
| Refined full PDE, $M_a=1$, first cycle | 0.00514961 to 0.04110368 |

After four cycles the refined $M_a=1$ PDE still spans 0.00456 to 0.03781. At
the two nearly zero-stress crossings, its strains were about 0.03003 and
0.00464. Replacing only the normal conditional state by fast equilibrium
reduced those values to about 0.00635 and 0.00627. The large loop is therefore
reversible finite-$M_a$ phase lag, not plastic hysteresis.

| $M_a$ | full-PDE strain range | conditional fast-$a$ range | maximum deviation |
|---:|---:|---:|---:|
| 1 | 0.0051496 to 0.0411037 | -0.0072324 to 0.0670660 | 0.0335102 |
| 2 | -0.0021915 to 0.0496306 | -0.0072237 to 0.0670722 | 0.0259908 |
| 5 | -0.0063658 to 0.0592780 | -0.0072208 to 0.0670820 | 0.0153630 |
| 10 | -0.0070016 to 0.0636568 | -0.0072211 to 0.0670884 | 0.0092257 |

No well-index transfer occurred in these one-well isolation checks. In the
actual three-well desktop run, plastic strain was below $5\times10^{-13}$ and
survival loss was about $1.2\times10^{-13}$. The refined four-cycle one-well
test had zero well activity. The small remaining fast-$a$ loading/unloading
difference at equal stress, about $7.5\times10^{-5}$, came from intrawell
configurational memory. It is reversible unless unloading and relaxation
demonstrate persistence. No residual plasticity is established here.

## Interpretation and architecture recommendation

The narrow desktop response is not a load-mapping error: the verified relaxed
$\kappa_{\mathrm{axial}}$ is used. It is not primarily an unresolved time step.
Two dominant causes were identified:

1. the broad 21-cell normal preview grid under-resolves equilibrium motion;
2. the current full PDE has $\omega\tau_a\approx1.28$, so finite $M_a$ creates
   a large reversible phase lag.

The preload Gibbs convention removes a zero-preload jump. Multi-cycle checks
show that initial transient is secondary.

The intended fast/slow architecture is plausible and its singular limit is
demonstrated, but adiabatic elimination is not accurate for the current default
finite-$M_a$ PDE at period 0.04. The production PDE is therefore not replaced.
A slow-$s$ solver should use $A_{\mathrm{eff}}$, retain moving opening
first-passage accounting, and be validated against full 2D calculations in a
resolved fast-$a$ regime before becoming canonical.

Hysteresis labels must remain distinct:

- finite-$M_a$ phase lag is reversible dynamic lag;
- intrawell $P_s$ memory can create reversible configurational lag;
- persistent well-index transfer after unload/relaxation is residual
  plasticity;
- survival loss is crack-opening first passage.

## Why model time is not seconds

For dimensional length coordinates and energy, overdamped mobility has units

$$
[M_{\mathrm{dim}}]
=\frac{\mathrm{length}^2}{\mathrm{energy}\,\mathrm{time}}.
$$

With $q_{\mathrm{dim}}=bq^*$, $G_{\mathrm{dim}}=\varepsilon_{\mathrm{LJ}}G^*$,
and $t_{\mathrm{dim}}=t_0t^*$,

$$
M^*=\frac{t_0M_{\mathrm{dim}}\varepsilon_{\mathrm{LJ}}}{b^2},
\qquad
t_0=\frac{M^*b^2}{M_{\mathrm{dim}}\varepsilon_{\mathrm{LJ}}}.
$$

The repository declares dimensionless $b$, LJ energy, and mobilities. It does
not supply a validated dimensional mobility or friction for this reduced
coordinate. Model time to seconds is therefore unresolved. It requires
atomistic relaxation data, measured dynamic modulus/phase, experimental
relaxation time, or another justified mobility/friction source. A specimen
characteristic length, area, or volume is not introduced for this microscopic
time nondimensionalization.
