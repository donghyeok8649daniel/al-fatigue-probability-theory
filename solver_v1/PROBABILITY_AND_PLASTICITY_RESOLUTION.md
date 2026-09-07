# Probability and plasticity resolution audit

This audit separates physical first passage, configurational well transfer,
and numerical residuals in the N=1 probability PDE. None of the quantities in
this document changes the verified force mapping
`f* = kappa_axial sigma/E`.

## Discrete physical bookkeeping

Only probability removed by the moving opening-basin mask contributes to
$A_{\mathrm{abs}}$. Numerical negative-density repair is accumulated in a
separate field. The canonical local fields are

$$
P_{\mathrm{init,local}}=A_{\mathrm{abs}},\qquad
S_{\mathrm{local}}=1-A_{\mathrm{abs}}.
$$

The independent mass check is

$$
R_{\mathrm{mass}}=M_{\mathrm{raw}}+A_{\mathrm{abs}}-1.
$$

Each reported first-passage flux is the opening mass removed during that record
interval divided by the interval duration. Consequently its discrete time
integral plus any initial mask removal reconstructs $A_{\mathrm{abs}}$; the
reported `flux_consistency_residual` verifies this identity without folding
$R_{\mathrm{mass}}$ into the physical probability.

## Configurational-well balance

For $s=bn+\xi$, $\xi\in[-b/2,b/2)$, the solver integrates $P_n$ over every
well. The Scharfetter--Gummel interface rates provide both signed net flux and
gross bidirectional activity at $s=(n+1/2)b$. The discrete check is

$$
P_n(t)-P_n(0)
=Q_{n-1/2}(t)-Q_{n+1/2}(t)-A_{\mathrm{open},n}(t)+R_n(t).
$$

$R_n$ is kept as `well_population_balance_residual`. A nonzero
`plastic_strain` is not classified as resolved transfer unless $P_n$, net/gross
boundary flux, and the strain all exceed their refinement-derived floors.
Residual plasticity additionally requires persistence after unloading and a
zero-load hold.

## September 2026 numerical audit

All values below are dimensionless model outputs. They are mechanism checks,
not calibrated pure-Al life or yield predictions.

For the exact current UI frequency of 25 cycles/model-time:

| Case | Stress range (MPa) | UI-grid absorbed mass | UI-grid max plastic strain | Classification |
|---|---:|---:|---:|---|
| Compression control | -900 to -100 | 0 | $1.86\times10^{-19}$ | Crack absent; plastic unresolved |
| Nonlinear mechanism test | -1100 to 2900 | $1.19\times10^{-13}$ | $4.88\times10^{-13}$ | Both unresolved |

The Case-B aligned-grid refinement gave:

| $(n_a,n_s)$ | Integrator / max $dt$ | Absorbed mass | Max plastic strain | Outside-well mass |
|---:|---:|---:|---:|---:|
| (15,24) | explicit / $10^{-3}$ | $7.89\times10^{-11}$ | $2.62\times10^{-12}$ | $1.01\times10^{-11}$ |
| (21,30) | explicit / $5\times10^{-4}$ | $5.20\times10^{-13}$ | $1.62\times10^{-12}$ | $6.27\times10^{-12}$ |
| (25,36) | explicit / $2.5\times10^{-4}$ | $1.89\times10^{-13}$ | $8.30\times10^{-13}$ | $3.20\times10^{-12}$ |

The observed local rare-event floor from this refinement is
$7.87\times10^{-11}$, larger than the finest absorbed signal. The plastic
strain floor is $1.79\times10^{-12}$, also larger than the finest signal.
Thus a displayed value near $10^{-16}$ or $10^{-11}$ in these runs is below
the demonstrated resolution and is not physical initiation or plasticity.

An aligned-grid constant-force audit showed the apparent plastic signal at
$f^*=3$ decreasing from $2.47\times10^{-11}$ on $(11,18)$ to
$4.42\times10^{-16}$ on $(19,30)$, while the discrete well balances remained
within numerical tolerance. This is grid-dependent well-boundary leakage, not
a resolved interwell transition.

A load/unload/zero-hold check likewise changed its final plastic strain from
$1.16\times10^{-12}$ on $(19,30)$ to $1.07\times10^{-13}$ on $(23,36)$.
Persistence is therefore not convergence-resolved and cannot be called
residual plasticity.

## Positive opening control and onset interpretation

At constant $f^*=5.1$ for 0.1 model time, opening absorption was
$7.72\times10^{-2}$, $8.83\times10^{-2}$, and $1.08\times10^{-1}$ on three
successively refined aligned grids. The signal is well above mass residual and
the observed grid variation, so the absorption mechanism is demonstrably
active. This corresponds to roughly 4.08 GPa through the present mapping and
is explicitly an extreme mechanism control.

The fixed-domain potential-of-mean-force profile loses its forward
configurational barrier near $f^*\approx2$, before the pristine normal opening
spinodal near $f^*=5.25$. However, the scanned finite-time full-2D runs did not
show a convergence-resolved interwell population transfer before opening loss
became resolved. Static barrier loss and dynamically resolved transfer are
therefore reported as separate indicators; no experimental Al yield stress is
claimed.

The one-cycle, one-model-time mechanism sweep was:

| Stress amplitude (MPa) | Peak $f^*$ | Min forward configurational barrier | Center opening barrier | Max outside-well mass | Absorbed mass |
|---:|---:|---:|---:|---:|---:|
| 50 | 0.0625 | 0.709 | 1.997 | $1.76\times10^{-15}$ | 0 |
| 150 | 0.1876 | 0.679 | 1.846 | $1.78\times10^{-15}$ | 0 |
| 1000 | 1.2506 | 0.424 | 1.094 | $3.28\times10^{-15}$ | $8.35\times10^{-25}$ |
| 1500 | 1.8759 | 0.0950 | 0.799 | $2.71\times10^{-14}$ | $2.72\times10^{-18}$ |
| 2000 | 2.5012 | 0 | 0.559 | $1.22\times10^{-12}$ | $3.63\times10^{-13}$ |
| 2500 | 3.1266 | 0 | 0.364 | $6.20\times10^{-11}$ | $6.08\times10^{-9}$ |
| 3000 | 3.7519 | 0 | 0.208 | $1.93\times10^{-9}$ | $1.52\times10^{-5}$ |
| 3500 | 4.3772 | 0 | 0.0898 | $2.37\times10^{-9}$ | $6.59\times10^{-3}$ |

These sweep values locate candidate mechanisms but do not by themselves certify
the very small well populations. The dedicated aligned-grid refinement above
classifies their interwell/plastic signal as unresolved. Thus the model's
static configurational barrier indicator changes first, while a dynamically
resolved configurational-transition onset was not established in the scanned
range before resolved opening loss.

## Statistical specimen aggregation

Three distinct notions are retained:

1. atomic geometry in the Bessel--LJ energy;
2. statistical correlation area $A_c$ supplied by external calibration;
3. effective stressed specimen area $A_{\mathrm{stressed}}$.

Under an unvalidated independent-equivalent-region approximation,

$$
N_{\mathrm{eff}}=A_{\mathrm{stressed}}/A_c,
$$

$$
P_{\mathrm{init,spec}}
=-\mathrm{expm1}\!\left[
N_{\mathrm{eff}}\mathrm{log1p}(-P_{\mathrm{init,local}})
\right].
$$

The area inputs are post-processing only. They do not scale strain, plastic
activity, the energy, or generalized force. If the local signal lacks a
convergence certificate, the physical specimen probability remains unavailable
regardless of $N_{\mathrm{eff}}$; only a clearly labeled mathematical
extrapolation is retained.
