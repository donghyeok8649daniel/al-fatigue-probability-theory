# Scoped aluminum line kinetics and fatigue-reference validation (v18)

Status: executed research benchmark, **not a production a/s mobility calibration**.
The scoped numerical checks below passed. Actual specimen yield, fatigue and
the production clock remain unvalidated; v17's material rejection is unchanged.

## 1. Three observables, not one rescaling

An ideal coherent interface fold, the stress at a specified macroscopic
plastic-strain criterion, and specimen fatigue crack initiation are different
observables. Changing a time unit cannot turn one into another. We preserve
the LJ pair, infinite Poisson/Bessel sums, all static fit vectors, and the
production Smoluchowski generator. No strength factor, fatigue-life fit,
empirical yield law, or specimen correlation area enters this study.

The first new kinetic benchmark is the measured position X of an existing
dislocation line. It is **not** the normal separation a or cell registry s.
For the low-speed source experiment,

    v = dX/dt = mu_tau * tau,
    [mu_tau] = m/(Pa s),
    B_line = b/mu_tau,       [B_line] = Pa s = J s/m^3.

We first fit mu_tau directly to measured velocity/stress pairs; this needs
no invented Burgers vector. B_line may be derived only with an explicit,
condition-labelled b. A finite straight segment of length ell has friction
B_line*ell and translational mobility 1/(B_line*ell) [m^2/(J s)]. Neither ell
nor this coordinate is the production two-row cell. There is no conversion
from these data to M_a,phys or a production t0.

## 2. Small-bow reference derived from the existing lattice elasticity

For the previously derived pinned line energy, let theta0 be its straight
character and T_line=(gamma+gamma'')(theta0)>0 the line stiffness [J/m].
To second order in bow slope, x in [0,L], fixed y(0)=y(L)=0:

    E[y] = (T_line/2) integral (y_x)^2 dx - tau(t) b integral y dx.
    Rayleigh[y_t] = (B_line/2) integral (y_t)^2 dx.
    B_line y_t = T_line y_xx + tau(t) b.

The stiffness, not gamma alone, follows from the anisotropic Schur kernel of
the unchanged LJ/Bessel bulk Hessian. R/r_core, pin span, line character,
temperature and source drag are explicit assumptions; none is fitted to yield
or fatigue. Core and finite nonlocal energy parts are still absent.

For y=sum_n y_n sin(n*pi*x/L),

    lambda_n = (T_line/B_line) (n*pi/L)^2,
    tau_n = B_line L^2/(T_line*pi^2*n^2),
    dy_n/dt + lambda_n y_n = 4 b tau(t)/(B_line*n*pi),  n odd,

with zero forcing of even modes. The static mean bow is
mean(y)=b*tau*L^2/(12*T_line). With exp(i omega t) convention, its normalized
transfer is

    G_line(omega) = (96/pi^4) sum_(n odd>=1)
                       [1/(n^4*(1+i*omega*tau_1/n^2))].

G_line(0)=1. Finite mode truncation has absolute transfer tail no greater
than (96/pi^4)/(3*N^3), where N is the highest included odd mode. This is a
bound, not a fitted accuracy. A spatial conservative discretization and time
refinement independently check this analytic response. A zero-load hold
decays mode n by exp(-t/tau_n); purely subcritical bowing leaves no permanent
registry change. In a periodic state work input equals nonnegative viscous
dissipation. This is reversible/anelastic motion, **not fatigue damage**.

Physical seconds in this separately labelled reference follow measured line
motion and declared geometry, not atomic mass or an invented a/s damping.
Its elastic/kinetic source mismatch must be displayed; a hypothetical pinned
span is not a measured source. Its time scale must never unlock the UI Hz mode.

### Dissipation identity

Multiplying the line equation by y_t, integrating and using fixed endpoints,

    d/dt [(T_line/2) integral y_x^2 dx]
        = tau b integral y_t dx - B_line integral y_t^2 dx.

For a periodic harmonic steady state with load amplitude p_hat=b*tau_hat,

    work/cycle = -pi p_hat Im[integral y_hat dx]
               = pi omega B_line integral |y_hat|^2 dx >= 0.

The finite-grid tridiagonal calculation checks both sides separately. An
observed loop from this equation is viscous dissipation of a reversible bow;
its area is not a crack-initiation probability or an irreversible slip count.

## 3. Validation contracts

- Velocity calibration: retain source pixels/axes, temperature, purity,
  dislocation character, selection rules, point scatter and held-out points.
  Pixel reading error is not experimental uncertainty. A fitted line is not
  an independent validation of all source microstructures.
- Yield: retain the 0.002 plastic-shear CRSS benchmark and source geometry
  requirements. First microscopic motion/bow-out is not that criterion.
- Fatigue: distinguish first detected crack, acoustic crack-growth marker,
  and final fracture. Report range versus amplitude, control mode, specimen,
  frequency, environment and censoring. Never fit mobility to fatigue life.
- Production: physical a/s mobility and t0 remain unavailable unless their
  own matched-coordinate calibration passes. Source laboratory seconds are
  real seconds but not automatically solver seconds.

## 4. Executed experimental velocity calibration

The audited Gorman author PDF has SHA256
`b4536b12a161babc565460facb83b4c64a685bced57cf0f48da9f7cb79df44d3`.
We extracted its native 2490 x 3305 page-5 image without resampling. Figure5c
contains a 23 deg C velocity-versus-resolved-shear panel. The transcribed
axes/pixels, reading allowance, source conditions and selection rule are in
`data/gorman1969_velocity_benchmark.json`. The source image is cached locally,
not republished as project data; extraction/replay is optional research tooling.

An affine two-axis map handles the small skew of the scanned axes. In particular
`10^6 dyn/cm^2 = 10^5 Pa`, not 1 MPa. All ten isolated 7--15-pixel rings in the
predeclared ROI are recovered by deterministic connected components. Marks
overlapping axes/other marks/the original fitted line are excluded rather than
invented. This is a selected subset, **not the complete source population**.

The x-then-y ordered indices1,4,7 were held out before fitting. Seven remaining
points fit v=mu_tau*tau by zero-intercept least squares. The source's own
mean(v)/mean(tau) construction and an all-point OLS refit are sensitivity
diagnostics, not extra held-out data.

| Quantity | Executed value | Interpretation |
|---|---:|---|
| Temperature | 296.15 K, source accuracy +/-5 K | not the 0 K static target |
| Resolved shear range | 0.394321--1.383478 MPa | existing leading line translation |
| mu_tau | 1.1627233044e-5 m/(Pa s) | seven-point scoped estimate |
| Training RMSE | 2.383294 m/s | substantial physical/selection scatter |
| Held-out RMSE | 3.501185 m/s | 23.62% of held-out mean velocity |
| Leave-one-out slopes | 1.104298e-5--1.218419e-5 | empirical sensitivity, not CI |
| Reading-only slope bound | 1.040818e-5--1.299157e-5 | conservative pixel perturbations only |
| All-ten-point OLS | 1.257187e-5 | does not replace the held-out fit |
| All-ten-point mean ratio | 1.261353e-5 | source-style construction on our subset |

The pixel halfwidths correspond to +/-0.0320923 MPa and +/-0.323312 m/s.
Held-out errors are much larger: image precision is not the limiting error.
Leading-line selection, unresolved line character, pinning and source stress
systematics remain. There is no claim of a unique exact Al mobility or an
independently validated specimen-wide drag law. The original source pools
edge/mixed dislocations; it does not identify the production a/s coordinates.

For an explicitly fixed **model 0 K b=2.8637824638 Angstrom**, the corresponding
conditional B_line=b/mu_tau is 2.4629956698e-5 Pa s. The b is not a newly
measured 23 deg C Burgers vector. Only mu_tau is fitted directly to the source.

## 5. Stress, frequency, grid and hold experiments

`run_line_timescale_validation.py` reads, but does not refit or promote, the
unchanged v17 nested candidate's bulk114/62/32 GPa tensor. Its corrected
FCC basis feeds the existing nonlocal LJ/Bessel Schur reference. The material
candidate still has `material_accepted=false`. This reuse constrains a line
elasticity comparison; it does not validate its rejected interface surface.

We choose explicit edge character, hypothetical L=0.2/1/5 micrometers,
R=L and r_core=2b. Those geometric choices are sensitivity inputs, not measured
source sizes. No value is adjusted to yield/fatigue. Angular32/64-mode
stiffness changes are below7e-21 J/m. A fully atomistic core or loop nucleation
barrier is not computed by this long-wave, outer-line calculation.

| Hypothetical L [um] | T_line [J/m] | Conditional tau_1 [s] | Lag at100 Hz [degrees] |
|---:|---:|---:|---:|
| 0.2 | 5.21577069e-10 | 1.91383902e-10 | -6.79998e-6 |
| 1 | 6.64934572e-10 | 3.75305552e-9 | -1.333482e-4 |
| 5 | 8.08292074e-10 | 7.71854768e-8 | -2.7424444e-3 |

The source-derived scalar drag with these hypotheses gives essentially relaxed
bowing at0.1/1/25/100 Hz. Its tiny attenuation can lie below the series tail;
the positive mode-truncation deficit must not be called physical dissipation.
The much higher frequencies chosen to make omega*tau_1=1 are **mathematical
rolloff tests**, not validated broadband Al kinetics. Source pulse velocities
alone do not establish frequency-independent drag at those frequencies, and
no inertia/core correction has been inferred. The 8-cycle high-frequency run
tests the declared overdamped equation, not an experimental fatigue life.
The dynamic perturbation amplitude0.01MPa is below the digitized experimental
range: extending the source's zero-intercept drag to that amplitude is an
explicit linear-reference extrapolation, not a validated low-stress mobility.

At omega*tau_1=1 the analytic infinite-series |G|=0.708055447,
phase=-44.261565 degrees; highest odd mode1487 has tail bound9.9912e-11.
The independent spatial refinement at L=1 um gives:

| Segments | Complex transfer absolute error | Phase error [degrees] |
|---:|---:|---:|
| 16 | 4.214754e-3 | -0.128729 |
| 32 | 1.054548e-3 | -0.032245 |
| 64 | 2.636903e-4 | -0.008065 |
| 128 | 6.592588e-5 | -0.002016 |

Errors decrease as dx^2. Across all60 spatial cases the maximum relative
work/dissipation residual is3.07e-15. A separate implicit time solve, using
128segments,8cycles then at least20tau_1 of zero-stress hold, yields:

| Steps/cycle | Measured amplitude ratio | Phase [deg] | Complex absolute error |
|---:|---:|---:|---:|
| 64 | 0.691799284 | -42.854832 | 2.365422e-2 |
| 128 | 0.699713386 | -43.561609 | 1.198039e-2 |
| 256 | 0.703807043 | -43.913221 | 6.038961e-3 |
| 512 | 0.705888886 | -44.088561 | 3.041534e-3 |

This is first-order dt convergence, not exact finite-dt agreement. Finest
relative complex error is0.430%; the final hold mean bow is4.10813e-19 m,
1.14463e-9 of the static mean. The finite exponential tail is retained, not
clipped. The model produces no permanent slip and defines no opening
probability (null, not zero).

The separate static stress sweep uses0.01/0.1/1/4 MPa shear. Linear versus
nonlinear anisotropic graph discrepancies grow with bow slope: for1um at
0.01/0.1/1/4MPa they are0.00153/0.01455/1.2154/9.1562%. At5um,4MPa the
graph is outside its outer-model fold; no static bow is falsely reported.
Outer-only folds38.491/9.814/2.386 MPa for0.2/1/5um demonstrate geometry
sensitivity, **not calibrated macroscopic yield strengths**. Defect density,
source activation/multiplication, core energy and the0.002 plastic-shear
criterion remain missing from a specimen yield comparison. The held-out
Krebs wire data already stored in v15 are not refitted using an assumed L=D/3
or an invented shear modulus.

## 6. Real fatigue data: endpoint and unit audit

`data/aluminum_fatigue_validation_v18.json` transcribes all seven Al groups
from Deschanel2017 Table1:99.95%Al, coarse elongated grains, room temperature
(no fabricated numericK), fully reversed tests. Sample replicates are pooled;
reported ranges are not confidence intervals or independently paired rows.
The text defines Delta as max minus min despite use of the word amplitude.
We therefore explicitly store both range and half-range amplitude.

| Control | Source full range | Amplitude | f [Hz, experiment] | Fracture cycles |
|---|---:|---:|---:|---:|
| Total axial strain | 0.005 | 0.0025 | 0.1 | 11600--15500 |
| Total axial strain | 0.0095 | 0.00475 | 0.1 | 3000--3800 |
| Total axial strain | 0.015 | 0.0075 | 0.1 | 580 |
| Axial stress | 50 MPa | 25 MPa | 0.1 | 17200 |
| Axial stress | 62 MPa | 31 MPa | 0.1 | 5300 |
| Total axial strain | 0.0076 | 0.0038 | 1 | 1360--5040 |
| Total axial strain | 0.015 | 0.0075 | 1 | 1450 |

Specimen dimensions and rates are in the machine-readable records. Constant
rate2*f*Delta agrees with the printed rates within their stated digit rounding;
we do not silently replace them by sinusoidal cycling. Source seconds=N/f
are valid experimental units, not a calibration of our model clock.

Interrupted microscopy reports PSBs by100cycles and crack signs by600 under
Delta epsilon=.0095; AE multiplets appear no earlier than1200cycles. These
are different observations, not identical first-passage endpoints. We retain
the AE/failure fractions but do not infer exact first-crack times or local
probability from them. In particular,17200 and5300 are final **fracture**
cycles, not local atomic absorption times or fitted fatigue lifetimes.

`fatigue_comparison_gate` refuses the present comparison for six independent
reasons: endpoint mismatch, missing a/s physical clock, specimen mapping,
loading protocol, microstructure, and probability convergence certification.
It checks declared prerequisites, not external truth. No production result
has been certified against these data. A positive declared gate alone would
still require an actual prediction to calculate a lifetime error.

## 7. Reproducibility and present limits

The fresh output directories are protected against overwrite. Optional source
image extraction needs PyMuPDF; numerical replay uses the committed small
transcription and NumPy/SciPy. No optional source-image dependency enters the
desktop startup path.

```powershell
py -3 -m solver_v1.run_line_kinetic_benchmark --out results/new_line_calibration
py -3 -m solver_v1.run_line_timescale_validation --kinetics results/new_line_calibration/calibration.json --out results/new_line_dynamics
py -3 -m pytest solver_v1/test_line_kinetic_validation.py -q
```

Use `--image <native-page5.png>` to replay all isolated source markers and
verify its raw SHA. This run did replay the image; a data-only replay must
not claim that it checked the source image again. Results are saved under
`results/strength_fatigue_kinetics_v18/`; scripts, fields and scope are distinct
from the production energy registry and physical-time calibration JSON.
Scoped Git attributes preserve the raw hashed v18 JSON bytes; they do not
change old scientific data or global checkout settings.

Another1MiB ranged request for the public300K NVT trajectory timed out before
HTTP response headers at20s; no new MD frames were obtained. The previously
audited1024frames/25.575ps are unchanged. A longer record alone would not
remove the documented underdamped/thermostat/mode-mapping validity checks.

**What remains unresolved:** material interface force/Hessian mismatch(v17),
nonlinear normal-relaxed defect cores and finite-source evolution; matched
source statistics and accumulation to specimen yield; fatigue initiation versus
growth/fracture mapping; valid coarse-coordinate slow kinetics for a/s; A_c
and spatial dependence. New literature data make these comparisons concrete;
they do not justify fitting forbidden strength/time scale factors.

The existing PDE remains LJ/reduced-hybrid, model-time and unaltered. The new
line benchmark is a separate deterministic continuum limit of the retained
lattice elastic reference, not a replacement empirical plasticity law. Its
reversible result is a mechanistic negative control, not failure to plot a
large enough strain. UI changes in this work are accessibility only: vertical
form/summary scrolling and a fixed solve action row; no new physics default.

## 8. Final verification and UI accessibility

Executed with Python3.13 on Windows; no numerical test was skipped:

| Suite | Passed | Skipped | Reported elapsed [s] |
|---|---:|---:|---:|
| New kinetics + desktop/localization targeted | 22 | 0 | 8.57 |
| Entire solver_v1 | 537 | 0 | 1708.52 |
| Entire app, final code | 34 | 0 | 172.38 |
| Desktop startup smoke | exit0 | -- | 2.49 |

Test logs/XML are local `.cache/research_validation_v18/` artifacts; aggregate
results are in `results/strength_fatigue_kinetics_v18/validation_manifest.json`.
Wall times include shared-machine conditions and are not performance guarantees.
Six CSV files and the fitted-velocity JSON reproduce byte-for-byte in a fresh
directory. Original-image markers were actually replayed; hashes also match
the staged JSON bytes. Generated SVG trailing whitespace was mechanically
normalized without changing non-whitespace content. Working/staged diff checks
pass. No source, result or test pass certifies specimen fatigue or the a/s clock.

Pre/Solve settings and the summary now scroll vertically. Run/convergence
buttons and progress stay outside the scrolling settings at the bottom.
Private widget bindtags preserve Matplotlib wheel zoom and prevent combobox
wheel changes; Page Up/Down and Tab-to-reveal work without solving again.
Both940x620 and940x480 actual Tk windows were exercised in ko/en, preserving
inputs/results/plot limits. Native-window captures separately confirm the
buttons are visible; they are not UI mockups. No language dictionaries or
numerical field meanings were changed. The new line-calibration JSON is
rejected by the production time loader, as required.

## Sources

- Gorman, Wood & Vreeland Jr (1969), *Mobility of Dislocations in Aluminum*,
  J. Appl. Phys. 40,833-841, DOI [10.1063/1.1657472](https://doi.org/10.1063/1.1657472).
  [Caltech author record](https://authors.library.caltech.edu/records/j5q43-k8847).
  Figure5c is 23 deg C, not an assumed 293 or 300 K. The source pools edge and
  mixed characters and selects leading, relatively unimpeded dislocations.
- Deschanel, Ben Rhouma & Weiss (2017), *Acoustic emission multiplets as early
  warnings of fatigue failure in metallic materials*, Scientific Reports7,13680,
  DOI [10.1038/s41598-017-13226-1](https://doi.org/10.1038/s41598-017-13226-1).
  This is an independent specimen-fatigue validation reference, not a lifetime
  calibration loss or a measurement of local atomic absorption probability.
