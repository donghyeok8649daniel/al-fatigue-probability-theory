# Finite-frequency phase resolution: v27

## Question and unchanged scope

Can a measured dissipative quadrature of the same conjugate periodic plane
gap be distinguished from thermal noise, a force-sign mistake, nonlinearity
and timestep sensitivity? This is an extension of the completed v26 response
test, not a refit of LJ, Bessel sums, the EAM embedding, or production mobility.
Al99 is an atomistic **reference**, not a replacement production potential.
The periodic 144-atoms-per-plane coordinate is still not a local isolated
interface. Its atomistic timestamps do not unlock the production PDE clock.

## Phase convention and work

For the exact generalized work and forces derived in v26,

    G_ext=-F q, F(t)=F_hat cos(omega t),
    q(t)=q_bar+Re[chi F_hat exp(i omega t)],
    chi=(c_cos-i c_sin)/F_hat.

The force sign belongs in the denominator. For positive real storage
response, passive lag has Im chi<0 and arg chi<0. A positive noisy estimate
is not repaired by clipping or taking its magnitude. At steady linear
response the mean power and work per cycle are

    <F q_dot>=-omega F_hat^2 Im chi/2,
    W_cycle=-pi F_hat^2 Im chi.

Finite-record work also contains transient and thermal exchange. We keep
sampled F*q_dot and integration-by-parts work, and compare each with the
whole-system internal-energy change. An energy residual is not dissipation.
The integer-cycle lock-in work and the sine integral of q are algebraically
the same functional; their agreement is NOT an independent physical test.
On a finite record integration by parts also includes [Fq]_start^end, which
need not vanish for a thermal trajectory. We save this endpoint separately.
Recorded atomic velocity power and internal-energy increments supply the
additional numerical work-balance checks, without claiming independence of
all these observables or converting thermal endpoint energy into dissipation.

## A priori frequency design

The completed 5 ns NVE equilibrium record from v25 is the only input to the
frequency selection, before any v27 forced trajectory is produced. It uses
the same source potential, initial restart, geometry and approximate 299 K
temperature as v26. Define

    F_sigma=kBT/sqrt(C0),
    Im chi(omega)=-omega S(omega)/(2 kBT),
    Var(c_sin) approximately 2 S(omega)/T_record.

S is the **two-sided** covariance spectrum without positive-bin doubling.
The finite-band estimate uses DPSS NW=3, two tapers and f±0.01 cycles/ps.
For planned single-record signal/noise 3,

    T_3=18 S/(F_sigma^2 |Im chi|^2).

Choose the lowest of the predeclared set
{0.05,0.2,0.5,1,2} cycles/ps with T_3<=360 ps. This yields normal 2 and
direct110 1 cycles/ps. These frequencies are atomistic response probes, not
fatigue loading conditions. They may expose inertial/memory behavior and do
not validate a constant overdamped mobility or a zero-frequency limit.

For each axis run ±F_sigma and ±F_sigma/2 at dt=2.5 fs, and ±F_sigma at
dt=1.25 fs. Twelve 400 ps trajectories, saved every25 fs; exclude first40 ps.
No new thermostat, damping, parameter fitting or force amplification is used.
All begin from the same restart, so signed runs are NOT independent replicas.
Positive/negative work can change temperature in NVE; first/last40 ps mean
temperatures are reported, not implicitly treated as a perfect fixed-T bath.

## Independent empirical null and its limitations

Take every complete nonoverlapping180/360 ps window from the unforced record,
and every represented plane. Apply the same offset/cos/sin fit with unit
virtual force. These are displacement quadratures in Angstrom, not a
nonzero susceptibility of an unforced crystal. For imposed force amplitude F,

    E_null=max_windows,planes |c_sin,null| / |F|.

A measured loss exceeds this **observed** null envelope only if
-Im chi>E_null. It is not a confidence bound, independent-replica estimate,
or guarantee against all future fluctuations. No sqrt(number of planes),
sqrt(number of blocks), or sqrt(2) sign-pair noise reduction is assumed.
Half-amplitude controls have twice the susceptibility noise envelope.

For 360 ps and full force, the old0.05 cycles/ps null envelopes are
2.96713e-5 (normal) and1.66043e-4 (slip) Angstrom^2/eV. The v26 expected losses
are much smaller. This independently supports retaining its unresolved
classification rather than interpreting positive estimates as negative drag.

The full/half driven-record sensitivity, ±force disagreement, amplitude
change, timestep change, and joint second/third harmonics are also retained.
Different timesteps generate different thermal histories: a timestep
difference is not automatically a pure deterministic discretization error.
Agreement within thermal envelopes is weaker than a quantified zero-error
extrapolation. No automatic production calibration flag is set.

## FDT comparison

Compare actual forced chi with both the independent covariance transform

    chi=[C0-i omega integral_0^infinity C(t)exp(-i omega t)dt]/kBT

and Im chi=-omega S/(2kBT). Preserve finite-time endpoint correction and
cutoffs2/5/10/20/40 ps. Preserve NW2/3, full/half records and spectral
half-bands0.005/0.01/0.02 cycles/ps. These are estimator-sensitivity ranges,
not statistical confidence intervals. Canonical FDT against finite NVE and
temperature drift remain explicit limitations, as do finite periodic size,
source-potential material accuracy, and coordinate/PMF mapping to production.

An additional retrospective quadrature audit keeps trapezoid/Simpson at
25/50 fs covariance sampling on the same40 ps lag interval. Simpson is NOT
automatically selected as truth. For normal response at2 cycles/ps,
coarsening to50 fs makes Simpson severely oscillation/alias sensitive and
even gives a positive imaginary estimate. This is a quadrature failure,
not negative friction. At25 fs the normal storage estimates are0.00266119
(trapezoid) and0.00263247 (Simpson) Angstrom^2/eV, a further roughly1%
sensitivity. In a smooth, decaying inertial covariance, Euler--Maclaurin
gives a leading trapezoid storage bias C0(omega*dt)^2/(12 kBT), apart from
the finite endpoint term and higher derivatives. It is not safe to assume
these higher terms vanish for fast phonons. Both raw rules and the failure
on coarsening are retained; finer-sampling validation of the time-domain
transform is not claimed. The reported dissipative FDT ranges instead use
the directly measured two-sided spectrum, with its own bandwidth/record
sensitivity. A time-domain cutoff or quadrature is never selected to force
agreement with the driven result.

## Reproduction

### Why a resolved finite-frequency phase is not an overdamped clock

For an equilibrium, reversible overdamped generator with reflecting
boundaries and the SAME coordinate and conjugate forcing, its spectral
representation has nonnegative weights:

    chi(omega)=beta sum_k w_k lambda_k/(lambda_k+i omega),
    w_k>=0, lambda_k>0, chi(0)=beta sum_k w_k.

Consequently Re chi(omega)<=chi(0) and |chi(omega)|<=chi(0). This is not
an assertion about a nonlinear absorbing driven ensemble far from equilibrium.
If the atomistic response in the selected band instead has an enhanced
storage response, inertia/memory cannot simply be replaced there by fitting
one positive overdamped M. The static covariance response and storage ratio
are therefore saved. Such a finite-band failure does NOT disprove a lower
frequency Markov reduction; it explicitly prevents importing a measured
atomistic high-frequency phase as the production kinetic calibration.

### Additional actual-engine benchmark

The v26 synthetic oscillator runner now accepts explicit optional frequency,
stiffness, drag and duration; its default output reproduces v26 exactly.
Tests at1/2 cycles per synthetic time use m=1,K=1000,drag=1, not Al parameters.
At dt=.0025 the phase errors are6.12e-7 and6.20e-6 rad, respectively,
decreasing at dt=.00125 to2.38e-7 and2.43e-6. These are much smaller than the
planned MD loss angle. Both complete three-dt engine runs and exact default
reproduction passed. They check the clock/forcing/fitting path, not thermal
statistics or quantitative Al friction.

### Commands

`python -m solver_v1.run_phase_resolution_study prepare --source ... --restart ... --out ...`
creates the protocol before new dynamics. `run --study ... --restart ...
--potential ...` executes the optional verified LAMMPS reference runner.
`analyze --study ... --source ... --out ...` rejects incomplete/mismatched
cases and saves null controls, individual/signed-paired response, FDT,
work balance and comparison tables. Fresh outputs are required; no existing
user or raw result is overwritten.

## Completed results and scientific disposition

All twelve400 ps runs completed, and all logged engine clocks agree with
step*dt within5.68e-14 ps. The analysis emitted36 individual-window fits,
18 signed-pair fits,1440 primary null windows and12 work audits. The plot
was visually inspected. Values below are signed-pair full360 ps fits at the
full force amplitude; negative phase means lag.

| Axis | Frequency [cycles/ps, reference MD] | dt [fs] | Re chi [A²/eV] | Im chi [A²/eV] | Phase [deg] | Loss / observed null envelope |
|---|---:|---:|---:|---:|---:|---:|
|normal|2|2.5|0.00264264749|-2.93428129e-5|-0.63616|1.2826|
|normal|2|1.25|0.00264287500|-4.08000543e-5|-0.88445|1.7834|
|direct110|1|2.5|0.0138111383|-3.29504081e-4|-1.36669|1.7484|
|direct110|1|1.25|0.0137958737|-3.52848104e-4|-1.46510|1.8723|

The observed full-record null envelopes are2.28782e-5 normal and1.88461e-4
slip A²/eV. Both full-force signed-pair losses exceed these envelopes at
both dt values. Every individual full-record phase is negative. Spectral
FDT sensitivity predicts Im chi in[-4.26782e-5,-3.24731e-5] normal and
[-3.46771e-4,-2.84710e-4] slip A²/eV. The driven results are compatible
with these ranges once the measured thermal envelope is retained; this is
not a tight relative-error certificate.

**Disposition: evidence for resolved finite-frequency lag against the
observed null controls, with limited precision; NOT a completed quantitative
phase-convergence or production kinetic calibration.** The actual dt change
in Im chi is about39.0% normal and7.1% slip. It cannot be separated from
different thermal histories here. Many180 ps blocks do not independently
clear their larger null envelopes. Half-force slip also remains below its
own full-record null envelope, despite a negative fitted phase of-1.184deg.
Thus the study does not certify all amplitudes/record lengths, and the
old0.05 cycles/ps v26 response remains unresolved. No favorable subset is
used to turn on a physical clock.

The normal/slip storage responses are12.08%/15.86% above their static
covariance susceptibilities. This exceeds the roughly1% normal covariance
quadrature sensitivity on the25 fs grid and matches the qualitative
inertial enhancement in equilibrium FDT. This band is NOT an appropriate
place to infer a single reversible overdamped M while preserving the same
static response. Lower-frequency and coordinate/PMF/material matching are
still required. Do not identify this reversible dynamic lag as plasticity.

### Harmonic and amplitude controls

Doubling force changes normalized storage by about0.21% normal and0.093%
slip. Raw second-harmonic amplitudes reach2.41% and9.98% of the fundamental,
respectively, but a nonzero harmonic is not automatically nonlinearity in
this thermal signal. `harmonics --study ... --source ... --out ...` performs
a separate fresh-output audit at2f/3f on the same unforced record. It
retains all individual harmonics and the force-even/force-odd combinations

    q_even,k=(q_hat,k(+F)+q_hat,k(-F))/2,
    q_odd,k=(q_hat,k(+F)-q_hat,k(-F))/2.

The expected quadratic even and cubic odd components ALL remain below
their observed null amplitude envelopes (largest ratio0.887). These
measurements do not establish a coherent nonlinear harmonic, nor prove
exact linearity below their noise floor. Individual harmonics that exceed
the observed envelope remain saved; only selecting the quieter pairs would
be misleading. See `harmonic_controls` for24 individual rows,12 pairs and
312 null windows.

### Work, heating and preservation

Whole-record work-energy residual maxima are0.009072eV (parts) and0.008294eV
(sampled velocity power), for the entire864-atom crystal. In the full-force
normal controls the work is0.747–1.140eV and the parts residual magnitude is
0.000142–0.002895eV. The source controls therefore contain substantial
externally supplied energy relative to the numerical work residual, unlike
v26's tiny loss. Endpoint terms and periodic fitted work remain separate.
Maximum COM velocity change is1.554e-14 A/ps. Mean run temperatures span
299.331–302.262K; first/last40 ps mean increases span0.469–4.542K. We did
not insert thermostat damping to suppress that heating. It limits fixed-T
FDT precision and the claim of a stationary driven ensemble.

Production LJ/Poisson/Bessel, static calibration, PDE, Ma/Ms, kinetic JSON,
UI and specimen area are unchanged. Physical seconds/Hz remain disabled.
The source EAM finite-frequency result is not a quantitative experimental-Al
friction validation, not a local-interface mobility, and not fatigue life.

### Executed validation

- 25 targeted tests passed (1.03s before final documentation).
- 683 full solver tests passed1790.49s with MD running concurrently.
- 34 app tests passed70.04s; desktop smoke passed; no skips.
- Actual-engine synthetic benchmarks: both frequencies, three dt each,
  decreasing errors; old default output reproduced exactly.
- Protocol re-generation was byte-identical. Incomplete study was explicitly
  refused without emitting a completed result. All12 actual MD cases and
  both final analyses completed. Per-run MD wall times sum7311.09s (two
  concurrent runs, not7311s of independent trajectory time).
- Final diff/commit/remote checks are recorded in the handoff/Git history.
