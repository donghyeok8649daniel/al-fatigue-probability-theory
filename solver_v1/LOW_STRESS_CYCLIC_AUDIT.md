# Low-stress cyclic audit: asking the fatigue question directly

This study follows the user's correction that multi-GPa ideal interface tests
do not establish ordinary metal fatigue. No parameter is refitted. The current
research surface still fails material-calibration gates. The production
TwoRowLJ/reduced hybrids, PDE, UI, and parameter files remain unchanged.

## 1. Experimental context is not a fitting target

Zhai, Briggs and Martin, *Acoustic microscopy of room-temperature fatigue
damage in aluminium single crystals*, Philosophical Magazine A 77 (1998),
957–980, DOI [10.1080/01418619808221222](https://doi.org/10.1080/01418619808221222),
report a resolved shear-stress amplitude of **4 MPa**, R=-1, room temperature,
air. Their [author institution record](https://scholars.uky.edu/en/publications/acoustic-microscopy-of-room-temperature-fatigue-damage-in-alumini/)
describes slip-band/surface observations and a change after 5 million cycles.
The accessible abstract does NOT determine all orientation, purity, frequency,
or specimen/defect parameters. We do not invent those missing inputs, reproduce
their lifetime, or identify an example 45-degree axis with their actual crystal.

The numerical protocols are:

- zero applied load, over every matched duration;
- declared resolved shear amplitude 4 MPa, fully reversed;
- nominal amplitudes 10,20,30,50 MPa at the explicit example e=(n+m)/sqrt(2),
  hence Tn=tau=sigma/2, R=-1;
- a separately labeled nonproportional probe, Tn and tau amplitudes 15 MPa,
  90-degree phase difference. This is not a fixed uniaxial loading direction.

Additional controls use normal-only 15 MPa, shear-only 15 MPa, and the same
30 MPa example with reversed shear direction. The 50 MPa example is an upper
stress probe, not a declaration of safe high-purity Al loading.

These stress levels are diagnostic choices, not universal Al fatigue limits.
The temperature convention 293.15 K is an explicit room-temperature test
setting, not an exact temperature extracted from the abstract. The 0 K energy
candidate is held fixed: thermal expansion/finite-T free-energy calibration
has not been supplied. kBT=8.617333262145e-5*T eV is not adjusted to cause slip.

## 2. Exact stress/work distinction

The interface diagnostic uses its existing conjugate work:

    G* = W_int[eV/cell] - f_a*(a*-h*) - f_s*s*,
    f_i* = T_i[MPa] 10^6 A_atomic_cell[m^2] L0[m] / E0[J],
    L0=4.05/sqrt(2) angstrom, E0=1 eV.

This is not a replacement of the production reduced calibration
f*=kappa_axial*sigma/E. It contains no A_c or mesh area. In the prior example,
Tn=tau=2 GPa corresponds to sigma=4 GPa, not 2 GPa nominal fatigue stress.
High-load results remain only idealized stability probes.

## 3. What is actually evolved, and what is deliberately NOT claimed

The new module is a **noncanonical reflecting-domain falsification experiment**.
It asks what follows IF one assigns the current energy per atomic interface
cell a collective-coordinate Gibbs distribution. That normalization is a
hypothesis, not a calibrated physical activation region.

    partial_t P = div[M(P grad G + kBT grad P)], M*=diag(1,.05).

It uses the unchanged production SG *energy-array numerical generator*, not
the production model constructor, opening saddle logic, or energy selector.
All normal/outer-registry boundaries in THIS diagnostic reflect. There is no
opening sink. Therefore `opening_probability=null`, not zero, in its output.
There is no claim of crack absence, first-passage validation, material
acceptance, or a newly production-valid PDE energy model. The purpose is to
reject unsupported fatigue interpretations before any such promotion.

The grid samples the analytic infinite LJ/environment sums directly. Proven
registry translations are reused; no finite-neighbor potential or spline is
substituted. Individual grid/timestep/domain refinements are separate from
the old energy-series tests. The direct110 scalar restriction omits transverse
registry motion and Shockley partial wells; conclusions do not cover all slip
paths. Static MPa cycles separately inspect both paths.

For direct110 only, reflection x->-x sends k*tau to k*tau-k*a1. It preserves
the triangular lattice and every radial pair sum and rotational/reflection
invariant moment norm. Thus W(a,s)=W(a,-s) on this path. Fine-grid preparation
can reuse this exact reflection after independent signed-point checks; it
does not average asymmetric values or assume Shockley112 has the same mirror.

Time and period are MODEL units only. No physical mobility or Hz is inferred.
Periods varied at fixed mobility diagnose finite-rate behavior; a period
change is not a material calibration. A finite number of model cycles cannot
establish experimental long-life behavior or absence of fatigue at millions
of cycles.

## 4. Initial ensemble and matched zero-load control

Primary initial density is Gibbs restricted to the central registry well.
The normal distribution is not artificially collapsed to a*. Releasing that
restriction at t=0 can cause thermal spreading even without cyclic loading.
Every driven state must therefore be compared against a ZERO-load run with
identical initial density, temperature, duration, grid and timestep.

A global reflecting-box Gibbs control is also tested. It has stationary well
populations and can have nonzero SG gross traffic. Its multiple initially
occupied wells do not constitute prior fatigue damage.

At each phase/cycle the stored density gives:

    n=floor(s/b+1/2), xi=s-bn,
    P_n=sum_{well n} cell_probability,
    <s>/h = <xi>/h + b<n>/h.

This is a LOCAL registry-shear decomposition, not a calibrated macroscopic
axial plastic strain. No chi multiplier is fitted. The absence of absorption
in this reflecting experiment gives

    Delta P_n = dt (J_left - J_right),
    Delta <n> = dt sum_interfaces J.

Backward Euler fluxes are evaluated with the same end-step density as the
generator. Raw mass, well/moment balance and minimum probability remain
visible; no negative-probability repair or renormalization is performed.

## 5. Gross SG traffic is not a continuum count of physical hops

An important additional audit: the SG net/forward/backward decomposition is
algebraically correct, but its gross rate has a grid-dependent interpretation.
With smooth P, spacing ds, D=M_s kBT and face probability density p_s,

    J_net = D/ds [B(dG/kBT) p_L - B(-dG/kBT) p_R] -> finite flux,
    J_gross = D/ds [B(dG/kBT) p_L + B(-dG/kBT) p_R]
            ~ 2 D p_s(face)/ds.

Thus gross nearest-neighbor traffic diverges as 1/ds: Brownian recrossings of
a sharp partition are not a finite committed interwell-event count. This does
NOT change the original SG formula or imply a bug in its net flux. Raw gross
traffic is retained with its grid spacing; net balance and well occupations
are the convergent directional diagnostics. Physical committed transitions
would require defined well cores/committors and an independently derived
reactive-flux observable. `abs(net)` is not that observable either.

## 6. Why per-cell energy is not automatically an activation barrier

W_int is a finite energy PER ATOMIC AREA of a spatially uniform displacement
of two infinite half-crystals. For N coherently translating atomic interface
cells, the total uniform energy/work is N times the per-cell quantity. Its
Boltzmann exponent is N G_cell/kBT, not automatically G_cell/kBT. A true local
slip nucleus instead has nonuniform displacement and an elastic/core boundary
cost, which the two rigid half-crystals do not represent.

Also, finite work of separation implies W_int tends to a finite constant as
a goes to infinity. The unrestricted zero-load Gibbs integral over the entire
normal half-line is then not normalizable. A reflecting-box Gibbs control is
not the full opening ensemble: a metastable conditional/QSD and an explicitly
validated escape treatment would be needed for that interpretation.

This observation does not set N, introduce A_c into the energy, or permit an
arbitrary multiplier chosen to suppress thermal hopping. A statistical
specimen correlation area is a different object. A one-cell diagnostic may
be useful precisely because it exposes thermal spreading that must NOT be
called irreversible macroscopic plasticity.

The analytic LJ/Bessel framework can remain the microscopic starting point
for a later spatial displacement theory. For example, eliminating harmonic
bulk variables z from a derived Hessian gives a Schur-complement kernel

    K_eff = K_uu - K_uz K_zz^+ K_zu,

with rigid translations/gauge treated explicitly. This is a possible route
to a nonlocal interface elastic cost, not an implemented dislocation model.
One must derive it from the same lattice and avoid counting the already
included interface harmonic interactions twice. No empirical yield law,
fatigue-life fit, or arbitrary stress concentration is inserted here.

## 7. Reproduction and acceptance

    python -m solver_v1.run_low_stress_cyclic_diagnostic --phase static
    python -m solver_v1.run_low_stress_cyclic_diagnostic --phase grid
    python -m solver_v1.run_low_stress_cyclic_diagnostic --phase cycles

Other grids, model periods, phase counts and domain sizes are explicit CLI
arguments. Raw execution records are under
`results/fcc111_active_interface/low_stress_v4/`. A successful numerical
experiment is not successful material fatigue validation. Conclusions must
distinguish driven excess over the zero-load control, reversible finite-rate
lag, thermal spreading, numerical resolution and untested opening loss.

## 8. Executed results, not a proposed future run

69 deterministic time-evolution protocols were executed. The grid energies
were independently evaluated at 31x75, 61x75, 61x147 and 91x219 cells, with
25/49/73 aligned cells per registry well. A separate 85x147 grid extends the
normal upper boundary from 2h to 2.5h. Registry domains contain 3,5,7 wells;
phase counts 64,128,256,512 independently test backward Euler time error.
Explicit/implicit comparison is additionally performed for a short 2-cycle
protocol, rather than claiming that every long implicit run has an explicit
duplicate. `execution_manifest.json` records commands, actual grid/parameter
hashes and run timings; total individual cyclic-run elapsed times sum to
1345.68 s (overlapping jobs, NOT total wall-clock duration).

The physical-stress inputs are MPa. All periods below are MODEL time, at
unchanged M*=(1,.05). The pristine eigenmode times are
tau_fast=.0436642952 and tau_slow=4.04431664. At periods .4,4,40,
omega*tau_slow=63.53,6.353,.6353, respectively. These pristine small-signal
diagnostics are not an exact finite-temperature nonlinear relaxation spectrum.

### 8.1 Static low-MPa cycles retrace the same intact well

Both direct110 and Shockley112 paths were followed through two complete
quasistatic cycles at 10,20,30,50 MPa nominal amplitude. Every sampled point
remained a stable local minimum. Maximum force residual was 1.8741e-14
eV/reduced coordinate; repeated-cycle coordinate differences were at most
2.3054e-15. At 30 MPa the direct110 normal strain spans
[-1.0158699e-4,1.0208683e-4], with s/b in
[-3.8476662e-4,3.8538788e-4]. No static well change occurred. This is not a
dynamic persistence test and does not establish an experimental elastic limit.

### 8.2 Actual cyclic transport is dominated by the zero-load baseline

The following is one consistently resolved comparison grid, 61x147,
128 steps/cycle, period40, 16 driven cycles followed by 8 zero-stress hold
intervals (drive time640, hold time320). P_out is registry population outside
n=0; it is **not crack probability**. The last column is a local mean registry
index, **not macroscopic plastic strain**.

| Protocol | P_out at drive end | Driven-minus-zero P_out | Mean n after hold |
|---|---:|---:|---:|
| Zero load | .02202894954 | 0 | -1.83e-14 |
| Resolved shear amplitude4 MPa | .02202907327 | 1.23725e-7 | -1.09314e-9 |
| Axial amplitude10 MPa, declared45-degree axis | .02202944369 | 4.94146e-7 | 1.00709e-7 |
| Axial amplitude20 MPa, same axis | .02203035104 | 1.40149e-6 | 4.05558e-7 |
| Axial amplitude30 MPa, same axis | .02203167161 | 2.72207e-6 | 9.14542e-7 |
| Axial amplitude50 MPa, same axis | .02203555259 | 6.60304e-6 | 2.54494e-6 |
| Tn/tau15 MPa, 90-degree phase | .02203161670 | 2.66716e-6 | -4.04412e-8 |

For 30 MPa, over 99.98% of the raw outside-well population is also present
without loading. Calling that entire population fatigue accumulation would
be wrong. Nevertheless, the small paired load-induced excess and directional
transport do not vanish on refinement; they should not be dismissed merely
because they are small. Normal-only15 MPa gives mean n=-8.95e-15; shear-only
15 MPa gives -4.17e-9 after hold. Reversing shear in the coupled30 MPa probe
changes 9.14542010e-7 to -9.14542029e-7 while preserving P_out. This supports
small directional rectification by coupled loading in this mathematical
system, not a new physical calibration or empirical plasticity law.

An initially global Gibbs ensemble remains stationary under zero load, to
numerical precision; it initially has 2/3 of its mass outside n=0 in the
three-well box. Thus outside-well population alone is not a history-independent
damage measure. Under the same30 MPa drive this alternate initial ensemble
leaves mean n=5.54726e-7 after hold, illustrating initial-state dependence.

### 8.3 Finite-rate loops and actual timestep convergence

At 30 MPa, the period4 last-cycle local registry-shear amplitude is
7.24357e-5 with phase -83.0624 degrees; period40 gives 4.60856e-4,
-42.2971 degrees on the same61x147,128-phase grid. This is the response of
the equations at fixed mathematical mobility, not a comparison of physical Hz.
The longer-period run also lasts ten times longer at fixed cycle count;
each duration has its own zero-load control.

For period40 and30 MPa, the independent time study gives:

| Steps/cycle | dt_model | Last-cycle <s>/h amplitude | Phase [degrees] | Mean n after hold |
|---:|---:|---:|---:|---:|
|64|.625|4.5587061e-4|-43.09045|9.2079225e-7|
|128|.3125|4.6085558e-4|-42.29707|9.1454201e-7|
|256|.15625|4.6343130e-4|-41.89760|9.1140640e-7|
|512|.078125|4.6474044e-4|-41.69714|9.0983617e-7|

Amplitude/phase changes decrease with refinement as expected for backward
Euler; the128-phase phase angle is NOT presented as fully time-converged.
The short explicit comparison's final density L1 discrepancy decreases from
2.65232e-6 to1.43903e-6 for64->128 phases (30 MPa, period.05,2 cycles).
This compares the same SG generator and analytic energy, not two different
material models. It does not certify the long-time phase error by itself.

At period40,30 MPa the first/16th-cycle external works are
1.32109e-6/1.55296e-6 eV/cell. Finite-rate dissipation and a nearly repeated
local response exist. The whole distribution is NOT fully periodic: thermal
spreading between wells continues, including during the zero-load hold.
An elliptical local loop is not by itself fatigue damage.

### 8.4 Resolution envelopes, persistence, and domains

For the refined91x219,256-phase30 MPa run, the driven-minus-zero P_out is
2.6984433e-6. Its conservative observed spatial/time/domain difference
envelope is1.6216931e-8, a signal/envelope ratio166.4. The held mean n is
9.1177955e-7 with observed envelope1.8402984e-8 (ratio49.5), dominated by
the3->5-well change. These are **empirical numerical difference envelopes**,
not confidence intervals, formal error bounds or physical certificates.
Raw P_out spatial differences are much larger than paired-excess differences;
the two floors must not be interchanged. `observed_resolution_envelopes.csv`
also states that50 MPa has no separate domain-refinement study; its smaller
partial envelope must not be mistaken for broader certification.

For30 MPa,61x147 equivalent spacing,128 phases, the held n values at3,5,7
wells are9.1454201e-7,9.3294502e-7,9.3309795e-7. The5->7 difference is
1.53e-10. In the7-well case, extending zero-stress hold from320 to1280 model
time changes n from9.3309795e-7 to9.3309274e-7, with final intrawell shear
1.40e-13. The local index contribution b<n>/h is then1.14280e-6. This is
short-to-intermediate-time persistent directional registry memory in the
declared model hypothesis; it is not validated residual Al plasticity. The
outermost-well mass at the long-hold endpoint is1.05e-5, explicitly monitored.
An unbounded periodic landscape with no hardening does not guarantee eternal
irreversibility; a finite reflecting box ultimately has a symmetric Gibbs state.

For period4, the32-model-time hold is insufficient: intrawell shear is still
-2.49e-7. Extending hold to320 reduces held n from9.73e-9 to1.68e-9.
That short-period endpoint must not be called a converged residual plasticity.

The expanded normal box changes30 MPa paired P_out by1.45e-12 and held n
by1.89e-12. Across all69 runs, maximum mass residual is9.6541e-12,
maximum per-step well-balance residual4.2262e-15, registry-moment balance
9.4560e-16 and decomposition residual2.0007e-16. No negative-mass repair
or probability renormalization was used. Cumulative net transfer is distinct
from gross SG traffic, whose mesh-dependent behavior was derived in section5.

### 8.5 Why this does not yet establish ordinary Al fatigue

The one-cell fixed-domain potential-of-mean-force sampled barrier converges
toward .18568 eV, about7.35 kBT at293.15 K. The4 MPa resolved shear work over
one full direct110 repeat is only .000507809 eV, .02010 kBT. Thermal spreading
and a weak load bias are therefore internally understandable; neither a
multi-GPa requirement nor an arbitrary mobility/energy rescaling is justified.

The zero-load mean normal strain is about .018402 relative to the0 K spacing
(at s=0 alone the conditional normal shift is about .010336). This is a
consequence of this one-cell anharmonic Gibbs hypothesis, **not a prediction
of bulk Al thermal expansion or residual opening**. Matched thermal baselines
are required before interpreting an unloaded nonzero normal strain.

The most defensible outcome is: small converging load-induced directional
registry memory and finite-rate loops exist in the tested mathematical
system; ordinary metal fatigue, a physical residual strain, and crack
initiation have NOT been validated. Missing ingredients are a justified
spatial/collective activation energy, local defect/slip-nucleus mechanics,
better independent material-energy agreement, and collective-coordinate
kinetic data. The reflecting test cannot answer crack/slip event ordering.
It contains no crack first passage: a missing crack result is not a zero one.
No specimen A_c or guessed stress-concentration multiplier is used to repair
these gaps. The next physical model must derive such mechanics from the same
lattice framework rather than fitting this diagnostic to a desired lifetime.

## 9. Data size, backups, and reproducibility

Phase histories are losslessly gzip-compressed CSV. All phase steps enter the
mass/flux audit and every cycle's moments; serialized histories keep full
first/last drive/hold cycles and quarter-cycle points elsewhere where stated.
Selected full-precision density snapshots include the initial, first,
penultimate/last driven and hold endpoint states; explicit cycle indices are
stored, never inferred from array position. Complete legacy snapshot archives
were preserved in ignored `.cache/low-stress-full-snapshots/` before sampling;
`output_compaction.json` and audit checksums document that reversible local
backup. Regeneration commands are in `execution_manifest.json`. Sampling and
compression only reduce stored output size, never grid/time accuracy.

## 10. Verification checkpoint

Actually executed with Python3.13.5, NumPy2.5.0, SciPy1.18.1:

- New targeted diagnostics: 11 passed,1.26 s.
- Full `solver_v1`: 223 passed,1192.85 s.
- `app`: 27 passed,2 skipped,155.84 s.
- Desktop startup smoke: passed,1.1416 s; historical TwoRowLJ
  a0=.7713438268704838 and kappa=86.29296488740997 remained unchanged.
- Working/staged diff whitespace checks: passed.
- All57 full-snapshot local backup SHA256 values: independently verified.

The first new temporary-file test attempt encountered a Windows OS-Temp
permission error (10 passed,1 setup error). It was rerun successfully using
a fresh explicitly named worktree-local `.cache` pytest temporary directory;
no existing user temporary directory was deleted or permission-modified.
The full suite used another fresh worktree-local directory. The diagnostic
results do not certify material calibration, first passage or physical time.
