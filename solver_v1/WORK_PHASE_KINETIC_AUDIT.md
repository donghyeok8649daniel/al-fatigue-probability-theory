# Native-work cross-check of the unresolved kinetic phase

## Scope

This is a new analysis of the **existing** v31 ten-record, 20 ns Al99 reference
MD campaign, not new MD and not production mobility calibration. Every input
trajectory is checked against the previously published SHA256 and restart
binding. Production LJ/Bessel energy, Smoluchowski drift, physical-time JSON,
temperature, mobilities and specimen aggregation are unchanged.

The question is narrower than "can we enable Hz": does saved-frame phase
estimation discard a measurable dissipation signal that native integration-step
work can recover? This is tested before spending more time extending trajectories.

## Exact work and response identity

The measured gap q is the mean upper-plane displacement minus mean lower-plane
displacement, projected on normal or direct110. Its conjugate force is F(t),
not force per atom. Opposing atomic forces are +/-F/Np. Use

    F(t) = F0 cos(omega t),  omega=2 pi f,
    q_response(t) = F0 [chi_R cos(omega t) - chi_I sin(omega t)].

For an arbitrary finite thermal record, without assuming periodic q,

    W = integral F dq
      = [F q]_t1^t2 + omega F0 integral q sin(omega t) dt.

Thus, over T containing an integer number of driving cycles, define

    L_q = 2/(F0 T) integral q sin(omega t) dt,
    L_W = 2 (W-[Fq]_ends)/(omega F0^2 T).

Both estimate -Im chi in the stated convention. Their equality is an
integration-by-parts identity even before a linear stochastic interpretation
is accepted. The endpoint term must NOT be called dissipation. A thermal
drift or initial/final displacement difference can contribute to it. Negative
finite-record loss is retained, never repaired using abs or clipping.

Units here: q in angstrom, F in eV/angstrom, W in eV, t in ps. L has units
angstrom^2/eV. These reference-MD ps are known; the production PDE clock is not.
The absolute time origin is preserved after slicing. Multiplying sliced time
by a different phase origin would be a phase error; no such error was found.

Native work uses the previously verified internal-step trapezoid (v30),
not the 25 fs saved-frame power. Coordinate quadrature is independently
evaluated at saved intervals 25, 50 and 100 fs. Whole/half/quarter records
after the first 100 ps are retained: 1900/950/475 ps, all complete cycles.
Coarsening stored frames is a sampling check, NOT an MD-timestep refinement.
Native and coordinate estimators are not independent physical samples.

## Actual results

`results/work_phase_audit` contains 168 signed driven block/stride comparisons,
168 null plane/block observations, 28 signed-pair block rows and planning data.
Null planes and adjacent blocks are NOT declared statistically independent.
Observed null envelopes are sensitivity scales, not confidence intervals.

| Initialization | Axis | Native-work loss [angstrom^2/eV] | Native loss / imaginary null |
|---|---|---:|---:|
| 35461 | normal | 7.3484807831e-6 | .777069 |
| 49277 | normal | -9.2221370640e-7 | -.107571 |
| 35461 | registry | 5.2067094196e-5 | 1.882055 |
| 49277 | registry | 5.7844978959e-5 | 1.081592 |

Across original-frame driven blocks, maximum absolute native/coordinate loss
difference is 4.04931739e-8 angstrom^2/eV. After signed pairing, the maximum
discrepancy / corresponding null is .000623525 (0.06235%). Native work therefore
confirms the old small/noisy phase; it does not expose a missing positive drag.
In particular the second normal response remains opposite-signed. The registry
full-record ratios are positive, but half/quarter uncertainty must not be
discarded in favor of the full-record ratio. Earlier complex uncertainty disks
are unchanged; their zero-drag inclusion is not replaced by these scalar ratios.

**Conclusion:** native-work storage/phase arithmetic is not the dominant
resolution problem in these records. This does not prove every MD systematic
error is absent. Earlier independent timestep controls remain relevant.

There IS a storage-resolution hazard on deliberate coarsening. Maximum loss
discrepancy over driven blocks increases from 4.0493e-8 at 25 fs to 9.2718e-7
at 50 fs and **1.52284e-4 at 100 fs**. The latter can exceed the entire desired
loss signal. Resolving the forcing frequency alone is insufficient: unresolved
thermal/high-frequency motion can contaminate its low-frequency quadrature.
Native work retains every integration step, so it exposes this error without
rerunning MD. Consequently 100 fs coordinate-only phase must not be substituted
for the validated 25 fs record merely to save storage or runtime. This observed
coarsening failure does not explain the unresolved sign in the ORIGINAL 25 fs
measurement, whose independent native-work result agrees closely.

## Cost of merely extending the same experiment

As a sensitivity calculation, suppose stationary null amplitude scales as
T^(-1/2), and request signal/null ratio r=3. Use the smallest independent
unforced FDT predicted loss across the previously declared full-record bands,
not an absolute-valued negative driven estimate. Then

    T_required = T_observed [r N_observed / L_predicted]^2.

Using the maximum full-record imaginary null across initializations gives:

| Axis | Predicted loss | Null | Conditional analyzed duration | Single-record wall estimate |
|---|---:|---:|---:|---:|
| normal | 4.36104673e-6 | 9.45666193e-6 | 80.406 ns | 39.82 h |
| registry | 3.96178888e-5 | 5.34813092e-5 | 31.161 ns | 15.43 h |

These are **not** required times proved by statistics, achieved precision,
confidence bounds or a scheduled campaign. Wall estimates scale the median
observed seconds/ps of these ten records; hardware contention and checkpoint
costs can change them. They apply to one record, not an entire signed/null
campaign. The older v28 spectral SNR estimates use another noise convention;
the different numbers must not be silently substituted for one another.
Long NVE forcing also requires checking heating; stationary scaling may fail.
No very large campaign was silently launched, and drive amplitude was not
increased to hide the cost or force a desired phase.

## Why this cannot alone supply production Hz

1. The measured normal dissipation still has unresolved sign; inverse drag
   does not have a finite certified upper mobility bound.
2. The reference finite-frequency plane-gap susceptibility is not automatically
   a constant diagonal local (a,s) mobility. Full-coordinate elimination can
   produce memory; scalar and full matrix inverses differ (v31).
3. Per-plane energy division requires BOTH mobility and thermal normalization
   changes. Importing a plane-average mobility into an atomic-cell PDE alone
   changes its generator, not merely its units.
4. An accepted matched Al material surface remains a separate gate. Better
   reference-MD phase precision cannot repair the failed static curvature fit.

The exact existing conversion M*=t0 E0 M_phys/L0^2 and f_Hz=f_model/t0 is
unchanged. A user-selectable guessed t0 would conceal these missing validations.
Consequently production seconds/Hz remain unavailable, explicitly, rather
than being populated from a finite-frequency diagnostic point estimate.

## Reproduction / tests

    python -m pytest solver_v1/test_work_phase_audit.py -q
    python -m solver_v1.run_work_phase_audit --study <v31 raw campaign> \
      --previous results/weak_replica_v31 --out <fresh report directory>

Analytic tests include signed forces, positive/negative/zero loss, arbitrary
phase origin, finite endpoint drift, quadrature agreement, invalid inputs,
and a hash-bound end-to-end synthetic report that never certifies a clock.
The runner refuses overwrite of an existing result directory. Full regression
completion is recorded in CURRENT_WORK_HANDOFF.md, not presumed here.
