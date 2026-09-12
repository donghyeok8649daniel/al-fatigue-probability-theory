# v29 — actual additional low-frequency conjugate-response measurement

Status: **12 new MD records completed (10 ns aggregate)**; analysis and
controls completed. This is a completed measurement, not a certified
production mobility calibration.

The user explicitly requested execution after v28 provided only reanalysis.
We preserve production LJ/Bessel, energy/temperature/mobilities, the kinetic
JSON and the disabled physical-clock gate. Al99 is a reference potential only.

## Predeclared protocol

- Same N6/864-atom FCC crystal,144 atoms per plane and same equilibrated restart.
- Same q=(mean upper plane displacement - mean lower plane displacement).e.
- Exact conjugate atomic forces +/-F e/144, evaluated at every internal step.
- Frequency0.2 cycles/ps for normal and direct110, not production fatigue Hz.
- Full probe4F_sigma, F_sigma=kBT/sqrt(C0), each sign1000 ps, dt2.5 fs.
- Amplitude control2F_sigma, each sign500 ps, dt2.5 fs.
- Full probe repeated at dt1.25 fs, each sign1000 ps.
- All records saved every25 fs; first100 ps excluded from phase estimation.
- Twelve records total10 ns; two subprocesses at once. No added thermostat,
  fitted friction, change to potential, or mobility tuning.

Half-force controls have shorter records deliberately for storage and harmonic
tests; they are NOT equal-precision loss measurements. Neither common restart
nor sign pairs/equivalent planes count as independent replicas.

Full source/restart/potential hashes and all cases are saved in protocol.json
before simulation. Existing complete cases may be resumed only after metadata
validation; partial/user files are not overwritten.

## Why this frequency and amplitude

v28's equilibrium planning gives SNR3 about827 ps normal/457 ps slip at4RMS
under linear stationary assumptions. The900 ps analysis interval targets this
measurement, not a chosen sign/result. Amplitude is experimental forcing, not
a kinetic coefficient. Its linearity must be tested rather than assumed.
At higher frequencies, the storage response was incompatible with the same
static-PMF overdamped model; this probe is at a lower frequency to test that
limitation directly. It does not by itself establish omega->0 convergence.

## Required analysis

1. Every signed record, full and half analysis windows; no favorable-window selection.
2. Exact engine time=step*dt and saved clock; conjugate work versus internal energy.
3. First/last100 ps and overall temperature, not hidden NVE heating.
4. Same unforced5 ns record, all complete900/450/400/200 ps null windows.
5. Complex response and inverse-drag disk bounds; no abs/clip or finite upper
   mobility when the empirical error set includes zero drag.
6. Second harmonic even-in-force and third harmonic odd-in-force, compared
   with matching unforced harmonic envelopes. Raw individual harmonics retained.
7. Force-amplitude and timestep differences, without separating deterministic
   error from chaotic thermal fluctuations unless evidence permits.
8. FDT spectral loss with NW2/3, full/half records, +/- .005/.01/.02 bands.

An empirical null envelope is NOT a confidence interval. Finite-band phase,
constant mobility, local-coordinate correspondence and real-Al kinetic
calibration are distinct gates. None automatically enables production Hz.

Run `python -m solver_v1.run_low_frequency_forcing_v29` with prepare/run/analyze
subcommands. Actual results, regression status and final Git state are recorded
in CURRENT_WORK_HANDOFF.md after completion.

## Actual results

Signed-pair susceptibilities use exp(+i omega t); negative phase means lag.
Units of susceptibility are Angstrom squared/eV. All below are reference MD
at 0.2 cycles/ps, NOT production physical-Hz predictions.

| Coordinate | dt (fs) | Re chi | Im chi | phase (deg) | loss / observed null |
|---|---:|---:|---:|---:|---:|
| Normal,4RMS | 2.5 | .00236788977 | -3.736226e-6 | -.0904053 | .9752 |
| Normal,4RMS | 1.25 | .00236481994 | -6.783835e-6 | -.1643610 | 1.7707 |
| Slip,4RMS | 2.5 | .01219256767 | -5.051132e-5 | -.2373634 | 1.8696 |
| Slip,4RMS | 1.25 | .01218669673 | -3.842335e-5 | -.1806469 | 1.4222 |
| Normal,2RMS | 2.5 | .00236061093 | -2.971428e-6 | -.0721212 | .2381 |
| Slip,2RMS | 2.5 | .01207293539 | -3.710131e-5 | -.1760750 | .5220 |

The finer full-force storage/static ratios are1.00296 normal/1.02231 slip,
much closer to static than v27's high-frequency probes. This does NOT certify
the exact same-PMF overdamped response. The loss changes about81.6% normal
and23.9% slip when dt is halved; these changes remain inside the summed
complex unforced envelopes, so deterministic timestep error and thermal
trajectory fluctuations cannot be separated here.

The independent equilibrium-spectrum sensitivity ranges for Im chi are
[-5.54253,-3.63620]e-6 normal and[-4.89928,-4.09960]e-5 slip. Both finer
driven loss intervals overlap these ranges. This is compatibility with FDT,
not a precision calibration or independent-replica confidence statement.

### Stronger forcing reveals nonlinear contamination

Normal even second harmonic is .000496761/.000484808 Angstrom at full force
for base/fine dt,6.414/6.260 times the corresponding observed null; half
force gives .000170679 Angstrom,1.807 times its (shorter-window) null.
Slip odd third harmonic is .000372382/.000267738 Angstrom,2.022/1.454 times
null; half force .000051430 Angstrom remains below null. Full/half slip
fundamental complex difference is1.058 times the summed null envelopes.
These are observed evidence of force-dependent/nonlinear response, not
formal significance levels. Increasing force did not give a clean linear
mobility experiment. Neither fitting an arbitrary damping nor selecting a
favorable sign/window repairs this.

### Work, heating and numerical controls

All12 engine clocks satisfy t=step*dt; max residual1.137e-13 ps. Saved force
error<=5.97e-12 eV/Angstrom; COM drift<=1.91e-14 Angstrom/ps. Whole-record
integration-by-parts work/internal-energy discrepancy<=.009527 eV. Direct
power sampled only every25 fs has discrepancy<=.061785 eV, and does NOT
converge merely by halving the internal timestep. This is an unresolved
power-quadrature control; it must not be called physical dissipation.
The integer-cycle phase/work identity is within2.51e-12 eV, but is algebraic
on the same q data, not independent validation. Endpoint work is saved
separately from fitted periodic work. First-to-last100 ps mean temperature
increases range .0682–5.0238 K; no thermostat was added to hide this heating.

### Mobility inference and decision

Point inverse-response mobilities normal base/fine are1.1770e11/6.4658e10,
slip2.3084e11/3.0317e11 m^2/(J s). These are finite-frequency scalar points.
The empirical complex control radius is the sum of fine null, dt difference,
amplitude difference and full-force half-window difference. Its inverse
image includes zero drag for BOTH coordinates. Mobility upper limits are
therefore null/unbounded, not finite calibrated constants. The control set
is deliberately transparent, not a statistical confidence interval.

The next defensible measurement is weaker forcing with longer records and
independent equilibrium restarts, explicit temperature dependence/control,
and internally accumulated work rather than coarse saved power. A subsequent
zero-frequency/common-local-coordinate test is still required. This is a
follow-up proposal, not an unperformed experiment reported as complete.
Production LJ/Bessel/PDE, M_a/M_s, static parameters, kinetic JSON and UI are
unchanged; physical seconds/Hz stay disabled. No actual Al fatigue/yield
validation is claimed.

Results: results/low_frequency_forcing_v29 (all signed/window responses,
nulls, FDT sensitivity, work checks, control comparison and inspected plot).
Tests actually executed: targeted27 passed3.79s; solver700 passed1993.89s;
app34 passed154.50s; desktop smoke passed. No tests removed or skipped.

## Methodological follow-up, not an Al calibration source

Vroylandt and Monmarché, *Position-dependent memory kernel in generalized
Langevin equations: theory and numerical estimation*, JCP156,244105(2022),
[DOI10.1063/5.0094566](https://doi.org/10.1063/5.0094566),
[author preprint](https://arxiv.org/abs/2201.02457), derives a position-dependent
memory/FDT framework and Volterra estimation route. This is a methodological
reference if memory prevents a constant-mobility reduction, **not a source of
Al mobility values**. No such kernel is substituted into production here.
