# Kinetics, load components and the GPa question

## Status and actual implementation

This audit starts from `09e0dc6ce85e0a81656a92995cb3c0816bb96e97`.
It does not change an energy parameter or promote a research model to the PDE.
Three different notions of dimension must be distinguished:

| Path | Load support | State / spatial dimension | Status |
|---|---|---|---|
| Desktop `solver_adapter` | One axial sigma/E, fixed chi=.2 | P(a,s,t), two state coordinates | Production reduced reference |
| `low_stress_cyclic_diagnostic` | Independent normal/shear amplitude, mean, phase | Same two state coordinates; reflecting boundaries | Research hypothesis test, no validated crack sink |
| `resolved_tensor_tractions` | Symmetric 3x3 Cauchy stress -> normal + two shears | Geometry projection only | No vector-registry PDE implied |
| Full FCC geometry | Two-dimensional registry plane, named scalar paths | Infinite 3D atomic geometry, static evaluations | Research reference |
| Discrete screw v6 | Anti-plane displacement on transverse rows/layers | 2D cross-section, one displacement component, uniform along line | Harmonic static research, not arbitrary 3D dynamics |

The old desktop direction entry was validated then discarded: it was not sent
to `UIAnalysisConfig` or the PDE. Its label falsely implied an orientation
projection. The entry is now explicitly disabled/not connected; scope and
result capability metadata describe the actual solver. No old force formula
is silently replaced. N>1/TT code elsewhere does not prove a calibrated vector
slip or spatial specimen solver.

## Seconds and Hz: available software, missing Al kinetic evidence

The conversion remains exactly

\[
M_i^*=t_0 E_0 M_{i,phys}/L_0^2,\qquad
t_{phys}=t_0t^*,\qquad f_{Hz}=f_{model}/t_0.
\]

The committed Al calibration still has **null** physical mobilities and clock.
No automatic fabricated seconds/Hz are enabled. Actual imported kinetic data
can now reach the UI through **Load kinetic calibration** instead of requiring
manual edits of the repository default. A file is bound to a model ID, static
parameter fingerprint and the same two-row cell coordinate definition.
Hybrid eV units and the Al-target length scale are checked. Unknown/stale
model, wrong coordinates or inconsistent t0/mobilities are rejected.

The old physical-mode path only converted labels/time while constructing
M=(1,.05), kT=.02 regardless of the calibration. That was unsafe for a different
measured mobility ratio or temperature. Now physical mode constructs the SAME
energy/PDE with the declared M* and kBT/E0. Model-time defaults and all static
energy/strain calibration remain unchanged. This is not tuning to get fatigue.

Loading a file does not rerun analysis, alter a prior result, or relabel its
axes. Switching basis converts the entered frequency to preserve the same
model period. Changing energy model disables an incompatible clock. A file
consistency check is not independent validation of its supplied provenance.

## Actual correlation-data import

For near-harmonic stationary overdamped collective coordinates,

\[
C(t)=\langle\delta q(t)\delta q(0)^T\rangle=e^{-Bt}C(0),\quad
B=M H_{phys},\quad C(0)=k_BT H_{phys}^{-1}.
\]

Thus the **ordered matrix** relation, not two separate raw-coordinate decay
fits, is

\[
B(t)=-\log[C(t)C(0)^{-1}]/t,\qquad M=B H_{phys}^{-1}.
\]

`kinetic_calibration_workflow.py` compares multiple positive lags, equilibrium
covariance, diagonality/positivity of M, and reconstructed C(t). The input
explicitly declares its measurement acceptance tolerance. Complex logarithms,
lag inconsistency, wrong fluctuation amplitude and appreciable off-diagonal
M are refused, not projected away and declared calibrated. Tests use exact
synthetic matrix correlations, not Monte Carlo data.

One clock chooses M_a*=1; measured M_s,phys/M_a,phys supplies M_s*. It is not
forced to .05. The rates-only helper remains, with its mode-assignment ambiguity.
Original trajectory block/window uncertainty, thermostat effects, inertial
transients and collective-coordinate normalization require scientific review.
Passing a matrix test does not by itself establish a Markovian Al reduction.

Input covariance CSV columns (include lag zero and >=2 positive lags):

```
time_seconds,C_aa_m2,C_as_m2,C_sa_m2,C_ss_m2
```

Metadata JSON must contain `model_id`, `length_scale_m`, `energy_scale_J`,
`temperature_K`, `source`, `relative_tolerance`, optional `notes`. Values must
come from the actual measured collective coordinates; no numerical Al example
is supplied because such data are missing. Run:

```
python -m solver_v1.import_kinetic_calibration metadata.json covariance.csv \
  --output measured_kinetics.json --diagnostics correlation_checks.json
```

Existing output files are not overwritten. Select the corresponding energy
model, load the resulting calibration, then select Physical time. Physical
mode changes temperature/mobility if the data require it; it is not an axis-only
trick. Tests verify the unchanged-generator special case exactly.

## Why published dislocation mobility cannot just fill M_a and M_s

Olmsted et al., *Atomistic simulations of dislocation mobility in Al, Ni and
Al/Mg alloys* (2005), DOI [10.1088/0965-0393/13/3/007](https://doi.org/10.1088/0965-0393/13/3/007),
[author manuscript](https://arxiv.org/abs/cond-mat/0412324), studies the position
and velocity of edge/screw dislocations, including low-speed phonon drag. It
is relevant to a future derived defect-coordinate dynamics, not evidence for
the mobility of a coherently displaced intact interface cell.

For that distinct line coordinate X, B_drag v=tau b. Since B_drag has Pa s
units, 1/B_drag has m^3/(J s), the mobility conjugate to energy PER LINE LENGTH.
The present a/s cell mobility has m^2/(J s). Even a finite-line conversion
1/(B_drag L_line) would require an independently specified physical line and
the mapping from X to the full slip/core field; it supplies no normal mobility.
No arbitrary length or A_c is inserted to make the units look compatible.

## Why the current stresses reach GPa

For the atomistic active interface only, the work conjugate is

\[
G_{int}=W_{int}-\frac{A_{atomic\,cell}L_0}{eV_J}
 [T_n(a-h)+\tau s],
\]

with tractions in Pa. Therefore

\[
\tau_{MPa}=\frac{eV_J W_s}{10^6 A_{atomic\,cell}L_0}.
\]

`unit_audit.csv` independently multiplies Pa*m^2*m/J_per_eV and tests 4,15,25,
50,100,2000 MPa. The direct SI work difference is zero at the saved precision;
round trips have ordinary floating error. For example 4 MPa supplies only
0.0005078089 eV over one reduced unit of slip. There is no factor-1000 error.
The production reduced stress mapping remains kappa*sigma/E, NOT this separate
research work expression.

Two earlier comparisons were physically different:

1. A perfect infinite rigid interface slipping coherently is a nucleation/ideal
   branch problem, not motion of existing dislocations at a specimen yield load.
2. Clamping a while varying s artificially prevents normal accommodation.

For fixed applied normal traction, follow a_n(s) solving W_a=f_n, W_aa>0.
The normal-relaxed potential and its EXACT derivatives on that branch are

\[
\Phi(s)=W(a_n(s),s)-f_n(a_n(s)-h),\quad
\Phi'=W_s,\quad
\Phi''=W_{ss}-W_{as}^2/W_{aa}.
\]

The first positive-to-negative Schur-curvature zero along the pristine branch
is a scalar-path shear spinodal. It is not every maximum later on a complete
registry curve, and not a fully relaxed vector-registry instability. The new
branch code follows continuation, requires force/stability residuals and
refines brackets. No force clipping or parameter change is used.

Actual zero-normal-traction results (33/65 point bracket studies agree):

| Energy, same interface conditions | direct_110 first normal-relaxed peak | shockley_112 first normal-relaxed peak |
|---|---:|---:|
| Unchanged analytic research candidate | 4768.897903 MPa | 2685.841510 MPa |
| NIST Mishin reference, comparison only | 4390.365017 MPa | 2345.062485 MPa |

The candidate direct110 fixed-gap first maximum is 7543.359226 MPa, falling
36.78% under normal relaxation. Same-condition source values also remain GPa.
The candidate Shockley first fixed-gap maximum is3012.198345 MPa, not the
much higher later registry peak. The source's coarse fixed-gap peak bracket
contained multiple curvature roots: one whole-interval root solve returned a
later2400.850 MPa maximum instead of the first2405.245 MPa maximum. Subdividing
the bracket and comparing33/65 sampling resolves the first root consistently.
Source/candidate differences in the relaxed peaks are ~8.6%/~14.5%; these do
not validate other branches or the full material model. In particular candidate
C11/C12/C44=87.816/64.833/49.271 GPa versus source114/62/32 GPa remain unfitted
errors (~-23.0%,+4.6%,+54.0%). No scalar multiplier can repair all three.

Roundy et al., *Ideal Shear Strengths of fcc Aluminum and Copper*, PRL82,2713
(1999), DOI [10.1103/PhysRevLett.82.2713](https://doi.org/10.1103/PhysRevLett.82.2713),
report a 0-K DFT ideal shear strength of 1.85 GPa for Al with five other strain
components relaxed. That is context for the ideal-strength scale, **not** a
matched target for these rigid half-crystals or an experimental yield value.
Clatterbuck et al., PRL91,135501 (2003),
DOI [10.1103/PhysRevLett.91.135501](https://doi.org/10.1103/PhysRevLett.91.135501),
also show why finite-wavevector instability may precede elastic instability.
One local positive Hessian or scalar spinodal is not complete material stability.

Prior nonlocal/discrete research does find MPa-scale forces for existing defect
pairs, but these are static driving/holding forces, not validated defect
velocities or fatigue life. Nonlinear vector cores, normal relaxation, finite
activation geometry, material fit and physical kinetics remain to be joined
consistently. The GPa issue cannot honestly be fixed by dividing stress by1000.

## Reproduction, plots and uncertainty

```
python -m solver_v1.run_stress_scale_audit
python -m solver_v1.report_stress_scale_audit
```

`results/kinetics_loading_audit/` contains all raw curves, first-branch peaks,
SI conversions, unchanged parameter hash, source hash, elastic errors and plot.
Normal branch residuals are below1.7e-14 eV/reduced-coordinate in this scan;
Schur residuals are below1e-12 in its curvature units. Refinement of33->65
checks root brackets, not all missing atomistic degrees of freedom or material
uncertainty. The earlier infinite-series validation remains applicable; no
energy derivative was changed here.

## The two previously skipped UI tests

They tested actual Tk windows: language/result/zoom preservation and specimen-
area postprocessing without a PDE rerun. They were not failing scientific
solver equations. Tcl initialization failed under the execution sandbox even
though the installed init.tcl was present. Approved native Tk execution
successfully ran both. Mixed runs also exposed intermittent native Tcl script
read errors; these are failures, not passed/skipped tests. The tests now only
skip a genuine headless Unix display condition; Windows Tcl errors propagate.
Final executed outcomes and timings are in the checkpoint validation record.
No machine-specific Tcl path is hardcoded and no system installation is edited.

Repeated mixed test runs also showed native script access failures on the
second/third separate Tk interpreter. Tests now share one module-scoped Tcl/Tk
interpreter and create/destroy an isolated Toplevel workspace per test.
`DesktopApp(root=...)` permits that lifecycle; normal startup still creates
one Tk root. Assertions, result isolation, and GUI rendering are not removed,
and no retry masks errors. This follows the single-interpreter recommendation
in the [Python Tkinter threading documentation](https://docs.python.org/3/library/tkinter.html#threading-model).
It removes repeated native initialization from the test design; it is not a
claim that the intermittent platform file-access cause has been proved in full.

The desktop smoke checks imports/conversion only. Actual window rendering
must be reported separately, never inferred from smoke success.

## Executed final validation

- Solver kinetic/stress targeted: 17 passed in 2.32 s.
- Final shared-interpreter Tk/i18n targeted: 10 passed in 2.37 s.
- Full solver_v1: 271 passed, zero skipped, in 1236.63 s.
- Full app: 31 passed, zero skipped, in 161.26 s, including actual Tk rendering.
- Desktop startup smoke: passed in 1.487 s.
- git diff --check: passed before staging; checked again for the final commit.

These are wall-clock test runtimes, not calibrated physical model times.
The earlier native Tcl failures and the eventual lifecycle correction are
recorded in `results/kinetics_loading_audit/validation_status.json` and the
handoff; they are not retrospectively counted as successful test runs.
Passing software checks does not remove the missing Al kinetic data or the
material-fit/defect limitations identified above.
