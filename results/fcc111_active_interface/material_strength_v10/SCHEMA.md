# material_strength_v10 — separate research, NOT an adopted material

Read `solver_v1/MATERIAL_TO_SPECIMEN_STRENGTH.md`. Static material calibration
and actual finite-source/specimen strength are different gates. No production
energy selector, probability generator, static historical parameters, kinetic
calibration or UI model was replaced.

## Folders and execution provenance

* `joint_fit/`: completed 25-point declared radial grid plus old decays and
  160 Powell calls, exact geometry/cohesion elimination and constrained
  coefficient solves. Powell exhausted its budget; no claim of a global or
  fully converged radial optimum. All actual trial parameters/status are saved.
  `calibration.json` and `profile_progress.json` are intentionally redundant
  final/checkpoint records, reasonably small. No surrogate trajectory replay.
* `opening_exchange/`: at that fixed radial choice, actual W_aa roots reveal
  negative W_a between fit points. Three coefficient re-solves add those
  locations as inequality constraints. Raw negative residuals and computed
  roundoff/tolerance errors remain visible; forces are never clipped. This is
  zero-registry path validation, not global/all-slip-path certification.
* `experimental_benchmark/`: short provenance and three readable points from
  a primary 99.99% Al wire experiment. The gamma=.2 point and 0.2%-criterion
  CRSS are unavailable, not interpolated or set to zero. The PDF and figures
  remain in ignored cache. No strength target is used in energy fitting.

The **four initial-attempt files in this directory root**:
`target_definition.json`, `profile_progress.json`,
`fit_and_heldout_residuals.csv`, `identifiability.json` are INCOMPLETE/SUPERSEDED.
They preceded the completed workflows above. The first fit-report attempt
failed because the bulk-only stage had zero D1/D3 sensitivity columns. It was
fixed by removing unobserved columns from that ablation, not by inventing data
or perturbing a singular column. The incomplete checkpoint explicitly remains
`completed=false`; do not call it an accepted calibration.

The first exchange launch then failed on a SOURCE plane with zero neighbors:
the empty gradient einsum returned NaN in the captured Windows/NumPy2.5 run.
The neighbor radius, spline values and derivative arrays were all empty;
the three gradient components were NaN while scalar/Hessian were zero. A
standalone repetition did not reproduce NumPy's allocation-dependent internal
failure, so its library internals are not claimed to be diagnosed. The code
now returns the mathematically exact zero jet for an empty SOURCE neighbor set.
Nonempty nonfinite sums still fail. This does not add a cutoff to LJ/Bessel or
change the published SOURCE cutoff. Thirty separated evaluations and a test
forbidding any empty einsum protect this exact identity.

## Interpretation of principal files

`target_definition.json` inside `joint_fit` contains all 25 source observations:
two exact, fourteen fitted, nine held out. Bulk values are 0 K rounded Mishin
elastic constants; interface values are calculated from the source file under
identical rigid-half geometry. Five/ten-percent and force discrepancy scales
are declared normalization, not experimental uncertainty. The reverse
barrier is dependent on saddle/fault energies and is not double-counted.

`fit_and_heldout_residuals.csv` compares candidate values **at source states**.
`stationary_states.csv` and `barrier_comparison.csv` instead solve each
candidate's own actual stationary states and downhill-connected basins.
The two types of comparison are not interchangeable. `ideal_folds_NOT_yield`
and `fold_independent_refinement` are static uniform-interface instabilities,
never measured specimen yield stress.

`opening_traction_extrema.csv` brackets W_aa=0 with 91/181 grids.
`extrema_independent_refinement.csv` uses 121/241 and tighter reciprocal
tolerance. The entire tested interval is 1<a/h<=12 at ux=uy=0. Sign-root
scanning alone does not prove the absence of every double root or feature
outside this interval. `opening_curve.csv` retains every sampled raw force.

`finite_q_stability.csv` is static Hessian eigenvalues in eV/L0^2, not frequency
or a proof over the whole Brillouin zone. Corrected crystallographic paths
and two independent direct cutoffs are reported. Canonical energy is infinite
reciprocal; `independent_direct_lattice.csv` explicitly retains the finite
validation error at radius/layer24,48,72. Do not call its algebraic tail zero.

`source_opening_extrema_audit.csv` records a SOURCE negative-traction lobe
(-785.6754MPa at a/h2.386631), verified with181/361 brackets, energy finite
differences and independent scalar-source values. Nonnegative traction in the
fit is therefore a **declared monotone-energy shape prior**, not a universal
physical law. `opening_shape_prior_ablation.json` removes that prior at fixed
decays and still shows large elastic errors. This does not certify the source
as real experimental cleavage, nor diagnose the source's cutoff/interpolation
versus inherent energy-family contribution to its oscillation. Preserve the
raw curve; do not clip either source or candidate forces.

`static_stress_and_unload.csv` records independent normal/x-shear/y-shear MPa
probes. `experimental_stress_range_static_check.csv` additionally probes the
MPa levels in the experimental plot. These are zero-K **local rigid-half
equilibria**, not finite wires or physical-time unloading. A static return to
the initial registry is not a dynamic fatigue validation.

`identifiability.json`: six coefficient directions after eliminating two
exact bulk constraints. In `joint_fit`, two additional log-decay columns are
included and step-refined. In `opening_exchange`, decays are held fixed and
the six-column result is explicitly conditional. Active inequalities and
model discrepancy preclude treating an unconstrained SVD as confidence limits.

Top-level `parameter_sets.csv`, `fit_quality_summary.csv`,
`strong_sensitivity_correlations.csv`, `strength_comparison_status.json` and
`material_decision.json` are generated by the reporter from completed runs.
The SVGs plot these calculations/short digitized points, not invented fields.
`verification.json` is the final actual test/execution ledger, not a material
readiness certificate.

## Invariants

L0=4.05/sqrt(2) angstrom; energy eV per interface cell; atomic interface area
sqrt(3)L0^2/2. Surface energies are J/m^2, physical tractions MPa, and static
curvatures eV/L0^2. A_c never enters these conversions. All original parameter
files are preserved. Mobility, seconds/Hz, finite-source strength, fatigue,
specimen spatial probability and UI promotion remain unavailable/unvalidated.
