# range_core_v12 — material-range and isolated-core research

These are STATIC LJ/Bessel calculations, not experimental yield, a finite
activation energy, calibrated fatigue or physical-Hz probability results.
No production PDE/UI parameter set was replaced.

## Material

- `material/target_definition.json`: unchanged matched0K source targets,
  normalization, fit/heldout roles and source-file SHA. NIST Al99/Mishin1999
  is a comparator only; its tabulated pair is not our analytic potential.
- `material/calibration.json`: actual39 predefined starts and292 coefficient
  profiles, both metrics and optimizer stopping status. One added rank2 radial
  range. The best evaluated point is retained even when an optimizer returns
  a worse point. This is not a global optimum or family-wide impossibility proof.
- `material/residuals.csv`: each individual residual with role/units/scale;
  no held-out curvature is used in the loss. No yield/source-size/kinetic target.
- `material/identifiability.json`: coefficients and log radial sensitivities;
  FD step refinement, local SVD/correlations, NOT confidence intervals. The
  tied-log-range tangent equals the old nested tangent only on equal ranges.
- `validated_summary.csv`, `stationary_states.csv`: independently re-solved
  perfect/fault/index-one stationary points. A local saddle is not a proof of
  a global transition path. Energies are per atomic interface area in J/m2.
- `finite_q_validation.csv`: independent direct-radius6/8/10 static force
  matrices along corrected cubic Gamma-X/L/K paths. Eigenvalues eV/L0^2 are
  NOT frequencies. The exact-C candidate is unstable despite uniform moduli.
- `static_load_unload.csv`: resolved normal/shear tractions in MPa, not nominal
  axial stresses without orientation. This static return is not a kinetic hold.

## Isolated core

`isolated_core/<case>/` preserves each actual optimizer run:

- `metadata.json`: material hash, SAME-potential elastic tensor, microscopic
  L0, free disk, transverse ring, fixed far-field condition and units.
- `checkpoint.csv`, `progress.json`: interrupted-run recovery/actual progress.
- `state.csv`: free-row logical indices, three actual displacements in L0,
  analytic energy gradients in eV/L0. NOT rescaled strain or probability.
- `summary.json`: energy per infinite-row repeat b*L0, actual free-force
  residual, solver stopping status, winding, numerical controls, elapsed time.
  Later runs include mandatory `final_stability_probe`. Earlier small-force
  centered states were discovered to be saddles; never infer stability from
  their optimizer or force flag. `*_mode_plus/minus` are explicit zero-load
  negative-mode controls, with eigenpair and finite-difference energy evidence.
- `radial_energy.csv`: the initially recorded partial site energies. Early
  cases predate the explicit linear-reference-work partition. Preserve them;
  use the consolidated `core_radial_energy.csv` for the audited partition.

The x atomic sum remains infinite Poisson/Bessel. Transverse environment
changes are neighborhood-refined; no total-tail certificate is inferred from
the reciprocal tolerance. The finite free disk has artificial image forces
and is NOT a measured pin spacing/source. All affected fixed-site embeddings
are included. Radius/ring changes are explicit numerical studies.

## Consolidated re-evaluation

`report_range_core_repair.py` reads completed states and re-evaluates analytic
forces/Hessians/partitions. It does not claim to rerun the fitting or relaxation.

- `core_summary.csv`: all actual conditions and convergence flags.
- `core_analytic_audit.csv`: same-row affine curvature vs full-bulk elasticity,
  exact summed linear-work derivative, saved-energy replay and selected lowest
  analytic Hessian eigenvalues.
- `core_morse_audit.csv`: actual analytic lowest-mode and independent two-sided
  energy-curvature checks at two displacement steps, preserving both initial
  negative-mode and final-state probes where available. A finite force stopping
  status alone does not imply a minimum.
- `core_domain_ring_comparison.csv`: inner fields compared on IDENTICAL row
  indices/radii. Domain enlargement is not compared via different total counts.
  x translations are compared modulo the actual row repeat; y/z are unscaled.
- `core_radial_energy.csv`: raw partial E, Taylor-linear reference work,
  quadratic-remainder partial E and `E_remainder-k ln(R/b)`. The sum of the
  linear term has zero free derivative and, for this boundary, zero total.
  This is an exact declared site partition, not a fitted force/energy correction.
- `core_vector_disregistry.csv`: actual three-component displacement jump
  between layers0/1. No scalar well index is assigned to vector partials.
- SVG figures visualize these calculated quantities, not a specimen colormap.
- `decision.json`: material, straight-core numerics, finite-source and physical
  yield/kinetics status kept separate.

The final report includes 21 completed actual relaxations and nine cases
passing the tested fixed-boundary Morse/force checks. The original saddle
unload coincides with a zero-load control to 5.32e-8 L0; stable-initial-state
25 -> 50 -> 0 MPa returns within 1.015e-7 L0. These are inner-row displacement
differences, not macroscopic residual plastic strain. Stable-domain/ring
convergence is incomplete. `verification.json` records actual regression runs
separately from scientific adoption, which remains false.

The core finite part is referenced to the already defined Burgers translation
b. It is not a guessed new core radius, nor a complete line-character-dependent
core energy. Do not put it into a finite-source activation probability without
deriving the missing finite geometry, mixed/edge character and kinetics.

All prior `yield_bridge_v11`, source datasets and historical parameters remain
unchanged. See `solver_v1/RANGE_AND_CORE_REPAIR.md` for the derivation and limits.
