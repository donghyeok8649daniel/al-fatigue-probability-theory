# yield_bridge_v11 — measured-yield route, not a yield calibration

This directory contains actually executed static research and source-data
reading. No production potential, PDE, mobility, area aggregation or UI change
is implied. Read `solver_v1/YIELD_STRENGTH_BRIDGE.md` for derivations/limitations.

## Authoritative outputs

- `experimental_benchmark/`: seven isolated Krebs2017 Figure2d markers at
  0.002 plastic shear, in the *published normalized* CRSS/G units. Missing
  modulus and MPa values are empty/null. Wire diameter is not pin spacing.
  Provenance includes source/figure/SI hashes, axis calibration, reading bands
  and selected pixel boxes. The PDFs/figures themselves remain ignored cache.
- `material_metric_refined/`: final actual fit using the independent-C metric
  and the active-face exact-bulk QP. `calibration.json` contains all radial
  profiles, starts, selected bounds, residuals and convergence/budget status.
  `progress.json` is the final completed progress snapshot, not another fit.
  `observations.json` binds unchanged0K targets/roles/historical source hashes.
  `residuals.csv` reports both old-H and new-C metrics; losses from different
  metrics are not directly interchangeable. No strength observable is fitted.
- `material_metric/`: **superseded initial calculation**, retained as an audit
  record. Its SLSQP exact-bulk implementation was sensitive to roundoff near
  coefficient bounds. The refined active-face run found a lower feasible
  exact-bulk/interface loss. Do not use the initial result as the final fit.
- `material_summary.csv`: actual bulk C, cohesion and interface root/endpoint
  values for historical and final candidates. Every `adopted` field is false.
- `material_stationary_states.csv`: stationary q=(a,x,y)/L0, energiesJ/m²,
  gradients eV/L0, Hessian eigenvalues eV/L0² and actual root/index validity.
  A valid root is not global material stability or a proved global MEP.
- `material_derivative_check.csv`: independent energy/gradient central FD,
  not canonical finite-difference forces. Steps and errors carry explicit units.
- `finite_source/elastic_metadata.csv`: exact model-derived cubic constants
  and leading-log line coefficients [J/m]. Source elasticity is a separately
  labeled0K comparator, never swapped into the LJ candidate without notice.
- `finite_source/angular_coefficients.csv`: line orientation, k and k+k''.
- `finite_source/refinement.csv`: angular comparisons against128samples;
  32/64/128/256segment graph solutions against the independently derived
  parametric branch. Angular fine-reference self differences are zero by
  definition, not an absolute error certificate. Shape error is max-bow/L,
  not a full pointwise curve norm; force residual is normalized by gamma_ref.
- `finite_source/source_shapes.csv`: actual analytic half-branch profiles,
  x/L and y/L, usingL=1micrometre, R=L, r_core=b as explicit hypotheses.
- `finite_source/hypothetical_source_thresholds.csv`: leading-log *outer-only*
  double-ended bow-out folds at predeclared source spans and cutoff ratios.
  NOT measured source lengths, fitted yield, single-arm sources or core barriers.
- `finite_source/stress_scenarios.csv`:125 combinations of model, hypothetical
  span and2/4/10/25/50MPa resolved shear.53 subcritical graph solutions;
  72 above-fold exclusions do not claim atomistic operation or multiplication.
- `finite_source/status.json`: completion, elapsed time and explicit limitations.
- TwoSVG figures plot the computed data. No experimental full figure is copied.
- `decision.json`: scientific decision. A finite-source leading-log reference
  was derived/checked; measured yield and material adoption both remain false.
- `verification.json`: actual tests and numerical audit, not forecast counts.

## Gauge / sensitivity caveats

All quantities use the existing per-atom/per-interface normalization. Parameter
order is(u,v,A,B,C,D3,D1,D2); u/v are the existing LJ power coefficients, not a
new pair law. Fixed-decay coefficient SVD does not include radial-parameter
uncertainty. The exact-bulk3-column Jacobian precedes active-A tangent removal;
the6-column comparison in `identifiability.json` removes only the original two
exact constraints for common comparison. Neither supplies material confidence.

The chosen positive-opening-force shape prior fromv10 is **not imposed** here:
the source itself has a verified negative traction lobe. Neither of the new
material candidates passes joint elasticity/interface/held-out validation.

`null` means unavailable/not justified, never zero physical yield. All source
sizes are explicit geometry hypotheses. No A_c, kinetic mobility, seconds,
Hz, empirical hardening, activation-volume fit or fatigue-life calibration is
introduced. Static elastic MPa does not require a kinetic time scale.
