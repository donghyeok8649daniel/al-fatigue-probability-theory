# v32 actual compatibility experiments

These are necessary static compatibility tests, not adopted material sets.
Read `solver_v1/RADIAL_CHANNEL_COMPATIBILITY_V32.md` for equations and units.

- `completed_trial`: 14 actual baseline/single-angular/dictionary LP profiles.
- `positive_power`: 9 actual alternative density-power profiles, same shapes.
- `density_channels`: 10 actual scalar-channel profiles. Linear extra density
  amplitudes are signed, quadratic amplitudes nonnegative; the stored
  `nonnegative` arrays in each profile are authoritative.
- `cross_square`: 10 actual PSD nonlinear-sector profiles.
- `coherent_vector`: 10 actual per-site mixed-vector profiles and analytical
  unit-coefficient finite-q tail bounds; moments are mixed before squaring.
- `vector_density_shape`: 8 positive coordination-polynomial vector profiles.
- `independent_shape`: 180 same-family profiles, completed at maxfev rather than
  optimizer convergence; eta12.899739. `independent_validation` is actual full
  matrix and excluded-q replay, not a new fit.
- `reference_density_shape`: four actual auxiliary-density fits (one simple
  exponential and three deterministic starts of the positive quadratic shape).
  It is NOT a material fit; source EAM density has an arbitrary auxiliary gauge.
- `quadratic_density_material`: one fixed scalar-shape material profile,
  eta19.349784, failed. New scalar density enters both embedding and screening.
- `quadratic_density_joint`: separate direct material-target shape search;
  120 evaluations, four explicitly failed LPs, maxfev (not convergence).
  Only five profiles retain positive LJ; best eligible eta19.349784. Lower
  objectives with a zero LJ coefficient are not eligible candidates.
  See `QUADRATIC_EXPONENTIAL_ENVIRONMENT_V32.md` for the exact transform/units.
- `quadratic_density_validation`: actual fresh-state, analytic-derivative,
  reciprocal/layer tolerance and finite-q radius/source comparisons. Failed
  material fit despite numerical derivative checks passing.
- `embedding_convexity_audit`: four fixed-shape coefficient-sign controls.
  Signed scalar quadratic leaves eta unchanged; signed angular quartic only
  reduces it to12.404692. No sign-relaxed energy is adopted or certified stable.
- `pair_gauge_shape_bound_audit`: three auxiliary pair/gauge value fits and
  excluded midpoint/derivative checks. Even the target-only ideal source
  density column does not make LJ equivalent to the source pair within the
  stated discrepancy. Both gauge fits activate the attractive-LJ lower bound.
  `pair_gauge_shape` is the preserved first output with a floating-nonzero-only
  positive_LJ flag; the bound-audited rerun is authoritative for interpretation.
- `quadratic_density_all_ranges`: broader eight-shape joint research search;
  completed320 profiles/4714.57s,318 feasible LPs,166 positive-LJ profiles,
  best eligible eta15.270143. Maxfev, not optimizer convergence.
  The four original validation
  states have already been inspected. The additional manifest declares four
  more states after search start but before final candidate/source inspection.
- `quadratic_density_all_ranges_validation`: actual new4/previously inspected4
  state validation. Fresh H error100.5632%, excluded-q249.6668%; not adopted.
- `spatial_density_gradient`: actual two-profile test of the true spatial
  density gradient, not the existing first moment. eta12.507100, no solution.

Each completed fixed-shape study saves its predeclared definition, actual
coefficient vectors, observations, excluded-q errors, SVD and independent LP
primal/dual diagnostics. No fresh off-grid validation gate has been passed.
The six completed channel studies contain61 profiles including baseline
replays. No eta<=1 candidate was found. The larger dictionaries are capacity
diagnostics, not recommended many-parameter energy laws.

The root `definition.json` is an explicitly preserved **failed first attempt**:
the observation constructor rejected role `excluded` before any LP ran. The
exporter now maps it to the valid `heldout` role while preserving exclusion
from the fit. Do not count that root definition as a completed study.

Raw precision/condition/radius meanings follow the v31 definitions. Source Al99
is target-only. Physical time, production a/s mobility, actual yield and
fatigue remain unvalidated. No production/UI/Ac inputs are altered.

Source bindings retain the historical loader context; its `energy_unit` is
not a new observable-unit declaration. Actual material outputs are eV/atom
(bulk), eV/interface atomic cell (interface), GPa (cubic constants), and
eV/L0^2 (displacement Hessians). Auxiliary radial/gauge studies instead use
Angstrom radii and eV/Angstrom^order pair derivatives. Every observation's
explicit `units` field and the derivation take precedence over loader context.
