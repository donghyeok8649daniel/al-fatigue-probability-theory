# v19: executed coordination-dependent environment calibration

Status: **completed research comparison; NOT accepted as an Al material**.
The original LJ/Bessel energies, production PDE/model registry, saved static
calibration and kinetic calibration are unchanged. No empirical yield,
fatigue-life, mobility, physical-time or specimen-area fit was performed.

Start commit: `53ced52f86e687bd3679bb85f09cd647cfc55913`.
See `solver_v1/COORDINATION_SCREENING_V19.md` for the explicit per-atom
derivation, source/units, staged fits, limitations and numerical results.

## Result routing

- Root: rational law, fixed inherited five shapes; fixed/free LJ pair CONTROL.
- `power_followup`: alternative x^p law on the same fixed shapes. It is not
  multiplied by the rational law. Previously inspected excluded data are
  retrospective; a separately declared40-observation set remains excluded.
- `fixed_range_validation`: source, rational fixed-pair and power free-pair
  candidates actually evaluated under independent conditions.
- `joint_refinement`: matched old five-shape vs new six-shape actual fits,
 104 development observations, five exact bulk constraints and sampled
 spectral inequalities. Old optimizer exhausted its20-evaluation budget;
 new optimizer stopped on ftol. All178 probes and stopping metadata saved.
- `joint_validation`: source and both joint candidates,297 static roots,
 independent40-jet observations, three paths, opening extrema, density/site
 direct checks, gradient/Hessian refinement and five dilations×243q×two radii.
- `final_report`: exact parameter replay, full equality-tangent sensitivity,
 corrected-unit residual tables, validation summary and actual static curves.
- `reproducibility.json`: separately rerun fixed-shape power fit equality.

`completion.json` and completed `calibration.json` are final. A
`checkpoint.json` or `progress.json` marked completed=false is an intentionally
preserved interim snapshot, not a contradictory final termination record.
Local `.gitattributes` keeps raw JSON/CSV bytes invariant for SHA-256 provenance.

## Read the negative result correctly

The extra shape changes loss1238.6110 ->1227.6349 in this bounded study, but
Haa27.7833 vs source19.6609 eV/L0² and50 MPa normal displacement0.000656211
vs0.000925366 Å still fail calibration. Bulk cohesion and all cubic elastic
constants are enforced targets, not independent material validation. Lower
loss does not prove better force/curvature predictions; both some excluded
errors and condition number worsen. A numerical test pass is not material
adoption or proof of actual yield/fatigue. These static states are not PDE
trajectories or physical-Hz predictions.

The six bulk-elastic CSV metadata rows in `joint_refinement/*/residuals.csv`
were corrected before commit to pair GPa predictions with GPa targets/names/
scales. Fitted vectors, losses and raw calibration JSON were not changed.
`final_report/*_audited_residuals.csv` independently regenerates these tables;
a regression recomputes all normalized residuals from displayed fields.

## Reproduce into NEW output directories

The runners intentionally refuse to overwrite existing result directories.
Use a verified Python3 interpreter and repository-relative paths:

```text
python -m solver_v1.run_coordination_screening --parent results/fcc111_active_interface/normal_response_v17/nested_range --out NEW_RATIONAL
python -m solver_v1.run_coordination_screening --parent results/fcc111_active_interface/normal_response_v17/nested_range --law power --out NEW_POWER
python -m solver_v1.run_coordination_shape_refinement --profile-directory NEW_POWER --out NEW_JOINT --max-nfev 20
python -m solver_v1.validate_coordination_screening --models NEW_JOINT/old_family NEW_JOINT/power_family --out NEW_VALIDATION
python -m solver_v1.report_coordination_screening --joint-directory NEW_JOINT --validation-directory NEW_VALIDATION --out NEW_REPORT
```

The first three commands execute optimization; the last two replay and
validate saved candidates, never silently refit. Different BLAS/compiler
environments can change final floating-point probes; compare residuals and
admissibility as well as hashes. Do not reinterpret an evaluation-budget stop
as optimizer convergence or a finite-q grid as a whole-zone theorem.
