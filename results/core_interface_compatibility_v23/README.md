# v23 — source-core / interface compatibility audit

Starting clean branch/fresh remote: `d1ac6695518d3738b8303687d8c6a7b58877203e`.
This is static research. No material, production PDE, kinetic or UI adoption.
The complete derivation and interpretation are in
`solver_v1/CORE_INTERFACE_COMPATIBILITY_V23.md`.

## What was actually executed

| Directory | What it establishes | What it does NOT establish |
|---|---|---|
| fixed_shape_audit |48 verified coefficient minimax LPs, saved full matrices/duals|Global nonlinear-family impossibility|
| box_profiles |4 same-objective spectral QPs with explicit original-scale boxes|Material adequacy merely from KKT|
| radial_minimax_corrected |160 bounded existing-shape profiles, subset/full replay|A converged positive material optimum|
| radial_full_validation |Full285 observations, sampled bulk tails, two source-core force configurations|A stable relaxed core or calibrated yield|
| rational_quartic_ablation |5 explicit prior rational-shape ablations; no needed improvement|A new adopted energy or supported core law|
| static_topology_corrected |Actual9 stationary states,3 saddle connections,15 mixed-traction states,33/65 opening brackets|Physical-time dynamics or fatigue|
| box_core_* / radial_joint_* |Actual nonlinear atomic attempts with individual summary/failure records|Automatic force/Morse/domain or material acceptance|
| core_resolution / core_replay |Frozen-force rings5/7/9, actual R8/ring5->7 and R8->10 re-relaxation, same-energy replay|Infinite-domain core convergence or an Al-shaped core|
| bulk_dispersion_audited_partition |66 static finite-q matrices, explicit scales, radius12/16 and tails|Phonon Hz, blind held-out certification or nonlinear core causation|
| shortwave_direct_energy / shortwave_direct_energy_refined |60 independent phase-class energy curvatures and amplitude refinement|A new finite-cutoff canonical energy or successful Al calibration|
| validation_ledger |Executed test/smoke summaries and byte hashes when finalized|Physical truth certified by a hash|

## Important correction

Positive Hxx at the SOURCE saddle coordinate does not imply no CANDIDATE
saddle. Full coupled Hessian, actual root location and downhill connectivity
are required. The v22 wider and v23 radial-joint surfaces both have verified
index-one interface saddles, including positive own Hxx. Their fault/saddle/
opening energies still differ substantially from the matched Al99 benchmark.
No old energy is altered by this correction of interpretation.

## Quantitative tradeoff

- Wider fixed-shape full115-jet minimax error:18.46445 original discrepancy
  units; even unconstrained signs give18.33840. Two offending jets alone are
  simultaneously feasible, so their signs alone do not prove incompatibility.
- New radial shape full115-jet minimax:13.92314. Its selected29-jet best
  vector has error9.98358 on selected data but40.23339 on the full data.
- Reprofiling the SAME core/interface objective at this new shape reduces
  source-core RMS to.0995689eV/L0(R8),.0969927(R16 excluded configuration),
  but raises other-interface squared loss to4520.30(v22 wider2276.47).
- Source / wider / radial-joint fault energies:
  .150479 / .079457 / .075309J/m².
  Saddle energies:.172002 / .119750 / .114757J/m².
  W(40h):1.741285 /1.840583 /2.524404J/m².
- These are matched0K rigid-interface energies, not independent experimental
  measurements at every geometry. The source is target-only Al99, with its
  original published interpolation/cutoff; LJ/Bessel remains the candidate.

All units, original target scales, seven equality anchors and declared signs
are preserved. The rational-quartic moment scale is not an area or length.
No A_c, desired strength or fatigue observable enters these fits.

Actual stable radial-joint cores exist at R8/ring5, R8/ring7 and R10/ring7.
Their force residuals are4.34e-10,3.30e-8,2.77e-12eV/L0; full fixed-boundary
Hessian minima are.172654,.172236,.130869eV/L0². The R2 inner field changes
2.19e-5L0 under neighborhood refinement and1.99e-3L0 under free-domain
enlargement. Sampled registry spans remain~.881-.886L0 vs source2.760L0;
this is not an infinite-domain certificate or a physical partial separation.

The new short-wave audit finds a separate material defect. At cubic
q=(.5,.5,0), a vector polarization has curvature375.872eV/L0² vs source
60.0863. Direct per-atom phase-class energy changes approach these values
under amplitude refinement. Small-q elasticity alone does not constrain
this response. No bulk matrix, source data or candidate force was clipped.

## Failure records and reproduction

Fresh directories are mandatory. Failures are not overwritten by retries:
`radial_minimax`, `static_topology`, initial core seed preflights, core budget
stops and the rejected reciprocal-trial descent remain explicit records.
The interrupted pre-final regression is recorded in `test_run_history.json`;
only the separate completed final XML is eligible for the final ledger.

```text
python -m solver_v1.run_core_interface_compatibility --out results/NEW_fixed_audit
python -m solver_v1.run_core_interface_box_profile --out results/NEW_box_profiles
python -m solver_v1.run_jet_shape_compatibility --out results/NEW_shape_search --max-fev 160 --max-seconds 1200
python -m solver_v1.validate_jet_shape_compatibility --search results/NEW_shape_search --out results/NEW_full_validation
python -m solver_v1.run_quartic_compatibility_ablation --out results/NEW_quartic_ablation
python -m solver_v1.run_v23_static_topology --out results/NEW_topology
```

The last command explicitly replays the documented stored candidate snapshots;
it does not refit them. `run_current_material_core` binds each atomic run to
its own parameter/state hash and records all actual solver options in metadata.
Final artifact bytes use local `.gitattributes` to avoid silent line-ending
changes. Do not rewrite scientific outputs after generating a final ledger.

Physical a/s mobilities, t0, actual seconds/Hz, finite-loop activation, actual
yield/fatigue and specimen correlation calibration remain unavailable.
