# v15 reproduction and adoption status

Run from this repository's root. `python` below means a verified Python3
interpreter with the repository dependencies, not an assumed launcher alias.
The executed Windows environment used Python3.13; no new global dependency
was installed. Use fresh output directories: the runners refuse overwriting
existing studies. No command below changes production defaults or launches a
physical-Hz fatigue simulation.

## 1. Actual final static calibrations

The following is an **optimization**, not a replay of a stored best vector.
The two deterministic starts, bounds, target source and independent stability
counterexamples are committed. The same source Al99 file is checksum-bound;
`source_and_targets()` obtains the target-only reference through the existing
reference loader. A first run may need network access for that source file.

```powershell
python -m solver_v1.run_tail_constrained_calibration --symmetry-channel --quartic-angular --exact-bulk --interface-development --starts results/fcc111_active_interface/interface_development_v15/strain_refinement_starts.json --additional-wavepoints results/fcc111_active_interface/interface_development_v15/strain_counterexamples.json --wavepoint-step .1 --stability-stretches .99 1.0037037037037038 1.01 --max-nfev 24 --out results/v15_replay_quartic
python -m solver_v1.validate_tail_calibration results/v15_replay_quartic --out results/v15_replay_quartic/validation_exact_tangent --grid-step .08
```

For the cross-invariant run add `--density-angular-cross`; for the final
cross-plus-Bloch run additionally add `--static-bloch-targets`, use18 maximum
function evaluations per start, and a different output directory. All have
the default three range bounds [.7,2,2] to [8,12,12], radius12, two-point
shape differences2e-4 and the SAME primal/KKT thresholds. Optimizer floating
point paths can depend on BLAS/versions; validate the returned vector and
each independent residual instead of requiring identical iteration counts.

Recorded completed studies:21 optimizer runs /1770 successful coefficient
profiles in `interface_development_v15/release_report/`, plus failed numerical
profiles and independent checks. This does not count every profile as an
independent material or prove a global minimum. The final stable quartic,
cross and cross-plus-Bloch heldout interface normalized RMS values are
78.29824,71.20608,78.90034 respectively: **none is adopted**.

The21st run is an actual independent reoptimization of the two-start
strain-stable cross case on the frozen final code. It executed112 profiles
in532.31s. Decays, coefficients, all predictions/residuals and the loss
1128.2343520341506 were exactly reproduced on this environment; timing fields
naturally differ. Independent243-q/radius12/16 and stationary/curve checks
were repeated. This is reproducibility of an inadequate material, not its
adoption. `validation_manifest.json` records exact software/test versions.

The existing report command regenerates summaries, not fits:

```powershell
python -m solver_v1.report_interface_development results/fcc111_active_interface/interface_development_v15 --out results/v15_replayed_report
```

`heldout_curvature_term_budget.csv` is an additional numerical diagnosis.
Reproduce its ungrouped columns with the returned `cache_type` from
`validate_tail_calibration.load_material`: form `cache_type(observations).
matrix(best['decays'])`, select `new_direct_039_Haa`, and multiply that row
elementwise by `best['coefficients']`. Observations come from
`development_observations(source, old_observations, states)`. Group the
explicit coefficient order u,v,A,B,C,D3,D1,D2,D_E,K3[,E_xI]. No finite
difference or additional optimization enters this budget.

The preceding20-run `final_report/` remains a historical report snapshot and
contains that term-budget CSV; `release_report/` includes the independent
21st optimization. Neither directory is silently overwritten by the runner.

## 2. Offline replay of the actual low-frequency MD audit

The full multi-gigabyte trajectory is not committed. A complete1024frame
plane-coordinate projection and its source metadata ARE saved. Copy them
into a new ignored input directory without overwriting any prior data:

```powershell
New-Item -ItemType Directory -Path .cache/v15_offline_input -ErrorAction Stop
Copy-Item -LiteralPath results/public_aluminum_validation_v15/evidence_report/projected_coordinates.npz -Destination .cache/v15_offline_input/plane_coordinates.npz -ErrorAction Stop
Copy-Item -LiteralPath results/public_aluminum_validation_v15/kinetic_plane_stationarity/source_metadata.json -Destination .cache/v15_offline_input/source_metadata.json -ErrorAction Stop
python -m solver_v1.run_zero_frequency_audit .cache/v15_offline_input --out results/v15_zero_frequency_replay
```

Projection SHA256 is
`2037b0104c33573abc27532517610ff84f308e6ba55ded08a0b0b2bf834a7eb6`.
All three finite-integral plans must remain in the report. The result is
below block/quadrature sensitivity, **not a calibrated mobility**. Three
failed longer-download attempts are separately recorded and are not used.

## 3. Direction-specific static and finite-source scenarios

```powershell
python -m solver_v1.run_tensor_traction_audit --models results/fcc111_active_interface/tail_calibration_v14/quartic_exact_bulk results/fcc111_active_interface/interface_development_v15/strain_stable_cross --out results/v15_tensor_replay
```

This actually solves297 stationary states for the source and two candidates.
All input tensor entries, local traction components, residuals and model
hashes are saved. It is not a time-dependent PDE or residual-plasticity test.

`finite_source_refined/` contains54 explicitly hypothetical source geometries/
tractions with angular, quadrature and stress-derivative refinement. The
outer-line model is derived in `FINITE_SOURCE_BARRIER_DERIVATION.md`; its
MPa critical stresses are not fitted experimental yield stresses, and its
finite energy is not converted to a rate without a matched core/kinetic model.

## 4. What was deliberately NOT changed

- LJ pair form, infinite Poisson/Bessel energy and per-atom environment sums.
- Previously stored reduced/full-FCC material parameter sets and references.
- Production `P(a,s)`/SG/opening and registry bookkeeping.
- Production scalar stress bridge, chi, mobilities, kT and energy selector.
- Uncalibrated physical-time JSON: no physical PDE seconds or Hz.
- Specimen A_c, independent-region certification and mesh/UI adoption gates.

The production UI still has no independent shear/direction input. Research
tensor projection and vector static calculations are not that feature. The
appropriate next scientific gate is accurate, stable interface/core forces
and a physically matched finite activated state, not a rescaled strain plot.
