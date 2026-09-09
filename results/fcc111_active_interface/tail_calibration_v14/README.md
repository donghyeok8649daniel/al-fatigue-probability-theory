# Actual Al static recalibration, v14

See `solver_v1/TAIL_CONTROLLED_AL_CALIBRATION.md` for derivation and limits.
All fitted targets are the same provenance-bound matched 0 K Al99 reference;
no fatigue/yield/mobility/area target has been added.

## Which files are results?

- `*/calibration.json`: completed actual deterministic optimization, with all
  profiles and every start's convergence/budget/numerical status. Completion of
  a runner does not mean every start converged or the material is accepted.
- `*/definition.json`: source hash, starts, observations and stage definition.
- `*/checkpoint.json`: intermediate recoverable profiles, NOT a completed fit.
  For completed runs these duplicate the final JSON's full profiles/optimizer
  logs. Their identical snapshots are preserved locally under the ignored
  `.cache/tail_calibration_v14/optimizer_checkpoints/`, not committed twice.
- `baseline`, `symmetry_channel`, `symmetry_exact_bulk`: early incomplete
  numerical development attempts. Their `attempt_status.json` explains why.
- `baseline_verified`: unchanged-family optimization with analytic spectral tails.
- `symmetry_exact_verified`: actual corrected Eg-channel exact-bulk fit.
- `cubic_exact_bulk`: separate density-cubic ablation; not selected.
- `quartic_exact_bulk`: selected **scoped static/bulk** comparison candidate.
- `rational_exact_bulk`: separate rational angular ablation; not selected.
- `quartic_range_continuation`: explicit even-range upper-bound continuation;
  its better training loss but worse held-out error does not justify adoption.

`quartic_exact_bulk/validation_refined` is the latest independent stationary,
finite-q/radius/analytic-tail and step-refined sensitivity audit. The earlier
`validation` is retained. Neither claims an all-wavevector mathematical proof.

`scoped_report` contains actual static MPa load/unload checks, source vs candidate
roots, curvature-refined opening extrema, plots and a portable scoped parameter
file. Its family-comparison table is a snapshot; if a later `final_report`
exists, use that for the completed-run inventory. The numerical states and
curves are the same candidate, not new fits.

The result is **bulk calibrated to the declared source; full interface NOT
accepted**. Do not relabel it either as entirely unattempted calibration or as
an experimental Al yield/fatigue model. The original historical material, core
states, production energy selector, PDE and kinetic JSON remain unchanged.

Units: eV, fixed L0=4.05/sqrt(2) Angstrom; explicit atomic-cell area converts
interface energies to J/m2 and independent tractions to MPa. Statistical A_c
is not used. Physical mobility, seconds and Hz have not been inferred.
