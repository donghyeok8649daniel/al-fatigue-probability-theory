# v22 — current nonlinear LJ/Bessel screw-core research

All materials here are **unadopted static research candidates**. `source_reference`
is a target-only Mishin Al99 comparator; it is not our production potential.
The unchanged starting material is bound by SHA256
`8736cb9d28f1430e3991b1046ef42544cefaa3c9a4f5d743e0d02931329f5c2e`.
No fatigue, experimental yield, finite-loop activation or production a/s
seconds/Hz calibration is established by this directory.

Read the derivation and final interpretation in
[CURRENT_MATERIAL_CORE_BRIDGE_V22.md](../../solver_v1/CURRENT_MATERIAL_CORE_BRIDGE_V22.md).

## Key outcome, not a material-adoption certificate

| Diagnostic | Original candidate | Wider joint trial |
|---|---:|---:|
| Same-source-core force RMS [eV/L0] | .28883360 | .15241668 |
| Other-interface normalized squared loss | 2215.203 | 2276.471 |
| R8/r5 reconstruction barrier [eV/straight-row repeat] | .0148345888 | .0097785373 |
| R10 sampled25–75% disregistry span [L0] | .87408147 | .87783192 |

Force mismatch drops47.23% and the verified reconstruction barrier34.08%, but
the relaxed core remains much sharper than the source (R10 span2.79084554).
Finite-opening/registry curvature sign errors remain. The wider fit also stops
at its outer evaluation budget, with an active saturation bound. All56 static
stress states in the final comparison were actually computed; the analytic
trials return to their initial cores after the tested signed50MPa protocol,
whereas the source has boundary-sensitive metastable translated states.
This is NOT experimental yield, fatigue, a thermal hold or physical-Hz data.

## Evidence routing

- `representation_validation`: independent direct-atom versus infinite Bessel
  values, gradients, Hessians and nonuniform per-atom assembly.
- `current_replay`, `source_replay`: replay bound state energies; no claim to
  have rerun minimization inside those reports.
- `stable_plus_R*`, `reconstruction_*`: original current candidate's actual
  stable branches, neighborhood/domain refinements and index-one saddles.
- `source_reference`: the same atomic geometry and boundary protocols using
  only the published source target; R8/R10/R16/R20 stress-return calculations,
  independent translated endpoints and explicit saddle descents.
- `force_compatibility`: exact fixed-shape coefficient nullspace lower bound;
  removing sign constraints is a counterfactual, not an accepted material.
- `constrained_force_profile`: sign/sampled-spectrum-constrained force-only
  profile and its worsening of the other interface observables.
- `core_informed_shape_probe`: initial actual bounded joint optimization.
- `core_informed_shape_continuation`: same-box continuation. `xtol` termination
  is distinguished from the larger-than-requested gradient optimality.
- `core_informed_wider_search`: explicitly enlarged original-centered search
  box, unchanged data/objective/energy terms. Evaluation-budget termination
  is recorded; the selected diagnostic snapshot is NOT adopted.
- `existing_screening_profile`: five actual profiles of the already implemented
  density factor x^z, retaining all other wider-search shapes and constraints.
  The zero-exponent profile replays the preceding fit; this is not a new law.
- `core_informed_validation`, `continued_probe_validation`,
  `wider_probe_validation`: complete target residuals, retrospective excluded
  states, local numerical sensitivity and separately bound candidate snapshots.
- `probe_core*`, `probe_stationary*`, `probe_stress*`: actual first-trial
  relaxed cores, independent minimizers and signed static load-return.
- `wider_core*`, `wider_stress*`: actual wider-search candidate calculations,
  including any failed attempts and recovery checks. Read the status, not the
  directory name, before classifying them.
- `core_energy_mechanism`: exact energy/force coefficient decomposition on a
  declared frozen interpolation, NOT a relaxed transition barrier.
- `farfield_load_audit`: actual work-conjugate atomistic stress. The nominal
  50 MPa shear is indeed about 50 MPa; absolute normal tail error is separate.
- `analytic_pair_tail_with_source`: rigorous infinite transverse LJ-pair
  force majorants and independent pair-gradient assembly on own/source states.
  This supersedes the own-state-only `analytic_pair_tail_audit` scope, not its data.
- `frozen_force_refinement`, `wider_frozen_refinement`: actual same-state
  neighborhood/tolerance changes versus material error. These changes are not
  a rigorous bound on all nonlinear environmental tails.
- `study_summary_final`: original/first-trial/source comparison including the
  independently reproduced R20 two-row shift; an intermediate report.
- `final_comparison`: final four-material comparison including the wider trial,
  actual56 stress states, sampled profiles, saddles, independent translations,
  separately labeled ring/domain changes and source-configuration multistarts.
- `wider_replay_final`: all four completed wider-candidate zero-load energies
  independently recalculated. This is replay, not new minimization.
- `validation_ledger`: actual JUnit counts/timing, newly executed desktop
  smoke, each completed/failed case and byte-level artifact inventory.

Some directories named `combined_analysis`, `final_analysis*`, `study_summary`
are earlier intermediate reports. They are preserved history, not the final
case inventory. Likewise `stable_R8r5` under the source is actually a saddle:
its negative Hessian was found and the subsequent descended states have
separate names. An optimizer return message alone never establishes stability.

## Units and what was actually loaded

Coordinates are reduced by L0=2.863782463805517e-10 m. Energies are eV per
straight infinite row's one-atom repeat, gradients eV/L0, Hessians eV/L0^2.
The prescribed far-field loads are **resolved shear sigma_xz**, not nominal
uniaxial stress. The standard signed protocol is
0 -> +5 -> +20 -> +50 -> 0 MPa and independently its negative branch.
Every material uses its own elastic tensor and exterior, without copying the
source elasticity into the analytic candidate.

Load-return here is a sequence of static stationary states, not a time hold,
fatigue cycle or kinetic response. A pre-existing winding-one straight screw
is supplied by the exterior; neither a source population nor nucleation rate
is fitted. Core translations depend on the finite free boundary. No arbitrary
line length is multiplied into an activation energy and no A_c appears.

## Preservation and reproduction

Each atomic calculation has `metadata.json` with the exact parameter binding,
geometry and numerical controls, then `summary.json` for an actually completed
calculation or `failure.json` plus checkpoint for a failed/stopped attempt.
`completed` means executed, not physically accepted; consult the independent
force and Hessian/Morse fields. Explicit recovery has a new directory and
retains the failed source. Exact commands and stages are described in the
derivation and runner modules; output directories must be fresh.

Byte hashes bind observations and parameter snapshots, not the physical truth
of a source. The local `.gitattributes` preserves their bytes across Windows
Git checkout. Do not normalize generated CSV/JSON newlines after hashing.
