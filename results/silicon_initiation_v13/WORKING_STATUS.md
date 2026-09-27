# Si first-crack initiation v13 - six-hour continuation

Start: 2026-09-28 04:19:39 KST / 2026-09-27 19:19:39 UTC.
Deadline: 2026-09-28 10:19:39 KST / 2026-09-28 01:19:39 UTC.
Baseline silicon-wafer-research c3a93f5d2e4df6e21a68d6080ed99cf6aee840bc.
Fresh fetch succeeded, local/origin matched, Si tree initially clean.

The user requested the governing equation explanation and six hours of continued
research. Begin with safely reusing the frozen 539/648 loading8 Hessian rows;
check the input/model/basis/gradient and three recomputed rows before adding 109.
Preserve v12 raw files. New results belong here. Matrix completion is distinct
from the six independent force/energy direction checks and material validation.

Then evaluate return8 and loading10 Hessians, compare grip ensembles, examine
large oxygen-containing structures and the prepared MPA comparison if time allows.
All stages are serialized. The authoritative running queue is in the ORIGINAL
worktree .cache/si-initiation-v13/queue_status.json. Do not restart old v12 queues.

Before a new stage inspect its output and queue processes to prevent duplicates.
Stop new heavy work in the final hour. Preserve actual completion/partial/failure,
independently replay results, bind source/raw hashes, record tests and Git status.
Final report must explain the current model and results to a new reader, not
merely list version history. First-crack A/B, finite-temperature F, mobility,
physical clock, true wafer initiation probabilities/lifetimes remain unvalidated.
No new crack seed or GPa-to-MPa scaling. No additional agents. Al/UI unchanged.

## Current state
08:33KST scientific outcomes finalized; no research process remains active.
Large-oxide retry ended23:25:06UTC at the unchanged PRE-ALLOCATION12GB guard,
exit1 /0predictions /0modelcalls. Idle admission12.53GB did not imply the same
available memory after imports/input preparation. Both failed directories are
preserved. No further retry or lower guard. The prepared finite precheck and
oxide raw replay did NOT execute because model allocation never began.
execution_accounting.json: final newforward1756/newAD2060; newDFT0/MD0.
execution_timeline.json preserves actual timings, loading8 unknown exit status
and both memory failures. COMPLETED_SUMMARY.md records results and limitations.
The portable report now supports --preview-outcomes; actual3outcome pages are
in report_inputs/final_sections.json.30page draft generated, pages24-26 visually
reviewed; final all-page QA follows actual source/manifest/Git verification.
Source/artifact freeze inputs are prepared. Manifest and Git/PDF verification
receipts are generated separately; do not infer their success from this note.
All former running-state paragraphs below are chronological records.

08:24KST zero-Hessian and independent replay COMPLETE, exited23:23:04UTC.
649newADrows +49newforward, six force/energy directions checked, exit0 both.
Minimum0.0150298619131eV/A^2, negative0; residual gradientmax2.44831e-6eV/A.
Independent spectrum error6.75016e-14eV/A^2; fineforceFD5.33338e-9eV/A^2.
Finite-prism current-length tangent63.58523GPa is NOT bulk Young modulus.
ZERO_LOAD_LOCAL_STABILITY_V13.md explains scope and unresolved initiation gates.
All former atomic controllers exited. Completed subtotal before oxide retry:
newforward1756 /newAD2060. No newDFT/MD; not physical time/material adoption.

Idle memory recovered to12.50GB, so one bounded larger-oxide retry started
23:24:51UTC with unchanged12GB admission/3GB runtime guards and old v12 runner.
Controller65288 (Windows launcher54992), child43756; ORIGINAL private
run_large_oxide_retry.py and large_oxide_retry_status.json are authoritative.
Scientific budget1200s/controller1500s, at least60min final-review reserve.
New folder large_oxide_321_1200_retry; old failed69frame attempt preserved.
If complete, finite-array/metric precheck then existing independent source replay
run automatically. Do not duplicate, edit running sources, or freeze manifests yet.
No retry outcome is assumed; read its actual result/status first.

08:10KST completed-stage accounting audited with new
audit_execution_accounting_v13.py. Draft in ORIGINAL cache accounting_draft.json:
1707 new forward evaluations /1411 new AD rows, matching the independent manual
sum. Running zero-Hessian counts remain null/excluded, not zero or planned totals.
Reused100MPa result,539old rows,1944grip rows and199MPA frames are separate.
No new potential or pytest run. After all atomic outcomes are settled, run
--final to save execution_accounting.json; REPRODUCE.md contains the command.
This is an execution ledger, not physical validation or final manifest.

08:04KST report portability completed. Use the VERSIONED generator
results/silicon_wafer_feasibility/build_current_model_report_v13.py, with explicit
--output/--font/--bold-font. report_inputs/baseline contains the original
explanation JSON and two figures; no private-cache input is required.
All28 extracted page texts AND page content streams match the previous draft
at the same display timestamp; portable page20 also visually reviewed.
report_portability_validation.json binds those sources and PDF hashes.
This is draft equivalence, not final all-page review or physics validation.
Final-addendum sections MUST equal the versioned report_inputs/final_sections.json.
That file is not written yet; finalize actual zero-Hessian outcome and accounting
first. Final source/artifact manifests, Git commit and final PDF are still pending.
REPRODUCE.md has the portable command. Retain the old cache builder as history,
but use the versioned generator for final delivery.
Zero-Hessian was369/649rows at08:04KST: compute75720 CPU3772.7s, RSS7.256GB,
available memory4.96GB. Actual process identity checked with psutil after the
normal-permission CIM query was denied. No duplicate atomic jobs or oxide retry.

07:35KST PDF draft28pages: added zero-load return(p18), nonlinear probe scope(p20),
fullQE outcomes(p23); all three rendered and visually reviewed. No empty pages,
replacement glyphs or missing fonts. report_draft_qa_28pages.json pins the draft.
Final all-page review remains pending. Builder automatically includes these
completed stages: final-addendum should add final zero-Hessian outcome and
actual accounting/failed/partial states, not duplicate existing result pages.
Zero-Hessian compute75720 is progressing, source energy replay0 and maximum
force difference2.22e-16eV/A. Only zero_hessian_status.json is now active.

07:32KST: fullQE COMPLETED and independent replay passed22:29:11UTC;
controller74344 exited.1294records(CP2K135+QE1159),199exactly reused,
1095new predictions+9controls=1104new calls. No new baseline/DFT/MD calls.
QE SiO2 RMSE0.158227->0.106350 over493frames, but196improved/297worsened.
QE pureSi1.273789->0.315089 over471frames,378improved/79worsened/14ties.
MPA remains a candidate; selected-subset and full-corpus element scores differ.

NEW zero-force returned-state649Hessian STARTED22:30:54UTC, controller56240,
wrapper64292. Status ORIGINAL .cache/si-initiation-v13/zero_hessian_status.json;
launcher run_zero_hessian_followup.py. Source return_zero_force/state_000/raw.npz,
ensembleforce/stress0, max5400s, deadline01:19:39UTC. At admission90min budget
and60min final-review reserve remained. Completed source relaxation is not a
curvature certificate; this job fills that gap. Do not duplicate it.
Named-state audit_dense_hessians_v13.py preserves the existing v12 numerical
replay, adds explicit --states selection and runs automatically after completion.
Syntax/CLI checked; actual new audit is pending. No additional pytest claimed.
Original queue and fullQE queue remain finished. Source manifest/report pending.

07:28KST: primary queue FINISHED22:18:07UTC; controller68516 exited.
Zero axial return converged:147new calls,617.234s, axial residual2.652e-6eV/A,
free force max3.330e-6eV/A. Independent replay passed. Same external0force,
different length: +2.23717257A / +7.43611% versus prior original zero-force gauge,
energy +8.15092314eV. Static history dependence only, not plasticity/crack proof.
Nonlinear72points+2source replays complete; NEW independent72point/36pair raw
replay passed, anharmonic_replay.json. New scope analysis separates transverse
force residuals: small line energy errors do not certify full basin harmonicity.
ZERO_LOAD_AND_NONLINEAR_SCOPE_V13.md records values/limits; plot visually checked.

Full-QE MPA follow-up STARTED22:23:21UTC, controller74344/wrapper41952.
Status ORIGINAL .cache/si-initiation-v13/full_qe_status.json, launcher
run_full_qe_followup.py. It runs complement1095 +9controls then independent
replay serially. Original199reused, no duplicate baseline calls. No other atomic
jobs should run concurrently. Primary queue remains finished, not restarted.
After full-QE finishes, inspect outcomes and remaining time before any new work.
Large-oxide retry not started: measured available memory11.63GB while idle,
still below existing12GB guard. Preserve failed folder and do not loosen guard.

07:11KST: the first-reader PDF draft now has25pages. New force-ensemble(p15),
local harmonic(p18), selected-MPA(p20) pages were rendered and visually reviewed.
Clarified overall free-atom RMS versus maximum per-atom RMS in p18 and reviewed
again. No empty pages, replacement characters or missing font glyphs. This is
NOT final full-page QA or delivery. Original-worktree
report_draft_qa_25pages.json pins this draft. The builder now automatically
includes these completed results; do not duplicate them in final-addendum.
Zero-force relaxation remains running; at100calls the axial residual was
2.22e-4eV/A, above the1e-5 success criterion. Do not call that checkpoint converged.

07:05KST follow-up: MPA199 and independent mpa_replay are COMPLETE, exit0.
199 selected frames +9 standard-calculator controls =208 actual new calls;
baseline reused,0 new baseline calls. Independent raw replay source arrays
match exactly. Group/element results below are descriptive, not adoption.
return_zero_force started22:02:07UTC, wrapper45716; nonlinear probes still queued.
Do not launch full-QE completion or retry large oxide concurrently.

Latest update 07:03KST: loading10 and all queued dense postprocessing completed.
grip_augmented completed and its NEW independent replay passed. Large oxide
failed BEFORE model allocation because available memory was below12GB; zero
frames evaluated. Failure folder and partial.json are preserved. MPA199 is now
running (wrapper41860/compute17668/controller68516); zero-force and nonlinear
stages remain queued. Do not repeat/interrupt them. Detailed new evidence below.
Only reconsider the failed oxide job after all atomic jobs exit, with at least
12GB actually available and sufficient time, using a NEW output directory.
Do not lower the memory guard merely to force execution.

Resume implementation tested: 17PASS/0.45s. Initial pytest setup failed on the
existing system TEMP permissions; a fresh workspace test directory passed.
Actual loading8 restart: old539-row hash/context/basis/gradient verified;
rows0/269/538 recomputed with exactly zero norm difference. Remaining109 completed,
plus all six independent force/energy directions. Minimum curvature
0.02947002388292215 eV/A^2, negative eigenvalues0. New forward49, AD112
(three validation rows +109 new rows); reused539. Numerical elapsed737.8675s
is not CPU time. Residual gradient0.000722eV/A remains; positive curvature
does not certify an exact stationary state, global minimum or first crack.
Queue controller v2 corrects zero-force replay inputs to state_000. The actual
atomic process was retained during controller replacement; no new duplicate run.
At04:41KST the finished adopted child had disappeared, but psutil supplied no
exit code so v2 stayed waiting. Verified child absence, complete outputs and
six direction checks, then replaced only that controller with v3. Original
exit code remains unknown. No atomic replay was repeated; return8 started
19:41:52UTC. v3 also stops both its Windows launcher and calculation children
on deadline. Private controller_recovery.json records the identities and hashes.
Authoritative queue now uses .cache/si-initiation-v13/queue_research_v3.py in
the original workspace. Do not restart v1/v2 or their recovery helper.
force100 copied byte-for-byte from v12 with reuse_provenance.json. Old summary
call counts are historical; v13 newforce100 calls/rows are zero.
Do not aggregate old copied counters as new work. Source/test binding is in
resume_implementation_validation.json; session_plan.json records scope/deadline.
Heartbeat si-7 ACTIVE every15 minutes until finalization, deadline above.
Old 100MPa649 Hessian complete; no duplicate calculation planned.

## Independent nonlinear replay prepared (04:55KST)
Added results/silicon_wafer_feasibility/replay_anharmonic_v13.py. It independently
reconstructs energy/force residuals from the raw arrays, including the existing
linear residual gradient, signed even/odd terms and the complete declared
displacement schedule. Source hashes, partial status and actual call counts
are checked. No model calls are made by this postprocessor.
Manufactured six-dimensional cubic/quartic data plus deliberate corruption,
missing final pairs, nonfinite values and double-counted calls: 10PASS/2.18s.
This is implementation validation, not completed Si nonlinear evaluation.
The scientific queue remains unchanged and return8 is running. Once the queued
anharmonic_probes stage has exited, run this NEW audit exactly once:

    python results/silicon_wafer_feasibility/replay_anharmonic_v13.py --results results/silicon_initiation_v13 --probes results/silicon_initiation_v13/anharmonic_probes --output results/silicon_initiation_v13/anharmonic_replay.json

Use the MACE research environment (NumPy/SciPy); no model allocation occurs.
Do not run it against actively written raw files. An actual partial result
remains partial. Tests/source binding: anharmonic_replay_validation.json.

## Geometry-only path initialization audit completed (05:09KST)
chord_geometry: all64,620 pairs between loading8 and return8, same grips exactly.
Free-atom labelled RMS1.054057A/max2.884706A. The straight chord approaches
1.710907414A at fraction0.493424705 (free pair253/286); this pair is2.329663/
2.358467A at its endpoints. Twenty independent scalar minimizations and
reverse/rigid-frame checks passed. Plot visually reviewed. New potential0/DFT0/MD0.
This warns about path initialization; no energy path, saddle, barrier, B state
or clock was found. CURRENT_MODEL_AND_PATH_GATE_V13.md records the governing
equation, scope, actual evidence and required next gates for a new reader.

## Sampling scope and full-QE follow-up prepared (05:28KST)
mpa_sampling_scope completed: exact199 selected sources pinned (CP2K135/QE64),
source-selection implementations agree. Baseline full/selected RMSE:
CP2K0.372760/0.518244 and QE0.691243/0.351422eV/A; QE compositions30/20.
This is CSV reanalysis, no new potential call. Official models share MPTrj;
their agreement is not an independent uncertainty bound. See MODEL_COMPARISON_SCOPE_V13.md.

NEW full-QE completion is PREPARED ONLY, not in the automatic v3 queue.
After the entire existing queue finishes, if at least3600s remain until01:19:39UTC,
and mpa_comparison plus mpa_replay are complete/pass, run max1200s:

    python results/silicon_wafer_feasibility/complete_mpa_qe_v13.py --cp2k results/silicon_initiation_v11/oxidized_surfaces_max320 --qe results/silicon_initiation_v11/oxide_mace_mtpu_1159_v2 --model LOCAL_MPA_MODEL --reuse results/silicon_initiation_v13/mpa_comparison --output results/silicon_initiation_v13/mpa_qe_full --max-seconds 1200 --deadline-utc 2026-09-28T01:19:39Z

LOCAL_MPA_MODEL is the original workspace .cache/si-initiation-v12/mace-mpa-0-medium.model.
Use the existing MACE environment/hidden process, and verify no other atomic job runs.
If complete, independently replay:

    python results/silicon_wafer_feasibility/replay_mpa_qe_v13.py --cp2k results/silicon_initiation_v11/oxidized_surfaces_max320 --qe results/silicon_initiation_v11/oxide_mace_mtpu_1159_v2 --reuse results/silicon_initiation_v13/mpa_comparison --result results/silicon_initiation_v13/mpa_qe_full --output results/silicon_initiation_v13/mpa_qe_full_replay

Expected:1294records, original199exactly reused, newQE1095plus9standard controls.
Do NOT report planned counts as actual. Syntax/CLI checks passed; candidate
inference and full raw replay remain pending. Preserve time for final review.

At05:39KST reproduced max(known_error, NaN) hiding a nonfinite comparison.
The new full-QE replay now rejects nonfinite scalar data on either side,
checks reported atom counts, and requires the exact nine new-frame controls.
Manufactured corruption/control tests:11PASS/0.30s, source-bound in
mpa_qe_replay_guard_validation.json. The preparation record keeps the historical
pre-edit validator hash; no physical result was evaluated by that old version.
Existing v12 code/results are unchanged. This is implementation validation only.

## Return8 completed; loading10 running (05:56KST)
return8 completed all648 Hessian rows and six independent force/energy directions.
Process exit0, numerical elapsed4090.7160592s (not CPU time). Newforward49,
newAD648; minimum curvature0.022821591769603626eV/A^2; negative eigenvalues0.
Maximum reduced gradient0.0007035344eV/A remains. Positive curvature of this
approximate stationary configuration is not an exact minimum/global stability
or first-crack certificate. The lowest curvature is smaller than loading8's,
but mode identity need not match; do not interpret it as a bulk modulus ratio.
The two completed8% states total98newforward/760newAD (including three resume
validation rows). This excludes loading10 in progress and reused force100.
loading10 started20:50:22UTC; wrapper50544/compute71224, controller68516.
05:56KST process identities and row65 progress confirmed. No duplicated jobs.
Matrix replay, localization and harmonic comparison remain queued after loading10.

## First-reader PDF draft prepared (06:04KST)
PDF skill create marker succeeded once; do not repeat it. Builder and authoring
status live in ORIGINAL .cache/si-initiation-v13/build_current_report.py and
report_authoring_status.json. Draft: ORIGINAL tmp/pdfs/si_v13_draft.pdf,21pages.
Not a final deliverable. Prior source hashes were checked before reusing the
reviewed model explanation. New completed8% curvature, chord geometry and
sampling-scope evidence are distinguished from pending atomic jobs.
First draft rendered all21pages; no empty pages or missing font glyphs.
Pages7/14/15/16 visually checked; after equation subscript formatting, page7
rerendered and checked. Final content still needs ALL-page visual review.
Final invocation requires --final-addendum PATH: JSON flags
reviewed_actual_outcomes/git_status_verified/manifest_verified must be true,
revision is the actual40-character source commit, and sections contains actual
follow-up outcomes, failures/partial results and new-versus-reused accounting.
No final addendum exists yet. Do not mistake prepared pages for completed science.

## Preserved inputs verified; final provenance tooling prepared (06:13KST)
Independently checked2,517 existing v12/v11/v9 source/artifact records against
working files AND c3a93f5 Git blobs. All matched. Exact receipt:
upstream_bundle_replay.json; recorded old63 tests were NOT rerun.
consolidate_validation_v13.py rechecked seven current source hashes and merged
the three already executed, disjoint17+10+11 tests into validation.json.
Recorded38 tests is not an extra38-test run; consolidation adds0test/model calls.
record_initiation_manifest_v13.py is prepared and CLI/import checked. It pins
new source/results, unchanged v12 artifacts and their original upstream chain.
Do not run --final until ALL scientific and record writers have stopped.
Then existing verify_initiation_bundle_v12.py accepts --results
results/silicon_initiation_v13 for working-file, INDEX and final-commit checks.
REPRODUCE.md records the exact workflow. No final manifests exist yet.
10% progressed to241/648 rows at06:13KST. Existing tracked Si files unchanged;
new v13 sources/results remain untracked. No old source or result was overwritten.

## Same-grip mode comparison completed (06:29KST)
same_grip_modes compares the two COMPLETE8% matrices only; no new model calls.
Lowest-direction squared overlap0.8835683142540757, acute angle19.95127385deg.
Holding the INITIAL lowest Cartesian direction fixed gives curvature
0.02947002388(initial) versus0.09774827402(returned)eV/A^2. Thus comparing
different minima alone does not establish uniform softening or a bulk modulus.
The endpoint chord has1.054057A/free-atom RMS. Initial-state quadratic
extrapolation predicts+222.234824712eV versus actual endpoint delta-9.719222938eV.
These are total finite-prism energies, not per-atom values or transition barriers.
First16 low-curvature modes capture42.5228%/27.5470% of the endpoint chord's
squared Cartesian displacement. This is geometric projection, not probability
or a physical reaction path. All first1/4/8/16 boundary gaps are retained.
replay_same_grip_modes_v13.py passed73 scalar/array comparisons using full
Cartesian projectors, SciPy principal angles and spectral Rayleigh/energy sums.
Maximum errors by units: angles3.5942e-8degree; energies1.9895e-13eV;
curvatures3.7192e-15eV/A^2; overlaps1.1102e-15. The replay's global maximum
mixes units and must not be reported as a single physical error. No new pytest.
Plot visually reviewed; current PDF draft22pages, added page15 reviewed.
Sources/raw/spectrum/protocol hashes pinned. Queued localization/harmonic
analyses remain separate; this did not duplicate them or edit any running runner.
