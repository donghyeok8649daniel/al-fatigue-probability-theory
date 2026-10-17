# Main manuscript workspace

The editable main manuscript is [Papeeria: 1st_Paper.tex](https://www.papeeria.com/p/86fa402b20a593c3dbca788ee535b601#/1st_Paper.tex).
The companion bibliography is [1.bib](https://www.papeeria.com/p/86fa402b20a593c3dbca788ee535b601#/1.bib).

This location was designated by the authors on 2026-10-18. This research integration chat collects workstream conclusions, reviews experimental evidence and author replies, and updates the manuscript. The six existing specialist workstreams retain their module and branch boundaries.

## Scope and evidence rules

- The manuscript describes the original nonlinear energy, original initial ensemble and absorbing opening rule. It does not claim calibrated aluminium fatigue life or experimentally identified crack initiation probabilities.
- Distinguish mathematical identities, results conditional on diffusion assumptions, software checks, and specimen evidence.
- Preserve missing physical scales, cell geometry, mobility and observation limits as unresolved; do not infer them from illustrative old draft values.
- Use reviewed workstream results. Receipt of a result is not acceptance. Do not repeat previous scientific calculations.
- Keep the previous Papeeria version before modifying it. Read the current remote text and detect intervening collaborative edits before writing.
- Add verified bibliography metadata and citations where used. A cited paper supports only the associated method or fact; it does not validate this material model.

## Synchronization

Papeeria is the authoring workspace. The authors explicitly approved public GitHub storage of the full manuscript on 2026-10-18. This directory contains the workspace link, a versioned 1st_Paper.tex snapshot and its 1.bib bibliography. Papeeria remains the editable authoring workspace. This is not a live connection or a promise of unattended access. Synchronization uses the authorized browser UI while login and editing remain available. No dedicated Papeeria plugin was found in the plugin search on 2026-10-18.

Research changes stay on research/20261018-* branches. Do not modify or merge main or probability-pde-solver-v1. The integration branch is research/20261018-integration, based on the recorded remote base c8ccff7553cbb2e66b0c37a76ee8a24808478824; other unpublished local research is not implicitly published.

For each substantive manuscript update record the source report, claim status, source-file hashes, remote save and compilation outcome, and the corresponding GitHub commit. Do not publish author-email addresses, private message bodies, or credentials in the public repository. An absence of new reviewed results requires no rewrite or notification.

## Current reviewed manuscript snapshot

- Revision: 20261018-current-model-rewrite-06
- [LaTeX manuscript](1st_Paper.tex), [bibliography](1.bib)
- [Archived source pair](versions/20261018-revision06/1st_Paper.tex)
- Papeeria compilation succeeded: ten PDF pages and five rendered references.
- Changed or reflowed pages 6-10 were visually reviewed. Unchanged page text on pages 1-5 was preserved.
- LaTeX SHA256: d9ba0705aae0e198b882abcffae37f1cbd4da5e53383ecea1eb6a634bdf26318
- Bibliography SHA256: fb9f1b343cf3dae56beefdc8d0ed85d1f69ef12fd9e1165dd4ff0019ba47b2f2
- Compiled PDF SHA256: 0b5d7cff77c654f3daa8c3e8c173fb05340c909807c53c55c7081d60a44b284f
- The raw XeTeX log retains one nonblocking inputenc-ignored warning. No undefined-reference/citation or overfull/underfull diagnostics were observed.

D11's conditional interior stopped-law construction is included after independent B8 mathematical review and core A14 measure/bookkeeping review. With the exact native energy, positive finite constant coefficients, a predeclared finite native horizon, deterministic L1 load, usual Brownian/initial filtration assumptions and an arbitrary open localization whose closure is compact in the energy-input domain, the interior equation has a unique continuous stopped law and a Borel initial-to-stopped-path kernel. This is an explicit standard cutoff/Picard/stopping construction, not a claim of new mathematics in the literature.

The same predeclared initial measure is preserved. Probability-law construction and linear extension to a finite raw initial measure are distinguished: total stopped-path mass is mu(E), and zero-localization mass is mu(E minus U). Outside-U initial mass is retained as constant proof-stopped paths, without conditioning, normalization or new physical absorption. Any initial absorbed mass stays in its separate original ledger. No extra state Jacobian, bin weight, time step or material cell volume is multiplied into the already supplied initial measure.

The conditional theorem does not certify that an actual original reflected/killed law has this restriction. Original background/initial/regulator source bindings and continuation after an earlier original opening remain unresolved. No actual localization region, wall, load, clock or physical scale is selected. Original global nonexplosion, a=0 accessibility, full reflection/corners, moving/table event equivalence and physical crack probability remain outside this result.

The D10 source map, core geometry workstream A13 and independent C implementation review now clarify continuous-background versus numerical-boundary scope. The positive-a energy kernel and unwrapped registry specify coordinate inputs. The finite SG box's outer numerical no-flux treatment does not identify a complete intended continuous background domain, lower-boundary treatment or reflection field/regulator. The inspected fixed-interval reductions and historical alternatives retain their narrower model scope. No numerical cutoff is adopted as a new continuous or physical wall.

This is a reviewed source-scope clarification, not a new existence theorem, proof of model failure or claim that no relevant evidence exists anywhere in the repository. Global well-posedness, continuum initial-law binding and exact/table opening-event agreement remain separate requirements. The conditional interior theorem above does not close the actual original-source gate; no localization exit is labelled physical opening and no actual domain/initial law is selected.

The independently reviewed E4/B7 symbolic nonvacuity result extends the retained E3 finite root-tube conditions. An exact identity shows that a positive curvature margin makes the lower force envelope positive. For fixed positive native parameters, declare an admissible finite R, then a positive force band below its cap, and then a finite A depending on that band. This establishes a nonempty symbolic family. It does not certify an arbitrary preselected operating band or one finite A for all bands approaching zero force. No actual endpoints, material loads, selector, left well or numerical certificate are supplied. The probability domain, initial ensemble and opening event are unchanged.

The core units workstream's A11/A12 source review and independent D9 probability-meaning review distinguish shared-clock compatibility from absolute clock identification. The family t0(c)=c*t0, Mphys(c)=Mphys/c preserves the same declared native law but represents different uncalibrated physical embeddings. Its physical-time load placement changes; no same-second first-opening CDF identity is inferred. Independently supplied energy/coordinate scales and matched absolute physical mobility can fix a candidate clock, with full-matrix compatibility still required. Agreement between two axis clocks without those independent inputs does not identify seconds.

Same-test measured frequency and an independently supplied and verified native full period, matched load waveform, cycle count, phase and initialization/time origins can specify a protocol-map candidate. That candidate does not independently validate the internal mobility or generator timing. Thermal/static anchors alone do not fix a kinetic clock, and a single random observed onset is not deterministically matched to a native path. Actual physical inputs, clock and observed-crack correspondence remain unresolved.

Revision03's A9 exact work-ownership correction and E3 conditional finite-tube boundary-location result remain intact. Work views are not added again to already counted energy-gradient work. Revision02's distinctions among physical coherent count, per-contribution normalization and passive unit changes remain, including the original r=1 reference and probability/cycle measures. Actual coordinate maps, scales, face geometry and constitutive closure remain unspecified. SB5 remains a separate declarative metadata result; actual evidence and an implemented spatial FVM validator are not supplied by its documentation review.

This is a working theoretical manuscript. Identities, reduced-diffusion assumptions, conditional theorems, software checks and experimental inputs remain distinct. The five verified bibliography entries are unchanged. New editorial source checks and rendering concern only this revision; no completed scientific calculation or prior approved-page QA was replayed.
