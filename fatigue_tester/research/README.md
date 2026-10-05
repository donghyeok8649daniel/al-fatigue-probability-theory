# Published research references for the physical tester

`reference_manifest.json` indexes committed research evidence by full source SHA,
path, SHA-256, byte size and pinned GitHub/raw URL. References are design inputs;
they do not authorize solver execution or make the tester a calibrated research model.

The physical tester remains under `fatigue_tester/`. Ten exact committed
reference documents from the probability and silicon branches are collected
under `published/<source-branch>/<original-path>`. Solver implementations and
desktop applications stay in their own branches. Uncommitted research work,
personal CAD files and local measurement logs are not included.

| Source | Commit inspected | Relevant boundary |
|---|---|---|
| `fatigue-tester` | `399ffb1f9d711439c4fbed395070dfc40050cc3c` | Existing hardware/BOM, generic firmware/HAL and telemetry |
| `probability-pde-solver-v1` | `c8ccff7553cbb2e66b0c37a76ee8a24808478824` | Project/surface-load interfaces, result meaning, uncalibrated model time, Al round-specimen reference |
| `silicon-wafer-research` | `8f9f754c9976b82851a8ad7e2a6a1b365e3bd5c2` | Si-specific material/load gates and historical device/CAD handoff |

These are the inspected publication snapshots, not a claim that all branches
still have these latest heads. A later update must record the newly read source.
`local_repository_path` identifies each collected file or existing tester file.
The SHA-256 and size bind the original source bytes; existing tester source
hashes describe the target base commit, so later local edits need their own
history. The collected files are preserved without rewriting their internal
links. If a link requires an omitted source-tree document, use the manifest's
pinned original GitHub URL.

| Collected reference | Purpose |
|---|---|
| [Project files](published/probability-pde-solver-v1/app/PROJECT_FILES.md) | Project-file versions and exchange boundaries |
| [Surface loads](published/probability-pde-solver-v1/app/SURFACE_LOAD_SETUP.md) | Load units and mesh identity |
| [Physical time](published/probability-pde-solver-v1/solver_v1/PHYSICAL_TIME_AND_MOBILITY.md) | Laboratory and model time distinction |
| [Result fields](published/probability-pde-solver-v1/solver_v1/RESULT_FIELDS.md) | Research result meanings |
| [Kinetics and fatigue validation](published/probability-pde-solver-v1/solver_v1/LINE_KINETICS_AND_FATIGUE_VALIDATION.md) | Validation scope |
| [Al round-specimen reference](published/probability-pde-solver-v1/solver_v1/data/aluminum_fatigue_validation_v18.json) | Published test conditions and attribution |
| [Material selection](published/silicon-wafer-research/app/MATERIAL_SELECTION.md) | Al and Si scope |
| [Si specimen loading](published/silicon-wafer-research/solver_v1/SILICON_SPECIMEN_LOADING_V6.md) | Orientation, cracks and distinct fixtures |
| [Historical CAD/device handoff](published/silicon-wafer-research/results/silicon_device_ensemble_v15/README.md) | Earlier device-constraint evidence |
| [Historical CAD source reference](published/silicon-wafer-research/results/silicon_device_ensemble_v15/CAD_REPOSITORY_REFERENCE.json) | Source and measurement provenance |

The numeric Al reference preserves its CC BY 4.0 attribution, DOI and full-range
convention. Its triangular-wave protocols and failure/crack observations are
not current machine ratings, independent replicate raw data, or calibrated PDE
predictions. Al round grips and Si wafer fixtures/loading methods remain separate.
The reference is the original repository's numerical transcription of
S. Deschanel, W. Ben Rhouma and J. Weiss (2017),
[Acoustic emission multiplets as early warnings of fatigue failure in metallic materials](https://doi.org/10.1038/s41598-017-13226-1),
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
This collection makes no changes to that snapshot; the file retains the
original transcription and range/amplitude conversion declarations.

Machine force/frequency/stroke/stiffness are unresolved. Existing BOM figures are
historical candidates. Physical seconds/Hz do not automatically map to the
production probability clock. Board-specific HAL and real force-loop validation
are incomplete. See `../docs/CAD_RESEARCH_HANDOFF_KO.md` for data handoff details.
