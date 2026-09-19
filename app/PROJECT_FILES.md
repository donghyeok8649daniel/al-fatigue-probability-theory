# Project files and desktop launcher

The top toolbar provides Open project, Save project and results (Ctrl+S), and
Save as. Ctrl+O opens a project. Save is available before analysis; after a run,
it also stores the numerical result. Finish/stop an active calculation first.

## .ftgsim version 2

This is a compressed ZIP container with ftgsim-manifest.json, state.json and
NumPy numerical arrays. Arrays retain shape, dtype, bits, NaNs and infinities;
loading uses allow_pickle=False. Member count/size and SHA-256 checks are
validated before applying state. Writes use a same-directory temporary archive
and atomic replacement, so serialization failure leaves an existing file intact.
Checksums detect corruption, not the truth of physical calibration claims.

The uncompressed archive budget is 16 GiB (2,048 members). Numerical arrays are
written directly to ZIP streams; loading validates each NPY header before
allocating its final array and streams bytes into it while computing SHA-256.
Large tetrahedral stress bases no longer require an additional full in-memory
NPY archive. Metadata is limited to 256 MiB and the manifest to 1 MiB. Existing
version-2 files remain readable; object arrays and executable pickle remain forbidden.

Saved state includes embedded geometry and surface mesh (mm), cylinder settings,
face selections, stored loads and balance operator face IDs, the tensor editor,
all numerical UI entry strings, energy-model ID, quality, time basis/calibration,
initial ensemble choice, local-only execution choice, selected plot and 2D plot limits. It includes the
last-run configuration separately from current edited inputs.
Older files without an initial ensemble choice retain the loaded-Gibbs default.

The optional `material` object stores `material_id` (`aluminum` or
`silicon_wafer`), a `doping_enabled` boolean, `dopant_species` (B/P/As/Sb), and
the raw `dopant_concentration_cm3` entry string. Draft concentration text may
be unfinished; active numerical configuration requires a finite positive value.
Hidden doping drafts survive switching to Al or disabling intentional doping.
The last-run configuration and result store their own material identity and
numeric active doping fields. Editing a Si draft never relabels an older Al result.
Files/configurations without material fields default to Al; the v2 schema is
unchanged. See [Material selection](MATERIAL_SELECTION.md) for the Si execution gate.

Results include all numerical histories, final density, numerical grid, strain
components, opening absorption, survival, interwell flux/occupancy diagnostics,
per-cycle summaries, model/units/provenance, resolution floors and certification
metadata. Python model objects are stored as descriptive parameter records,
never executable serialized objects. Loading restores plots without rerunning
PDE or turning a saved uncertified result into a certified result.

AI credentials and chat history are not project data. Arbitrary Python classes
and object arrays are rejected. Project opening is replay, not automatic solve.
A saved physical calibration still passes the existing validation gate.

The prior .ftgsim v1 files belong to the historical FEM/finite-chain application.
They must be opened with that original AlFatigue executable. Version 2 explicitly
identifies aft.probability-project/2 and refuses to reinterpret legacy physics.
No conversion between the two numerical theories is implied.

## Windows integration

Run app/install_windows.ps1 with -Pythonw pointing to the installed pythonw.exe.
It builds a small AlFatigueProbability.exe launcher under LocalAppData and makes
an Al Fatigue UI desktop shortcut. The launcher uses the installed interpreter
and current worktree; it is not a standalone bundled solver distribution.
launcher.paths is a machine-local file and is not committed. A --check invocation
runs the UI startup smoke. A .ftgsim command-line argument opens that project;
opening a project can create a separate window so it is not silently discarded
by the ordinary single-instance launcher. Existing default file associations and
Windows UserChoice are retained; Open With also offers the new launcher.

Close/reopen an older UI process to get new code. Existing running windows and
unsaved inputs are not force-closed by installation.
