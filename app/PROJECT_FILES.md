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

Saved state includes embedded geometry and surface mesh (mm), cylinder settings,
face selections, stored loads and balance operator face IDs, the tensor editor,
all numerical UI entry strings, energy-model ID, quality, time basis/calibration,
local-only execution choice, selected plot and 2D plot limits. It includes the
last-run configuration separately from current edited inputs.

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
