# `.ftgsim` file format version 1

## Status and scope

`.ftgsim` is the implemented portable project/result format for this
repository. Public searches performed when selecting the name found no
recognized pre-existing `.ftgsim` format; this is a practical collision check,
not a legal registry guarantee.

The file does not broaden the physics. Its active material scope remains pure
single-crystal aluminum, uniaxial normal tension, and crack initiation. A bundle
may carry 3D specimen geometry while `analysis_dimension` remains 1. The
`loading_axis` is the Cartesian tensile direction, distinct from
`crystal_loading_direction_hkl`. Poisson diameter change is kinematic only and
does not introduce transverse applied stress.
Displaying an extruded 2D/3D bar does not add shear or multiaxial mechanics.

## Container

An `.ftgsim` file is ZIP/ZIP64 and contains ordinary UTF-8 JSON, CSV, PNG, SVG,
TXT, or Markdown. It never contains Python pickle, embedded scripts, solver
executables, or dynamically imported plugins.

Required root members:

```text
ftgsim-manifest.json
setup.json
geometry.json
display.json
```

Optional checksummed results currently use:

```text
results/nodes.csv
results/elements.csv
results/metadata.csv
results/summary.json
results/probability_elements.csv
results/initiation_elements.csv
```

An imported portable mesh may additionally be embedded as exactly one of:

```text
geometry/source.obj
geometry/source.stl
geometry/source.ply
geometry/source.vtk
```

The reader checksum-validates the member before extracting it and never treats
mesh content as executable code.

## Manifest

Version 1 uses `format`, semantic `schema_version`, `bundle_kind`, UUID,
creation time, generator/provenance, physical scope, member paths, and SHA-256
checksums. Readers reject unsupported major versions. The version-1 tensile
reader also rejects unknown setup inputs so a parameter is never silently
ignored.

## Separation of concerns

- `setup.json`: material, scalar normal loading, waveform and solver settings;
- `geometry.json`: mesh dimension, loading axis and bar/mesh references;
- `display.json`: view, selected scalar field and deformation scale;
- `results/`: generated data only.

Changing display state cannot change the simulation definition. The FEM
cross-sectional area is not the atomic representative area $A_0$, correlation
area $A_c$, or a count of independent probability samples.

When `results/initiation_elements.csv` exists, it is an optional scalar result
channel containing survival, cumulative initiation, outflux and hazard. Its
setup status must still state whether physical parameters are calibrated; file
presence alone is not evidence of an aluminum lifetime law.

For the single-crystal model, `geometry.json` records both the normalized
loading-axis vector and integer `crystal_loading_direction_hkl`. It also records
whether the scalar axial modulus was supplied directly or projected from all
three cubic constants. Orientation metadata never implies shear-fatigue or
multiaxial support.

## Security limits

The reader validates all paths before extraction, rejects absolute paths,
`..`, backslashes, symbolic links and unsupported member types, limits member
count and decompressed sizes, and verifies SHA-256 checksums. Only checksummed
CSV/JSON members below `results/` can be extracted. Bundle content is never
executed.

## Commands

Open an existing project:

```powershell
py -3 -m app.desktop_ui model.ftgsim
```

Run the headless FEM smoke case and save a project:

```powershell
py -3 -m simulations.fem_tension_app --headless-smoke `
  --output-dir results/data/ftgsim_example `
  --preview-dir results/figures/ftgsim_example `
  --save-project examples/tensile_demo.ftgsim
```

The GUI's `Save .ftgsim` button writes beside the selected output directory.
No Windows registry keys or administrator-level file associations are changed.
