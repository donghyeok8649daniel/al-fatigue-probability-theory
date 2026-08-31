# CAD-style mesh viewport

The common geometry viewport can inspect a model before analysis, a completed
mesh, or geometry shown while separate solver results are being updated. The
viewer never changes the active constitutive scope: all fatigue input remains
the loading-axis scalar normal stress.

## Supported input

| Format | Current support |
|---|---|
| OBJ | vertices, polygon faces and lines |
| STL | binary and ASCII triangular surfaces |
| PLY | ASCII vertices and faces |
| legacy VTK | ASCII `POLYDATA` points, polygons and lines |
| STEP/IGES | not supported; no CAD-kernel backend is installed |

STL, OBJ, PLY and VTK inputs already contain discretized geometry. The current
viewer does not claim to generate a volume mesh from an unmeshed CAD solid.
That requires a separately selected meshing backend in a later stage.

## Navigation

- left drag in 3D: orbit;
- middle drag: pan;
- right drag in 3D: zoom;
- mouse wheel: cursor-centered zoom in 1D/2D and centered zoom in 3D;
- `Reset`: restore extents and the default 3D camera;
- `1D` / `2D` / `3D`: switch display projection without changing mechanics.

The red arrow is the scalar normal-loading axis. It is not a tensor or shear
criterion.

When an imported geometry is open, `Save .ftgsim` embeds that OBJ/STL/PLY/VTK
member with a SHA-256 checksum. Opening the project extracts it to the selected
project output area without silently overwriting a different existing file,
then restores the interactive viewport.

## Commands

Open directly in the integrated application:

```powershell
py -3 -m simulations.fem_tension_app examples/cube_3d.obj
```

Use the standalone viewport:

```powershell
py -3 -m simulations.mesh_viewer examples/cube_3d.obj
```

Generate a headless preview:

```powershell
py -3 -m simulations.mesh_viewer examples/cube_3d.obj `
  --save-preview results/figures/mesh_viewer/cube_3d.png
```
