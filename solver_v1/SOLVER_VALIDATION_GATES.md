# Solver readiness and later specimen workflow

## Current decision

The user requested stage-by-stage physical checks and a solver report BEFORE
approving a specimen UI redesign. The full-FCC Al active-interface solver is
not ready. This is a scientific readiness failure, not a request to manufacture
plasticity or physical time.

The current audit establishes independent bulk elasticity and consistent
energy/length/area units, but its candidate surfaces fail interface validation.
The existing N=1 deterministic reference PDE remains available; it must not be
confused with a fully calibrated spatial aluminum fatigue solver.

## Gates

1. Geometry and energy: infinite lattice, counting, analytic derivatives,
   independent direct sums and numerical tails.
2. Material calibration: independent observables, stability, identifiability,
   source/temperature/relaxation consistency, correct dimensional normalization.
3. Interface: physically credible GSF/cleavage and coupled saddle topology,
   not just positive curvature at one bulk minimum.
4. Probability evolution: generator compatibility, survivor moments, net/gross
   registry transport, opening absorption, independent numerical resolution.
5. Kinetics: actual coordinate mobility data before seconds/Hz. Model-time
   reference mode can remain available but is not a physical-frequency prediction.
6. Spatial specimen model: derived mechanics/local-state coupling, boundary
   conditions, orientation and mesh convergence. One local density is not
   automatically a field over a CAD mesh.
7. Report evidence and obtain the user's confirmation before UI implementation.

Completed numerical tests do not automatically satisfy physical gates.
Read results/aluminum_full_fcc_calibration/audited_v2/readiness.json.
Steps 3 onward must not be silently unlocked by a small bulk least-squares loss.

## Intended UI after validation and confirmation

The requested workflow is retained:

    Geometry import / editable cylinder specimen default
       -> Meshing
       -> Preprocessing (material, orientation, load, constraints, numerics)
       -> Solving
       -> Postprocessing on the mesh

Before implementation, agree on the supported import format and distinguish
CAD solids from triangulated surface meshes. Preserve the original file and
units. A default cylinder's dimensions must be explicitly editable; they are
specimen geometry, not hidden constitutive calibration parameters.

The mesh postprocessor should expose selectable physical fields, for example:

- applied/local stress and strain components;
- local opening initiation and survival;
- well occupancy, net registry transport, gross configurational activity;
- physical resolution masks, kept separate from numerical residual fields.

Spatial fields require an actual mechanical solution and a defined mapping to
each local probability state. Do not simply paint a single global probability
on a mesh and claim a spatial crack-risk prediction. If loading/state are
uniform, a uniform field may be the correct result and should be labeled.

Cell IDs, point/cell association, units, model ID, time basis, physical floor,
and mesh/solver provenance must accompany each exported field. Numerical
residuals must not appear under physical crack-probability titles.

## Lightweight implementation principles for the later UI

- Separate geometry, mesh, preprocessor, solver worker, result store and viewer.
- Keep PDE work off the UI event loop; provide progress and cancellation.
- Reuse identical local histories and cache postprocessing where mathematically
  justified. Do not silently lower solver accuracy.
- Changing a field, colormap, language, view, or specimen correlation area must
  not rerun the PDE.
- Keep Preview distinct from verified analysis. Missing/under-resolved values
  get an explicit mask, not an invented physical value.
- Preserve Korean/English central localization and view/zoom state.
- Do not add heavy CAD/meshing dependencies before selecting the minimal
  supported workflow and confirming the scientific coupling.

## Scales must remain separate

A_atomic_cell belongs to atomic energy normalization.
A_c belongs only to statistical specimen probability aggregation.
Mesh areas/volumes belong to geometry, quadrature and structural mechanics.
They are not interchangeable, and refinement must not change the material law.

Physical kinetic mobility and spatial correlation remain unavailable unless
independent data are supplied. No static energy or cylinder dimension creates
a physical time scale.
