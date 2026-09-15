# Single load input and interactive surface diagnostics

The axial stress mean/amplitude widgets live only in Face Load. Their StringVars
are also the values read by DesktopApp._config; Pre has material/time/settings.
Presets modify these same entries. Saving a face selection is not required to
run the independent local axial PDE. The setup surface is not a solved spatial
stress field. Nonzero shear or custom tensor expressions explicitly block this
scalar solver instead of being silently ignored. They remain preparatory inputs.

Left click picks the closest intersected triangle and its connected coplanar
patch. Ctrl-click toggles patches. Right drag rotates, middle drag pans and wheel
zooms. STL/OBJ do not retain CAD topology, so curved surfaces are faceted patches;
this is not a complete Autodesk solid-model selection engine. Group selection
remains available for cylinder top/bottom/lateral faces.

Local refinement marks selected triangle edges. Adjacent polygons include the
same shared edge midpoints and are triangulated with an interior centroid.
Consequently no T-junctions are introduced. Geometry, area and enclosed volume
are preserved; imported curved geometry is not reconstructed. The 20,000-face
limit applies. Remeshing invalidates saved face assignments.

## Surface extrapolation

The separate mesh map assumes the same local opening probability p(t) at every
surface patch. It never presents this as a spatial stress/probability solution:

    log S_i(t) = (A_i/A_c) log1p(-p(t))
    P_i(t) = -expm1(log S_i(t))
    product_i S_i = exp[(sum_i A_i/A_c) log1p(-p)]

Color variation comes only from patch area. Subdivision changes the individual
patch probability but preserves aggregate survival. All values are explicitly
uncertified mathematical extrapolations, even when p is below numerical
resolution. This map does not alter or create a physical certification field.
It uses the full displayed surface area, not the manually entered stressed area;
this whole-surface-equivalence assumption is explicit in the map. No local PDE
rerun occurs when viewing, rotating, refining, changing time index or language.

## Geometric volume versus characteristic correlation volume

For a closed, consistently oriented triangular boundary,

    V_geom = abs(sum_tri v0 dot (v1 cross v2))/6.

Coordinates are translated near the origin to reduce cancellation. Open meshes
or inconsistent orientation are rejected. Edge closure is not a proof that an
import has no self-intersections/cavities; the UI calls this a volume estimate.

A bulk initiation law could use an independent volume-region hypothesis:

    N_eff,V = V_eff/V_c
    log S_spec = (V_eff/V_c) log(1-p_volume).

It needs a bulk-region local probability, a correlation volume V_c and a spatial
definition of V_eff. The present local event is normal opening of an interface;
its p cannot automatically be relabeled p_volume. V_geom equals V_eff only under
an explicitly uniform whole-body loading/hazard assumption. A surface mesh does
not contain volume cells, and assigning arbitrary volumes to its triangles would
not supply a volume probability field. Therefore this change computes V_geom but
does not invent V_c or replace A_c. If a supported bulk-initiation law is derived,
area and volume hazard contributions can be combined in log survival with their
own statistical assumptions; they must not count the same events twice.

Neither A_c nor V_c affects the energy, f*=kappa sigma/E, or strain.

## Actual local PDE check

`python -m app.run_load_input_audit` evaluates the unchanged TwoRowLJ reference
at 25 cycles/model-time, one cycle, E=69 GPa, on three grid/integrator settings.
The generated `results/ui_load_audit/runs.json` contains actual runs. Combined
grid/time changes are not an independent directional refinement certificate.

At the finest 61x90 grid:

| MPa mean +/- amplitude | opening p | max abs plastic strain |
|---|---:|---:|
| 50 +/- 100 | 0 | 1.90006e-18 |
| -500 +/- 400 | 0 | 6.22355e-19 |
| 900 +/- 2000 | 5.06131e-15 | 1.21977e-12 |

Extreme-load opening changes from 1.19004e-13 (preview) to 1.06951e-14 to
5.06131e-15. This is unresolved under the observed resolution change.
Extreme-load plastic strain changes by ~26% between the last two settings;
no converged residual plasticity claim follows. No unload/hold was run here.
All nine runs satisfy p=cumulative absorbed mass exactly and survival=1-p;
maximum strain-decomposition residual is 6.33235e-18. Flux consistency residual
is zero under the solver's stored discrete flux convention (not independent
continuous-flux quadrature). Stress mapping gives kappa=86.29296488740997 and
150 MPa / 69000 MPa = .002173913, f*=.187593402: no MPa/GPa factor error.

This checks numerical/input consistency, not material calibration. The active
default remains TwoRowLJ, not a validated Al bulk/interface model. Physical
collective-coordinate mobilities and production seconds/Hz remain unavailable.
