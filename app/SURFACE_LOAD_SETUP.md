# Surface load setup (research preparation, not a spatial solver)

The original LJ/Bessel energy, PDE, calibration, local probability and physical
clock are unchanged. Arbitrary surface tensors and balancing tractions are NOT
silently reduced to the scalar axial production PDE. Unsupported configurations
are refused at Solve. No structural displacement/stress field is calculated.

## User workflow

Generate/import the surface mesh, rotate using the right mouse button (middle
button pans, wheel zooms), and left-click a visible connected planar patch.
Ctrl-click toggles patches. This uses triangles, not native CAD curved-face IDs.
The top/bottom shortcuts remain optional global-Z conveniences, not a geometry
assumption for imported specimens.

Enter global XYZ stress components in the 3×3 matrix. All entries have MPa
units. Mirrored entries must contain identical expressions; this conservative
rule guarantees symmetric Cauchy stress throughout the history. Expressions
support the existing restricted sin/cos grammar, not arbitrary Python.
The historical default is XX. Cylinder-axis tension requires ZZ (axis Z).
Frequency now belongs to Load. Each added load snapshots its expressions,
numeric parameters, model frequency, selected triangles and area. Loads add on
overlapping triangles. Delete and re-add to revise a stored load.

## Resultants and balance

For outward unit normal n_i, centroid r_i and area A_i:

    t_i(t) = sigma_i(t) n_i
    F(t) = sum_i A_i t_i(t)
    T(t) = sum_i (r_i-r_ref) cross [A_i t_i(t)]

MPa × mm² = N; torque is N mm. r_ref is the mesh vertex mean. A closed,
consistently oriented surface is required. Surface validation does not certify
self-intersections, nested cavities or a multi-body CAD assembly.

Equal tensile stress on opposite faces gives opposite forces because normals
are opposite. Equal numerical stress on arbitrary unequal areas does NOT imply
equal force. A supported specimen can be balanced by reactions; these are not
inferred here. The optional correction is an explicitly added applied traction,
not an elastic support reaction.

On user-selected correction faces, minimize sum_i A_i |delta t_i|² subject to
six force/torque constraints. Set y_i=sqrt(A_i) delta t_i. The columns of B are
sqrt(A_i)[I; cross(r_i-r_ref)/L], with L the bounding-box diagonal. Then

    y = -B^+ [F; T/L].

Require rank(B)=6. The UI asks once before activating this operator. It is
evaluated against the complete instantaneous load, so it is not a t=0-only
correction. The dialog reports t=0 resultants explicitly. Repeated clicking
does not reconfirm an active correction. Adding/removing loads or remeshing
clears approval; imported correction faces are a proposal requiring approval.
The matrix entries are never overwritten by this correction.

## Declarative format and AI boundary

`aft.surface-loads/1` JSON is the initial setup exchange format. It includes
units, frozen load data, correction faces, and SHA-256 binding to actual vertex
and triangle arrays. Import on a different mesh is refused. Import does not
execute scripts/PDE, change energy calibration, or approve corrections.
It is a surface-load format, not yet a complete geometry/solver language.

An external AI chat provider is NOT connected in this change. The proposed
future boundary is AI -> validated setup draft -> user review -> import.
No credentials are assumed and no model/research data are uploaded. Provider
configuration, data-sharing consent and full solve-schema validation remain
separate work. A deterministic expression parser is not presented as AI.

## Limits

No FVM/FEM volume mechanics, support reactions, multi-axis production generator,
material calibration or physical seconds/Hz are added here. Physical time mode
is not accepted for storing surface histories in this version. Mesh changes
invalidate stored triangle assignments rather than silently transferring them.
Local refinement after loads therefore requires reassigning loads.
