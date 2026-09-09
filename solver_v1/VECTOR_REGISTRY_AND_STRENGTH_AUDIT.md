# Full vector registry, constrained paths, and ideal interface strength (v9)

Status: **static research/reference**, not an accepted Al fatigue solver.
The v8 vector-row calculations exposed a surviving small-cell fault near
`u=(b/2,sqrt(3)b/6)`. They also showed that zero scalar x winding does not
establish disappearance of a vector defect. Before introducing a finite
source or changing material parameters, this audit checks the underlying
TWO-dimensional registry landscape against the same-condition Al benchmark.

No LJ, exponential density, angular/embedding coefficient, mobility,
temperature, production energy selector or probability equation is changed.
The old TwoRowLJ a0/kappa and reduced fits remain reference results.

## 1. Geometry and the quantity being calculated

Use the already verified +ABC geometry and orthonormal plane frame. The
interface state is now explicitly

    q = (a, u_x, u_y),     u = u_x e1 + u_y e2,
    delta_a = a-h,
    r_Rk^2 = (kh+delta_a)^2 + |R+k tau+u|^2.

All intrabulk spacings remain h, and each half-crystal is rigid. This is not
the fully atom-relaxed dislocation-core geometry of vector_core_v8. The
interface primitive area, coordinate and energy units are

    L0 = 4.05/sqrt(2) angstrom = 2.863782463805517e-10 m,
    b = 1, h = sqrt(2/3), E0 = 1 eV,
    A_atomic_cell = sqrt(3)b^2 L0^2/2 = 7.102490842787125e-20 m^2,
    gamma [J/m^2] = W [eV/interface cell] E0/A_atomic_cell.

`A_c` does not occur. A per-cell uniform-interface barrier is NOT a finite
dislocation-loop activation energy. No arbitrary coherent area or line length
is introduced to turn it into one. Static optimizer steps are not time.

`direct_110` means u=(s,0), with full repeat b. A Shockley translation can be
represented by tau=(b/2,sqrt(3)b/6), of length b/sqrt(3). The full direct
translation a1 can be decomposed into tau followed by a1-tau. This is a
specified crystallographic registry construction, not an inferred specimen
orientation, an empirical slip law or a claim that a straight path is the
global minimum-energy transition path.

## 2. Same energy, full derivatives

Cross-interface pair multiplicity is still k:

    Delta W_pair = sum_(k>=1) k [w(kh+delta_a,k tau+u)-w(kh,k tau)].

For a lower-half atom at depth r,

    Delta rho_r = sum_(k>=r+1) [rho_plane(active)-rho_plane(bulk)],
    Delta W_emb = 2 sum_r [F(rho_bulk+Delta rho_r)-F(rho_bulk)].

Inversion about an interface midpoint exchanges the two half-crystals. Its
scalar density is unchanged and rank-l STF moments acquire (-1)^l. Thus the
two halves have equal per-atom embedding/invariant energies, also for a
general u, justifying the factor 2 here. The angular energy remains

    Delta Q_l,r = sum_(k>=r+1) [Q_l,plane(active)-Q_l,plane(bulk)],
    W_l = 2 D_l sum_r ||Delta Q_l,r||^2,   l=1,2,3.

ALL neighbor changes are summed before F or the per-atom squared norm. The
rank-2 background is zero only because the reference bulk is cubic. The
existing constructor enforces that restriction. There is no global-crystal
embedding and no added/double-counted local GSF or elastic penalty.

The canonical plane sum is unchanged:

    w(d,delta) = (1/A*) sum_G w_hat(d,G) exp(iG.delta),
    w_hat_p(d,G) = 2pi/Gamma(p) (|G|/2d)^(p-1) K_(p-1)(d|G|),
    rho_hat(d,G) = 2pi kappa exp(-dQ)(1+dQ)/Q^3,
    Q = sqrt(kappa^2+|G|^2).

The LJ G=0 layer difference is still evaluated by the exact Hurwitz-zeta
identity in `_power_mean_difference`. Nonzero modes and environment layers
are tolerance-controlled, not finite-neighbor replacements.

For any Fourier channel c(d,G), all ten independent value/derivative entries
are evaluated together:

    c, c_a, iGx c, iGy c,
    c_aa, iGx c_a, iGy c_a, -Gx^2 c, -Gx Gy c, -Gy^2 c.

The existing exact exponential-moment recurrences supply c, c_a, c_aa for
each STF rank. A real orthonormal STF projection preserves the Frobenius
norm. Complete reciprocal shells are retained; truncation uses the absolute
derivative-coefficient envelope, not cancellation of odd shells. Coefficients
are cached at their ACTUAL separation d; they are not frozen at a reference
radius. The last-shell/last-layer values are empirical diagnostics, not a
rigorous bound on the entire nonlinear energy tail.

With z_r=rho_bulk+Delta rho_r and i,j in {a,x,y},

    (W_emb)_i = 2 sum_r F'(z_r) rho_r,i,
    (W_emb)_ij = 2 sum_r [F''(z_r) rho_r,i rho_r,j + F'(z_r) rho_r,ij],
    (W_l)_i = 4D_l sum_r Q_r : Q_r,i,
    (W_l)_ij = 4D_l sum_r [Q_r,i:Q_r,j + Q_r:Q_r,ij].

`FullRegistryInterface` returns W, its 3-vector gradient and symmetric 3x3
Hessian. Old scalar projections are independently reproduced. Finite
differences are validation only. No source interpolation enters this model.

## 3. Independent Al benchmark and provenance

The reference is Mishin et al., *Interatomic potentials for monoatomic metals
from experimental data and ab initio calculations*, PRB 59,3393 (1999),
[DOI](https://doi.org/10.1103/PhysRevB.59.3393), using the
[NIST Al99 setfl release](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/).
The NIST entry/conversion provenance was rechecked on 2026-09-09. The source
file checksum remains

    60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284.

`MishinVectorInterfaceReference` computes per-atom source densities and
energies on exactly the same rigid half-crystals, 0 K, lattice 4.05 angstrom.
Its source-defined cutoff and cubic interpolation belong to the BENCHMARK,
not to our LJ/Bessel theory. It is neither experimental truth nor a selectable
production potential. The old source's scalar-path results are reproduced,
including the physical-to-reduced first/second derivative factors L0/L0^2.

The tiny source normal residual at exactly 4.05 angstrom is retained. Its
local zero-normal-force interface spacing is 0.816498342 L0, compared with
h=0.816496581 L0. Energies of stationary states are measured relative to the
model's own local perfect stationary state, without changing the bulk units.

Candidate parameter SHA256 remains

    9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655.

It is `angular_monotone_opening` from matched_v3, with NO refit. The
`parameter_metadata.path=direct_110` field records that preserved scalar
constructor; it is NOT a restriction of the new vector evaluator or its grid.
The new manifest's coordinate order is explicitly [a,ux,uy].

## 4. A constrained maximum is not necessarily a transition saddle

At the old direct midpoint (h,b/2,0), W_x=0 but

    candidate: W_a=-3.46542, W_y=-1.82985 [eV/L0],
    source:    W_a=-2.28785, W_y=-1.40263 [eV/L0].

It is therefore NOT a stationary point of the full state. It needs normal
and transverse holding tractions. Likewise a high-symmetry on-top state has
TWO negative Hessian eigenvalues after normal relaxation: it is not an
index-one slip saddle.

We solve bounded local stationary equations with the analytic Hessian. A
reported minimum must have positive full Hessian, and a transition saddle
exactly one negative eigenvalue. Small force residual, nonzero eigenvalue
margins and no contact with a numerical search-window bound are required.
This explicitly rejects the source's flat separated region: zero force at
its cutoff is not evidence of an intact minimum/saddle.

Perturbing the unique negative mode of the Shockley saddle in both signs,
then minimizing all three coordinates, reaches the pristine and intrinsic-
fault minima. This verifies LOCAL basin connectivity. It does not certify
the globally lowest saddle among every possible atomistic/finite-defect path.

| Normal-relaxed quantity [J/m2] | unchanged candidate | matched source |
|---|---:|---:|
| connected forward saddle above pristine | 0.173237970 | 0.172002365 |
| intrinsic fault above pristine | 0.126800169 | 0.150479477 |
| reverse saddle above intrinsic fault | 0.046437801 | 0.021522888 |
| old fixed-a direct midpoint (constrained) | 0.679976011 | 0.603038353 |

The forward saddle agrees within about0.72%, but the fault energy is about
15.7% too low and the reverse barrier is about2.16 times the benchmark.
Matching one forward barrier therefore does NOT establish correct defect
stability or recovery. These newly relaxed stationary quantities were not
fit targets. Existing bulk C11/C12/C44 errors (about87.82/64.83/49.27 versus
114/62/32 GPa) also remain; this audit did not recalibrate them.

The candidate's partial saddle is at (a,ux,uy) approximately
(0.842845,0.309064,0.178438), with eigenvalues
(-1.96009,6.93162,20.57499)eV/L0^2. The source saddle is at
(0.849292,0.371406,0.214431), with (-1.86212,4.54159,12.81703).
The fault minimum remains near u=tau for both. This supports interpreting
v8's small-cell vector fault as a registry-fault state rather than accepting
its scalar b/2 partition jump as newly accumulated plastic strain. It does
NOT prove that the periodic small-cell fault is an isolated Al dislocation.

## 5. Correct conjugate work and the first uniform-interface fold

For specified normal and two resolved shear tractions,

    G(q)=W(q)-A_atomic_cell L0/E0 [Tn(a-h)+tau_x ux+tau_y uy],
    f_i*=T_i A_atomic_cell L0/E0.

Here T in Pa, E0 in J/eV, q dimensionless. No arbitrary chi, correlation area,
mass, mobility or time conversion appears. The original reduced-model
f*=kappa sigma/E remains unchanged. These are declared interface tractions;
a nominal specimen stress requires a specified crystallographic projection.

For pure x shear, let z=(a,y). At each prescribed x follow the local branch

    W_z=0,  H_zz>0,  z'(x)=-H_zz^-1 H_zx,
    Wbar_x=W_x,
    Wbar_xx=H_xx-H_xz H_zz^-1 H_zx.

The first positive-to-negative zero of Wbar_xx, with H_zz still positive,
is the local stress-controlled fold. It has a zero eigenvalue of the FULL
3x3 Hessian and a small normal/transverse force residual. It is not just the
first failed root solve. We bracket from the pristine branch, refine the
root, and repeat with independent31/61 and41/81 point scans.

`z(x)` is a **locally transverse-equilibrated fixed-x branch**, not a proven
global minimum-energy path. Away from stationary points its gradient need
not be parallel to that curve. This distinction matters for reaction rates.

The constraint-ablation experiment changes NO parameter:

| First required x-traction maximum [GPa] | candidate | source |
|---|---:|---:|
| normal a and transverse y fixed | 7.54336 | 6.46164 |
| normal a free, transverse y fixed | 4.76890 | 4.39037 |
| normal a and transverse y free | 3.12402 | 2.76750 |

Only the last row is the unrestrained three-coordinate uniform-interface
fold. The others require the recorded holding tractions. Removing artificial
path constraints materially lowers the predicted IDEAL strength. Nevertheless
even the independent Al benchmark is in GPa for this flawless rigid-half
experiment. Thus converting GPa to an experimental MPa yield by a units patch
or manual factor would be wrong. Existing/mobile defects, finite sources,
core energetics and kinetic data constitute a different physical problem.

The free candidate fold is at (0.830353057,0.134502071,0.046509103)L0;
the source at (0.822867359,0.129867226,0.043408744)L0. The41/81 root scans
agree in traction to better than2e-8MPa. That is numerical repeatability,
NOT a material uncertainty estimate; the model/source difference is about
12.9%. The other two Hessian eigenvalues stay positive at both folds.

## 6. Actually executed lower-stress mixed scenarios

32 static states were solved across both models:

- pure x shear: 0,25,50,150,0,-50,0MPa;
- mixed (Tn,tau_x,tau_y): (0,0,0),(25,25,15),(50,50,25),
  (150,100,-50),(0,0,0)MPa;
- compression/shear: (0,0,0),(-50,50,15),(-150,100,25),(0,0,0)MPa.

All remain on the positive-Hessian pristine branch. At150MPa pure shear the
candidate has ux=.00385251, uy=.0000439106 L0, not an interwell jump.
Maximum solved force residual is7.65e-15eV/L0. Actual registry displacement
after static unload is at most1.76e-15L0 across both models: numerical
zero, not tiny claimed plasticity. This is athermal static unload, NOT a
finite-temperature zero-stress hold, probability flux test or fatigue test.
No crack probability is calculated from these static states.

## 7. Independent numerical verification

17 new tests check full gradients/Hessians, scalar-path recovery, threefold
covariance, lattice translations, per-site component sums, direct real-space
validation, reciprocal/depth tolerance refinement, saddle connectivity,
Schur curvature, source units/derivatives, static load/unload and invalid or
flat-source-state rejection. Existing tests are not removed or weakened.

Primary run:122.69s,882 actual grid points,22 fixed/stationary/descended
states, two resolutions of both relaxed branches. Constraint/stencil
rechecks:29.46s. Wall times include concurrent regression work and are not
kinetic quantities.

At a=1.09h,u=(.237,.087), direct radial/layer24/48/72 energy errors versus
the infinite reciprocal model are8.85e-6,1.16e-6,3.48e-7eV/cell. The finite
real-space algebraic tail is not represented as machine-precision agreement.
At fixed-spacing intrinsic fault the corresponding errors are9.73e-7,
5.87e-8,1.13e-8eV/cell. The direct sum is validation only.

Changing analytic reciprocal/neighborhood tolerances2e-11 to2e-13 at four
stationary states changes energy by at most4.64e-13eV/cell, gradient by
2.81e-12eV/L0 and Hessian by4.50e-12eV/L0^2. First/second derivative
central differences at step1e-5 have candidate maximum absolute errors
6.09e-9eV/L0 and7.35e-8eV/L0^2 over the sampled off-path states.

The tabulated source needs its own stencil audit. A coarse stencil crosses
spline knots, producing an8.39e-4eV/L0^2 error at one state. That failure is
retained in the CSV; refinement to1.25e-6 lowers it to2.48e-10. At smaller
steps roundoff eventually stops improvement for both evaluators. Neither
the source's coarser interpolation/stencil error nor the direct finite-tail
error is concealed by reporting only the best derivative check.

See `results/fcc111_active_interface/vector_registry_v9/SCHEMA.md` for
machine-readable outputs and `verification.json` for actual final regression
counts/timings. No old calibration/result file is overwritten.

The first full suite had321 passes and one new SOURCE finite-difference NaN
failure. Isolated runs,9600 derivative evaluations and3000 allocation-stress
checks did not reproduce it. Finite-output fail-fast checks were added, not
clipping, looser tolerances or skips. The final full suite passed322 tests,
app31, targeted17, smoke and diff checks. All940 saved special/stress/fold/grid
states were re-evaluated with the guards and matched the stored numerical
values exactly. The original intermittent anomaly's ROOT CAUSE REMAINS
UNIDENTIFIED; a later pass is not falsely described as a proved root-cause
fix. It is explicitly retained in verification.json and the handoff. This
is another reason not to promote the research/source comparison into a
production material certification.

## 8. Interpretation and next gate

The missing registry direction materially biased the constrained path. The
new reference eliminates that particular limitation and can apply two shear
components consistently with the same analytic energy. It is not a new
3D spatial probability PDE. Normal opening, local registry slip and finite
dislocation generation/motion remain distinct physical geometries.

Next material work should constrain the independent elastic tensor AND the
fault energy, saddle position/curvature and reverse barrier, with held-out
vector states. A good USF value alone is insufficient. Any next analytic
embedding extension requires a documented compatibility/rank audit, not
parameters chosen to force slip. Defect/source calculations must then use
that independently assessed surface and actual finite geometry, without
arbitrary pins or activation length. v8's pair annihilation is not itself a
yield measurement or proof of failure of every possible defect state.

Physical M_a, M_s and t0 remain unavailable. Seconds/Hz, fatigue lifetimes,
correlation area/independence and spatial specimen calibration are NOT
provided by this result. Production energy registration and the later
CAD/meshing UI gate remain closed. The scientific result is a better-defined
vector static mechanism and a quantified remaining material discrepancy,
not "plasticity successfully implemented".

Reproduce with an available Python3 interpreter:

    python -m pytest solver_v1/test_vector_interface_reference.py -q
    python -m solver_v1.run_vector_registry_audit --output <new-result-directory>
    python -m solver_v1.run_vector_registry_rechecks

The recheck command intentionally uses the default v9 output location for
its independently computed companion files. The primary command refuses to overwrite a completed
result; use a fresh directory for an independent execution.
