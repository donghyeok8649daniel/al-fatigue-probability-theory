# Full current site energy on infinite Bessel rows — v22 research

Study started2026-09-10 19:49UTC (2026-09-11 04:49KST).
Starting local/fresh remote a8e4130934ab8bb67f2a657a489d269d10b20b73,
clean probability-pde-solver-v1. User allows about four hours of continued
research. No production/PDE/UI/kinetic/default-material parameters are changed.
Additional explicitly unadopted joint core-force trial coefficients are saved
separately; they are not a replacement of the bound v20 reference.

## Question and staged work

The v21 line uses only the same energy's elastic tensor. It therefore cannot
test the Bessel energy's nonlinear core resistance. The older isolated-core
code accepts the earlier17 channels and quadratic invariants, while the v20
candidate additionally contains the radial Eg quadrupole, quartic rank3,
nonlinear even saturation and an optional density/rank1 product. Passing only
its base surface would OMIT energy terms. No existing production result is
asserted to have done that; this is the representation gap to close first.

1. Derive and implement ALL current per-atom terms on infinite atomic rows.
2. Compare energy/first/second derivatives against the existing full interface
   model and independent direct atomic row values; preserve old core tests.
3. Bind the unchanged positive-LJ v20 candidate and its actual elastic tensor.
   Diagnose isolated screw minima, independent force/Hessian/domain/ring errors.
   Candidate material remains unadopted; this is not actual Al yield validation.
4. If stable/converged, test signed static stress and translation resistance;
   otherwise preserve the failed state and identify the missing physics/error.
5. Run targeted/full regressions, save complete/incomplete status, commit/push
   only after tests. Keep at least the final regression window in the time budget.

## Channel and energy contract

Each row stays infinite along e1. The actual neighbor displacement is

    r_l,ij² = (l b + dx_ij)² + dy_ij² + dz_ij².

Sum l over all integers with the existing Poisson/Bessel formulas including
G=0. Sum every retained neighbor-row contribution at each site before applying
the nonlinear law. Use22 channels: pair,rho,Q1(3),Q2_even(5),Q3(7),Q2_odd(5).
The last five are the already existing radial quadrupole, not a new potential.
They use the exact STF basis already used by the other moments. Define

    x_i=(rho_bulk+delta rho_i)/rho_ref,
    QE_i=(Q2_odd,i-eta Q2_even,i)/N,
    f(I;k)=I/(1+k I), f'=1/(1+kI)², f''=-2k/(1+kI)³.

The fixed eta,N,reference curvature are copied from the SAME interface model.
The "Eg" label denotes its cubic-reference tuning; all five STF quadrupole
components remain away from that reference. "Odd-range Q2" is still an even-
parity rank2 tensor: the radial range comes from the rank3 kernel.
For scalar embedding F and current coefficients c,

    e_i = .5 delta Phi_i + F(rho_i)-F(rho_bulk)
        + D3 ||Q3_i||² + K3 ||Q3_i||⁴
        + D2 f(||Q2_even,i||²;alpha/N2)
        + DE f(||QE_i||²;alpha)
        + D1 g(x_i)||Q1_i||².

No global density/plane-by-plane embedding. No reference gauges change during
strain/core motion. The cubic reference Q=0 does not suppress noncubic moments.
With summed channel z, for I=||Bz||²:

    I_z=2 B^T Bz, I_zz=2 B^T B,
    (f(I))_zz=f'' I_z I_z^T+f' I_zz.

The product g(x)I additionally needs BOTH mixed terms
g' (x_z I_z^T+I_z x_z^T) and g'' I x_z x_z^T. The atomic Hessian follows

    E_UU = sum_i [ z_U^T e_zz z_U + sum_c e_zc z_c,UU ].

Default historical site response remains unchanged. A separate subclass supplies
the full current law. It is a static research interface, not a PDE registration.
Neighbor-row ring and free-domain truncation remain separate from the infinite
atomic-row sum; neither is silently called infinite-domain convergence.

## Original candidate core and physical status

The initial current-material R2/r3 centered state is a stationary saddle:
force residual1.19e-7eV/L0, lambda_min=-2.17874eV/L0^2. Both signed negative-
mode descents lead to equal-energy fixed-boundary minima, E=.68878052884eV
per straight-row repeat. R4/6/8/10 branches and independent ring5/7 refinement
have been executed. This is not a material-adoption, Peierls or actual-yield gate.

## Matched source-only comparison (added after representation verification)

The existing checksum-bound NIST Al99 source is now an independent comparator
on the SAME vector-row geometry and static boundary protocol. It is **not**
substituted for LJ. Its published finite cutoff defines that external source;
the canonical research candidate remains an infinite Poisson/Bessel sum.
The source kernel differentiates its own radial tables in dimensional Angstrom,
then multiplies derivatives by L0_A and L0_A^2 for our reduced coordinates.
Both models count one-half pair energy and apply F to each atom's density.

For the source's actual reference environment,

    V_ik(R)=partial_i partial_k [phi(R)+2F'_bulk f(R)],
    D_iJ=sum_R f_i(R) R_J,
    C_iJkL=[.5 sum_R V_ik(R) R_J R_L + F''_bulk D_iJ D_kL]/Omega.

This is the affine energy tensor at the declared geometry, including any
prestress; it is not a silently substituted experimental modulus. Its
anti-plane block must equal the separately assembled row Hessian. The source
neighborhood guard requires alpha(ring+1)-2 max|u_transverse| > r_cut/L0,
where alpha is the smallest singular value of the transverse lattice map.
This encloses all source cutoff neighbors even after displacement. No analogous
finite-support claim is made for the candidate's algebraic LJ tail.

Primary provenance: Mishin et al., PRB59,3393(1999),
[NIST source and conversion notes](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/).
The pre-existing verified source SHA is60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284.
This is static atomistic potential validation, not experimental yield or fatigue.

## Translation and loading experiment definitions

The Burgers/line direction is x=e1. A primitive in-plane translation
a2=(b/2,sqrt(3)b/2,0) shifts the core glide position by d=sqrt(3)b/2 along y.
The glide direction differs from the Burgers direction. For integer j shift k,

    u_i^seed = u_(i-k a2)^old - u_affine(r_i-k a2) + u_affine(r_i).

Only free rows use this initial guess. Exterior Dirichlet data stay fixed;
all final coordinates minimize the unchanged total energy. No final position
is imposed. A surviving translated minimum is evidence of a metastable state
on THAT boundary, not an infinite-domain Peierls stress. A low-force result
must additionally pass its Hessian/Morse check. Same-boundary saddle differences
remain energy per straight-row repeat, not finite-loop activation energies.

All MPa schedules in this study are **resolved shear tau_xz**, not nominal
uniaxial tensile stress. With Burgers vector b e_x and line tangent e_x,
the Peach--Koehler driving force per length is

    f_line=(sigma dot b_vec) cross e_x
          = b(tau_xz e_y-tau_xy e_z).

The imposed linear far field has tau_xy=0, hence positive tau_xz drives+y,
consistent with the sign of the matched source-core translations. Converting
nominal axial stress requires a specified loading axis and its Schmid
projection; no crystallographic loading orientation or factor is invented here.
The winding-one boundary introduces a **pre-existing straight screw**, not
homogeneous defect nucleation or a measured source population. This distinction
is essential when comparing ideal strength, glide onset and specimen yield.
No absorbing opening boundary is present in this static core experiment;
opening probability and slip-versus-crack event ordering are not calculated.

The triangle circulation sums wrapped x-phase differences after subtracting
affine shear. It reports projected integer winding only; partial splitting,
core width and full vector Burgers classification are not inferred from it.
Adjacent-plane profiles retain the actual logical-row ABC offset. They are
not an interpolated or rescaled continuum slip field.

## What makes this a dislocation model without replacing Bessel

Instead of a single rigid interface coordinate, let each transverse atomic
row i=(j,l) have a three-component displacement U_i, constant along the straight
line x=e1. With its actual ABC reference vector r_ij, define the row-channel
kernel K by the infinite atom sum and z_i=sum_j[K(r_ij+U_j-U_i)-K(r_ij)].
The reference bulk density is added back before F. The displacement-dependent
energy is the sum of the per-site e_i above, and its exact force is

    partial E/partial U_l = sum_(i,j) (delta_lj-delta_li) K_ij,U^T e_i,z.

For free-degree difference matrices B_ij and site derivatives D_i=z_i,U,

    H = sum_i D_i^T e_i,zz D_i
      + sum_(i,j) B_ij^T [sum_c e_i,zc K_ij,c,UU] B_ij.

This is both the matrix-free Hessian and the new explicitly assembled matrix;
their equality is independently tested. No harmonic core substitution, local
gamma Hessian double counting, force clipping or empirical plastic law appears.
The imposed screw topology is the exterior phase circuit

    sum_circuit wrap(delta(2 pi U_x/b)) / (2 pi) = +/-1.

Spatially nonuniform registry, its core and the surrounding elastic field
therefore emerge from the SAME atomic energy. The far-field elastic tensor is
the same candidate's bulk Hessian. This is a straight q_x=0 reference, not a
finite curved line, emitted loop, 3D dynamic PDE or macroscopic plastic-strain
history. The physical atomic lattice is not an arbitrary continuum mesh.

## Representation and finite-domain results

The independent real-atom evaluator uses no reciprocal coefficient or Bessel
derivative. At7 declared displacement/separation points the128-image reference
agrees in all22 channels; the complete small-core energy differs by
6.66e-16eV/row repeat, free atomic gradient by2.73e-14eV/L0. Whole nonuniform-core
central differences converge quadratically: energy derivative error1.81e-6 ->
7.11e-9 and Hessian-vector error1.06e-4 ->4.14e-7 as the step decreases from
2e-4 to1.25e-5 L0. These test small arithmetic errors, not material accuracy.

On the unchanged candidate, R4 L-BFGS and Newton-CG energies differ by4.6e-14eV;
their inner2L0 displacements differ by7.10e-8L0. For the same branch/ring5:

| Free-domain comparison | max inner2L0 displacement change / L0 |
|---|---:|
| R4 -> R6 | 0.00664853 |
| R6 -> R8 | 0.00373818 |
| R8 -> R10 | 0.00208094 |
| R8 ring5 -> ring7 | 3.16359e-7 |
| R10 ring5 -> ring7 | 4.12060e-7 |

Domain effects dominate the tested transverse-neighborhood effect. Neither is
a proof over the infinite domain or every channel tail. The original Al-target
v20 parameter SHA8736cb9d28f1430e3991b1046ef42544cefaa3c9a4f5d743e0d02931329f5c2e
is unchanged. Its C11/C12/C44=114/62/32GPa versus the source's actual
113.79497/61.55332/31.59481GPa at nominal4.05Angstrom: each has its own tensor.

Absolute stored core energies at DIFFERENT free radii or neighbor rings are
not automatically energy-convergence errors. The sum includes every exterior
site affected by a free row; enlarging the ring changes that energy-site set
and includes more fixed far-field energy. For example the wider candidate's
R8 energy is .8758608333eV at ring5 and .9219939748 at ring7, despite the small
local-coordinate change. Compare identical-boundary/operator energy differences
for a saddle or minimizer check. The combined report explicitly labels when
an absolute-energy difference does NOT use the same represented operator.
No state-independent exterior contribution is hidden by changing the energy
reference after a run; full infinite-domain elastic/core matching is still open.

The sampled adjacent-row25–75% span is0.874L0 for the candidate R10/r7,
versus2.760/2.791/2.830L0 for source R8/10/16. Crossings are interpolated only
between explicitly listed atomic-row samples; brackets are saved. This is a
diagnostic span, NOT a fitted continuum core radius or partial separation.
The much sharper candidate core motivates independent force validation.

## Core symmetry saddles and low-stress return

About center(d/6,h/2), the zero-load symmetry maps(j,l)->(-j,1-l), with residual
displacement w=U-U_V transforming as(w_x,-w_y,-w_z). Energy and involution are
tested; midpoint Newton roots use the ORIGINAL energy and indefinite Hessian.
They are not relaxation dynamics. Full finite-domain spectra give Morse index1.
Both signs of its negative mode at R8 descend to the two actual stable endpoints.

| Model | R / L0 | ring | reconstruction barrier eV / straight-row repeat |
|---|---:|---:|---:|
| Unchanged candidate | 8 | 5 | 0.0148345888 |
| Unchanged candidate | 10 | 5 | 0.0145539015 |
| Unchanged candidate | 8 | 7 | 0.0148341800 |
| Mishin target | 8 | 5 | 0.0003362987 |
| Mishin target | 10 | 5 | 0.0006445512 |
| Mishin target | 12 | 5 | 0.0008267211 |
| Mishin target | 16 | 5 | 0.0010162721 |

The source barrier has large relative boundary sensitivity. It is not yet an
infinite-domain number. Source coarse energy-curvature stencils(.002,.001L0)
are inadequate near its finely tabulated spline features:13-step refinement
shows gradient-based curvature agreement to2.84e-11eV/L0^2. Very small energy
differences instead become roundoff-limited. This distinction is preserved.

Executed static schedules are0 -> +5 -> +20 -> +50 ->0MPa and a separate
0 -> -5 -> -20 -> -50 ->0MPa. All8 states for each of current R8 and source
R8/R10 pass finite-boundary force/positive-Hessian checks. Candidate return
errors in inner2L0 are7.36e-9 and2.76e-10L0. No persistent registry change is
resolved by that static protocol. Source negative loading switches to its
twofold partner, with return error<7e-10L0 to that partner; the switch occurs
by-20MPa at R8 but only by-50MPa at R10 on this coarse schedule. That domain
dependence and the distinction between reconstruction and full lattice glide
forbid calling those stresses an Al Peierls/yield measurement. The triangle
containing projected circulation can jump even on a smooth stable response;
it is not a continuous core-position observable.

## Fixed-shape source-core force compatibility

At each identical frozen source atomic configuration, evaluate the candidate
and source energy gradients in eV/L0. Mechanical force is their negative;
the RMS difference is identical. Do not splice source moduli into a relaxed
candidate: freezing source exterior positions defines a target configuration,
not candidate constitutive equilibrium. The smaller interior-force evaluator
copies ALL source free-atom coordinates needed outside its inner disk and
retains every affected site; it is tested against the full-disk gradient.

For fixed radial shapes theta, let c=(u,v,A,B,C,D3,D1,D2,DE,K3), where the
pair is u/r^12-v/r^6. Frozen force design D(theta) is exactly coefficient-linear:

    g_candidate = D c, E c = y_exact,
    c = c_particular + S Z z, E S Z=0.

S is a positive coefficient scaling, Z spans the null space. The exact7 bulk/
pristine-interface rows have rank7, leaving3 coefficient directions. The
unconstrained LS residual is the component of g_source-D c_particular outside
range(D S Z); no allowed coefficient choice can do better at THAT theta.
For inner2L0 (18 rows /54 components):

| Calculation | core gradient RMS eV/L0 | other interface squared residual |
|---|---:|---:|
| Unchanged candidate | 0.28883360 | 2215.2032 |
| Signs relaxed, exact7 retained | 0.11323126 | not an admissible material |
| Declared signs + sampled bulk spectrum/tails | 0.11524311 | 17100.5123 |

The signs-relaxed minimizer needs K3<0. With the declared sign constraints the
optimum has K3=0 and positive LJ, but other interface loss increases7.72 times.
Its exact residual1.35e-13 and KKT residual2.8e-14 show this is not a failed QP.
R10/perturbed source states excluded from that loss retain the same mismatch.
The source's own equilibrium gradient is~2e-10eV/L0; candidate R8/R10 RMS
changes only0.00220eV/L0. This is not roundoff or simply the source free boundary.
This is a FIXED-SHAPE incompatibility bound, not proof against the full family.

## Limited joint shape correction, not an adopted material

No new energy term is added. Vary the five existing log radial/saturation
shapes in a declared +/-0.2 neighborhood, keeping exponent0 and the inherited
global bounds. At every shape solve the equality/sign/sampled-spectrum QP.
For the102 remaining interface development observations and54 core components:

    L = L_interface / 2215.203159 + L_core / 4.504941858.

These block scales make fractional-error tradeoffs visible; they are NOT
reported physical uncertainties or a unique material weighting. No yield,
fatigue, desired barrier or desired core width enters the targets.
38 actual profiles took492s. The8-evaluation shape optimizer hit its budget;
the lowest positive trial is not its accepted endpoint. It gives source-core
RMS0.166348 (excluded R10:0.165500, R16:0.164122; perturbed R16:0.170744/
0.170550), but interface loss2473.438 is11.66% worse than baseline. It is NOT
adopted. Existing reference parameter files remain byte-for-byte unchanged.
Saved profile differences at the last accepted point have singular values
1.91376,.882827,.717567,.520339,.242411 (condition7.895); these describe a
scaled numerical local Jacobian, not structural identifiability/statistical CI.
The strongest column correlation is about-.799 (scalar and rank1 log range).

The trial was separately re-relaxed on actual atomic rows. Its first R4/r5
stationary state has force5.73e-9 but minimum curvature-1.229405: a saddle, NOT
a minimum. A negative-mode descent gives a stable core; subsequent R8/r5 and
R10/r7 minima pass full force/Hessian checks. At R10/r7, independent Newton-CG
minimization and full-matrix stationary Newton give E=1.03240616887445 and
1.03240616887443eV, with inner2L0 field difference3.19e-8L0. The latter has
force2.84e-13 and minH=.148568717. Root solving is not physical dynamics.

Crucially, its sampled25–75% span is.87750L0 (R10/r7), barely changed from
.87408L0 in the original candidate, versus2.83L0 in the source. Frozen-force
error improvement therefore does NOT establish a corrected relaxed core.
Its actual eight-state R8/r5 signed load-return protocol remains stable and
returns within1.89e-9/+ and1.35e-9/- L0 of its initial core. This is no resolved
static persistent motion on that boundary, not a calibrated yield prediction.
Additional explicitly source-configuration multistarts preserve the target's
own material and exterior; copying coordinates is not adopting source forces.
The wide source seed returns to the original narrow candidate within1.48e-10L0,
and to the trial's narrow core within2.80e-9L0 (inner2L0). The original first
Newton line search exhausted its reciprocal-mode limit. An independent
L-BFGS attempt hit its declared1000s budget; both failures/checkpoints remain.
An explicitly bound checkpoint recovery then attained force4.36e-11/minH.2162692.
No exhausted series, iteration budget or optimizer message was relabeled as a
stable state without a new force/Morse check.

The actual frozen interpolation from original narrow coordinates to source
wide coordinates has endpoint energy change+.21011882eV(original),+.09774920
(trial),-.15029062(source), under EACH model's own exterior. Directional work
and every coefficient energy sum match independently. In the declared gauge,
the original positive costs areD3 .09973749, radial Eg .08150329 andK3 .03846605
eV, whereas the LJ pair change is-.00617674eV. Term energies are not separately
gauge-invariant observables; they locate a mismatch in this specific parameter
decomposition, not a theorem that all environmental terms must be weakened.
The interpolated path is NOT an optimized transition path or activation barrier.

### Continued and explicitly wider deterministic shape search

The first continuation retained the original log-shape box of radius .2,
the SAME source force rows, seven exact anchors, coefficient signs, sampled
spectral constraints and block-normalized objective. It completed 27 profiles
with an `xtol` termination, but optimality 2.26e-4 exceeded the requested
2e-6 and the saturation upper bound remained active. Its core-force RMS was
.166719 eV/L0 and other-interface squared loss 2470.133. This is a local
numerical termination, not a unique physically identified material.

A separate explicitly recorded search enlarged only that original-centered
log box from .2 to .6; it added NO energy term, strength target or mobility.
Seventy actual profiles took 945.07s. The selected profile is now the accepted
optimizer endpoint, but `max_nfev=15` stopped the outer optimization and its
optimality is 5.98e-4; the saturation upper bound is again active. Save this
as an incomplete optimization with a usable diagnostic candidate, not a
converged Al parameterization. The five shapes plus fixed screening exponent:

    (2.703305047819125, 7.633706329130604, 11.044328645438164,
     172.74231833476296, 6.911784604905995, 0).

In the unchanged coefficient order (u,v,A,B,C,D3,D1,D2,D_Eg,K3):

    (.014208235724468011, .03800279268021126,
     2.095246728098389, 1.2139518504113378, 2.290246746698926,
     14.697901569610833, 2.5621632452682843, 2.1542273322160246,
     1.5545373800178843, 1006.4071116002428).

The seven scaled exact residuals are below 5.95e-14. The inner coefficient
QP satisfies its declared sign/sampled-spectral constraints, but a finite
sample is NOT a proof of all-wavevector bulk stability. The actual source
force RMS falls .28883360 -> .15241668 eV/L0 (47.23%); the other-interface
loss is 2276.471 versus 2215.203 originally (+2.77%). Excluded R10/R16 states
give .15133584/.149745 eV/L0, with perturbed R16 states .155481/.154861.
These were independently evaluated, not inferred from the fitted R8 error;
they are retrospective exclusions, not blind external validation.

Actual saved one-sided finite-difference profiles yield singular values
(1.854919,.915460,.553356,.475618,.247517), numerical rank five and condition
number 7.494. This is a LOCAL numerical sensitivity, not structural uniqueness,
a confidence interval or evidence that the active bound is physical.

The remaining errors are material: a declared registry-saddle curvature is
predicted +.182166 instead of -.186425 eV/L0^2; one excluded opening curvature
is -1.234287 instead of +.903279 eV/L0^2. Near opening 1.17h, normal force is
1.121366 instead of 1.639823 eV/L0. Those sign and magnitude failures must not
be hidden behind the improved joint objective. `wider_probe_validation`
contains ALL observable rows and the separate unadopted snapshot. Actual
relaxed-core checks and their failed attempts are recorded separately; no
snapshot flag is rewritten to pretend they were performed inside the fit.

The existing rank1 factor g(x)=x^z was also checked, not replaced or silently
omitted. At fixed five wider-search shapes, five predeclared z values were
independently coefficient-profiled with the SAME objective and constraints:

| z | joint loss | core-force RMS [eV/L0] | other-interface loss |
|---|---:|---:|---:|
| 0 | 1.306122 | .152417 | 2276.471 |
| -.5 | 1.307321 | .150988 | 2290.634 |
| +.5 | 1.308835 | .153185 | 2276.244 |
| -1 | 1.313410 | .148725 | 2322.136 |
| +1 | 1.314770 | .153423 | 2287.451 |

The z=0 profile was required to replay the preceding independent fit. A
negative exponent slightly reduces frozen-force error at the price of worse
interface error; none of these five probes improves the declared joint loss.
This rules out a simple fixed-shape screening switch as the demonstrated
solution, not all joint six-shape fits or all analytic environmental families.
No new snapshot or relaxed-core calculation is claimed for these extra probes.

### Actual wider-candidate core and stress checks

The wider candidate was then solved, not merely scored at frozen source
coordinates. The near-core R8/r5 minimum has force4.89e-10eV/L0,
lambda_min=.197706701eV/L0^2. Starting instead from the broad source core
first exhausted the Newton trial's reciprocal limit, then a separate L-BFGS
run hit its480s budget. Both failures remain. A declared recovery of that
checkpoint reached force3.75e-11 and the SAME minimum: energy difference
1.10e-14eV and inner2L0 coordinate difference8.13e-10L0. The recovered state
was independently force/Morse checked; the earlier optimizer did not pass.

R8 ring5->7 changes inner2L0 coordinates9.27e-7L0. At fixed ring7,
R8->R10 changes them .00185908L0. The R8/r7 force is6.27e-11 and its Hessian
positive; R10/r7 force3.68e-7 is below its declared5e-7 criterion and its
minimum eigenvalue .1395258 is positive. This separates neighborhood and
domain effects, not an infinite-domain material certificate. The four
saved zero-load state energies were independently replayed with zero error
at saved floating precision in `wider_replay_final`.

The sampled quarter-to-three-quarter disregistry span at R10 is .87783192L0:
almost unchanged from the original .87408147, and still far below the source
R10 value2.79084554. This is NOT a fitted continuum core radius. The original
material's too-sharp profile has not been cured by its47.23% frozen-force
error reduction.

An actual index-one reconstruction saddle has energy .88563937058125eV,
force3.52e-10eV/L0 and barrier .00977853732535eV per straight-row repeat at
R8/r5. This is34.08% below the original same-geometry .01483458876477.
Both actual negative-mode descents reach the two stable endpoints with force
about4.08e-10 and matching endpoint coordinates within6.04e-10L0. The new
barrier is still about29 times the SOURCE's R8 reconstruction barrier, but
these values are NOT full-translation Peierls stresses or finite-loop barriers.

All eight actual wider-candidate states in the signed5/20/50MPa protocol
are stable. The +50->0 and -50->0 states match the original zero-load state
within7.98e-9 and1.03e-9L0, respectively. There is no resolved persistent
static shift on this tested boundary. This neither proves no thermal glide
nor establishes a dynamic hold. The total final comparison contains56 actual
stress states: original, first trial, wider trial (R8), and the source at
R8/R10/R16/R20. The source-only shifted endpoints remain distinct evidence,
not motion that can be attributed to our analytic candidate.

## Finite-boundary translation connection, not an infinite-line yield result

The source R16/r5 +50MPa-to-zero endpoint is independently reproduced by
initializing a one-a2-translated core at zero load and relaxing without fixing
its final position. The two energies differ1.78e-15eV and inner2L0 atomic
coordinates2.98e-9L0. This verifies the SAME finite-boundary metastable state,
not just a jump in the projected triangle location.

For the original and this independently matched endpoint, an analytic full-
matrix stationary solve finds one negative mode and no unresolved modes:

    E_original = 1.390228544224666 eV / straight-row repeat
    E_shifted  = 1.391986573541316
    E_saddle   = 1.392316917715047
    forward barrier = .00208837349038
    reverse barrier = .00033034417373.

Both actual negative-mode descents return to those two stable endpoints.
Saddle force is2.50e-10eV/L0, lambda_min=-.111790169eV/L0^2; the independent
matrix-free eigenpair residual is4.61e-14. Thirteen actual finite-difference
stencils give best gradient-curvature error2.33e-11, while energy curvature
is limited to1.37e-6 by spline/roundoff tradeoffs. All stencils are retained.
The unequal endpoint energies reveal the fixed exterior's image bias. They
must NOT be averaged away or fitted to an assumed Peierls stress. Domain
enlargement is a separate test. No finite-loop rate is obtained from this energy.

For R20, a separately initialized two-a2-translated zero-load source state
matches the +50MPa-to-zero state to 9.33e-15eV and 1.63e-9L0 (inner2L0).
This further checks actual atomic coordinates, not only the projected triangle.
The registry-equivalent coordinate comparison wraps x differences modulo b,
because every infinite atomic row admits an integer atom relabeling. It is
a local-environment comparison, not a measured accumulated specimen slip.
The R16 and R20 protocols retain DIFFERENT translations; the applied exterior
and its image forces still affect the distance travelled. These results
establish a finite-boundary metastable shift in the SOURCE comparator, not
domain-independent Al yield, a Peierls threshold or a specimen residual strain.

## Actual stress and length/energy units

For the declared affine field, with x the line direction,

    F = I + e_x tensor (gamma_y e_y + gamma_z e_z), det F=1,
    P_iJ = (1/Omega) partial E_atom/partial F_iJ,
    P_iJ = (1/Omega) e_z dot sum_R K_,i R_J, J=y,z,
    Omega = b*d*h, [P]=eV/L0^3.

For these two columns Cauchy stress equals Piola stress, since F_yJ and F_zJ
are unchanged unit rows. Convert by(eV_J/L0_m^3)/1e6 to MPa. The applied
strain solves the SAME model's anti-plane elastic2x2 block against(0,tau).
This is nominal linear traction control; actual nonlinear stress is audited,
not silently prescribed to equal its label.

At nominal+50MPa the actual sigma_xz is49.98906634(original),
49.98836746(joint trial),50.00338664(source)MPa. Ring7/11/15 changes in this
shear response are negligible at the shown precision. Energy derivatives
independently match the work conjugacy. Thus a factor1000 strength discrepancy
is NOT an MPa/GPa conversion error in these calculations.

Absolute normal prestress is different: truncated transverse LJ tails leave
sigma_zz=-1.49149/-.416568/-.170452MPa for original ring7/11/15, and
-7.35007/-2.05992/-.842887MPa for the trial. These are **unconverged tail
diagnostics**, not converged physical normal pressures or a new applied load.
The source's own nominal4.05Angstrom state instead has-.272675716MPa, invariant
after its full published cutoff is enclosed. Neither prestress is clipped.
The +50MPa normal reaction increment is-.940541(original),-.996669(trial),
-.036750(source)MPa; it is separately converged over those rings. Local relaxed
core ring comparisons and absolute homogeneous stress tails are different
observables; a small former error does not certify the latter.

## Relation to the probability theory and actual strength

### Infinite transverse LJ force majorant

The free-disk boundary and the pair neighbor-ring truncation are independent.
For the latter a rigorous, deliberately conservative infinite remainder can
be bounded without replacing the Bessel potential by a finite cutoff. Let

    T = [[d,d/3],[0,h]], alpha=sigma_min(T), K=retained integer ring,
    U_perp,max=max_i |(U_y,U_z)_i|, D_perp=2 U_perp,max,
    eta=D_perp/alpha, q0=K+1-eta>0.

The screw exterior has Uy=Uz=0, so this is a global transverse-displacement
bound, not an assumption about unknown exterior coordinates. Integer square
shell k has8k rows and each reference transverse separation is >=alpha*k.
Choose each pair's equivalent x displacement in[-b/2,b/2], exploiting the
EXACT infinite-row periodicity; then |delta U|<=D=hypot(b/2,D_perp).
For phi=c r^-m,

    ||Hess phi|| <= |c|m(m+1) r^(-m-2).

A positive unimodal row sum at any x phase is bounded by twice its maximum
plus its real-line integral divided by b. For nu=(m+2)/2 this gives

    sum_n [r_perp^2+(n b+x)^2]^-nu
      <= 2 r_perp^(-2nu) + I_nu r_perp^(1-2nu),
    I_nu=sqrt(pi) Gamma(nu-1/2)/(b Gamma(nu)).

The mean-value theorem with r_perp>=alpha*k-D_perp bounds the missing pair
force difference from its reference. The missing reference force sums to
zero by inversion on the symmetric omitted rings. Complete all shells with

    T_p=8 alpha^-p [zeta(p-1,q0)+eta zeta(p,q0)],
    B_m=D |c|m(m+1)[2 T_(m+2)+I_((m+2)/2) T_(m+1)],
    ||g_pair,infinite-g_pair,K|| <= B_12(u)+B_6(v).

These are exact upper SERIES sums (Hurwitz zeta), not fitted cutoff constants.
Independent one-center pair-force assembly agrees with the full per-site
coefficient gradient; actual ring3/5/7/9/13 changes lie below the majorants.
At the original relaxed core, ring5/7/13 upper bounds are1.218e-3/2.359e-4/
1.097e-5eV/L0, while the observed7->9 pair-force change is3.24e-7eV/L0.
The bound is intentionally loose. Frozen source-core states are checked
separately. **This certifies only the LJ pair remainder**, NOT all nonlinear
environmental tails, the infinite free domain, material accuracy or yield.
No computed force is corrected/clipped by the bound.

### Probability and specimen closure remain separate

The present step supplies the atomistic energy E({U_i}) and its derivatives,
not a phenomenological plastic-flow law. For a future FINITE physical atomic
domain with independently justified mobility, the same deterministic framework
would be

    partial_t P(U,t) = sum_i div_Ui [ M_i(P grad_Ui E + kBT grad_Ui P) ].

Boundary conditions and probability loss still distinguish registry transport
from opening. Reducing to a measured core coordinate X(U) requires its actual
conditional measure/metric and the collective mobility; simply assigning the
single-cell M_s to a dislocation line is not justified. A straight-periodic
calculation's E per repeat is not a finite total activation energy. The number
of repeats or line length cannot be chosen to create an Arrhenius factor.
Finite kink/loop geometry and its energy must be derived first.

Dislocation motion is also not macroscopic yield: the v21 signed swept-area
observable requires actual source populations, specimen normalization and
accumulated motion. A metastable core shift at a few MPa in an atomistic source
potential supplies one missing mechanism check, not experimental Al0.2%-offset
yield, fatigue or residual strain. The minimal scalar a/s production path and
its uncalibrated kinetics remain unchanged.

There is also no permissible additive shortcut E_total=E_v21_line+E_v22_core.
The atomistic core calculation already contains the same material's elastic
field. A future coupling must match and subtract overlapping contributions
under the SAME boundary, energy, geometry and normalization. A sampled
finite part E_R-K log(R/b) is not an adopted core constant without convergence;
an arbitrary core radius or statistical area cannot replace that matching.

### Next mathematically compatible step: a finite-period line basis

This is a derivation/implementation route, NOT an implemented kink solver.
For N atoms per row superperiod B=Nb, retain a vector U_(i,mu) for each
transverse row i and basis atom mu=0,...,N-1. For an atom (i,mu), a neighbor
(j,nu,k) has vector

    r_(i mu,j nu,k) = R_j-R_i + [(nu-mu)b+kB] e_x
                     + U_(j,nu)-U_(i,mu).

For every basis pair sum k from minus to plus infinity with the same Poisson/
Bessel identities, now at period B. Sum ALL basis/row contributions into
each atom's environmental moments BEFORE applying its nonlinear site law.
The period-cell energy is (1/2) sum_(i mu,j nu,k)' phi(r) plus the per-atom
embedding/angular energies, with a same-boundary reference subtraction.
The prime removes only the true self atom, not all neighbors on its row.

In the uniform-line limit U_(i,mu)=U_i, the exact regrouping identity is

    sum_(nu=0)^(N-1) sum_k f(x+(nu-mu)b+kNb) = sum_n f(x+nb).

Equivalently the reciprocal phase sum over the basis is zero unless its
supercell mode is a multiple of N. Thus this route recovers the current
infinite-row theory exactly; it does not require a finite atomic cutoff.
Each uniform period-cell energy is N times the current per-repeat energy,
not a newly adjustable activation energy.

Coincident transverse rows need separate analytic self handling. For 0<x<B,

    sum_k |x+kB|^(-2p)
      = B^(-2p) [zeta(2p,x/B)+zeta(2p,1-x/B)],
    sum_k exp(-kappa |x+kB|)
      = [exp(-kappa x)+exp(-kappa(B-x))]/[1-exp(-kappa B)].

At x=0 excluding self, these become 2*zeta(2p)/B^(2p) and
2/[exp(kappa B)-1], respectively. The positive-transverse-distance Bessel
implementation must NOT be naively evaluated at zero distance. Distinct atoms
may not be allowed to coincide; no repulsive force is clipped to permit that.

A finite kink-pair excess energy would require an actual nonuniform line
stationary point, two verified endpoints, and period/position/boundary
convergence. Merely multiplying the present straight-row barrier by N does
not supply it. Collective-coordinate kinetics and conversion to seconds stay
additional independent requirements after the energy problem is solved.

## Units, reproducibility and remaining gates

L0=2.863782463805517e-10m, b*=1, E0=1eV. Core E is eV per straight periodic
row repeat bL0, gradient eV/L0, Hessian eV/L0^2. Energies cannot be assigned an
arbitrary line length to manufacture a finite activation barrier. A finite
loop/kink, spatially varying line, real source population and interactions
must be derived/validated before a specimen yield calculation is closed.
M_a,phys/M_s,phys/t0 remain unavailable; no physical a/s seconds/Hz are enabled.
No local probability, A_c, specimen extrapolation or production solver changes.

More explicitly, the numerical coordinate is U*=U_phys/L0, so the saved
gradient number is partial E[eV]/partial U*, interpreted in the force unit
eV/L0; g_phys=(eV_J/L0_m)*g_saved. Likewise H_phys=(eV_J/L0_m^2)*H_saved.
The line energy unit is eV/(b*L0) after dividing the periodic-repeat energy
by its actual repeat, not an adjustable conversion factor.

Results and scripts are under `results/current_material_core_v22/` and the
similarly named runners. `current_replay`/`source_replay` independently replay
saved energies; `combined_analysis` compares actual fields/stress return and
saddles. `force_compatibility`, `constrained_force_profile`, and
`core_informed_shape_probe` distinguish exact lower bound, constrained fit and
limited nonlinear search. `core_informed_validation` contains excluded-state
checks and an explicitly unadopted research snapshot, not a default material.
Incomplete optimizer convergence, negative Hessians and all trial data remain
visible. Final regression timing/status is saved with the study handoff.

## Decision after the actual experiments

| Question | Evidence / status |
|---|---|
| Does the full current energy generate a nonlinear vector core without replacing Bessel? | Yes, in the declared straight-line finite-free-domain research geometry; independent atom sums and derivatives agree. |
| Can source-core force mismatch be reduced without a yield target/new energy term? | Yes,47.23% here, with exact bulk/initial-tangent anchors preserved. |
| Is the core now quantitatively source-like? | No: the relaxed profile remains much sharper; reconstruction resistance is reduced but still too large in this matched comparison. |
| Does it match the whole interface material dataset? | No: finite-opening and registry curvature sign errors remain. The outer fit also has a budget termination and active bound. |
| Does the analytic candidate show persistent low-stress motion in these static branches? | No resolved persistent shift in the tested R8 signed50MPa return. This is not a thermal or specimen-yield verdict. |
| Does the source comparator move? | Yes, independently reproduced shifted finite-boundary states; distance remains boundary-sensitive. |
| Is actual yield/fatigue or a/s physical time calibrated? | No. A finite line/source geometry, source populations/interactions, valid material fit and collective kinetics remain separate missing requirements. |
| Production or specimen UI adoption? | No gate opened; all historical default models and parameter files retained. |

The next calibration work must compare the actual source forces and interface
jets jointly under transparent uncertainties/constraints, not lower a barrier
by hand. More broad deterministic searches may still find a better member of
the existing family; the local failures do not prove its global impossibility.
Any proposed extra analytic radial/environment term needs a demonstrated
compatibility/rank benefit plus new excluded states before adoption. In
parallel, finite-period line-basis derivation above is a way to preserve the
infinite-sum theory when eventually deriving kink/loop activation, not an
already implemented or kinetically calibrated solver.

## Final executed regression and case inventory

The final validation ledger was generated on2026-09-10T23:30:17UTC, after all
atomic studies finished. It records107 atomic case entries, including4 failed
attempts preserved with their independent recovery links, and0 incomplete
atomic cases. An atomic run's completion is not a material-adoption statement;
in particular the outer parameter optimizer's budget termination remains.

| Actually executed check | Result | Elapsed |
|---|---:|---:|
| Targeted new/existing row-core tests |77 passed,0 skipped|48.56s|
| Full solver_v1 suite |638 passed,0 skipped|1195.90s|
| app suite, including real GUI checks |34 passed,0 skipped|60.64s|
| Newly executed desktop smoke |exit0|0.90666s|

`validation_ledger/executed_tests.json` binds the actual JUnit artifacts and
records smoke stdout/stderr. `case_status.csv` retains stability, saddle,
failure and recovery distinctions. `artifact_hashes.csv` binds the saved
scientific artifacts by bytes; it does not certify their physical validity.
The results directory's attributes preserve those bytes across Git checkout.
Existing TwoRowLJ a0=.7713438268704838 and kappa=86.29296488740997 are still
reported by the actual desktop smoke. These passing checks do not open the
Al-material, specimen-yield, finite-loop or physical-a/s-time gates.
