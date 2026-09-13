# v32 — a minimal analytic density-shape diagnostic

This is an isolated **research** extension. The LJ pair, triangular/ABC geometry,
per-site environmental sums, conservative probability generator, stored material
sets and production mobility/clock are unchanged. Source Al99 remains target-only.
The result is not a published or calibrated aluminum EAM potential.

## 1. A specific primitive obstruction, not a general no-go theorem

For positive amplitudes and decays,

    rho(r) = sum_j A_j exp(-k_j r),
    w_j = A_j exp(-k_j r)/rho(r),
    (log rho)' = -sum_j w_j k_j,
    (log rho)'' = sum_j w_j(k_j-<k>)^2 >= 0.

The single exponential has zero logarithmic curvature. The actual source
Al99 auxiliary density at r=2.863782463805517 Angstrom has logarithmic curvature
-2.358450363316104 Angstrom^-2. Thus adding positive pure exponentials cannot
reproduce that radial shape. The source auxiliary density is gauge-dependent,
not an experimentally measured electron density. This identity does **not**
prove that the full energy, including angular terms, cannot be represented by
the original family. It motivates a bounded, separately verified hypothesis.

## 2. Small positive analytic extension

Let y=r-c and

    f(r)=C exp(-k r) [(r-c)^2+w^2], C,k,w>0,
    f'=C exp(-k r) [2y-k(y^2+w^2)],
    f''=C exp(-k r) [2-4ky+k^2(y^2+w^2)].

It is positive at every finite r, smooth, and exponentially decaying. Its log
curvature is 2(w^2-y^2)/(y^2+w^2)^2 and may be negative. This adds a center and
width to the exponential shape, not an arbitrary high-order embedding fit.
The r=0 point is excluded from lattice neighbors, as in the original theory.

If r is in Angstrom, k is Angstrom^-1, c,w are Angstrom and C has units of
auxiliary density/Angstrom^2. C is not an embedding energy coefficient.
For r_phys=L0 r*, use

    k*=k_phys L0, c*=c_phys/L0, w*=w_phys/L0, C*=C_phys L0^2.

For material calculations C is normalized once so the full perfect FCC density
is one. It must not be renormalized during deformation. The independent scalar
shape then has only three parameters (k*,c*,w*). The old scalar exponential
decay retained in the parent constructor is an implementation placeholder,
not an extra identifiable parameter of the new scalar energy.

## 3. Exact transform, no finite-neighbor replacement

With the same 2D Fourier convention as the FCC plane kernel, set
Q=sqrt(k^2+|G|^2). The Yukawa transform is

    integral exp(-k sqrt(d^2+r^2))/sqrt(d^2+r^2) exp(-iG.r) d^2r
       = 2pi exp(-dQ)/Q.

Differentiating with respect to k gives the exponential transform

    T_exp(d,G;k)=2pi k exp(-dQ)(Q^-3+d Q^-2).

Exponential decay justifies differentiation under the integral at k>0.
Since partial_k exp(-kr)=-r exp(-kr),

    T_quad=C [partial_k^2+2c partial_k+c^2+w^2] T_exp.

Poisson summation gives the infinite plane density

    rho_plane(d,delta)=(1/A_atomic_cell) sum_G T_quad(d,G) exp(iG.delta).

ABC phases, all reciprocal vectors including G=0, and opposite-vector pairing
are retained. Registry derivatives multiply each term by iG_j or -G_i G_j.
Normal derivatives commute with the k derivatives. The code evaluates a small
polynomial recurrence of k^h d^j Q^-p exp(-dQ), not finite differences. For a
term before the common exponential,

    partial_k -> h k^(h-1)d^j Q^-p
                 -k^(h+1)d^(j+1)Q^(-p-1)
                 -p k^(h+1)d^j Q^(-p-2),
    partial_d -> j k^h d^(j-1)Q^-p-k^h d^j Q^(-p+1).

The Fourier identity is exact; finite reciprocal/layer evaluation is
semi-analytic. Three small shell/layer envelopes are an **empirical** stopping
criterion, not a rigorous global remainder bound. Independent direct sums and
tolerance refinement are required. The same-plane self-excluded normalization
is an infinite triangular radial series with an analytic Voronoi envelope and
incomplete-gamma remainder, not an unexplained fixed neighbor cutoff.

For the transform derivative of order n, unit conversion is

    T_star^(n)(d_phys/L0,G_phys L0)=L0^(n-2) T_phys^(n)(d_phys,G_phys).

## 4. Same per-atom interface and bulk energy

For every site, sum all neighbor contributions first, x_i=sum_j f(r_ij), then
apply the existing embedding

    F_i=-A sqrt(x_i)+B(x_i-1)+C_emb(x_i-1)^2.

The notation C_emb is distinct from the radial normalization C. Reference
subtraction includes F(0) for cohesion. At the active interface sum
F(x_i(a,u))-F(x_bulk) over both half-crystals; no global crystal embedding and
no sum of separately embedded planes is used.

The new scalar density is also used in the existing per-site coordination
screening x_i^z ||Q1_i||^2. Replacing it in F but leaving the old scalar density
in screening would define a different model and is not what is tested here.
Angular moment shapes remain separate fixed research kernels; all LJ and other
angular contributions are unchanged. F', F'' and the screening product rule
give analytic energy gradients and Hessians.

In pristine centrosymmetric bulk define

    r_h(q)=sum_R [1-cos(q.R)] Hessian_R f,
    r_g(q)=sum_R sin(q.R) grad_R f.

The three embedding Bloch columns are

    H_A=-r_h+(1/4)r_g r_g^T,
    H_B=2r_h,
    H_C=2r_g r_g^T.

Rank1 screening uses x_bulk=1,Q1_bulk=0, so its pristine harmonic column is
unchanged at fixed angular shape. These expressions are independently tested
against the actual six-phase sinusoidally displaced per-atom embedding energy.
The direct bulk validation operator carries analytical omitted-neighbor bounds;
the state-dependent interface uses infinite Poisson sums.

Hydrostatic scalar curvature currently follows the inherited five-point density
strain derivative, with step refinement and an independent analytic direct
radial sum check. It is labeled numerical/semi-analytic, not a new analytic
closed-form elastic derivation. Interface forces/Hessians are analytic.

## 5. Executed staged evidence

Results live in `results/radial_channels_v32`.

1. `reference_density_shape`: 17 fitting radii and 16 midpoint validation radii.
   Only values enter the loss; derivatives do not. Three fixed starts converge
   to the same quadratic minimum. Squared normalized loss decreases from
   116.860461 (single exponential) to 10.786467. Best midpoint normalized error
   is still 1.529582, not a complete within-tolerance match. Scales are explicit
   model discrepancies, not experimental uncertainty.
2. Auxiliary shape: k=3.12148296056 Angstrom^-1, c=1.95976189237 Angstrom,
   w=.256088514351 Angstrom. C=619.016499265 is in the source auxiliary gauge,
   not a physical eV parameter. The four-coordinate value-fit Jacobian singular
   values are approximately 316.667,40.0860,17.7034,1.96883. They do not establish
   identifiability of the material energy.
3. `quadratic_density_material`: use that scalar shape, hold the old angular
   shapes, and refit the energy amplitudes against the **material** targets.
   Seven exact anchors, 115 inspected interface curvatures and four fit-q
   matrices give eta=19.349784, excluded-q maximum relative error=230.8101%.
   Positive LJ coefficients and tested eigenvalues alone do not rescue this
   failed fit. The source-density improvement did not imply energy improvement.
4. `quadratic_density_joint`: a separately declared deterministic three-shape
   search fits those material targets directly, not the auxiliary density.
   Physical search bounds are k in [1,6] Angstrom^-1, c in [0,3] Angstrom,
   w in [.05,2] Angstrom. These are exploratory limits, not calibrated priors.
   Angular shapes are held and amplitudes profiled with the same constraints.
   Completion/budget status and final numbers are in its summary; an interrupted
   checkpoint is not completion. A bound or maximum-evaluation exit is not a
   proof of global inadequacy.

The joint search actually completed 120 evaluations in 317.05 seconds and
stopped at maxfev, not optimizer convergence. Four coefficient LPs failed and
were recorded explicitly; their infinite objective marks an uncertified
evaluation, not proved physical instability or proved mathematical
infeasibility. This driver retained the failure message but not the raw LP
status; that limits the diagnosis of these four evaluations. No finite
clipped penalty replaces them. SciPy's line search emitted warnings when interpolating
infinite values; no convergence claim is made. Only five evaluated profiles
retained both positive LJ coefficients. The lowest objective around16.0416
drives a base LJ coefficient to zero and is ineligible. The best eligible
profile remains the initial eta19.349784. No production bound was changed.

A subsequent independent reconstruction/replay of evaluations12,47,95,115
returned LP status2, HiGHS model_status Infeasible in all four cases. This is
saved in `quadratic_density_joint/failure_replay.json`. It diagnoses the fixed
shape coefficient constraints; no independent Farkas certificate or global
nonlinear-family impossibility theorem is claimed.

`quadratic_density_validation` independently evaluates the four predeclared
fresh interface states and all eight Bloch matrices. Maximum fresh Hessian
relative error is95.05495%; maximum excluded-q error230.81006%. In contrast,
analytic-gradient/FD disagreement is2.2621e-10 eV/L0 and analytic-Hessian/FD
disagreement8.0874e-10 eV/L0^2. The Hessian change on tightening reciprocal/
layer tolerance is5.1514e-14 eV/L0^2. The largest independent finite-q radius
change is8.1112e-5 with omitted-neighbor bound1.9558e-4 eV/L0^2. Thus these
large source discrepancies are not explained by the measured differentiation
or truncation errors. Positive eigenvalues at the eight tested q points do
not prove whole-Brillouin-zone stability.

The validator initially called the bulk context as if it were an interface
evaluator and stopped before producing files. It now explicitly constructs
the target-only vector interface with the same L0. The successful results
above are the actual corrected run, not the failed invocation.

Tests cover Hankel integration, real/reciprocal full jets, units, bulk density
normalization, normal/shear elasticity, hydrostatic refinement, Bloch response,
periodicity, interface finite-difference verification and full calibration-
matrix/actual-energy replay at every nonbulk observation. Tests passing means
the implementation is consistent, not that the material calibration passed.

### Completed joint variation of scalar and angular ranges

`quadratic_density_all_ranges` actually completed320 evaluations in4714.57s:
318 feasible LP profiles,166 positive-LJ profiles and2 explicitly infeasible
LPs. Powell stopped at maxfev, **not convergence**. The best eligible training
eta is15.27014270, versus the unchanged original-family12.50738498. Its scalar
kernel is C=123375.1210, k*=11.95719813, center*=.689670396, width*=.082373474;
C is the fixed bulk-density gauge, not an independently identified material
constant. The coefficient and angular-shape vectors are in the completed JSON.
No lower-LJ-bound solution is adopted.

The actual independent validation contains four newly declared states and
four previously inspected excluded states, labeled separately. Maximum new
Hessian relative error is100.56319%; excluded-q error249.66680%. In contrast,
gradient/Hessian difference errors are3.4359e-10/3.6565e-10 in their eV/L0 and
eV/L0^2 units. Tightening reciprocal/layer tolerance changes H by2.1316e-14;
the largest q-radius change is9.0594e-6 against a2.1844e-5 bound. Minimum
eigenvalue across the tested q matrices is1.438731, not a whole-zone proof.
The large material errors are therefore not explained by these numerical
checks. Varying the previous angular ranges with the new scalar shape does
not make the found candidate usable. This remains a finite search over a
declared domain, not a global impossibility theorem.

### Nested exponential endpoint control (completed v33)

The earlier finite width bound [.05,2] Angstrom does not include the exact
pure-exponential limit. To avoid confusing that restricted search with the
whole quadratic family, a further control uses the SAME decay k and two
kernels independently normalized at the same pristine bulk density:

    f_lambda=(1-lambda) C_exp exp(-kr)
             +lambda C_quad exp(-kr)[(r-c)^2+w^2],
    0<=lambda<=1.

These are mixed at the per-neighbor density level before nonlinear site
energy, NOT a mixture of fitted energies. For lambda>0 this is exactly the
same quadratic family with

    C_eff=lambda C_quad,
    w_eff^2=w^2+(1-lambda) C_exp/(lambda C_quad).

At lambda=0 the original exponential material matrix is used exactly, not a
huge-width floating approximation. No additional radial functional family,
physical correlation length, density gauge or production parameter is added.
The current parent k and angular shape are held fixed. Centers0 and.68432638,
width.08942317 and lambda=.0001,.001,.01,.1,.3,1 are declared before evaluation.
Here lambda is the mixture coordinate; eta elsewhere is the minimax residual.

`run_nested_density_control_v33` actually completed all13 LPs in81.584s.
The exact exponential baseline prediction replay error is0; maximum density
mixture/normalization error is2.00e-15. The smallest positive-LJ minimax error
remains the baseline12.50738498. At lambda=.0001 the two centers give
12.50807594 and12.50905852, respectively; larger sampled mixtures do not
improve it. Some profiles set an LJ coefficient to0 and remain ineligible.
Maximum equality residual is1.10e-12 and dual gap3.58e-13. Thus the tested
nested directions do not resolve the material conflict. This is a bounded
continuation control, not a proof over all centers, decays or angular shapes.
Results are in `results/rank_four_environment_v33/nested_density_control`.
No inspected state is renamed a new blind validation state.

## 6. The three requested physical gates remain separate

### Independent auxiliary pair/gauge check

The actual source EAM has the exact linear gauge transformation

    phi_new=phi_source-2 B rho_source,
    F_new(x)=F_source(x)+B x.

The factor2 follows from the pair half-count. Thus replacing the source pair
by LJ while retaining the source density/embedding up to this gauge would
require phi_source(r)=phi_LJ(r)+2 B rho_source(r). This is a narrowly defined
representability test, not a theorem about all possible new embeddings or
angular environments.

`pair_gauge_shape_bound_audit` fits only source auxiliary pair values at17
radii2.1 through6.1 Angstrom, with16 excluded midpoints and derivatives.
Its discrepancy scale is max(.05|phi_source|,.005 eV), explicitly not an
experimental uncertainty. No target pair function is put into production.

| Radial diagnostic | Squared normalized value loss | Excluded max value error [eV] | Excluded max second derivative error [eV/Angstrom^2] |
|---|---:|---:|---:|
| LJ only | 788.955969 | .103124 | 26.641491 |
| LJ + ideal source-density gauge term | 608.689214 | .078275 | 22.704898 |
| LJ + quadratic-density gauge term | 645.987160 | .082570 | 22.731664 |

The two gauge fits put the attractive LJ coefficient on its optimizer lower
bound (floating values about5e-21 and9e-21 are **not** meaningful retention of
an attractive pair). Optimizer active masks are recorded, not an invented
coefficient floor. `pair_gauge_shape` preserves the first output, whose
positive_LJ flag only tested floating nonzero; use the bound-audited rerun.
Numerical coefficients and losses are unchanged. This auxiliary calculation
shows why merely matching the density and assuming the embedding's linear
gauge absorbs the source pair difference is insufficient. It neither changes
the LJ base nor establishes the impossibility of the full analytic family.

### Material, kinetics and specimen gates

An acceptable common static Al landscape is required before treating core or
finite-source thresholds as material predictions. No failed fit is fed into a
new yield/fatigue run for the appearance of progress. Physical local a/s
mobility additionally needs a justified local PMF/thermal normalization and
resolved low-frequency dissipation. Neither the source auxiliary density nor
its fitted decay supplies time. Existing production seconds/Hz remain disabled;
A_c remains only an external statistical specimen-aggregation parameter.
