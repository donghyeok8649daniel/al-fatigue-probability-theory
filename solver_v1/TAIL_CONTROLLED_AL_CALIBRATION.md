# Tail-controlled Al material calibration (v14)

This is an actual deterministic recalibration, not just replaying v13's
rejected coefficients. Production energies/PDE/UI, original parameter files,
physical-time status and statistical correlation area are unchanged. The
target is the provenance-bound, matched **0 K Mishin Al99 static reference**;
it is not an experimental room-temperature yield or fatigue-life fit.

The source remains Mishin et al., Phys. Rev. B59,3393(1999),
[DOI](https://doi.org/10.1103/PhysRevB.59.3393), with
[NIST release notes](https://www.ctcms.nist.gov/potentials/testing/Download/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/2/Al99_releaseNotes_1.pdf).
SHA256 of the target file is
`60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`.
The tabulated potential only supplies comparison data; it never replaces LJ.

## 1. Why calibration previously failed

The previous lowest-loss family could exploit a negative finite-q Hessian.
A constraint imposed at radius8L0 had a small negative mode at radius20L0,
44 times the last refinement change. V13 also showed a separate conflict:
its C'=(C11-C12)/2 was15.14GPa versus26GPa, although C44 was30.66 versus32.
No mobility change, empirical yield cutoff or force clipping can repair those
static energy errors. An actual constrained optimization is required.

## 2. Infinite-lattice tail bound inside fitting

Perfect FCC Voronoi volume is V=b^3/sqrt(2) and its covering radius is
c=b/sqrt(2). For retained neighbor radius R>2c, any omitted site's Voronoi
cell satisfies |x|>R-c and |R_site|>=|x|-c. Thus, for a positive radial f
decreasing beyond R-2c,

    sum_(|R_site|>R) f(|R_site|)
       <= (4pi/V) integral_(R-2c)^infinity (u+c)^2 f(u) du.

For f=u^-p, p>3, this integral is

    (4pi/V) [r^(3-p)/(p-3)+2c r^(2-p)/(p-2)+c^2 r^(1-p)/(p-1)],
    r=R-2c.

For u^m exp(-k u), it is the sum of three upper incomplete-Gamma integrals.
The code checks k(R-2c)>m rather than assuming monotonicity. Phase envelopes
retain the acoustic limit:

    |sin(q.R)| <= min(1,|q|r),
    1-cos(q.R) <= min(2,|q|^2 r^2/2).

A rank-l STF bond derivative has operator norm at most
w(r)[l r^(l-1)+k r^l]. Orthogonal STF projection cannot increase that norm.
If L_R is a retained collective derivative and its tail norm is at most dL,

    ||2 L_infinity L_infinity^dagger - 2 L_R L_R^dagger||
       <= 2[2||L_R||dL+dL^2].

Apply the same square-difference bound to scalar F'' density variations.
Every coefficient column H_j(q) obtains its own error e_j(q). For fixed
ranges a sufficient tested-wavevector condition is

    H_R(q;c) - sum_j |c_j| e_j(q) I >= 0.

All sign combinations of the signed coefficient sector give linear
halfspaces. Negative polarizations supply successive separating constraints.
This is an analytically bounded omitted-neighbor error, NOT an arbitrary
positive margin or the last observed radius difference. Remaining Bessel
normalization roundoff is independently checked; it is not claimed zero.
This condition only covers the tested wavevectors: full-zone searches and
independent radius refinement are separate validation.

The canonical energy stays the infinite Poisson/Bessel evaluator. The direct
force-constant operator is an independent validation/admissibility calculation,
not a replacement finite-neighbor energy.

## 3. Deterministic optimizer and declared objective

Retain the same 5%-per-cubic-constant and10%-interface discrepancy scales,
with.25GPa for zero-force residuals. These are modeling tolerances, not
experimental uncertainties. The old held-out observations stay excluded.
The unchanged family has8 linear amplitudes and3 nonlinear log ranges, with
force/cohesion exact. A separate staged run also holds all3 C targets exact;
that is a staged inverse problem, not a claim of zero uncertainty in Al data.

At fixed ranges c=c0+Zz removes exact equalities. SVD of the normalized design
whitens the objective into Euclidean projection onto linear inequalities.
An active-face solve and nonnegative dual check provide primal/KKT diagnostics.
Redundant sign constraints can have nonunique multipliers: an unconstrained
minimum-norm dual is not an admissibility test. Actual optimizer failures are
saved and do not count as converged starts. Every range evaluation recomputes
the energy jets and spectral columns. No target is clipped or silently removed.

## 4. Smallest symmetry-resolving angular test

For the existing full per-atom STF quadrupole Q_k, cubic symmetry splits its
strain derivative into Eg (tetragonal) and T2g (off-diagonal shear) responses.
The normal [111] strain isolates T2g, since its cubic diagonal is isotropic.
In the plane frame its STF normal response is

    d_alpha Q_k = B_k diag(-1/3,-1/3,2/3).

B_k is SIGNED and evaluated from the existing infinite layer derivatives,
not recovered ambiguously from a squared elastic constant. The current
engineering shear derivative is also computed as a tensor, not just its norm.
Use the two radial ranges already present in v12/v13:

    eta=B_odd/B_even,
    Q_E=(Q_odd-eta Q_even)/N,
    N^2=2||d_gamma(Q_odd-eta Q_even)||^2,
    E_extra=sum_atoms D_E ||Q_E||^2,  D_E>=0.

The normalization only fixes the amplitude gauge. It gives H_hydro=H_normal=0,
H_shear=1 per unit D_E, hence a positive C' increment without changing C44.
There is ONE added fitted amplitude, not another fitted length or orientation.
The full invariant is rotationally covariant; the cubic reference only selects
a radial basis. Environment sums precede per-atom norms and half-crystal sums.
If the ranges coincide the combination vanishes: the unidentifiable gauge is
explicitly rejected, not divided by zero. D_E=0 exactly recovers the old model.

This use of radial moment channels is compatible with the general invariant
construction discussed by [Shapeev(2016)](https://doi.org/10.1137/15M1054183),
but this particular analytic construction/normalization is derived here. It
does not inherit a published Al potential's validation, nor claims of novelty.

Analytic plane jets are combined BEFORE taking the norm. Their derivatives
and direct real-space checks are independent. A separate baseline cache avoids
colliding the two rank2 invariants in the legacy evaluator. Finite-q tails of
the linear combination obey (dQ_odd+|eta|dQ_even)/N.

## 5. Separate minimal shape ablations, not automatic additions

The Eg term enables stable exact-bulk fits but leaves interface error. Two
single-amplitude tests are kept separate; they are not both automatically added.

1. Density-shape test H(x)=(x-1)^3 in the per-atom embedding. Its value, first
   and second derivatives vanish at x=1, but H(0)=-1 MUST enter cohesion.
   Its amplitude is signed. Total-energy stability is tested independently;
   global convexity is not claimed. F',F'' remain exact polynomials. This
   option is not justified merely by a marginal improvement in fit loss.
2. Nonlinear angular test J3=sum_atoms ||Q3_atom||^4 with K3>=0. Because the
   existing odd moment is zero in any affine centrosymmetric bulk, this term
   has zero bulk jets. Since Q3=O(delta q) at a pristine cut, it also leaves
   the perfect-interface Hessian unchanged. It tests the rigid ratio between
   harmonic interface stiffness and finite-fault energy in the old quadratic
   angular term. It is NOT the square of a global sum over atoms or layers.

For I=Q:Q, I_i=2Q:Q_i and I_ij=2(Q_i:Q_j+Q:Q_ij),

    (I^2)_i=2 I I_i,
    (I^2)_ij=2(I_i I_j + I I_ij).

Only after summing each atom's full infinite environment do we form these
expressions and sum over sites. Neither ablation changes LJ, the exponential
density kernel, the probability equation, stress units or kinetic mobility.

## 6. Acceptance is scoped

A stable bulk calibration, an interface fit, whole-material validation and
experimental yield are different accomplishments. Saved tables retain actual
dimensional errors, held-out shape, roots/Morse indices, spectral tails, range
bound hits and Jacobian conditioning. A fitted coefficient vector is not
automatically promoted. The physical energy unit is eV and L0=4.05/sqrt(2)
Angstrom; A_atomic_cell converts eV/cell toJ/m2. A_c never enters these fits.
Static force constants are eV/L0^2, not phonon Hz. Collective-coordinate
mobilities and seconds/Hz remain unavailable without independent kinetic data.

Machine-readable actual optimizations and validation belong to
`results/fcc111_active_interface/tail_calibration_v14/`. Numerical failures and
incomplete earlier attempts must remain explicitly marked, not reported PASS.

## 7. Executed calibration, ablations and bound continuation

These are ACTUAL coefficient/range optimizations, not imported or manually
adjusted coefficients. `run_tail_constrained_calibration` recomputes the
canonical infinite-series observation matrix at every range evaluation.
`validate_tail_calibration` and `report_tail_calibration` only replay/evaluate
completed fits and perform independent checks; they do not refit.

| Family / completed run | Objective | Held-out normalized RMS | C11/C12/C44 [GPa] |
|---|---:|---:|---|
| Unchanged 8-amplitude family, tail-constrained | 31.6274 | 17.3376 | 98.3205 / 68.0351 / 30.6616 |
| Added Eg channel, bulk exact | 52.3792 | 17.1191 | 114 / 62 / 32 |
| Eg + cubic-density ablation, bulk exact | 52.3156 | 17.1033 | 114 / 62 / 32 |
| Eg + per-site quartic rank3, bulk exact | 14.2312 | 18.3918 | 114 / 62 / 32 |
| Same quartic, rational shape ablation | 13.4472 | 18.3871 | 114 / 62 / 32 |
| Same quartic, even-range bound 12 -> 24 | 9.85820 | 34.0228 | 114 / 62 / 32 |

The first row permits elastic residuals whereas the other rows enforce those
targets exactly. Its objective is NOT the same constrained problem. The
original held-out points remain excluded in every optimization. They are now
development/model-selection validation, not an untouched external test set.
The bulk targets are the already documented rounded zero-K Mishin values
(114/62/32GPa and3.36eV/atom), not newly recomputed high-precision elastic
constants from the finite tabulation. Interface states are actually recomputed
from the hash-bound Al99 file at the same nominal4.05Angstrom lattice. Exact
constraint satisfaction means matching these declared targets, not that their
physical uncertainty vanishes.
No whole-family global optimum or impossibility theorem follows from these
finite deterministic starts and bounds. Individual convergence/budget/numerical
stops and all evaluated profiles remain in each completed `calibration.json`.
`successful_starts` in the inventory counts SciPy success flags, not a global
optimum or small projected gradient. In particular the poor quartic second
start stops by xtol with optimality1176.86 and loss207.12; it is not selected.
The selected quartic basin is reproduced by two different starts. Its range
bound and remaining sensitivity/discrepancy still prevent a unique-fit claim.

The cubic-density term buys almost no improvement and is **not selected**.
The rational shape is the analytic per-site function

    H(I) = I^2/(1+alpha I),  I=Q3:Q3, alpha>=0,
    H'(I)=I(2+alpha I)/(1+alpha I)^2,
    H''(I)=2/(1+alpha I)^3.

It nests the quartic at alpha=0, retains zero pristine harmonic contribution,
and is nonnegative/increasing/convex. Three deterministic starts find alpha
about 973, but the held-out improvement is negligible. It is **not selected**
over the simpler quartic. Expanding the even-range bound improves training loss
but again hits the bound and worsens held-out RMS to 34.02. Thus neither an
arbitrary search bound nor lowest training loss is a physical adoption rule.
The range12 quartic is retained for **scoped bulk/static comparison**, not full
material adoption. The continuation results do not establish radial convergence.

## 8. Scoped parameter set and physical normalization

The selected research family is exactly

    E = (1/2) sum_ij [u r_ij^-12 - v r_ij^-6]
      + sum_i [-A sqrt(x_i)+B(x_i-1)+C(x_i-1)^2
               +D1 ||Q1_i||^2 +D2 ||Q2_i||^2 +D3 ||Q3_i||^2
               +D_E ||Q_E,i||^2 +K3 ||Q3_i||^4].

Every sum is per atom with the existing infinite exponential environment;
Q1/Q3 share the odd range, Q2 has the even range, and Q_E is the derived radial
combination. This is an LJ-based analytic angular environmental research
extension, not the baseline square-root-only EAM and not a published MEAM.
All distances/moments above use fixed reduced L0 coordinates. The amplitudes
are energies in eV, x is dimensionless. The LJ functional form is still exactly

    sigma_LJ=(u/v)^(1/6), epsilon_LJ=v^2/(4u).

Ranges (scalar, odd, even):

    (2.7829540002031936, 5.646840689773203, 11.99999999743871).

Coefficients (u,v,A,B,C,D3,D1,D2,D_E,K3), eV in the stated moment gauge:

    (.3265081090541026, 1.5095995103677002, 10.747234706529515,
     16.53828047329518, .2218083326229816, 4.327308845811441,
     8.487797648795736, .7980971835516364, .2647734651412821,
     139933.99432203436).

The large K3 is a coefficient of a small fourth-power moment, NOT an atomic
bond energy of 140keV. Its value is gauge-dependent. Analytic jets at the actual
fitted amplitude are independently checked with tighter lattice tolerance and
finite differences; small-amplitude tests alone would not establish this.
The scoped JSON stores the coefficient convention, source hash and fitted-file
hash, permitted uses and explicit prohibition of production/kinetic adoption.

L0=2.863782463805517e-10 m, a_lat=4.05 Angstrom,
h/L0=sqrt(2/3), and A_atomic_cell=sqrt(3)L0^2/2. Energy/cell is converted to
J/m2 by eV_to_J/A_atomic_cell. Cohesion is energy/ATOM relative to separated
atoms, including F(0)=-B (and -cubic amplitude only in that separate ablation).
No statistical correlation area, element area or fitted activation area enters.

The physical static traction conversion used in the research scenarios is

    f_reduced = A_atomic_cell L0 T_phys / eV_to_J,
    G_int = W_int - f_reduced dot (delta a, ux, uy).

This does not change the production reduced mapping f*=kappa sigma/E. The
three prescribed tractions are interface-normal and two in-plane components;
they are not a claim of a fully spatial 3D PDE or a specified specimen orientation.

## 9. Independent outcome: what is calibrated and what fails

Actual bulk force is about -2.7e-15 eV per reduced strain, cohesion is
3.3600000000eV/atom, and C11/C12/C44=114/62/32GPa. This **is a completed
scoped bulk calibration** to the declared 0K source targets; do not describe
everything as still uncalibrated. A source EAM fit is not independent
experimental Al validation.

The candidate's OWN full (a,ux,uy) stationary points were solved, not merely
sampled at source coordinates. Perfect/fault/saddle have Morse indices0/0/1
and force norms below2e-15eV/L0. Their comparison is:

| Quantity | Al99 source | Candidate | Interpretation |
|---|---:|---:|---|
| Relaxed intrinsic-fault energy [J/m2] | .15047948 | .14922925 | -0.831% |
| Relaxed adjacent saddle [J/m2] | .17200237 | .17912665 | +4.142% |
| Perfect-interface Haa [eV/L0²] | 19.66074 | 24.09191 | +22.54% |
| Perfect-interface Hxx [eV/L0²] | 4.39355 | 3.99449 | -9.08% |
| Fault Haa [eV/L0²] | 14.14055 | 30.17498 | +113.39%; fails |
| Saddle Haa [eV/L0²] | 12.71734 | 29.29311 | +130.34%; fails |
| Rigid work at a=40h [J/m2] | 1.741285 | 1.732383 | -0.51%; finite-separation point |

The final row is not an exact infinity limit; LJ retains an algebraic tail.
The fixed-a direct110 curve has a much larger error than the relaxed Shockley
minima/saddle (see actual curve/held-out CSVs). They cannot be swapped or
described as one uniformly validated gamma surface.

Opening also fails shape validation: the candidate has a traction minimum
-1278.06MPa at a/h=2.578812; the matched source has its own smaller negative
lobe -785.675MPa at2.386631. Do not hide the candidate error by force clipping,
or assert all negative traction is a numerical bug when the source has one.
Both129 and257 independent curvature-root brackets found the same extrema.
The first fixed-registry normal traction peak is9.11682GPa at1.294730h,
versus source12.96957GPa at1.185846h. These are ideal homogeneous cut tractions,
NOT measured yield strength or a coupled finite-defect instability.

155 independent irreducible-wedge q points at radii12 and16L0 have positive
minimum eigenvalues above the analytic tails. The smallest sampled robust
margins are .00620815 and .00642097eV/L0² (near acoustic q). Three continuous
searches of lambda_min/|q|² converge near the L point: H_min=14.49429 with
tail<=.00073904eV/L0². This fixes the previous observed negative-mode problem;
it is **not a proof at every q**, finite strain, or every defect configuration.

The force/cohesion-eliminated local sensitivity has11 directions, singular
values approximately

    442.513,280.105,112.097,25.2355,22.3063,15.2407,6.53836,
    3.73701,3.29563,1.61245,.226276,

with condition1955.64. Halving log-difference step2e-4 ->1e-4 changes a
Jacobian entry by at most5.12e-5. This local numerical rank is not global
identifiability or confidence; the range bound and model discrepancy remain.
Column correlations/right singular vectors are saved, not replaced by a
claim of a unique physical Al parameterization.

## 10. Actual stress checks and readiness

32 source/candidate static states cover pure shear, mixed normal+two-shear,
compression+shear and unloading, with stresses from -150 to150MPa. Candidate
force residual<=1.64e-14eV/L0; unloading registry displacement<=1.61e-16L0.
The interface stays locally stable and returns to its pristine registry.
This is reversible static interface response, not a dynamic zero-stress hold,
new interwell plasticity, experimental0.2% yield or fatigue. No PDE was run
using the unaccepted material. Existing core results were not relabeled as
results of these newly calibrated coefficients.

The current classification is **bulk-calibrated, interface-incomplete analytic
research candidate**. Good relaxed fault energies do not outweigh wrong vector
curvatures/direct-path/opening shape. The missing mechanism is not repaired by
inventing mobility, slip-source length, characteristic area, empirical yield
or arbitrarily adding many polynomial coefficients. Next material work should
use these full-vector held-out failures to test the smallest derived angular
shape alternative and require independent data, not optimize a desired yield.
Any new targets taken from this validation must be relabeled development data.

Production TwoRowLJ/reduced hybrids, static historical parameter files, PDE,
kinetic calibration JSON, UI defaults and A_c semantics are unchanged. Actual
M_a_phys/M_s_phys/t0 and physical seconds/Hz remain unavailable.

## 11. Executed regression checks

New targeted tests:22 passed in67.94s. Full `solver_v1`:437 passed in792.45s.
`app`:31 passed in108.85s, no skips. Actual desktop startup smoke passed in
2.17s and reported the unchanged TwoRowLJ a0/kappa. The Windows py launcher
was unavailable; the verified Python3.13 interpreter ran equivalent commands.

All current calculations completed before this handoff. Final JSON/CSV/SVG
parsing, diff inspection and Git remote verification are separate release checks;
see CURRENT_WORK_HANDOFF.md and the final commit record. The scientific status
is not changed to full-material PASS merely because regression tests pass.
