# Matched-condition calibration and physical scenario audit

The previous bulk-only candidates were not accepted because a small bulk
residual did not establish valid fault or surface physics. In particular one
had negative work of separation. This study attempts a repair without
replacing the analytic LJ pair, the infinite Poisson/Bessel representation,
or the production probability theory.

Status: STATIC RESEARCH. No candidate is admitted to the probability PDE on
the basis of a least-squares loss alone. See the generated matched_v3 tables
for actual runs and individual errors, not just best-fit parameters.

## 1. Fix the comparison before adding physics

The earlier Lu DFT targets used atomic relaxation and a different lattice
constant. Our rigid half-crystals were not the identical physical experiment.
We now evaluate the published Mishin Al potential on exactly the same rigid
FCC(111) geometry, a_lat=4.05 angstrom, at 0 K. This is a mapped ATOMISTIC
MODEL benchmark, not experimental truth or a new claim of a published Al fit.

Source: [Mishin et al., PRB 59, 3393 (1999)](https://doi.org/10.1103/PhysRevB.59.3393),
[NIST record](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/),
[NIST source file](https://www.ctcms.nist.gov/potentials/Download/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/2/Al99.eam.alloy).

NIST credits the setfl conversion to Chandler Becker (2008) from files
provided by Yuri Mishin. This repository uses that file only to produce
independent target observations. It is NOT an energy option in our PDE.

SHA256 of the checked release:

    60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284

The source's published cutoff defines that reference potential. It does not
truncate our canonical LJ or exponential infinite sums. We separately compare
cubic and linear interpolation of the source. The source evaluator performs
per-atom density aggregation before F and excludes self interactions.

Reproduction, using an available Python 3 interpreter:

    python -m solver_v1.reference_eam_targets --download
    python -m solver_v1.run_matched_interface_study --phase calibration
    python -m solver_v1.run_matched_interface_study --phase scenarios

Download is explicit, checksum-verified, cached outside Git, and refuses to
overwrite a different file. Calibration actually optimizes again; scenarios
read the saved fit vectors transparently. Neither command changes defaults.

## 2. Matching units and reference observations

    L0=4.05/sqrt(2) angstrom, E0=1 eV,
    h/L0=sqrt(2/3), A_atomic,target=7.1024908428 angstrom^2,
    gamma[J/m2]=E[eV/interface cell] E0/A_atomic,
    T[Pa]=W_q* E0/(A_atomic L0).

For an imperfect fit, the bulk is relaxed isotropically; a and b change
together. The potential, L0, and density-reference gauge do not change. Use
the actual candidate's area and volume for J/m2 and GPa. The stored physical
lattice error is not silently removed by redefining units.

Mapped source observations calculated in this study:

| observation | reference |
|---|---:|
| Cohesion | 3.359999988 eV/atom |
| Rigid (111) work of separation | 1.741285309 J/m2 |
| Rigid direct110 maximum | approximately 0.60303835 J/m2 |
| Rigid Shockley maximum | approximately 0.18971119 J/m2 |
| Rigid intrinsic SF | 0.156630432 J/m2 |

The two maxima are curve diagnostics, not independently prescribed fitting
locations. The fit instead uses specific displacement states: direct half
registry, Shockley half partial, Shockley endpoint, 20% interface opening,
and full separation. Other curve points and stress branches are validation.
Endpoints used in fitting are not included in the held-out curve RMSE.

The existing rounded Mishin bulk targets are retained: lattice 4.05 angstrom,
cohesion 3.36 eV/atom, C11=114, C12=62, C44=32 GPa. They are 0 K atomistic
targets, not mixed with room-temperature measurements. Source interpolation
changes the reported interface quantities by less than about 3e-7 J/m2 in
the tested curves; this is tiny compared with the model's fault mismatch.

## 3. Gauge-fixed deterministic objective

At a fixed scalar exponential decay, define

    u=4 epsilon_LJ sigma_LJ^12, v=4 epsilon_LJ sigma_LJ^6,
    x=rho/rho_ref,
    F(x)=-A sqrt(x)+B(x-1)+C(x-1)^2.

rho_ref is fixed at target FCC. Isolated atoms have F(0)=-B+C. Cohesion must
include this term, not just negate the bulk energy.

All ten observations are coefficient-linear at fixed decay:

    y=Q(kappa)[u,v,A,B,C]^T,
    L=sum_i [(y_i-y_target,i)/scale_i]^2.

Scales are explicitly declared MODEL DISCREPANCY allowances: .10 eV per unit
normal strain, 2% cohesion, 5% each independent elastic curvature, and 10%
each interface energy. They are not measured standard deviations; no
statistical confidence interval is inferred from them.

Primary profiles use u,v,A,B,C>=0 and explicitly record boundary solutions.
The square family fixes B=C=0; the linear family fixes C=0. A separate signed
B,C exploratory solve tests whether failure is just a convex-sector artifact.
An optimum approaching v=0 is not converted into a positive LJ parameter by
clipping and is not accepted as a physical parameter set.

Density-decay grids, deterministic starts, all profile values and local solver
statuses are saved. A failed/worse local refinement does not erase a better
grid candidate. Search bounds are not physical parameter confidence ranges.

## 4. Scalar-family outcome and what was repaired

The actual joint study found:

| family | sum of normalized residuals squared | conclusion |
|---|---:|---|
| square-root only | 300.3876 | interface and bulk mismatch |
| square-root + linear | 121.38024 | separation repaired, ISF still wrong |
| additional convex C term | 121.38024 | C=0 selected; no improvement |
| positive two-exponential density | 121.38024 | returns to effectively single decay |
| positive squared envelope, extended range | 109.52176 | other targets improve, ISF still wrong |

The linear candidate's actual zero-pressure lattice is 4.027424998 angstrom,
cohesion 3.357395388 eV/atom, and C11/C12/C44 approximately
108.90/70.63/51.58 GPa. Its separation work is approximately 1.80754 J/m2,
an improvement over the negative or extremely large bulk-only predictions.
Its intrinsic SF is approximately -0.00036268 J/m2, however, versus the
matched +0.15663043 reference. It is not accepted as an Al slip landscape.

The independent normalized log sensitivity for (u,v,A,B,decay) has condition
approximately 50.23 in this joint observation set. This is better local
information than the old underidentified fit, but cannot cure model bias.
It is not a confidence interval or proof of global uniqueness.

## 5. Observable compatibility, not repeated parameter guessing

At each fixed decay a linear program maximizes the intrinsic-fault energy,
subject to all OTHER targets being within three declared discrepancy scales.
Zero coefficients are even allowed, enlarging the feasible set. For example,
at decay 2.8 the upper bound is about -8.41e-5 eV/cell, whereas the source
requires +0.06943468 eV/cell. Higher sampled decays can make the other
constraints infeasible already. All sampled outcomes are saved.

This certifies incompatibility at those fixed decays in the stated coefficient
sector, NOT over a continuous unbounded parameter domain. The signed-C study
and enlarged shape search are reported separately. We do not claim a global
impossibility theorem from a grid or an optimizer exit code.

Two positive exponential components add their DENSITIES before F. A second
minimal positive envelope is q(1-d q)^2, q=exp[-kappa(r/L0-1)]. Its exact
expansion uses the same Poisson sums at kappa,2kappa,3kappa. These tests retain
analytical structure and do not introduce a finite-neighbor production
potential. They did not resolve the intrinsic-fault problem in the tested
ranges and are not promoted merely because they have more parameters.

## 6. A separate angular research hypothesis

The first-plane symmetry cancellation and the failed radial profiles motivate
testing environmental directionality, not changing the LJ pair. The minimal
odd third-moment invariant is derived explicitly in
ANALYTIC_ANGULAR_ENVIRONMENT.md. At fixed radial ranges it adds one linear
coefficient, D, with exact zero contribution to homogeneous centrosymmetric
FCC bulk. It is not an empirical slip law.

Initial matched fits lower the objective to approximately 44.37 when radial
ranges are tied and 28.29 with an independent angular range. They make the
ISF positive, but still compromise other observations. This is improvement,
NOT automatic acceptance. Positive bulk and surface energies alone do not
validate a complete curve, full Hessian, spinodal, or first-passage dynamics.

## 7. Stress scenarios and interpretation

Each implemented candidate and the source reference is actually evaluated in:

- pristine equilibrium and small compression;
- normal tension from 0.05 to 10 GPa, including loss-of-branch checks;
- compression control down to -2 GPa;
- pure resolved shear from 0.05 to 4 GPa;
- an EXPLICIT example axis halfway between plane normal and slip direction;
- load/unload to zero on the tracked quasistatic branch;
- both direct110 and Shockley112 registry paths, partial opening and separation.

The multi-GPa cases are idealized mechanism/stability probes, not ordinary
calibrated Al fatigue conditions. Orientation is not inferred. Normal traction
and resolved shear are stated separately. For an explicitly supplied unit
loading vector e, normal n, and in-plane direction m,

    Tn=sigma(e.n)^2, tau=sigma(e.n)(e.m),
    G=W-A_atomic[Tn delta_a_phys+tau s_phys].

This conjugate work is for the separate interface reference and does not
overwrite the verified reduced-model f*=kappa sigma/E convention. A_c does
not occur in the calculation.

A stable root must have a small force residual and positive two-coordinate
Hessian. Load continuation uses substeps and stops on failure, without
jumping to a distant well and calling it plastic flow. An unresolved root is
NOT declared a spinodal. Positive scalar-path Hessian also does not establish
stability in the omitted transverse registry direction.

Quasistatic unloading is NOT a dynamic zero-stress hold or evidence of
persistent plasticity. Static barriers do not establish dynamic event ordering.
Physical M_a, M_s, t0, seconds, and Hz remain unavailable. No PDE run on a
failed interface surface is presented as calibrated fatigue dynamics.

## 8. Output and acceptance policy

All outputs are under results/fcc111_active_interface/matched_v3, separate from
the preserved audited_v2 and old calibration records. Fit tables retain units,
source, scales and individual residuals. Full curve/stress and derivative
comparisons are saved; figures show actual calculated data.

The original reduced models, their static fit files, opening probabilities,
plasticity bookkeeping, specimen aggregation, and UI remain unchanged.
Solver readiness requires physically defensible scenarios and independent
convergence, not visual strain amplitude or a successful optimizer message.
Static-interface adoption and any UI workflow redesign remain separate gates.

## 9. Nested results after the first scenario failure

The I3 candidate's local normal tangent was 188 GPa versus 126.45 GPa in the
matched source. This prompted the explicitly derived vector companion I1,
not a mobility adjustment. Eleven observations now include the source local
normal curvature; that observable is no longer held out. The I1 coefficient
is allowed to be signed and its TOTAL energy stability must be checked.

| nested study | objective | important limitation |
|---|---:|---|
| I1+I3, soft reference geometry | 33.87411 (11 observations) | bulk C44 54.37 vs 32 GPa; LJ epsilon 8.81 eV |
| I1+I3+C, exact geometry/cohesion | 33.10337 | attractive LJ coefficient numerically zero; rejected |
| I1+I2+I3, exact geometry/cohesion | 29.70446 | again attractive LJ coefficient zero; rejected |
| I1+I2+I3+C, signed B/D1/D2 sector | 29.19092 | finite LJ, but large elastic mismatch and opening overshoot |

Objectives with ten and eleven observations or with exact equalities are NOT
directly interchangeable scores. C is the previously declared convex scalar
embedding extension, not an extra empirical strain-evolution parameter.

Why I2? Both odd moments vanish under all affine centrosymmetric strains,
so they cannot repair bulk anisotropic elastic mismatch. The traceless second
moment vanishes in cubic FCC but responds to anisotropic strain. Its analytic
derivation and independent real-space verification are in the angular document.
It shares the angular decay; no new radial-range parameter was added.

The final row gives, at its actual cubic equilibrium:

    a_lat=4.05 angstrom, cohesion=3.36 eV/atom,
    epsilon_LJ=.41127127248 eV, sigma_LJ=2.30447390855 angstrom,
    scalar decay=1.93889283776, angular decay=4.87723238255 [per L0],
    C11=86.43583, C12=64.68436, C44=50.67033 GPa,
    intrinsic fault=.13410878077 J/m2 (reference .15663043247),
    local normal tangent=147.82995 GPa (reference126.44968),
    local registry tangent=30.38239 GPa (reference28.25751).

Held-out direct/partial GSF curve RMSE is .04689350/.01349077 J/m2, each over
46 points with training fractions excluded. The full signed log-absolute
Jacobian condition is 4685.63; the exact-geometry/cohesion tangent condition
is 1833.91. These are local inverse-curvature diagnostics, not experimental
confidence intervals or proof that the extra parameters have unique meaning.

The broad opening test found an ADDITIONAL failure: at a/h=3, W=1.88623 J/m2
and traction=-.65538 GPa although the asymptotic separation is about1.7502 J/m2.
It overshoots then descends, unlike the matched reference. Positive local
Hessian and positive separation energy had not excluded this behavior.
Accordingly the candidate is not adopted. A separate deterministic
`monotone_opening_calibration` imposes W_a>=0 at declared a/h values using the
same coefficient-linear analytic derivatives. It does not clip forces or
change the fit targets/weights. Sampled inequalities must be independently
refined; they are not a proof of monotonicity between sample points. Oscillatory
disjoining pressure is not ruled out for every imaginable material, but here
it conflicts with the selected matched atomistic reference.

## 10. Reproduction of all additional checks

    python -m solver_v1.odd_moment_calibration
    python -m solver_v1.run_odd_moment_quadratic_probe
    python -m solver_v1.run_constrained_odd_calibration
    python -m solver_v1.even_moment_calibration
    python -m solver_v1.even_moment_calibration --convex --signed --free-linear
    python -m solver_v1.run_matched_interface_study --phase odd-scenarios
    python -m solver_v1.run_matched_interface_study --phase even-scenarios
    python -m solver_v1.run_extended_identifiability
    python -m solver_v1.run_static_bulk_stability
    python -m solver_v1.run_coupled_interface_grid --family even_moment
    python -m solver_v1.monotone_opening_calibration

The first five/calibration commands actually optimize. Scenario commands
explicitly read saved fit vectors and execute fresh force/energy calculations.
Registry-symmetry and layer-decomposition checks are reproducible with
`python -m solver_v1.run_registry_symmetry_audit`.

The exact radial first-plane cancellation is an important negative theoretical
result: the matched source obtains about .06941008 eV of its .06943468 eV
intrinsic fault from its pair radial shape, while the tested smooth LJ/scalar
hybrid gets a tiny negative farther-plane contribution. This is NOT permission
to replace our LJ with that tabulated pair. It identifies what a new analytic
environmental description must reproduce without breaking the base theory.

Static stress continuation and finite-q stiffness tests do not establish
plastic flow, first-passage probabilities, residual plasticity, or time/Hz.
Nor do these finite family/parameter searches prove that EVERY possible
analytic LJ-based many-body model is impossible. They locate specific
compatibility failures and preserve all evidence for the next research step.

## 11. Final constrained candidate and honest readiness decision

The coarse opening inequality grid missed a negative interval near a/h=3.5.
The added regression reproduces this failure using a synthetic analytic curve,
without inserting that curve into material physics. The actual refinement
solves W_aa=0, adds the worst force minimum as another coefficient inequality,
and repeats. Four exchanges at fixed radial ranges reduced the minimum force
to -1.03e-14 eV per reduced coordinate, below the measured arithmetic/reciprocal
refinement scale. Raw signs are retained. 80/160 root-bracket grids locate the
same stationary forces on 1.001<=a/h<=12; this is not a global interval proof.

The physically constrained objective is29.47407. No fitting weight, mobility,
temperature, LJ functional form, external-work convention, or target value
was altered to obtain this result. Final parameters and all constraint
iterations are in monotone_opening_refined.json; earlier failed fits remain.

| quantity | matched reference | final constrained research candidate |
|---|---:|---:|
| Lattice constant [angstrom] | 4.05 | 4.05 |
| Cohesion [eV/atom] | 3.36 | 3.36 |
| C11 [GPa] | 114 | 87.81617 |
| C12 [GPa] | 62 | 64.83263 |
| C44 [GPa] | 32 | 49.27105 |
| Intrinsic fault [J/m2] | .15663043 | .13100554 |
| Work of separation [J/m2] | 1.74128531 | 1.76301375 |
| Local normal tangent [GPa] | 126.44968 | 147.29486 |
| Local registry tangent [GPa] | 28.25751 | 31.80526 |

The LJ coefficients remain finite: epsilon=.15445751473 eV,
sigma=2.45246438877 angstrom. Signed log-absolute local Jacobian condition
is1093.33. Its two-bulk-equality tangent condition580.08 excludes the additional
active opening inequalities and is not a complete constrained confidence region.

Held-out direct/partial GSF RMSE is .05270159/.01601694 J/m2 over46 non-training
points each. Compression, tension, resolved shear, explicit mixed orientation,
and static unload were actually re-run on this candidate. Static unload
returns on the intact branch; it is not a residual-plasticity experiment.
The sampled finite-q minimum remains positive: .08074197 at12 L0 and
.08074287 at20 L0. No claim is made for the entire Brillouin zone.

This repairs significant thermodynamic/path inconsistencies and provides a
mathematically explicit static research surface. It STILL fails the declared
independent Al elastic targets and has held-out curve error. A nearly zero
traction plateau is not established by the source data. It is therefore not
a quantitatively Al-calibrated potential, not a finished fatigue solver, and
not admitted to the production probability PDE or a new mesh UI workflow.

Final validation and source/parameter hashes are saved under matched_v3.
New targeted tests35 passed; full solver212 passed; app27 passed,2 skipped;
desktop smoke passed. These numerical outcomes are separate from the material
non-adoption decision. Source-cache-dependent tests were actually run with the
verified source in this execution; on an offline fresh clone they explicitly
skip until the source is provided, rather than inventing values.

For the final extra checks:

    python -m solver_v1.run_monotone_opening_refinement
    python -m solver_v1.run_matched_interface_study --phase monotone-scenarios
    python -m solver_v1.run_static_bulk_stability --extra-monotone-radius 20
    python -m solver_v1.run_coupled_interface_grid --family monotone_opening
    python -m solver_v1.run_matched_research_summary
