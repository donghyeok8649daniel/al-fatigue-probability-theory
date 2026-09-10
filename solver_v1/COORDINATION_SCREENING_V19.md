# Coordination-normalized rank-one response: v19 research protocol

This is a static calibration experiment, NOT an adopted Al potential, a
production PDE change, a yield fit, or a kinetic calibration. The v17
reference force/Hessian mismatch remains the starting evidence. Its independent
range parameter gave negligible improvement. All old parameters/results stay.

## Derivation before fitting

For each atom, retain the SAME infinite exponential scalar and vector sums:

    x_i = sum_j w_scalar(r_ij) / rho_ref,
    Q_i = sum_j w_1(r_ij) r_vec_ij,  I_i = Q_i dot Q_i.

The existing fixed gauges give x_bulk=1 and Q_bulk=0 in perfect cubic FCC.
The v17 site term is D1 I. Test ONE dimensionless shape z in [0,1]:

    g_z(x) = 1 / (1-z+z*x),
    E1_i = D1 g_z(x_i) I_i.

At z=0 this is exactly the old term. At z=1 it is the density-normalized
dipole invariant I/x. Intermediate z continuously mixes a constant reference
normalization with the actual local coordination density. This is an explicit
analytic modeling hypothesis, not a published or calibrated Al screening law.
No coefficient depends on applied stress or desired plasticity. There is no
new atomic cutoff, area, source length, or energy amplitude. The LJ pair,
scalar embedding, ranks2/3 and all Poisson/Bessel kernels are untouched.

For every represented atom x>0, hence the denominator is positive. As an
isolated atom's neighbors recede, I/x tends to zero when 2*k1>k_scalar
(polynomial distance factors do not defeat exponential decay). This condition
must be checked for the tested microscopic ranges; positivity alone does not
prove all-geometry coercivity. z<1 has a finite denominator as x tends to zero.

All environment sums are taken BEFORE the site function, then site energies
are summed. For the symmetric active interface:

    W1_int = 2 D1 sum_depth g_z(x_depth) I_depth.

Do not divide the already-summed interface invariant by a global density.
No subtraction of a large infinite bulk constant is necessary: Q_bulk=0.

For coordinate indices p,q, d=1-z+z*x:

    g=1/d, g_x=-z/d^2, g_xx=2*z^2/d^3,
    I_p=2 Q dot Q_p,
    I_pq=2 Q_p dot Q_q+2 Q dot Q_pq,
    (g I)_p=g I_p+g_x x_p I,
    (g I)_pq=g I_pq+g_x(x_p I_q+x_q I_p+x_pq I)
                 +g_xx x_p x_q I.

These derivatives, not finite differences, define the research evaluator.
Independent finite differences and direct per-site real-space sums test them.

In an affine Bravais crystal Q=0 by inversion at ANY strain, so this term
does not change homogeneous cohesion or elasticity. At the reference x=1,
g=1; the full pristine harmonic finite-q column is unchanged for fixed D1.
At a fixed-material cubic dilation x_bulk need not equal1:

    H1_z(q) = g_z(x_bulk) H1_old(q).

There are no omitted derivative cross terms at harmonic order because Q=0.
Density tail uncertainty must be propagated through g, and direct displaced
site energies must use the same factor. This is not a proof of whole-zone
stability or nonlinear core validity. At an intact interface the original
Haa is also unchanged for fixed coefficients; any fitted improvement must
come from resolving the finite-deformation/coefficient tradeoff, not hiding
a tangent rescale.

## Declared fit and validation

First use the completed v17 nested shape as a fixed numerical CONTROL,
not as an independently measured set of microscopic ranges. Compare z=0,
.25,.5,.75,.9,.99,1 and a bounded scalar refinement on the SAME104 v17 fit
observations. Five bulk equalities remain exact. Separately keep/remove the
inherited positive-LJ pair CONTROL to expose its effect. Do not silently
insert a positive coefficient floor if an optimum loses the LJ pair.

The48 previously inspected v17 excluded observations remain excluded but
are now retrospective validation, not fresh blind data. Declare a new set
of off-grid/off-path states in the runner before fitting; they do not enter
the objective or selection. All scales retain their old units and meaning
(model-discrepancy scales, not experimental standard errors). Source is the
unchanged0K rigid-interface Al99 reference, not an experimental yield target:
[Mishin et al.1999 / NIST original potential](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/).
The tabulated source is target-only. It does not replace our analytic pair.

Save every scalar profile, stopping reason, coefficient/SVD, individual
residual and source/parent hash. Fit data and source are not overwritten.
Actual execution results are added below only after the calculations run.

Length unit L0=4.05/sqrt(2) Angstrom; energy unit1eV; interface cell area
sqrt(3)*L0^2/2 converts to J/m^2 and traction. A_c is never used. Production
M_a,phys/M_s,phys/t0 remain unavailable; seconds/Hz remain disabled there.
No material/core/PDE/UI gate is opened by this protocol alone.

## Follow-up declared after the rational control

The first actually executed rational profiles give fixed-pair optimum
z=.093812, loss1255.701 versus1256.822 at z=0, and a free-pair optimum at
z=0. Thus stronger low-density weighting does not cure the error. The
free-pair one-sided objective derivative at0 is positive. This motivates
testing the opposite sign of density response, NOT forcing plasticity:

    g_p(x)=x^p, -1<=p<=1,
    g_x=p*x^(p-1), g_xx=p*(p-1)*x^(p-2).

This one-shape log-linear interpolation connects the already checked I/x,
the old I, and the simplest density/vector mixed invariant xI. It is a
separate ALTERNATIVE to the rational factor, never their product or two extra
parameters together. All previous chain-rule/harmonic statements apply.
For negative p the isolated-neighbor condition is 2*k1+p*k_scalar>0.
The bounds choose the three minimal invariant conventions, not a universal
material parameter bound. A boundary optimum is not an identified exponent.

Keep the same104 development observations and exact bulk constraints.
The40 rational-probe validation observations are now inspected retrospective
checks; declare40 different state/jet observations before this fit. No old
held-out observation is silently promoted into the objective. Record the
baseline on each exact comparison set. This follow-up is a fixed-range
functional diagnostic, NOT a global refit or proof of family impossibility.

## Joint refinement control declared before execution

The power profiles also produce only negligible improvement. A fixed-range
test alone cannot locate the best nonlinear material surface. Therefore run
a matched deterministic full-shape least-squares comparison: the five inherited
log shapes with power fixed0 versus those same five shapes plus linear power.
Both release the inherited pair CONTROL, retain the same104 fit observations,
five exact bulk rows, and all inherited fixed-material spectral checks. Starts
are the saved fixed-shape free-pair optima; max20 objective iterations per
family, numerical2-point shape differences1e-4, all probes/checkpoints saved.
The inherited windows are k_scalar[.7,8], k_odd/k_even/k1[2,12],
alpha_even[1e-5,1e6]; p[-1,1]. These are numerical search constraints, not
measurements. Record an active bound or evaluation-budget stop without claiming
global or statistical identification. No arbitrarily positive LJ floor is added.

All80 v19 excluded observations have now been inspected and are retrospective.
`validate_coordination_screening.independent_states()` predeclares another10
off-grid states/40jets before joint refinement finishes; they do not enter
any objective. Its sampled-q/dilation tests, site-direct calculations, analytic
derivative refinement and full normal/two-shear MPa scenarios are independent
checks. Static unload return is not a kinetic hold or residual plasticity.

## Executed calibration and physical assessment

All results below were actually computed on 2026-09-10. The reproducible
directory is `results/fcc111_active_interface/coordination_screening_v19/`.
The rational profiles, power profiles and joint runs contain respectively
55, 27 and 178 coefficient solves (260 total, excluding exact replay and
sensitivity probes). The power fixed-range study was rerun into a separate
scratch directory: all profile dictionaries and summary/residual CSV bytes
were identical. The joint optimizer is an actual refit, not saved-vector replay.

| Fit | Additional shape | Normalized squared loss | Normal Haa [eV/L0²] |
|---|---:|---:|---:|
| Rational, inherited pair CONTROL | z=0.0938120 | 1255.70144 | 27.56433 |
| Rational, free pair | z=0 | 1254.65751 | 27.54652 |
| Power, inherited pair CONTROL | p=-0.0835107 | 1255.61056 | 27.56930 |
| Power, free pair | p=0.252580 | 1253.79097 | 27.57020 |
| Joint five-range/shape control, p=0 | fixed 0 | 1238.61101 | 28.05579 |
| Joint five-range/shape + power | p=-1 boundary | 1227.63489 | 27.78327 |
| Matched Al99 source at reference geometry | not fitted | not applicable | 19.66091 |

The joint five inherited coordinates are k_scalar, k_odd, k_even,
alpha_even and k_rank1. The power candidate has these values
`[2.7365634817, 6.2554200402, 12, 426.3611734048, 7.4963989807]`.
Its ten coefficients in the existing `[u,v,A,B,C,D3,D1,D2,D_E,K3]` gauge are
`[0.007755080406, 0.050440701277, 2.005226493261, 1.440242107158,
2.477524510270, 36.892601027718, 2.814640515881, 1.954090139209,
1.611977393846, 33582.79464013082]` (full precision in calibration.json).
Since u=4*epsilon*sigma^12 and v=4*epsilon*sigma^6, the positive LJ pair is
epsilon=0.08201927679 eV, sigma/L0=0.7319266013. Large raw coefficients
cannot be compared without their invariant normalization; no physical
uniqueness or independently measured density range is claimed.

Both joint fits preserve force equilibrium, cohesion3.36 eV/atom and
C11/C12/C44=114/62/32 GPa, with normalized exact-row residual below8.5e-14.
These are enforced bulk targets, not held-out validation. The joint old
control stopped at20 function evaluations (90 actual profiles including
shape differences): **not optimizer-converged**. The power run satisfied
ftol after16 function evaluations/88 actual profiles, with optimality0.710:
this is a small-objective-change stop, not a global optimum or a vanishing
gradient certificate. Saved best feasible probes differ slightly from the
last accepted iterates; both are retained. k_even reaches its search bound12
in both and p reaches -1. Neither boundary is a measured material constraint.
Two of the90 old-family trial profiles have a zero-pair closure and are
excluded from the eligible positive-LJ candidates; all saved final candidates
have a resolved positive LJ pair. All260 profiles' normalized bulk-row
residuals stay below2.49e-13, including the recorded inadmissible closures.
The legacy JSON key `unconstrained_pair_closure_loss` means release of the
inherited fixed-pair CONTROL only; original nonnegative-coefficient and
spectral constraints still apply. It does NOT mean unconstrained regression.

The extra shape lowers the bounded-study loss by0.8862% versus the joint
old control. Its normal curvature is still41.3122% high. Its registry
Hxx=4.444126843 is1.1502% high, but agreement at this one tangent does not
validate the rest of the registry/opening surface. The fixed-range optimum
p>0 switching to a joint optimum p=-1 illustrates a shape/coefficient
tradeoff; it is not physical identification of a density exponent.

### Exact-constraint identifiability

The completed reporter re-evaluated every stored observation from the exact
coefficient vectors (maximum replay error0) and differentiated all five
inherited shapes, not only the added power. Eliminate the five bulk equations:

    M_E(k)c=y_E,
    D=diag(max(1,abs(c))), Z=null(M_E D / scale_E),
    c_,j = D z_,j,
    (M_E D / scale_E) z_,j = -M_E,j c / scale_E,
    y_,j = M_,j c + M c_,j.

Retain the five independent coefficient-tangent columns M D Z plus the
five/six shape columns. Shapes1..5 use log coordinates; the power uses a
linear coordinate. Compare2e-4/1e-4 derivative steps, with a second-order
one-sided power derivative at its bound. The old/power equality-tangent
Jacobians have numerical ranks10/10 and11/11 and condition numbers228.658
and887.527 in this explicitly normalized convention. Singular values are:

    old:   [939.261,225.013,92.0451,59.4802,24.3390,14.0457,
            10.5689,8.66454,6.60179,4.10772]
    power: [927.541,229.750,99.4694,60.5523,24.6730,17.0290,
            12.1975,8.62724,7.15740,3.70527,1.04508].

Maximum Jacobian step changes are2.48e-4 and2.49e-4, respectively. The
extra shape worsens local conditioning. `column_cosines` describes similarity
of normalized tangent effects, NOT a statistical parameter correlation or
confidence interval. Inequality tangent cones/active bounds are not imposed
in this diagnostic; positive numerical rank is not global identification,
physical parameter uniqueness, successful material fit or statistical precision.

### Truly excluded states and actual MPa scenarios

The validator's additional40 jets were excluded from all fits. The source
and each joint candidate were also actually solved under11 tensor scenarios
with9 signed loading/unloading states each:297 local equilibria. They include
three cubic normal axes, three shears, biaxial/mixed loading and explicit
interface normal/shear tractions. All297 local roots and the pristine/fault/
index-one-saddle checks passed; maximum force residual1.762e-14 eV/L0.

| Static measurement | Al99 source | Joint old control | Joint power candidate |
|---|---:|---:|---:|
| Opening change under50 MPa pure interface normal traction [Å] | 0.000925366 | 0.000649758 | 0.000656211 |
| Slip-x under4 MPa pure interface shear1 [Å] | 0.000330998 | 0.000318936 | 0.000327231 |
| Normal-relaxed fault energy [J/m²] | 0.1504795 | 0.1747662 | 0.1739073 |
| Coupled index-one slip saddle energy [J/m²] | 0.1720024 | 0.1977860 | 0.1973711 |
| First fixed-registry maximum traction [MPa] | 12969.568 | 10216.947 | 10194.649 |
| Opening energy at a/h=40 [J/m²], not an exact infinite limit | 1.7412853 | 1.9102495 | 1.9016837 |

Thus the new candidate underpredicts50 MPa normal displacement by29.09%,
despite correctly converting energy, coordinate, atomic area and traction.
The4 MPa shear response is closer (-1.14%). Static return to zero force is
not a time-dependent hold or residual-plasticity test. GPa ideal coherent
traction is **not experimental yield strength** and was not fitted as yield.
Neither the stress scenarios nor the barrier table establish dynamic event order.

Normalized RMS over ten excluded states per field:

| Field | Joint old | Joint power | Power within declared scale |
|---|---:|---:|---:|
| Energy | 1.01729 | 1.01279 | 6/10 |
| Normal force | 5.09213 | 5.16679 | 1/10 |
| Normal curvature Haa | 5.95695 | 5.90727 | 1/10 |
| Registry curvature Hxx | 8.17946 | 8.44712 | 2/10 |

The force and registry-curvature validation actually worsen. A smaller
training objective is not a material acceptance criterion. All128 earlier
excluded observations remain outside the loss but are retrospective after
inspection, not fresh blind validation or independent experimental data.

### Numerical checks and reporting-unit correction

For each joint candidate the independent checks include five fixed-material
cubic dilations(.98,.995,1,1.007,1.02),243 wavevectors and radii12/16;
all sampled eigenvalue-minus-tail margins are positive, minimum0.00560971
eV/L0². This is NOT an all-wavevector or all-strain proof. Direct displaced
site energies verify the dilated harmonic operator for three polarizations.
The new per-site term agrees with a refined independent real-space sum to
3.33e-16 in its unit-coefficient energy. Gradient/Hessian finite differences
are refined at4e-5,2e-5,1e-5; worst final Hessian error3.26e-7 eV/L0².
These errors cannot explain the8.12 eV/L0² normal-curvature discrepancy.

The first opening-traction root agrees between65/129 independent composite
brackets. The source has additional extrema; its second root moves slightly
under bracketing refinement. These recorded roots are not a proof that every
source-spline or coupled instability was found. The new fitted curves each
have one detected maximum in the scanned interval; forces are never clipped.

Final file inspection found a **new joint-run CSV export metadata bug**:
GPa predictions for C11/C12/C44 were paired with the pre-transform mode
names/targets/units. The optimizer itself used the correct transformed
targets and units. The exporter and its sensitivity/report observation
metadata are corrected; a regression recomputes each displayed residual.
Only those six CSV rows' metadata were repaired, not any fitted vector,
physical result or calibration JSON. `final_report/*_audited_residuals.csv`
independently replays the completed vector with the corrected metadata.

## Decision

**Do not adopt this extension into production.** It preserves and verifies
the analytic structure, but the tested local density/vector interaction does
not resolve the normal-force/curvature calibration defect. The fixed-pristine
Hessian identity explains why it cannot independently repair the tangent:
the fit must trade other coefficient contributions against that tangent.
These bounded searches do not prove that the entire analytic family is
impossible, and the old control's evaluation limit remains an open optimizer
qualification. A further material hypothesis must be justified against the
joint force/curvature and directional data, not a desired plot or yield value.
No new field is propagated into production PDE, core dynamics, specimen UI,
mobility, physical seconds/Hz or A_c. The numerical/code tests can pass while
the material validation fails; both facts are retained.

## Actual final regression

Targeted31 passed in48.19s (18new v19 tests plus existing related tests).
After the unit-export correction, full solver555 passed in975.44s,
app34 passed in65.27s; zero failures/errors/skips. The earlier pre-correction
full run also passed553 in614.20s, but does not replace the final rerun.
The final run overlapped part of the sensitivity report; timings are observed
wall times, not an isolated performance benchmark. Desktop smoke exited0 in
2.257s with unchanged historical a0=0.7713438268704838 and
kappa=86.29296488740997. Working and staged `git diff --check` passed.
The existing data/calibration/production PDE/UI/kinetic files are unchanged.
These regressions establish implementation compatibility, NOT acceptance of
the new material hypothesis.
