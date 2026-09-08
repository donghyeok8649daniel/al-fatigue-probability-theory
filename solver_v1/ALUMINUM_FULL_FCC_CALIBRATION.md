# Full-FCC Al static calibration: independent-mode audit

## Status

This is a static research calibration, NOT a validated Al fatigue solver.
The initial six-observation experiment was structurally underidentified.
Its numerical rank-five claim is rejected. Current results are in
results/aluminum_full_fcc_calibration/audited_v2 and
results/fcc111_active_interface/audited_v2.

The production TwoRowLJ, reduced hybrid parameters, probability PDE, stress
mapping, mobility, and specimen aggregation are unchanged. No active-interface
surface is admitted to production. Negative scientific findings are retained.

## 1. Reference targets and conditions

The bulk fit uses the 0 K Mishin et al. (1999) EAM reference:

| quantity | value |
|---|---:|
| FCC lattice constant | 4.05 angstrom |
| atomization/cohesive energy | 3.36 eV/atom |
| C11 | 114 GPa |
| C12 | 62 GPa |
| C44 | 32 GPa |

Source: [Mishin et al., PRB 59, 3393](https://doi.org/10.1103/PhysRevB.59.3393);
[NIST potential record](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/).
The model-reference values are not experimental exactness. The source ledger
and original aluminum_reference_targets.csv remain unchanged.

Held-out contextual comparisons:

- Lu et al. LDA-DFT: direct USF 0.250, Shockley USF 0.224, ISF 0.164 J/m2.
- Mishin EAM (111) surface 0.870 J/m2, separation reference 1.740 J/m2.
- Wei et al. DFT (111) surface 1.06 J/m2, separation reference 2.12 J/m2.

[Lu et al.](https://doi.org/10.1103/PhysRevB.62.3099) used atomic and volume
relaxation and a 3.94 angstrom DFT lattice. The
[author preprint, sections 2 and table I](https://arxiv.org/pdf/cond-mat/9903440)
makes those conditions explicit. Our rigid halves at the EAM-scale lattice
are NOT a like-for-like realization of these relaxed GSF targets.
[Wei et al.](https://doi.org/10.1039/C6RA08958E) supply a relaxed-surface DFT
context, not an unrelaxed half-crystal benchmark.

Finite-temperature and room-temperature targets do not enter the 0 K loss.
Vacancy formation remains unrepresentable in this uniform state space.
GSF is not homogeneous shear; cleavage is not uniform plane expansion.

## 2. Physical normalization

The fixed coordinate unit is L0=b_target=4.05/sqrt(2) angstrom.
Numerical energies are eV per bulk atom; E0=1 eV=1.602176634e-19 J.

    h_target/L0 = sqrt(2/3)
    A_atomic,target = sqrt(3) L0^2/2 = 7.1024908428 angstrom^2
    V_target = A_atomic,target h_target = 4.05^3/4 angstrom^3

Thus a strain-energy curvature in eV/atom divided by V_atom, followed by the
eV-to-J conversion, is an elastic modulus. These are not two-row-cell energies.
A_c, specimen stressed area, and element area never enter the calibration.

For imperfect target fits, a pressure-free reference requires isotropic
relaxation of BOTH b and h, not a normal root at fixed b:

    b_model = lambda
    h_model = sqrt(2/3) lambda
    a_lat,eq = lambda a_lat,target
    A_atomic,eq = lambda^2 A_atomic,target
    V_eq = lambda^3 V_target

The potential parameters and density-reference gauge remain fixed during this
strain. Renormalizing the environment at each strained state would change the
energy function and invalidate the derivatives.

## 3. Why the legacy rank claim failed

For [111] normal strain alpha and equal in-plane biaxial strain beta:

    epsilon = alpha n tensor n + beta(I - n tensor n)
    X = C11 + 2 C12

    Haa = V(X + 4 C44)/3
    Hbb = 4V(X + C44)/3
    Hab = 2V(X - 2 C44)/3
    Hbb = 2 Haa + Hab

The three curvatures contain TWO independent cubic combinations, not three.
At a cubic configuration the beta force is twice the alpha force. Consequently:

- legacy stage 1 (two forces + cohesion): structural rank at most 2;
- legacy stage 2 (those + three curvatures): rank at most 4;
- five effective parameters cannot be identified by those observations.

Finite differences produced a small fifth nonzero singular value and apparent
condition numbers 3e7--5e7. This numerical value was not new material data.
The old test asserting rank 5 has been replaced by identity/rank-limit tests.

The first legacy candidate, when independent shear is measured, gives
C11 about 95.22, C12 about 71.39, C44 about 32 GPa. It did not fit all cubic
elastic constants. Historical vectors/CSVs remain with explicit superseded
status; the archival replay cannot overwrite audited_v2 results.

## 4. Missing independent mode

For hydrostatic deformation and engineering simple shear:

    F_hydro = (1+eta) I
    F_shear = I + gamma m tensor n
    m = [1,-1,0]/sqrt(2), n = [1,1,1]/sqrt(3)
    gamma = s/h111

The existing homogeneous stack exactly represents this affine shear because
every adjacent plane acquires the same additional displacement s m. It is NOT
the local active-interface slip coordinate used for a GSF calculation.

At zero pressure:

    H_etaeta = 3V(C11 + 2 C12)
    H_aa     = V(C11 + 2 C12 + 4 C44)/3
    H_gammagamma = V(C11 - C12 + C44)/3

These three combinations are independent and recover all cubic elastic
constants. Cubic positive moduli establish local homogeneous strain stability,
not global FCC-vs-other-phase stability or a full phonon spectrum.

Normal force, normal curvature and shear curvature use the existing analytic
reciprocal derivatives. Isotropic density curvature uses a five-point strain
stencil of the SAME infinite lattice density. Pair isotropic derivatives use
exact homogeneity of r^-12 and r^-6. Stencil and sensitivity steps are refined.

## 5. Model and gauge

The pair remains exactly LJ. The scalar environment remains exponential.

    W_bulk = 4 epsilon (sigma^12 S6 - sigma^6 S3) + F(x)
    x = rho_full / rho_ref
    F0(x) = -A sqrt(x)
    F1(x) = -A sqrt(x) + B(x-1)

Each full-FCC atom uses the total neighbor density BEFORE applying F.
rho0=1 and r_e=L0 fix the redundant raw density scale. rho_ref is the density
at the ideal target FCC geometry for that decay parameter and is held fixed
during strain and interface opening.

F1(0)=-B. Cohesion is therefore -B-W_bulk, not simply -W_bulk.
This atomized-reference correction is explicitly tested.

Linear B is gauge-equivalent to pair redistribution in an unrestricted EAM.
Here the pair functional form is deliberately fixed to LJ, so B is retained
only under that declared hybrid convention. A unique physical pair/embedding
decomposition is not claimed.

## 6. Deterministic staged calibration

Run the actual computation:

    python -m solver_v1.run_full_fcc_calibration_audit

This is not a replay of saved fitted numbers. At fixed density decay d:

    u = 4 epsilon sigma^12
    v = 4 epsilon sigma^6
    y = Q(d) [u,v,A,B]^T

    sigma = (u/v)^(1/6)
    epsilon = v^2/(4u)

The square-root family omits B. Positivity-constrained linear least squares
at each fixed d gives a deterministic variable-projection profile.
Zero coefficients are not clipped into small positive physical parameters.
49 logarithmic d points from 0.35 to 8 and bracketed local-minimum refinements
are evaluated. This interval is a declared search range, not a proof of global
optimality over every possible decay.

Stages: equilibrium/cohesion (underdetermined), then normal/hydrostatic
elasticity, then independent shear. Registry/opening targets stay held out.
All profile points, boundary solutions and stages are saved.

Five independent residuals use:
- force scale 0.10 eV per unit normal strain;
- cohesion scale 0.0672 eV/atom (2%);
- each independent curvature scale 5% of its target.

These are transparent model-discrepancy scales, NOT experimental uncertainty
estimates. The saved objective is sum r_i^2, not half that sum.

Three additional deterministic starts run bounded least squares in the
previous search box (epsilon<=2 eV, sigma/L0<=1.35, A,B<=12 eV, etc.).
Those are search restrictions, not measured physical parameter bounds.
The unrestricted coefficient profile is reported alongside them.

## 7. Fit results and identifiability

| candidate | epsilon eV | sigma/L0 | decay | A eV | B eV | sum r^2 |
|---|---:|---:|---:|---:|---:|---:|
| square profile | 0.00002846 | 1.73784 | 1.56771 | 3.86276 | 0 | 10.924 |
| linear unrestricted | 18.05740 | 0.699155 | 5.176403 | 3.702510 | 55.33476 | about 7e-13 |
| linear historical box | 0.84661 | 0.806966 | 1.36765 | 10.16648 | 12 (bound) | 0.92938 |

The square profile did not fit the independent elastic targets within the
declared 5% scales. This attempted/profiled failure is not a global theorem
that every square-root model is impossible.

The unrestricted linear candidate fits the bulk observables but has strong
energy cancellation:

    pair energy at target = -54.9923 eV/atom
    embedding relative to isolated atom = +51.6323 eV/atom
    cohesion = 3.36 eV/atom
    cancellation ratio = 31.73

It is not justified as a physically identified aluminum parameter set.
The independent log-Jacobian singular values are approximately:

    6797, 1335, 133.7, 3.978, 3.411

with condition number about 1993. Sensitivity-step refinements now retain the
independent smallest modes; unlike the old rank-five result they are not merely
duplicate-observation noise. These local sensitivities do not constitute
confidence intervals or proof of global uniqueness.

The historical-box starts agree on a different compromise with B at its
declared upper bound. This bound dependence is visible, not hidden by tuning.

## 8. Actual relaxed bulk reference

| candidate | a_lat [angstrom] | cohesion eV/atom | C11 | C12 | C44 [GPa] |
|---|---:|---:|---:|---:|---:|
| square profile | 4.04913 | 3.36140 | 91.95 | 60.07 | 49.51 |
| linear unrestricted | 4.05000 | 3.36000 | 114.00 | 62.00 | 32.00 |
| linear historical box | 4.04693 | 3.35830 | 110.46 | 64.21 | 37.31 |
| initial legacy fit | 4.05000 | 3.36000 | 95.22 | 71.39 | 32.00 |

The modulus conversion uses each relaxed atomic volume. Target-geometry
residuals and relaxed-equilibrium values are separate columns in the CSV.
A zero-pressure reference cannot be inferred just from a small fitted loss.

## 9. Active-interface physical validation

The new candidates are probed only as STATIC references. GSF is evaluated
along both explicit paths, and opening is compared with the separated-limit
calculation. Energy per interface cell is converted using the candidate's
actual microscopic area:

    gamma[J/m2] = W[eV/cell] E0 / A_atomic[m2]
    traction[Pa] = W_a E0 / (L0 A_atomic)

Preliminary magnitudes persist under the dimensional/refinement audit:
- unrestricted bulk-exact candidate: separation work about 12.76 J/m2;
- historical-box candidate: separation work about -12.3 J/m2;
- legacy candidate: about 1.146 J/m2.

The negative separation work is a decisive failure: separated surfaces are
energetically favored over that candidate's locally stable bulk. It is rejected.
The unrestricted candidate's bulk success does not validate its large
decohesion energy or its mismatched GSF. Relaxed DFT and rigid-half values
are contextual rather than exact fit comparisons, so no family-wide
impossibility is inferred from them.

## 10. Numerical status and release gate

The audit records separate:
- five-point strain and log-sensitivity step refinements;
- reciprocal/layer/neighborhood tolerances;
- independent direct radii and explicit layer counts;
- active-interface analytic gradient/Hessian checks;
- large-opening asymptote and fixed-registry loaded opening saddles;
- original reduced reference/static parameter regression tests.

Only fixed-s=0 opening saddle/traction limits are currently calculated.
The first normal-curvature zero is followed from the intact branch; one
unimodal traction maximization over a wide interval can miss a narrow local
metastable peak. In particular, the negative-separation-work candidate still
has a small positive local barrier. This does not restore thermodynamic
credibility. Constrained minima with negative registry curvature are NOT
stable two-coordinate intact states; that curvature is saved explicitly.
They are not a full coupled slip/opening spinodal or a dynamic event-ordering
prediction. The fully coupled search is deferred while the physical landscape
fails validation.

Readiness remains FALSE for:
- physically validated full-FCC Al solver;
- active-interface production PDE;
- specimen-mesh probability UI.

No extra polynomial embedding term, mobility tuning, or fatigue fitting was
introduced to manufacture readiness. Next work requires matched interface
targets/relaxation and controlled family assessment. The existing deterministic
probability theory remains intact.

M_a,phys, M_s,phys and t0 remain unavailable. Physical seconds/Hz are disabled.
Atomic calibration does not calibrate specimen correlation area or spatial
field coupling.
