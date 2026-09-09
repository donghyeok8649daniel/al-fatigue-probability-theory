# Independent environment ranges and the core-matching contract

Research stage `range_core_v12`; no production/PDE or kinetic change.

## Material hypothesis before fitting

The previous actual fits expose a simultaneous bulk/interface conflict. A
single angular exponential range was imposed on all STF ranks1/2/3. Odd ranks
vanish in any affine centrosymmetric bulk; rank2 does not. There is no symmetry
identity that equates their radial ranges. Test the smallest relaxation of
that imposed equality, adding ONE microscopic shape parameter:

\[
 Q^{(l)}_i=\sum_{j\ne i} w_l(r_{ij}^*)\operatorname{STF}[(r_{ij}^*)^{\otimes l}],
 \quad w_l(r^*)={e^{-k_l(r^*-1)}\over \rho_{ref}(k_l)},
 \quad k_1=k_3=k_{odd},\ k_2=k_{even}.
\]

`rho_ref(k)` is the existing fixed perfect-FCC sum at b/L0=1. It is evaluated
once for each range, not refitted to a deformed configuration. Density scaling
remains a gauge; D_l has eV units in this convention. k_l is dimensionless,
equivalent to physical decay k_l/L0. No new mesoscale or probability area is
introduced. Scalar environmental density and the LJ pair remain unchanged:

\[
 E={1\over2}\sum_{i\ne j}\phi_{LJ}(r_{ij})+
   \sum_i\{F(\rho_i)+\sum_{l=1}^3D_l Q^{(l)}_i:Q^{(l)}_i\}.
\]

Every environment is summed before its norm or embedding is applied. Each
rank retains the same analytically differentiated2D Poisson exponential
transform/Bessel lattice sum. Only its decay argument differs. Forces and
Hessians still satisfy

\[
 \partial_\alpha(D_l Q:Q)=2D_l Q:Q_\alpha,\qquad
 \partial_{\alpha\beta}(D_l Q:Q)
 =2D_l(Q_\alpha:Q_\beta+Q:Q_{\alpha\beta}).
\]

At k_even=k_odd the original energy, gradients, Hessian and bulk calibration
matrix are recovered. At fixed three radial ranges the same eight coefficient
columns remain EXACTLY linear. The rank2 column alone changes, including its
bulk contribution. This is not addition of an empirical elastic/yield law.
Angular-environment models are precedents, not transferable calibration;
see [Baskes1992](https://doi.org/10.1103/PhysRevB.46.2727). The present energy
does not implement that published MEAM potential.

The v11 target roles, units and declared5%-C discrepancy are retained. Neither
yield nor source size nor kinetics enters the fit. First test the common-range
limit and sensitivity rank; then perform deterministic nested fits, held-out
state checks and stability/finite-q diagnostics. An improved loss is not
material adoption. New signed coefficients require total stability checks.

## Core/source contract

The v11 gamma(theta)=k_line(theta)ln(R/r_core) retains only outer elasticity.
The next admissible matching object is a core/finite-part energy derived from
the SAME atomistic environment, with explicit subtraction convention. A
change of matching radius must cancel between inner and outer energies:

\[
 E_{inner}(R_2)-E_{inner}(R_1)
 \longrightarrow k_{line}\ln(R_2/R_1).
\]

This identity, its force counterpart and domain refinement must be checked
before eliminating an arbitrary cutoff. The current infinite-row core has
no variation along its line. It cannot by itself certify a curved finite
source, core character dependence or emitted-loop plastic strain. A real
source geometry is independent input, not a parameter chosen to fit yield.
The importance of resolving the core under consistent GSF/relaxation
conditions is illustrated by [Lu et al.2000](https://doi.org/10.1103/PhysRevB.62.3099).
No published Peierls stress is inserted as our yield cutoff.

## Isolated straight-core boundary and exact finite-neighborhood forces

The previous opposite-signed periodic dipole can annihilate. That is not an
isolated dislocation experiment. `isolated_screw_core.py` instead declares a
single straight x-parallel screw, with fixed anisotropic Volterra displacement
outside a finite free disk. The three internal displacement components relax.
The core center is a declared inter-row initial/boundary location, not a
fitted core radius, empirical pinning population or tuned yield threshold.

For the SAME candidate cubic tensor rotated with the audited +ABC frame,
anti-plane symmetry gives

\[
 A=C_{xyxy},\ B=C_{xzxy},\ D=C_{xzxz},\quad
 A u_{,yy}+2B u_{,yz}+D u_{,zz}=0,\qquad
 p={-B+i\sqrt{AD-B^2}\over D},
\]
\[
 u_x={b\over2\pi}\arg[y-y_c+p(z-z_c)],\quad
 u_y=u_z=0\quad\hbox{on the far boundary}.
\]

This fixes one Burgers winding rather than fixing two nearby opposite cores.
Changing a row displacement by b leaves that infinite atomic row equivalent;
no numerical clipping is used to preserve the winding. Under applied shear,
add the homogeneous displacement with

\[
 \begin{bmatrix}A&B\\B&D\end{bmatrix}
 \begin{bmatrix}\gamma_y\\\gamma_z\end{bmatrix}
 =\begin{bmatrix}0\\\tau_{xz}\end{bmatrix}.
\]

This is a static Dirichlet loading experiment. The boundary can exert image
forces. A small nanometre disk is not an experimental micrometre pin spacing.

Row indices have transverse positions
\(y_{jl}=\sqrt3 b(j+l/3)/2,\ z_{jl}=hl\) and reference x offset
\(b(j+l)/2\). The x sum is still the infinite 1D Poisson/Bessel series,
analytically differentiated at the ACTUAL three-component row separation.
Let T_c be its pair, scalar-density and STF channels. For a refined symmetric
transverse neighborhood N, the site changes are

\[
 \Delta t_{ic}=\sum_{r\in N}[T_c(R_r+u_{i+r}-u_i)-T_c(R_r)],\quad
 \rho_i=\rho_{bulk}+\Delta t_{i\rho},\quad Q_{il}=\Delta t_{il},
\]
\[
 \Delta e_i=\tfrac12\Delta t_{i,pair}
  +F(\rho_i)-F(\rho_{bulk})+\sum_l D_l|Q_{il}|^2.
\]

Perfect cubic bulk moments vanish; this is not silently extended to noncubic
backgrounds. Same-row contributions are unchanged and are retained in the
exact bulk density. The finite-neighborhood approximation freezes omitted
ENVIRONMENT CHANGES, not the canonical infinite atomic-row pair potential.
It needs transverse-ring refinement; reciprocal convergence alone is not
all-neighbor convergence. The old periodic energy-tail cancellation is not
claimed for this nonperiodic boundary problem.

If I is the set of free rows, include energy sites E=I+N and all their
neighbors. Fixed sites in E can have changed density and MUST contribute F.
Summing just the free-site energies would omit a genuine force. With
\(w_{ic}=\partial e_i/\partial t_{ic}\), each directed environment bond gives

\[
 g_{i\to j,\alpha}=\sum_c w_{ic}\partial_\alpha T_c(R_{ij}+u_j-u_i),
\quad \partial_{u_i}E\mathrel{-}=g_{i\to j},\quad
 \partial_{u_j}E\mathrel{+}=g_{i\to j}.
\]

Here \(w_{pair}=1/2,\ w_\rho=F'(\rho_i),\ w_{Q_l}=2D_l Q_{il}\).
The analytic Hessian differentiates BOTH w and T. Thus it retains F'' and
all full-environment angular products. Tests compare energy differences,
analytic gradients, Hessian-vector products, symmetry and force sums.

For tensor units eV/L0^3 and b in L0, the circular-annulus energy coefficient
per straight-row repeat b*L0 is

\[
 k_{row}={b^3\sqrt{AD-B^2}\over4\pi}\quad[\mathrm{eV/row\ repeat}],
 \qquad E_{annulus}=k_{row}\ln(R_2/R_1).
\]

Independent angular quadrature verifies this continuum coefficient. The
actual relaxed atomistic site energies still must approach the same annular
slope under disk AND environment refinement. A converged free force alone
does not remove core-radius ambiguity or certify a Peierls/yield stress.
Returning a static unload branch is not a physical-time zero-stress hold.

### Do not equate a raw partial site-energy sum with continuum strain energy

The actual new core check initially found a large annular-energy mismatch
even though the independently calculated row affine Hessian agrees with the
same full-bulk tensor. The local energy partition contains a Taylor-linear
reference term. Let

\[
 t_i^{(1)}=\sum_{r\in N}g_r^0\cdot(u_{i+r}-u_i),\qquad
 g_r^0=\tfrac12\nabla T_{pair}(R_r)+F'(\rho_{bulk})\nabla T_\rho(R_r),
\]
\[
 e_i^{(\ge2)}=\Delta e_i-t_i^{(1)}.
\]

In the perfect lattice \(g_{-r}^0=-g_r^0\). With the complete affected-site
closure E, \(\partial_{u_j}\sum_{i\in E}t_i^{(1)}=0\) for every free row j.
For our screw far boundary only u_x is prescribed and \(g_{r,x}^0=0\) by the
reference row phases. The total linear work is zero up to roundoff. Thus this
is an **exact site-energy redistribution**, not a change of potential,
force, Hessian, minimizer or total fixed-boundary energy. A partial radial sum
does not cancel the same linear boundary work. Comparing its raw value to
quadratic continuum elasticity would be inconsistent.

We now retain and export BOTH raw partial sums and the explicitly
Taylor-remainder partition. An independent finite-difference test verifies
that the total subtraction is constant under free variations. No fitted
energy offset or force repair is used. The straight-core finite-part candidate
in the declared reference convention is

\[
 E_{c,b}(R)=\sum_{|r_i-r_c|<R} e_i^{(\ge2)}-k_{row}\ln(R/b).
\]

Only a plateau under radius, free-domain AND transverse-environment refinement
can support extracting this finite part. Changing the logarithm's reference
from b to r_ref shifts the finite part by exactly
\(k_{row}\ln(r_{ref}/b)\), cancelling the corresponding outer-energy change.
This does not yet provide the core CHARACTER dependence or nonlocal finite
terms needed by the v11 curved source. It is not a finite activation energy.

## Mandatory Morse check and unloaded control

Actual checks found that the centered Volterra seed can relax to a symmetric
STATIONARY SADDLE: force residuals near1e-7 eV/L0 coexist with negative full
core eigenvalues. A force stopping criterion cannot certify a minimum. This
is not repaired by tuning mobility, the potential or the force tolerance.

`core_stability.lowest_core_mode` evaluates the analytic Hessian-vector
operator, its lowest eigenpair and independent two-sided ENERGY curvature
at two refined perturbation sizes. New core runs save this audit even when
it fails. For a resolved negative mode v, execute separate deterministic
initializations u+eta*v and u-eta*v, then relax the SAME full energy with the
SAME far boundary. The sign of v is fixed by its largest component for
reproducibility. eta is a numerical starting perturbation, not a thermal
noise amplitude, empirical barrier change or physical plasticity parameter.

Both lower states and the original saddle are retained. A minimum requires
small actual force AND positive lowest curvature above its eigenpair
residual, plus independent domain/environment checks before physical use.
An applied stress may merely dislodge the unstable initial saddle. Therefore
any residual compared to that original state must be checked against these
ZERO-LOAD relaxed controls, before it can be called stress-induced plasticity.
The corrected physical comparison starts from a validated stable core, not
from an arbitrary symmetric Volterra placement. Static minimization still
does not supply physical kinetic relaxation times or a fatigue lifetime.

## Actual material fit and held-out checks

The loss remains the declared v11 matched-condition loss; the elastic scales
are 5% of each independent source C. The source is the 0 K Al99/Mishin1999
atomistic comparator, NOT experimental fatigue or kinetic data. Its tabulated
potential generates targets only. The LJ pair form is still our pair model.

| Quantity | Source | New separate-range candidate | Status |
|---|---:|---:|---|
| C11 [GPa] | 114 | 87.51655 | -23.23%, fails |
| C12 [GPa] | 62 | 69.99781 | +12.90%, fails |
| C44 [GPa] | 32 | 32.80588 | +2.52% |
| Relaxed fault energy [J/m2] | 0.150479 | 0.151279 | stationary branch checked |
| Index-one stationary energy [J/m2] | 0.172002 | 0.168804 | relocated saddle checked |
| Opening endpoint at 40h [J/m2] | 1.741285 | 1.827637 | finite endpoint, not exact infinity |
| Hxx at the SOURCE saddle [eV/L0^2] | -0.186425 | +0.323297 | held-out curvature has wrong sign |

The Hxx comparison evaluates the SAME source state; it is distinct from the
candidate's relocated index-one saddle. The target roles/scales/provenance
and every residual are in `material/target_definition.json` and `residuals.csv`.
The old coefficient file is unchanged. Neither source size nor empirical
yield strength enters the loss.

Actual deterministic work: 39 fixed starts and 292 coefficient profiles,
438.55 s. The common-range loss 46.7495 improves to 38.9096 (16.77%); the
held-out normalized RMS changes from 12.2948 to 9.4344. This is insufficient
for adoption. The added microscopic range did not solve all material errors.
The 5%-C Powell search exhausted 180 evaluations; no global optimum is claimed.
The exact-C search returned a worse point, so its best previously evaluated
incumbent was correctly retained instead of using the last optimizer vector.

For reproducibility the new decay vector is
`(scalar, odd, even) = (2.57910137477, 5.99653483935, 4.96232918974)`.
The coefficient order `(u,v,A,B,C,D3,D1,D2)` is
`(0.1472597760,0.7173508523,6.1574921554,8.7079955822,1.6195061171,58.15543954,5.711563466,65.73586990)`.
Here u/v are the existing LJ coefficient convention, not new pair functions.
The unchanged hard reference constraints give h/L0=sqrt(2/3) and cohesive
energy 3.36 eV/atom. These exact constraints do not certify other elastic modes.

After removing the two hard reference constraints, the full local
coefficient/log-range sensitivity has rank 9 in the declared coordinates,
condition number 196.70 and singular values
`(272.7721,77.4650,64.2812,57.9784,17.1588,7.4450,3.4825,2.4482,1.3867)`.
Radial FD step 4e-4 to 2e-4 changes entries by at most 1.85e-4. This is local
numerical sensitivity, not global identifiability or uncertainty intervals.
The tied-log-angular tangent is the old nested tangent only when the ranges
are equal; this candidate is not on that submanifold.

The exact-C incumbent reproduces 114/62/32 GPa but has loss 223.0087 and an
independent finite-q instability. At corrected cubic q=(.375,.375,0), its
minimum eigenvalues for direct radii 6/8/10 are approximately
-25.7835/-25.8213/-25.8290 eV/L0^2. A positive uniform elastic tensor and
perfect-interface Hessian missed this failure. The separate-range candidate
is positive at the sampled q points, with minimum 0.008111 eV/L0^2 near q=0;
finite sampling is not a proof for the entire Brillouin zone. **Neither
candidate is an adopted Al material.** The data do not prove impossibility
of every member of the analytic family.

## Actual core tests and the false-residual diagnosis

The core runs use the historical matched angular candidate except the
explicitly named `range_R4_r4_zero`. The far tensor is always calculated from
that SAME potential, not replaced with source elastic constants. Length is
L0=(4.05/sqrt(2)) Angstrom, b/L0=1; static traction is converted through
1 eV=1.602176634e-19 J. No statistical correlation area is used.

For the historical candidate k_row=0.2775806192 eV per row repeat. The
independent infinite-row affine Hessian agrees with its full-bulk tensor to
5.12e-10 eV/L0^3 at ring4. Summed Taylor-linear-work free gradients are near
1e-16 eV/L0, and their total work is roundoff zero. The reported radial
redistribution therefore changes no core energy/force. The remaining
domain/annulus convergence requirement is not waived by this identity.

| Case, zero traction | Minimum Hessian eigenvalue [eV/L0^2] | Finding |
|---|---:|---|
| Centered historical R4/r4 | -0.833110 | stationary saddle |
| Centered historical R6/r4 | -0.742302 | stationary saddle |
| Centered historical R6/r8 | -0.754056 | saddle persists on ring refinement |
| Centered historical R8/r6 | -0.713219 | saddle persists on domain refinement |
| Centered new-range R4/r4 | -1.028191 | stationary saddle |
| Historical R6/r4 mode + | +0.397762 | lower fixed-boundary minimum |
| Historical R6/r4 mode - | +0.397762 | equal-energy lower minimum |
| Historical R6/r6 stable + | +0.404567 | same stable branch, larger ring |

The two R6/r4 lower energies are 1.032538443572 eV/row, below the centered
energy by 0.009221774 eV. They agree within 6.4e-14 eV. Maximum free force is
at most 2.38e-7 eV/L0. Two-sided energy curvature at displacement steps
.002/.001 L0 agrees with the analytic negative/positive modes as refined;
see `core_morse_audit.csv`. A raw force-converged flag on an earlier case is
not a stable-core certificate.

Use the same 18 inner rows (radius <2 L0) for every displacement comparison.
Starting at the R4/r4 SADDLE, the old 25 -> 50 -> 0 MPa sequence leaves a
0.1257949 L0 change. The unloaded state is within 5.3196e-8 L0 of an
independently relaxed ZERO-LOAD negative-mode control. The apparent residual
is explained by relaxation of the unstable starting state, not demonstrated
stress-induced plastic flow.

Repeat from the stable zero-load control:

| Resolved shear traction [MPa] | Max inner displacement change [L0] | minH [eV/L0^2] |
|---|---:|---:|
| 0, stable reference | 0 | 0.641965 |
| 25 | 0.00370464 | 0.654572 |
| 50 | 0.00737592 | 0.664007 |
| 0, after static unloading | 1.01433e-7 | 0.641964 |

Forces are at most 2.05e-7 eV/L0. No resolved persistent displacement is
established on this branch. These are local core displacement differences,
not macroscopic plastic strain; static unloading is not a physical-time hold.
The comparison is not an experimental Al yield prediction at 25 or 50 MPa.

Stable-branch domain R4 -> R6 changes the inner field by 0.013723 L0; at R6,
ring4 -> ring6 changes it by 0.0012741 L0. Thus domain/core matching is still
incomplete. A plausible finite-part curve for the larger UNSTABLE centered
state cannot replace a converged stable branch. All 21 computations are
preserved; nine pass the tested fixed-boundary force/Morse checks. Earlier
states without a Morse test remain explicitly `stability_unchecked`.

## Verification, scope and next scientific gate

Actual commands use an installed Python 3 interpreter (no py launcher):

```
python -m pytest solver_v1/test_core_stability.py solver_v1/test_isolated_screw_core.py solver_v1/test_range_resolved_material.py -q
# 23 passed, 16.52 s
python -m pytest solver_v1 -q
# 394 passed, 631.95 s
python -m pytest app -q
# 31 passed, 148.56 s; no skips
python -m app.desktop_ui --smoke
# PASS, 4.66 s; historical a0/kappa unchanged
```

Final full solver tests used process-local BLAS/OMP thread counts of 1,
without changing tolerances or physics. Independent checks include direct
atomic versus reciprocal rank channels, gradient/Hessian finite differences,
equal-range exact regression, affected-site force balance, and actual negative
finite-q modes. Synthetic Morse tests are labeled synthetic, never MD data.

This step fixes a real stability/initial-state interpretation failure and
improves an explicitly unadopted analytic calibration. It does NOT close:

- simultaneous Al bulk/GSF/vector-curvature/finite-q material validation;
- stable infinite-domain screw core and edge/mixed character energies;
- finite source with measured geometry and emitted plastic-strain criterion;
- collective-coordinate mobility, seconds/Hz, spatial A_c or fatigue life.

The production LJ/Bessel PDE, calibrated-reference files, UI and time status
remain unchanged. Next follow the stable core branch under domain/tail
refinement; independently resolve the material-fit conflict. Do not force
either result to match an empirical yield stress or promote this research
surface merely because all code regression tests pass.
