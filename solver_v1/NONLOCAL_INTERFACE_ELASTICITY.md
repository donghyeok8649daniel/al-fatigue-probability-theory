# Nonlocal interface elasticity and the low-stress defect question

Status: a **static long-wavelength research reference**, not a production
fatigue solver, atomically resolved dislocation core, or accepted Al fit.
The saved `angular_monotone_opening` LJ/Bessel surface is unchanged. This
work introduces no fitted stress threshold, kinetic mobility, or activation
area. It follows the outstanding normalization issue in
[LOW_STRESS_CYCLIC_AUDIT.md](LOW_STRESS_CYCLIC_AUDIT.md).

## 1. The missing energy is spatial, not a mobility adjustment

The previous W_int(a,s) translates an entire infinite half crystal rigidly,
and reports energy per atomic interface cell. A local patch of slip is a
different configuration: surrounding intact material must deform. Multiplying
one-cell W by an invented area/count does not derive that deformation or the
critical configuration. The statistical specimen parameter A_c has no role.

For a smooth interface displacement jump d(x1,x2), consider the long-wave
energy reduction of the SAME atomistic potential:

    E_bulk = (1/2) sum_(halves) integral C_ijkl u_i,j u_k,l dV,
    C_ijkl = (1/V_atom) partial^2 E_atom / partial eta_ij partial eta_kl
             at the stress-free reference (with the strain convention fixed).

The local misfit term is

    E_misfit = integral_interface gamma(d) dA,
    gamma(a,s) = W_int(a,s) / A_atomic_cell.

This is an explicitly declared **continuum bulk + local interface**
approximation. It is not an exact atomistic partition at wavelengths of one
or two lattice spacings. Its use does not replace the infinite LJ pair sum,
the density/Bessel sums, or the per-atom embedding/angular terms.

The general connection between elasticity and an atomistically computed GSF
surface is established in Peierls--Nabarro research; for Al, see
[Lu, Kioussis, Bulatov & Kaxiras, PRB 62, 3099 (2000)](https://doi.org/10.1103/PhysRevB.62.3099),
and their [author preprint](https://arxiv.org/abs/cond-mat/9903440).
Those authors compared DFT and EAM GSF inputs and found that fine core
structure affects predicted dislocation properties. That is a reason to
validate our core, not permission to adopt their numbers or claim the present
research potential is calibrated. No empirical plastic evolution is imported.

## 2. A discovered +ABC/cubic-axis correspondence bug

The existing plane geometry has

    e1=(1,-1,0)/sqrt(2), e2=(1,1,-2)/sqrt(6), e3=(1,1,1)/sqrt(3),
    tau=(a1+a2)/3, z_l=l h, Delta_l=+l tau.

These e vectors define the geometric construction frame. Attaching the usual
cubic [100]/[010]/[001] components to that frame WITHOUT checking the sign of
tau is incorrect. For example, the first +tau upper atom, written using those
vectors, is

    (tau_x e1 + tau_y e2 + h e3)/a_lat = (2/3,1/6,1/6),

which is not a conventional FCC half-integer coordinate in those cubic axes.
The **generated lattice is still FCC**, in a different orientation. A proper
pi rotation about e3 supplies consistent cubic components:

    B_stack = rows(-e1,-e2,e3), det B_stack=+1,
    r_cubic = r_plane B_stack,
    r_cubic/(a_lat/2) in Z^3, with even component sum.

For the same first upper atom this gives (0,1/2,1/2). The tests verify all
sampled positive/negative layers and in-plane translations. No atom, ABC
shift, scalar slip direction, lattice energy, or fitted parameter was changed.

The new `plane_basis_in_stacked_cubic_axes()` is used for tensor rotation and
`StaticBulkHessian.crystallographic_wavevector()`. The latter previously
assigned cubic Gamma-X/L/K labels in the wrong orientation. Historical
`matched_v3/finite_q_static_stability.csv` vectors remain historical samples;
their old path labels are NOT silently certified. Fresh correctly labeled
samples are `nonlocal_v5/finite_q_corrected_paths.csv`.

`RegistryPath.direction_3d` remains a vector in the original GEOMETRIC frame,
not a newly relabeled cubic orientation. Local direct_110 and shockley_112
paths and their energy values are unchanged. A laboratory orientation must
always specify which crystal axes it uses.

This error was detected by an independent acoustic-limit test: the erroneous
cubic tensor gave relative matrix errors 0.372--0.403 for in-plane/oblique
waves, with the sign of an off-diagonal coupling reversed. After correction,
at |q|L0=0.005 and direct validation radius32L0, all four directions agree
within 1.86e-5. Errors decrease with q and are independently radius-refined.
Normal-wave-only checks would have missed the orientation error.

## 3. Half-space elimination without an added material coefficient

Let z=x3 be the normal, t a unit in-plane Fourier direction, and k>0. Write

    u(x,z) = u_hat(z) exp(i k t.x),
    A_ik = C_ij kl t_j t_l,
    B_ik = C_i3 kl t_l,
    D_ik = C_i3 k3,
    t_hat = D u_hat' + i k B u_hat.

Force balance gives

    d/dz [u_hat; t_hat/k] = k N [u_hat; t_hat/k],

    N = [ -i D^-1 B,                  D^-1;
          A-B^T D^-1 B,       -i B^T D^-1 ].

Three modes decay into z>0 and three into z<0. For invariant subspace columns
[U;T], define surface impedances (OUTWARD half-space traction):

    Z_plus  = -T_plus U_plus^-1,
    Z_minus = +T_minus U_minus^-1.

Ordered complex Schur subspaces are used instead of individually inverting
degenerate isotropic eigenvectors. Elastic units are scaled before this
linear algebra. Bulk Kelvin-matrix positivity, decay-mode count, Hermitian
residual, Riccati residual, and surface positivity are checked. Recorded
roundoff symmetrization does not repair a negative elastic eigenvalue.

For d=u_plus-u_minus, minimize the two half-space energies over their common
displacement. Equal/opposite interface tractions imply

    u_plus = Z_plus^-1 f, u_minus = -Z_minus^-1 f,
    K(k) = k (Z_plus^-1 + Z_minus^-1)^-1,
    E_el = (1/2) integral d^H K d d^2k/(2pi)^2.

Thus the nonlocal cost is DERIVED, not a tuned gradient penalty. At k=0,
K=0: uniform translations do not strain either half crystal. The local
uniform registry/opening stiffness belongs to gamma and is counted ONCE.
Adding the full local interface Hessian again to K would double count it.

This zero-mode observation is necessary, not sufficient to prove the local
gamma plus continuum split is atomically exact. Atomic surface relaxation and
short-wave elastic corrections must be matched before such a claim.

## 4. Independent checks and the screw special case

For isotropic lambda,mu with nu=lambda/[2(lambda+mu)], the identical-half-space
jump kernel in the (wave,tangential-transverse,normal) frame is

    K/k = diag(mu/[2(1-nu)], mu/2, mu/[2(1-nu)]).

The general Schur code recovers this including isotropic degeneracy. An
independent linear finite-element integration of bulk elastic energy over
z in [0,16/k], followed by a finite-matrix Schur complement, converges to the
same anisotropic surface impedance at the expected second-order rate. These
finite elements are a validation calculation, NOT an atomistic neighbor cutoff.

For the FCC mirror-symmetric screw case, displacement is along local e1,
variation along e2, and all fields are constant along the screw line e1.
The scalar gradient energy is

    (1/2) integral [A w_x^2 + 2 B w_x w_z + D w_z^2] dx dz.

The upper decay exponent is

    p = (sqrt(A D-B^2)|k| + i B k)/D,
    Z_plus = Z_minus = mu_bar = sqrt(A D-B^2),
    K_s(k) = mu_bar |k|/2.

With Delta=C11-C12,

    A=(Delta+4C44)/6, D=(Delta+C44)/3,
    B=-(Delta-2C44)/(3sqrt(2)) in the current +ABC frame,
    mu_bar=sqrt(C44 (C11-C12)/2).

Changing the stacking/frame sign changes B, not mu_bar. The full tensor test
is nevertheless essential, especially for general polarizations/directions.

## 5. Static nonlocal registry equation, units, and invariance

For the one-dimensional screw reference, with the normal gap clamped at h,

    E/line = (1/2) integral s (mu_bar/2)|D_x|s dx
             + integral gamma(h,s) dx - integral tau_res(x)s dx,
    delta E/delta s = (mu_bar/2)|D_x|s + gamma_s - tau_res,
    Hessian[v] = (mu_bar/2)|D_x|v + gamma_ss(s) v.

Here E/line is J/m, NOT joules. In code x,s use fixed L0, gamma uses eV/L0^2,
mu_bar and traction use eV/L0^3, and E/line uses eV/L0. The conversions are

    gamma_code = gamma_SI L0^2 / eV_J = W_cell[eV]/A_cell[L0^2],
    C_code = C_SI L0^3/eV_J,
    tau_code = tau_SI L0^3/eV_J,
    (E/line)_SI = (E/line)_code eV_J/L0.

The Fourier convention on length L is s_m=(1/L) integral s exp(-ik_m x)dx,
k_m=2pi m/L, so E_el/line=(L/2) sum K_m |s_m|^2. q=0 carries no elastic
cost. At zero applied traction, s -> s+b changes neither energy term. Mesh
refinement changes quadrature, not the material law or a coherent activation
area. L is the DECLARED periodic computational domain, not a hidden material
characteristic length. It is independently enlarged below.

The gamma Fourier series is a controlled spectral replay of the unchanged
analytic potential, not a new constitutive fit or a spline. Samples64/128
and independent off-grid analytic energy/gradient/Hessian checks are saved.
At128, maximum errors are 3.28e-15,7.27e-13,1.38e-10 in the respective
reduced energy/force/curvature units. No such table is registered in production.

Linearizing inside a well gives

    s_hat(k) = tau_hat(k)/(gamma_ss(0)+mu_bar |k|/2).

The implemented nonlinear Newton-CG solve additionally verifies the full
force residual and positive intrawell curvature. A force-free saddle is
REJECTED. Iteration number is not time; unloading here is static, not a
dynamic residual-plasticity hold test.

## 6. Actual unchanged candidate and long-wave checks

Source parameter file:
`matched_v3/monotone_opening/fitted_candidates.json`, SHA256
`9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655`.

- L0=2.863782463805517e-10 m; A_atomic_cell=7.102490842787125e-20 m^2.
- C11=87.81616955, C12=64.83263414, C44=49.27104858 GPa.
- mu_bar=23.795197936 GPa (3.48818093568 eV/L0^3).
- Pressure residual=-2.45e-14 GPa; same static fit, no optimization.
- gamma_ss(0)=5.71023805274 eV/L0^4.
- Fixed-gap direct_110 ideal shear traction=7543.359 MPa. This is a
  derivative maximum, NOT experimental yield, a coupled spinodal, or fatigue.

The comparison with the original static target C11/C12/C44=114/62/32 GPa
STILL fails quantitatively. A new positive nonlocal kernel does not fix that
material discrepancy. The source constants imply mu_bar=28.844410204 GPa,
retained only as a separate material-sensitivity comparator.

Direct finite-radius static K(q) is used ONLY for independent validation of
the analytic infinite potential: radii12/20/32 L0, qL0=.04/.02/.01/.005,
four directions. The radius32, q=.005 relative error is at most1.86e-5.
Fresh correctly oriented Gamma-X/L/K samples have positive eigenvalues
(minimum over the reported samples/radii0.0669003 eV/L0^2). This is not a
proof of stability over the entire Brillouin zone or finite-amplitude states.

Small sinusoidal tractions4/15/25 MPa were actually solved at wavelengths
4/8/16/32/128/512 b. The nonlocal/local linear amplitude ratios are
0.67578/0.80653/0.89290/0.94342/0.98523/0.99627. Short wavelengths are
explicitly outside an established atomistic accuracy regime. The long-wave
responses remain intrawell and statically unload to zero within the force
tolerance. This is not cyclic fatigue or a dynamic hold calculation.

## 7. Pre-existing screw dipoles: an explicitly restricted variational probe

To distinguish creation in a perfect crystal from movement of existing
defects, specify a pair of opposite straight screw defects. In a periodic
box, a periodized arctangent trial has separation d and width w:

    s_0 = b d/L,
    s_m = [2b/(L k_m)] sin(k_m d/2) exp(-|k_m|w), m!=0.

Its elastic energy is independently available as an exact infinite Fourier
sum:

    E_el/line = mu_bar b^2/(4pi)
       log[1 + sin^2(pi d/L)/sinh^2(2pi w/L)].

The identity follows from sum r^n/n=-log(1-r), so it checks Fourier
normalization, +/-mode counting, and finite-mesh error independently.
The external work is EXACTLY -tau_res b d per unit line length. At fixed
separation, the width is minimized using the unchanged gamma. Width is a
variational coordinate, not a fitted material parameter or an adopted core.

Define only the restricted-family force:

    tau_balance = (1/b) partial_d(E_el+E_misfit)/line,
    partial_d G/line = b(tau_balance-tau_res).

Positive tau_res-tau_balance favors expansion WITHIN THIS FAMILY. It is not
a measured velocity, event count, fatigue lifetime, or proven nucleation
saddle. At d>>w,

    tau_balance,periodic ~ mu_bar b/(2L) cot(pi d/L),
    tau_balance,isolated ~ mu_bar b/(2pi d).

This is the analytically derived defect-interaction scale; it explains why
motion of an existing widely separated defect pair need not require the
multi-GPa ideal resistance of a pristine rigid interface. The logarithmic
size dependence is also consistent with the classical dislocation problem;
see [Nabarro (1947), original-paper repository](https://wiredspace.wits.ac.za/items/71216603-7d8e-485b-8cc7-02e522c19877/full).
No experimental yield stress or phenomenological hardening was inserted.

### Executed separation and stress scenarios

| d [nm] | d/b | trial balance, L/d=8 [MPa] | isolated far-field [MPa] |
|---:|---:|---:|---:|
|9.1641|32|112.21742|118.34761|
|36.6564|128|28.05040|29.58690|
|146.6257|512|7.01254|7.39673|
|293.2513|1024|3.50627|3.69836|

For each separation, explicit resolved shear0/4/5/10/15/25/50 MPa was applied
to the energy derivative; see `dipole_stress_scenarios.csv`. These are stated
initial defect geometries, not inferred Al dislocation populations. At zero
stress, opposite defects attract. Reversing stress reverses the work term;
there is no built-in irreversible fatigue accumulation.

Material uncertainty matters: at293.25 nm the target-elastic comparator gives
an isolated balance4.483 MPa rather than3.698 MPa. Thus the prediction at
4 MPa CHANGES SIGN between these elastic descriptions. Claiming a calibrated
4 MPa Al transition would be incorrect. At5 MPa both far-field comparators
favor expansion for this SPECIFIED defect pair, still not physical fatigue.

### Numerical convergence is NOT core validity

At d/b128, L/d8, dx/b=.25/.125/.0625/.03125:

    tau_balance = 28.63133374,28.05266482,28.05040423,28.05040422 MPa.

The final spatial refinement change is1.15e-8 MPa. Independently, at dx/b=.0625
and L/d4/8/16/32/64:

    tau_balance = 23.23775308,28.05040423,29.20596108,29.49204783,29.56339690 MPa.

The remaining periodic-image effect is quantified by the analytic cotangent
formula; L/d64 is about0.08% below the isolated leading value, not exactly
an isolated boundary condition. Values from L/d8 alone are not domain-converged.
Changing the specified core width .15/.28/.5/1 b at d/b1024 varies the
L/d8 balance by about6.0e-5 MPa; far-field drive is far less core-sensitive
than an activation barrier.

However, the optimized trial width is0.27789 b (~0.0796 nm), narrower than
atomic row spacing. Also its FULL Euler residual remains about1.2 GPa:
stationarity with respect to width alone is NOT a solved dislocation core.
The energy-weighted spectrum has q90 b roughly0.67--1.23, not uniformly
small. A very fine numerical grid cannot repair this atomistic limitation.
The continuum/arctangent core is therefore **not physically certified**.

## 8. What remains to reach a physical probability model

The present kernel resolves the previously missing *long-wave elastic cost*
without modifying LJ/Bessel or introducing A_c. It does not solve all local
activation physics. Next requirements are:

1. Derive the discrete half-crystal force-constant/Green operator, retaining
   the analytic infinite lattice (including shifted reciprocal phases) and
   the scalar/angular many-body density derivatives. For an explicitly
   defined jump constraint B and gauge-fixed harmonic Hessian H, the exact
   constrained linear relation is K_int=(B H^+ B^T)^-1 where defined. Its
   zero-mode, translation constraints, and relaxed uniform interface must
   be established, not guessed by subtracting the rigid W_int Hessian.
2. Validate the resulting core against a compatible atomistic reference,
   allow vector registry/partial splitting and normal relaxation as needed.
   The clamped-gap scalar direct_110 path is not the minimum-energy 2D GSF path.
3. Define a finite dislocation loop/slip patch and its spatial energy measure.
   One infinite straight-line energy in J/m cannot be used in exp(-E/kBT).
   A chosen line length is not an acceptable hidden activation parameter.
4. For genuine probability evolution, retain the existing deterministic
   Smoluchowski structure over the justified spatial collective coordinates.
   Marginalizing an interacting field produces conditional neighboring forces,
   not automatically the old independent single-cell generator. Do not add
   Monte Carlo or a phenomenological plastic flow rule.
5. Independently improve/validate the material surface and obtain actual
   collective-coordinate kinetic information. No static energy, elastic
   constant, defect spacing, atomic mass, or A_c supplies physical seconds/Hz.

Opening probability, residual plasticity, and fatigue-life predictions were
NOT computed here. No production energy selector, PDE generator, kinetic JSON,
specimen aggregation, or UI default was changed. The UI redesign gate remains
closed pending solver evidence and user confirmation.

## 9. Reproduction and verification records

    python -m solver_v1.run_nonlocal_interface_reference --prepare --scenarios
    python -m solver_v1.report_nonlocal_interface_reference

The first command performs actual Bessel sampling, elastic checks and static
scenarios. The second only summarizes/plots saved data; it is not a re-fit or
a substitute for execution. Results and provenance are in
`results/fcc111_active_interface/nonlocal_v5/`.
The final test/commit checkpoint belongs to CURRENT_WORK_HANDOFF.md and
`validation_status.json`; unfinished/interrupted tests must not be counted.
