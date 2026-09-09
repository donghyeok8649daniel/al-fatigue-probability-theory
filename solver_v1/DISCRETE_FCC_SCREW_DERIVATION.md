# Infinite-row FCC screw reduction: what was actually missing

Status: **separate static research calculation**, not a production probability
generator or a calibrated aluminum dislocation model. Results are in
`results/fcc111_active_interface/discrete_screw_v6/`. The unchanged parameter
source is the `angular_monotone_opening` candidate used in `nonlocal_v5`.
Its elastic-target errors remain; this calculation is not a new fit.

The user's question was why the remaining low-stress mechanism was unresolved.
Two distinct defects in the preceding calculation needed a quantitative fix:

1. A one-width arctangent trial family did not satisfy the full Euler equation.
   More spatial samples of that same family could not remove its force error.
2. The core was narrower than an atomic spacing. A continuum `mu*|q|/2` kernel
   was only a long-wave approximation, not the atomic finite-wavevector energy.

We now solve **all profile modes** for the existing continuum reference, and
independently derive the **actual infinite-row discrete harmonic energy**.
They must not be combined by silently double-counting a local GSF Hessian.
This fixes the identified mathematical approximations at their respective
levels; it does not certify a nonlinear atomic core or a fatigue probability.

## 1. Geometry: infinite line, discrete cross-section

Take the same FCC(111) lattice and the explicit direct-110 line direction e1.
Let `d=sqrt(3)*b/2` be the transverse spacing between rows in one plane.
The two-dimensional primitive vectors and ABC translation are

\[
 a_1=(b,0),\quad a_2=(b/2,d),\quad \tau=(b/2,d/3).
\]

For line index j, layer l, and atom index n along the line,

\[
 R_{njl}=(nb+\delta_{jl},\ d(j+l/3),\ hl),\qquad
 \delta_{jl}=b(j+l)/2\pmod b.
\]

Here h is the **actual cubic bulk equilibrium spacing**, not an independently
changed density normalization. For a perfect FCC reference h/b=sqrt(2/3).
Every complete row has delta=0 or b/2, so it is invariant under x -> -x.

The allowed displacement is `u_(j,l) e1`, constant along the infinite line.
This is an exact anti-plane subspace, not the full vector core state space.
The (j,l)=(0,0) line has zero relative displacement between any of its atoms:
its energy is a constant, not an omitted force. The self atom is not counted.

The exact generated +ABC orientation uses
`geometry.plane_basis_in_stacked_cubic_axes()` for comparison with cubic C.
The earlier geometric/lab basis must not be reused as a cubic tensor frame.
No FCC coordinates, LJ parameters, or historical reduced energies were changed.

## 2. LJ row stiffness: reuse the verified Bessel identity

For transverse radius A=sqrt(y^2+z^2),

\[
 S_p(A,\delta)=\sum_{n\in\mathbb Z}[A^2+(nb+\delta)^2]^{-p}
             =C_{p0}(A)+\sum_{m\ge1}C_{pm}(A)\cos(g_m\delta),
 \quad g_m=2\pi m/b,
\]
\[
 C_{pm}(A)=\frac{4\sqrt\pi}{b\Gamma(p)}
       \left(\frac{\pi m}{bA}\right)^{p-1/2}K_{p-1/2}(g_m A).
\]

Therefore the exact row curvature is

\[
 \phi_{\mathrm{row},xx}=-4\epsilon\sum_{m\ge1}g_m^2
 [\sigma^{12}C_{6m}-\sigma^6C_{3m}]\cos(g_m\delta).
\]

`row_lj_slip_curvature` reuses the mode coefficient from `lattice_bessel.py`.
There is no new pair potential or finite-neighbor cutoff. The reciprocal zero
mode is independent of line translation and drops out of every slip derivative.

## 3. Environmental density and analytic moment derivatives

For f(r)=C exp(-kappa r), the same one-dimensional transform is

\[
 H(A,G)=\int_{-\infty}^{\infty}e^{-\kappa\sqrt{A^2+x^2}}e^{-iGx}dx
       =\frac{2A\kappa}{Q}K_1(AQ),\quad Q=\sqrt{\kappa^2+G^2}.
\]

The curvature of the row density is

\[
 \rho_{\mathrm{row},xx}=-\frac{2C}{b}\sum_{m\ge1}g_m^2
 H(A,g_m)\cos(g_m\delta).
\]

For STF angular moments we need the first derivative of
`exp(-kappa*r) STF(R^rank)` with respect to the displacement along x.
Multiplication by x^p has transform i^p partial_G^p H; a derivative in x
multiplies that by iG. A compact analytic recurrence avoids finite differences:

\[
 H_j:=\partial_{G^2}^j H
 =2A\kappa(-A/2)^j Q^{-j-1}K_{j+1}(AQ),
\]
\[
 H_G=2G H_1,\quad H_{GG}=2H_1+4G^2H_2,\quad
 H_{GGG}=12GH_2+8G^3H_3.
\]

If a rank-r raw tensor component contains p x-indices, u y-indices and v
z-indices (p+u+v=r), its row derivative is

\[
 B^{\mathrm{raw}}_{p,u,v}
 =\frac{2C_\mathrm{ang}}{b}\sum_{m\ge1}
 \operatorname{Re}\{i g_m i^p H^{(p)}(A,g_m)e^{ig_m\delta}\}y^u z^v.
\]

Apply STF **after this sum** to obtain B_1, B_2, B_3. Normalization
`C_ang` and `kappa_ang` are read from the unchanged angular surface, not
replaced with the scalar-density parameters. This matters for the current fit.

The implementation independently checks these formulas against differentiated
large real-space rows, including non-special phases. At the perfect phases,
the scalar first density derivative is exactly zero for every complete row.

## 4. Full per-atom EAM/angular harmonic energy

At perfect cubic bulk the scalar density is rho_bulk, and Q1=Q2=Q3=0.
Let the row-relative anti-plane displacement be Delta u_R=u_(i+R)-u_i.
The harmonic scalar EAM expansion contains

\[
 \Phi_R=\phi_{\mathrm{row},xx}(R)
       +2F'(\rho_\mathrm{bulk})\rho_{\mathrm{row},xx}(R).
\]

The factor 2 follows from pair counting versus per-site embedding, not an
arbitrary convention. The additional scalar term would be
`F''(rho_bulk) |sum rho_row,x Delta u_R|^2/2`; it vanishes **in this specific
subspace** because each rho_row,x is zero. It must not be dropped in a general
three-component or line-dependent deformation.

For each angular rank r,

\[
 \delta Q_r(i)=\sum_R B_r(R)\Delta u_R,
 \quad E_{r}^{(2)}=D_r\sum_i\|\delta Q_r(i)\|^2.
\]

All rows are summed before taking the norm. D1 can be negative in the saved
research candidate; it is neither clipped nor replaced by its absolute value
in the physical Hessian. Absolute values are used only in tail diagnostics.
The Q2 term assumes cubic reference; noncubic Q_bulk requires extra terms and
is rejected here. This is not conventional calibrated MEAM.

## 5. Exact anti-plane bulk symbol and ABC Bloch phase

Use the logical-row Bloch convention

\[
 u_{jl}=u_0 e^{i(q_y d j+\theta l)},\qquad
 q_z=(\theta-q_y d/3)/h.
\]

The symbol is

\[
 K(q_y,\theta)=\sum_R\Phi_R(1-\cos\varphi_R)
   +2\sum_{r=1}^3D_r\left\|\sum_R B_r(R)(e^{i\varphi_R}-1)\right\|^2,
 \quad \varphi_R=q_y d j+\theta l.
\]

K is in eV/L0^2 per atom. `2 sin^2(phi/2)` and `expm1(i phi)` avoid loss of
accuracy near the translational mode K(0,0)=0. We do not assign an atomic mass
or convert a static eigenvalue to Hz. The full 3D real-space Hessian independently
confirms that the e1 polarization decouples for q_x=0 in this geometry.

### Convergence, not a hidden finite cutoff

Every row-slip derivative has zero reciprocal-zero term. Its leading
transverse decay is approximately exp(-2*pi*A/b), with algebraic prefactors.
Thus, after infinite analytic summation along each line, the remaining
discrete cross-section sum is exponentially convergent. This is fundamentally
different from cutting off the original LJ neighbors at a finite radius.

Absolute row envelopes are accumulated before Bloch cancellation. In each
angular channel let S_r=sum_all ||B_r|| and T_r=sum_omitted ||B_r||. Since
||L_r||<=2S_r, the omitted symbol satisfies

\[
 |\Delta K|\le2\sum_{\rm omitted}|\Phi_R|
               +16\sum_r|D_r| S_r T_r.
\]

This combines quantities in **stiffness units**, not unlike pair-curvature and
moment-gradient units. It bounds omission within the extra validation rings;
the infinite remainder beyond them is still a convergence estimate. The code
reports used rings, validation rings, last-shell envelope and the sum of extra
validation rings. It calls this an **observed absolute-shell remainder estimate**,
not a rigorous infinite-series error bound. We additionally compare independent
2/3/4/5/6/8/10-ring evaluations and a large 3D direct validation sum.

## 6. Rigid-cut Hessian is an independent counting check

If all upper layers translate by s and all lower layers remain fixed, cross-
interface separation l occurs l times. The exact curvature is

\[
 H_{\mathrm{rigid}}=\sum_{l>0,j}l\Phi_{jl}
 +4\sum_rD_r\sum_{d_0\ge0}\left\|\sum_{l>d_0,j}B_r(j,l)\right\|^2.
\]

The 4 is two symmetric half-crystals times the second derivative of a squared
moment. This reconstructs the separately implemented full active-interface
W_ss at s=0, including all scalar and angular terms. It is a much stronger
check than comparing only a bulk elastic coefficient.

For the actual candidate:

- original W_int,ss = 4.945211215219163 eV/L0^2 per interface cell;
- row reconstruction = 4.945211215219173;
- absolute discrepancy about 9.77e-15 (1.3e-13 across earlier looser controls);
- pair plus redistributed scalar embedding = 2.521193077910488;
- angular contribution = 2.424018137308685.

## 7. Infinite-stack relaxation and the ACTUAL constrained variable

For a declared linear jump B u, minimizing the bulk harmonic energy gives

\[
 K_{\rm jump}=(B H^+ B^\dagger)^{-1}.
\]

We implement two explicitly **different** constraints. Neither is silently
called the unique atomistic definition of a spatial disregistry field.

1. `logical_rows`: s_j=u_(j,1)-u_(j,0). The two rows have physical transverse
   positions differing by d/3. In Bloch form B=e^(i theta)-1.
2. `spectral_same_y`: translate the upper-layer Fourier interpolation by -d/3
   before taking the difference, B=e^(i(theta-q_y*d/3))-1. This is a defined
   interpolation of discrete row data, not a pair of adjacent atoms.

For offset eta=0 or d/3 respectively,

\[
 C_{\rm jump}(q_y)=\frac1{2\pi}\int_{-\pi}^{\pi}
 \frac{4\sin^2[(\theta-\eta q_y)/2]}{K(q_y,\theta)}d\theta,
 \qquad K_{\rm jump}=C_{\rm jump}^{-1}.
\]

A midpoint quadrature avoids the removable q=theta=0 pole. Automatic doubling
requires **two consecutive** converged changes. Small q can require many
points (the executed q L0=.002 study reached 65536). One fixed coarse angular
grid would miss the narrow acoustic contribution. Unstable sampled bulk modes
or exhausted refinement are errors, not silently repaired inverse matrices.

An independent finite-layer matrix assembled directly from row coefficients
and per-site angular operators agrees with this Green/Schur result.

At q=0 both constraints coincide and

\[
 K_{\rm jump}(0)=4.795625552378889,
\]

about 3.02% lower than the rigid-cut stiffness. This is legitimate surrounding-
layer relaxation, not a mobility fit or softening inserted by hand.

## 8. Even the low-q matching needs the constraint definition

Let the acoustic tensor in the plane frame be A=C_0101, B=C_0201,
D=C_0202, V the volume per atom, and define

\[
 c=VD/h^2,\quad \alpha=d/3-hB/D,\quad
 \beta=h\sqrt{AD-B^2}/D.
\]

Near the acoustic pole,

\[
 K(q,\theta)\simeq c[(\theta-\alpha q)^2+\beta^2q^2].
\]

Writing theta=alpha*q+v separates the even v^2 and constant terms in the
compliance integral. The nonanalytic linear contribution is

\[
 C_{\rm jump}(q)-C_{\rm jump}(0)
 =\frac{(\alpha-\eta)^2-\beta^2}{2c\beta}|q|+O(q^2),
\]
\[
 K_{\rm jump}(q)-K_0
 =K_0^2\frac{\beta^2-(\alpha-\eta)^2}{2c\beta}|q|+O(q^2).
\]

The independent numerical small-q integral verifies this formula. An atomic
gap constraint is not the same as a mathematical displacement discontinuity
between two cut continuum half-spaces; its slope need not equal
`A_atomic_cell*mu_bar/2`. The mismatch is not fixed by fitting a new coefficient.

In particular, at q L0=1 the calculated stiffness is about 8.83393 for logical
rows and 8.08444 for spectral same-y interpolation. At the transverse zone edge
they are 26.19934 and 22.34258, respectively. The old continuum-plus-rigid-local
comparator differs substantially in that regime. Full tables preserve the units.

### Do not double count local energy

K_jump includes its local K0 already. Adding the existing gamma'' on top would
double count local stiffness. Subtracting K0 and then adding a nonlinear rigid
GSF curve is **also not automatically consistent**: the correct local nonlinear
curve would have to be relaxed under the same constraint as the Schur kernel.
We do neither. A future fully nonlinear discrete row/cross-section energy must
be validated before a nonlinear atomistic core or probability path is declared.

## 9. Removing the one-width ansatz in the continuum reference

This is a separate controlled calculation with the *previous* continuum
functional (see `NONLOCAL_INTERFACE_ELASTICITY.md`):

\[
 E[s]/\text{line}=\frac12\int s\frac{\bar\mu}{2}|D|s\,dx
                 +\int\gamma(s)dx.
\]

Fix registry content Q=int s dx=b*d_content, describing a specified existing
opposite screw pair rather than creating a new pinning potential. Minimize over
**all sampled profile values**, not a width parameter. The Euler condition is

\[
 (\bar\mu/2)|D|s+\gamma'(s)=\lambda,
 \quad \int sdx=Q,
\]

where lambda is the constant holding traction. Projection removes only the
constraint reaction; it must not be mistaken for a force-balance error.

Fourier square-root preconditioning changes optimizer coordinates, not energy.
L-BFGS objective termination alone is insufficient. A projected Newton-CG
polish and direct maximum-force test check stationarity independently.
`force_converged` is distinct from the optimizer's `success` flag. These are
static constrained stationary solutions, not automatically proven global minima
or nucleation saddles. Translations and the fixed-content mode require gauges.

The variational identity

\[
 \frac{dE_{\min}}{d d_{\rm content}}=b\lambda
\]

is independently checked by re-solving nearby contents. A uniform applied
traction tau contributes -tau*Q; hence tau-lambda is a **static growth drive**
for content. With Q held fixed, varying tau does not simulate a trajectory or
change the constrained profile. We label the stress table accordingly.

The content-equivalent separation and the distance between s=b/2 crossings
are not identical because intrawell far-field offsets and tails contribute
to the mean. Both are saved. Energy retains eV/L0 of dislocation **line**,
never eV of a finite activation event.

### A second protocol error: matching content does not match core separation

At nonzero holding traction the intact intrawell far-field shift is
approximately `lambda/gamma''(0)`. Its integral grows with the domain length.
Thus fixing Q while enlarging the domain changes the actual core separation.
For content 128 b^2, L/content=4,8,16,32 produced actual crossing separations
127.6933,127.2580,126.4450,124.8193 b. Those are **not identical defects**.
The corresponding traction trend must not be mislabeled an outer-boundary
convergence study for a fixed physical pair.

`--matched-domains` instead solves a bracketed outer problem: adjust the
constraint Q until the *actual* crossing separation equals 128 b in every
domain. Material parameters and the functional are held fixed. Intermediate
force convergence and the existence of exactly two crossings are mandatory.
At dx/b=.0625, L/d=4,8,16,32 the resulting tractions are approximately
23.237335,28.049974,29.205530,29.491616 MPa. These agree with the derived
finite-period far field; the remaining periodic correction is not hidden.
The actual L/d=64 solve gives 29.56296548 MPa, against the isolated leading
asymptote 29.58690135 MPa. The last domain still has a .0809% periodic effect.
At L/d=8, dx/b=.0625 -> .03125 changes the matched-separation traction by
only 1.95e-10 MPa. Unlike the fixed-content table, this comparison holds the
same actual defect geometry. There were 38 independently solved constrained
profiles in the matched-domain study (850.30 seconds for this execution).

For fixed content 128 b^2 and L/content=8, dx/b=.25,.125,.0625,.03125 gives
holding tractions 9.842586,26.003904,28.231490,28.231613 MPa. All four have
small force residuals, yet the coarse values are wrong: numerical lattice
pinning in a continuum discretization is not physical atomic Peierls pinning.
The final spatial refinement change is .000123 MPa. This is a numerical
convergence result for the stated continuum constraint, not an atomistic
material error bar or a plastic-signal certificate.

## 10. What this does and does not resolve

Verified components:

- analytic infinite-row pair/density/STF derivatives;
- same-energy 3D finite-q comparison with radial refinement;
- exact reconstruction of the pre-existing rigid interface Hessian;
- infinite-layer Schur compliance, phase and finite-layer checks;
- actual automatic layer-phase and transverse-ring refinement;
- static low-MPa harmonic responses for different spatial wavelengths;
- full-profile rather than one-width continuum stationarity;
- force, mesh, domain, conjugate-work and registry-content checks.

The old 1.2 GPa *internal* trial-profile residual is removed by solving the
actual continuum variational equations, not by changing LJ, temperature,
mobility, chi, or an empirical slip law. The atomic short-wave omission is now
quantified with an independently verified lattice Hessian rather than hidden
behind numerical mesh refinement.

Still unresolved, and therefore **not production-valid fatigue**:

1. A nonlinear discrete atomic core with a consistent spatial disregistry;
   two-component partial splitting and normal relaxation may be important.
2. Material fit: the saved C11/C12/C44 remain about 87.816/64.833/49.271 GPa,
   not the 114/62/32 GPa reference set. No new Al calibration was performed.
3. A finite nucleus/loop energy. Straight-line energy is J/m; an arbitrary
   line length or statistical correlation area cannot convert it into a
   physical activation barrier.
4. Collective-coordinate physical mobility, seconds and Hz remain unavailable.
5. No new crack opening probability, accumulated fatigue damage, experimental
   yield stress, or residual aluminum plasticity is inferred from these static
   calculations. The UI/production model and gate remain unchanged.

## 11. Reproduction and provenance

Run the actual calculations, not only the plotter:

```text
python -m solver_v1.run_discrete_screw_reference --kernel --profiles --matched-domains --report
python -m pytest solver_v1/test_discrete_fcc_screw.py solver_v1/test_nonlocal_registry_reference.py -q
```

Inputs are the unchanged saved analytic candidate and the checksummed
`nonlocal_v5/prepared_misfit.npz`. The latter's Fourier energy/force/Hessian
have already been checked against off-grid analytic interface evaluations.
Raw tables, quadrature histories, timing, profile data, and separate physical
versus numerical status are saved under `discrete_screw_v6`.

The methodological continuum/dislocation context remains the primary sources
listed in `NONLOCAL_INTERFACE_ELASTICITY.md` (Nabarro 1947; Lu et al. 2000).
All numerical values above were computed from this repository's candidate,
not copied from those papers or claimed as experimental Al predictions.

## 12. Next nonlinear calculation: explicit contract, NOT completed code

To close the remaining atomic-core gap, the same row sum permits a full
nonlinear discrete anti-plane energy before any continuum approximation.
This is the next proposed calculation, not a result produced by this module.
For cross-section sites i and neighboring rows R define

\[
 \eta_{iR}=\delta_R+u_{i+R}-u_i,
\quad
 \Delta\rho_i=\sum_R[\rho_{\rm row}(A_R,\eta_{iR})
                         -\rho_{\rm row}(A_R,\delta_R)],
\]
\[
 \Delta Q_{ri}=\sum_R[Q_{r,\rm row}(A_R,\eta_{iR})
                         -Q_{r,\rm row}(A_R,\delta_R)].
\]

With one line-repeat b per cross-section site, the finite state-dependent
energy difference is

\[
 \mathcal E_b[u]=\frac12\sum_{i,R}
 [\phi_{\rm row}(A_R,\eta_{iR})-\phi_{\rm row}(A_R,\delta_R)]
 +\sum_i[F(\rho_b+\Delta\rho_i)-F(\rho_b)]
 +\sum_{i,r}D_r[\|Q_{rb}+\Delta Q_{ri}\|^2-\|Q_{rb}\|^2].
\]

Q_rb=0 only for the stated perfect cubic reference. The scalar F is per site;
**F'_bulk cannot be frozen outside the harmonic expansion**. Infinite-row
reciprocal-zero terms cancel in each displacement difference, so the transverse
state-dependent interaction remains exponentially convergent. A nonlinear
implementation must recover the independently verified harmonic symbol on
linearization and the existing rigid GSF curve when rigid cut motion is imposed.

The line energy is `mathcal E_b/(b*L0)` in J/m after converting eV to J. It is
not a finite thermal activation energy. Specifying a finite loop/patch requires
additional spatial geometry rather than a numerical energy multiplier.

There is an additional kinematic requirement: translating a whole infinite
atomic row by b relabels identical atomic positions. Energy is periodic in its
row phase, while unwrapped registry transport records history. A dislocation
calculation must specify the Burgers/winding content, boundary conditions and
unwrapped slip consistently. A discontinuous integer row relabeling alone is
not proof of a resolved dislocation core. Do not replace these constraints by
an arbitrary spring or a tuned barrier. Full vector registry and normal opening
remain further extensions requiring their own stability and convergence tests.
