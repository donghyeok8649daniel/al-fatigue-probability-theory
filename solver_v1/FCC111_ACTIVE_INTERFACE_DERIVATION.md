# FCC(111) active-interface static reference

**Status: research only, not a validated Al constitutive surface.** Initial
parameter examples below are legacy underidentified candidates. Current
independent-mode calibration and physical gates are in
ALUMINUM_FULL_FCC_CALIBRATION.md and the audited_v2 result directories.
The production PDE and desktop defaults are unchanged.

## 1. Geometry and collective coordinates

Select the interface between planes 0 and 1. The lower half is `l<=0`; the
upper half is `l>=1`. With calibrated bulk spacing `h`, define

\[
 z_l=lh\quad(l\le0),\qquad
 z_l=lh+\delta a\quad(l\ge1),\qquad a=h+\delta a,
\]

and

\[
 u_l=0\quad(l\le0),\qquad u_l=s e_s\quad(l\ge1).
\]

Only one interface opens and slips. This differs fundamentally from the
homogeneous reference `Delta_l=l(tau+s e_s)`, which deforms every adjacent
plane pair.

The implemented paths are the explicit `direct_110` full-Burgers path and the
held-out `shockley_112` segment described in `FCC111_FULL_STACK_DERIVATION.md`.

## 2. Finite pair-energy difference

Lower-lower and upper-upper distances are invariant under rigid half-crystal
motion and cancel exactly. At layer separation `k`, exactly `k` lower/upper
plane pairs cross the cut per primitive interface cell. Hence

\[
 \Delta E_{pair}(a,s)=\sum_{k=1}^{\infty}k
 [w_{plane}(kh+\delta a,k\tau+s e_s)-w_{plane}(kh,k\tau)].
\]

This finite difference is evaluated with the existing 2D reciprocal kernels.
For the reciprocal-zero power term, with `r=2p-2` and `x=delta_a/h`, the
infinite mean difference is exact:

\[
 \sum_{k=1}^{\infty}k[(k+x)^{-r}-k^{-r}]
 =\zeta(r-1,1+x)-x\zeta(r,1+x)-\zeta(r-1).
\]

Its opening derivatives follow analytically:

\[
 D'(x)=-r[\zeta(r,1+x)-x\zeta(r+1,1+x)],
\]

\[
 D''(x)=r(r+1)[\zeta(r+1,1+x)-x\zeta(r+2,1+x)].
\]

Only exponentially decaying reciprocal corrugation is truncated, under an
explicit derivative-aware tolerance. No infinite bulk constant is numerically
subtracted.

## 3. Per-atom EAM density bookkeeping

Let `rho_bulk` be the full-FCC environment density. An atom at depth `r` from
either side of the interface has density change

\[
 \Delta\rho_r(a,s)=\sum_{k=r+1}^{\infty}
 [\rho_{plane}(kh+\delta a,k\tau+s e_s)-
  \rho_{plane}(kh,k\tau)].
\]

The correct embedding difference per interface cell is

\[
 \Delta E_{emb}=2\sum_{r=0}^{\infty}
 [F(\rho_{bulk}+\Delta\rho_r)-F(\rho_{bulk})].
\]

The factor two represents equal-depth atoms on the two sides. This is a sum
of per-atom embedding energies. It is not one global `F` of the entire crystal
density, and it is not the old reduced `2F(rho)` convention transplanted to a
multi-plane crystal. Far from the interface `Delta rho_r` decays exponentially,
so the neighborhood terminates using derivative-aware small-layer criteria.
These are empirical tail proxies, not rigorous total-error bounds. Independent
energy/derivative refinement is required. The legacy field named
estimated_tail_absolute combines internal reduced-coordinate convergence
indicators (including density) and is NOT an error measured in eV. Refine
reciprocal modes, layer count and direct radius independently.

For each coordinate `i,j` the analytic chain rule is

\[
 (\Delta E_{emb})_i=2\sum_r F'(\rho_r)(\Delta\rho_r)_i,
\]

\[
 (\Delta E_{emb})_{ij}=2\sum_r
 [F''(\rho_r)(\Delta\rho_r)_i(\Delta\rho_r)_j
 +F'(\rho_r)(\Delta\rho_r)_{ij}].
\]

## 4. Reference and opening limits

Define

\[
 W_{int}(a,s)=\Delta E_{pair}+\Delta E_{emb}.
\]

At the perfect state, all value differences cancel term by term:

\[
 W_{int}(h,0)=0.
\]

For the calibrated bulk candidates the residual gradients are below the bulk
fit/root tolerance and the interface Hessian is positive.

As `a` tends to infinity, all cross-interface pair and density contributions
vanish, leaving two free half-crystals. Therefore

\[
 W_{sep}=\lim_{a\to\infty}W_{int}(a,0)
\]

is the model work of separation per primitive cell and
`W_sep/A_atomic_cell` is compared with `2 gamma_111`. `A_c` never enters.

## 5. GSF limit

At fixed perfect opening,

\[
 \gamma(s)={W_{int}(h,s)-W_{int}(h,0)\over A_{atomic\,cell}}.
\]

The direct path has period `b`. The Shockley endpoint is `b/sqrt(3)` along
`[2,-1,-1]/sqrt(6)`; the entire primitive translation has period `sqrt(3)b`.
The initial underidentified candidates predict a Shockley endpoint essentially
degenerate with perfect FCC, so it misses the nonzero Al intrinsic stacking-
fault energy.

## 6. Coupled landscape and Hessian

The implementation returns analytic/semi-analytic
`W_a,W_s,W_aa,W_as,W_ss`. At the first bulk-equivalent candidate,

\[
 H_{int}=\begin{bmatrix}
 16.39502854 & -2.95\times10^{-15}\\
 -2.95\times10^{-15} & 2.89468722
 \end{bmatrix}\ \mathrm{eV/cell}
\]

in reduced coordinates. Candidate two gives `(16.36561647, 2.86957073)` on
the diagonal. Symmetry makes the pristine mixed curvature nearly zero, but the
saved 2D surface shows strong coupling away from the minimum: opening lowers
registry corrugation and shifted registry changes normal traction.

## 7. Numerical validation

At `a=1.05h,s=0.17b`, reciprocal versus independent radius-50/layer-40 real
space validation gives, for candidate one:

* energy absolute error `5.56e-7 eV/cell`;
* normal derivative error `1.37e-5`;
* registry derivative error `2.88e-8`;
* Hessian errors from `1.47e-9` to `1.48e-6`.

The reciprocal calculation uses about 10 interface layers. Analytic gradients
and Hessians also agree with centered finite differences. Machine-readable
curves, surface grids, convergence tables, and plots are in
`results/fcc111_active_interface`.

## 8. Held-out validation and barrier interpretation

The two initial, underidentified bulk candidates predict:

| quantity | candidate range | Al reference |
|---|---:|---:|
| direct `<110>` USF | 0.438--0.449 J/m2 | 0.250 |
| Shockley USF | 0.0764--0.0768 J/m2 | 0.224 |
| intrinsic SF | -0.00027 to -0.00024 J/m2 | 0.164 J/m2 |
| separation work | 1.089--1.145 J/m2 | 2.12 J/m2 DFT |

Thus a stable bulk fit is not sufficient. These candidates do not validate an
Al slip/crack landscape, but two failures do not prove global family
insufficiency. The Lu DFT values include atomic/volume relaxation at a different
lattice constant; rigid-half comparisons are contextual, not identical
geometric targets. Separation work is an asymptotic energy difference, not a
finite-load saddle barrier. Dynamic event ordering is not inferred.

## 9. External work and the status of chi

For interface area `A_atomic_cell`, normal opening `delta a_phys`, and slip
`s_phys`, the natural work is

\[
 G_{int}=W_{int}-A_{atomic\,cell}
 [T_n\delta a_{phys}+\tau_{res}s_{phys}].
\]

For uniaxial Cauchy stress `sigma e tensor e`, interface normal `n`, and slip
direction `m`,

\[
 T_n=\sigma(e\cdot n)^2,\qquad
 \tau_{res}=\sigma(e\cdot n)(e\cdot m).
\]

Therefore the static interface model does not need an arbitrary `chi` inside
its energy. Orientation projections remain necessary to map nominal axial
stress to normal and resolved shear traction. If a scalar axial strain is later
reported, its projection is derived from the same vectors; it must not be
retuned to force a desired mechanism.

## 10. Relationship to production theory

`TwoRowLJ` remains the verified reduced reference and the current PDE path is
unchanged. The active interface is labeled a static research reference. It is
not connected to the probability PDE because its Al interface validation has
not passed. The initial fit omitted independent shear information; the new
independent bulk fit still does not establish interface or kinetic validity.

Physical `M_a`, `M_s`, `t0`, seconds, and Hz remain unavailable. No static
barrier result supplies an overdamped kinetic scale.

## 11. Independent elastic refit and dimensional audit

The new audit solves isotropic, zero-pressure FCC equilibrium before
constructing either half. Both the in-plane lattice and normal spacing change
together. Parameters retain the fixed unit L0=b_target and the density gauge
is held fixed. In particular h/b=sqrt(2/3) at every perfect bulk reference.

Use the actual candidate's A_atomic for energy/area conversion:

    gamma[J/m2] = W_int[eV/cell] * E0 / A_atomic[m2]
    T_n[Pa] = dW_int/da* * E0 / (L0 A_atomic)

The following independent checks are recorded in audited_v2:

- reciprocal/layer tolerances 2e-8, 2e-10, 2e-12;
- direct (radius,layers)=(24,24),(48,24),(24,48),(48,48),(80,80);
- GSF along both explicit paths with real displacement units;
- opening out to 100h and an independent separation-limit evaluator;
- normal force, full local Hessian, and analytic-derivative finite differences;
- fixed-s=0 opening saddles at 25%,50%,75%,95% of peak traction.

The normal loaded diagnostic uses G=W_int-f_n(a-h) at ZERO resolved shear.
It is constrained to s=0 and does not select an arbitrary axial specimen
orientation. Its traction maximum is a constrained opening limit, NOT the
full coupled slip/opening spinodal. Registry curvature is saved to expose
possible loss of transverse stability.

Current independent-bulk-fit candidates fail the physical readiness gate:
one gives very large separation work and the historical-box candidate gives
negative separation work. The latter is energetically unstable to cleavage
despite a positive small-strain Hessian. Neither is admitted to the PDE.

The fully coupled static saddle map is deferred until an acceptable local
energy surface exists. Physical-time first passage additionally requires
independent collective-coordinate kinetic data. A negative
physical result must not be concealed by adjusting mobilities or plotting
scales.
