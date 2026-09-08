# Full FCC(111) infinite-plane-stack reference

## 1. Scope and model status

This document defines a new **static research/reference** geometry. It does
not replace `TwoRowLJ`, the reduced analytic LJ--EAM hybrid, or the current
desktop production choices. The existing `W(a,s)` is a reduced two-row /
projected lattice energy. The new object is a per-atom, homogeneous-spacing,
infinite FCC(111) plane stack. It is fuller crystallographically, but it has
not been refitted or kinetically calibrated.

The four distinct objects are

\[
 \phi(r),\qquad w_{\rm plane}(d,\boldsymbol\delta),\qquad
 W_\infty(a,s),\qquad
 G(a,s;f)=W_\infty(a,s)-f[(a-a_0)+\chi s].
\]

The old normal-chain/two-row energy is not added to `W_infinity`; doing so
would count the same atomic interactions twice.

## 2. FCC(111) orthonormal basis

For conventional FCC lattice constant `a_lat`, choose

\[
 \mathbf e_1={1\over\sqrt2}(1,-1,0),\quad
 \mathbf e_2={1\over\sqrt6}(1,1,-2),\quad
 \mathbf e_3={1\over\sqrt3}(1,1,1).
\]

Direct evaluation gives
`e_i dot e_j = delta_ij`. The in-plane nearest-neighbour spacing and
interplane spacing are

\[
 b={a_{\rm lat}\over\sqrt2},\qquad
 h_{111}={a_{\rm lat}\over\sqrt3}.
\]

In `(e1,e2)` components the triangular primitive vectors are

\[
 \mathbf a_1=b(1,0),\qquad
 \mathbf a_2=b(1/2,\sqrt3/2).
\]

Consequently

\[
 |\mathbf a_1|=|\mathbf a_2|=b,\quad
 \mathbf a_1\!\cdot\!\mathbf a_2=b^2/2,
 \quad A_{\rm atomic\,cell}={\sqrt3\over2}b^2
 ={\sqrt3\over4}a_{\rm lat}^2.
\]

`A_atomic_cell` is microscopic crystallographic normalization. It is not the
statistical correlation area `A_c`, and it is not a finite-element area.

## 3. Direct and reciprocal triangular lattices

Atoms in one plane occupy

\[
 \mathbf R_{mn}=m\mathbf a_1+n\mathbf a_2,qquad
 |\mathbf R_{mn}|^2=b^2(m^2+mn+n^2).
\]

Solving `a_i dot b_j = 2 pi delta_ij` gives

\[
 \mathbf b_1={2\pi\over b}(1,-1/\sqrt3),\qquad
 \mathbf b_2={2\pi\over b}(0,2/\sqrt3).
\]

Thus `G_hk=h b1+k b2` and

\[
 |\mathbf G_{hk}|^2={16\pi^2\over3b^2}(h^2-hk+k^2).
\]

The implementation groups complete equal-magnitude shells using the integer
quadratic form `h^2-h*k+k^2`; the first shell degeneracies are 6, 6, 6, 12
for indices 1, 3, 4, 7.

## 4. Exact ABC stacking phase

One gauge for the B-plane shift is

\[
 \boldsymbol\tau={\mathbf a_1+\mathbf a_2\over3}.
\]

Because `3 tau=a1+a2` is a lattice vector, `l tau mod Lambda` generates
`A,B,C,A,...`. Adding any lattice vector to `tau`, or choosing another
primitive-cell representative, is a gauge change and leaves every lattice
sum invariant.

## 5. Full registry and scalar paths

The general registry is a two-component displacement. Relative to plane zero,

\[
 \boldsymbol\Delta_l=\boldsymbol\tau_l+\mathbf u_l-\mathbf u_0.
\]

The implemented homogeneous scalar reference applies the same extra shift to
every adjacent pair,

\[
 \boldsymbol\Delta_l=l(\boldsymbol\tau+s\mathbf e_s)\pmod\Lambda.
\]

This is a homogeneous shear path, not an isolated single stacking-fault
interface. The explicit options are:

* `direct_110`: `e_s=e1`, period `b`. This matches the repository's existing
  direct `<110>` target: `s=b/2` maps to an `a_lat/4<101>` saddle and `s=b`
  to a full `a_lat/2<101>` translation.
* `shockley_112`: direction `[2,-1,-1]/sqrt(6)`, full lattice period
  `sqrt(3)b`, Shockley increment `b/sqrt(3)`. This is held out and is not the
  current scalar calibration path.

No unspecified crystallographic orientation is inferred.

## 6. Full atomic distance and self exclusion

For homogeneous normal spacing `z_l=l a`, the reference-atom displacement is

\[
 \mathbf r_{mnl}=\mathbf R_{mn}+\boldsymbol\Delta_l+l a\mathbf e_3,
 \qquad r_{mnl}^2=(la)^2+|\mathbf R_{mn}+\boldsymbol\Delta_l|^2.
\]

Only `(m,n,l)=(0,0,0)` is excluded. A future nonuniform stack replaces
`|l|a` by `d_ij=sum_{k=i}^{j-1} a_k`; that extension is not needed here.

## 7. LJ pair energy and plane kernel

With

\[
 \phi_{\rm LJ}(r)=4\epsilon_{\rm LJ}
 [\sigma_{\rm LJ}^{12}r^{-12}-\sigma_{\rm LJ}^{6}r^{-6}],
\]

define

\[
 S_p^{2D}(d,\boldsymbol\delta)=
 \sum_{\mathbf R\in\Lambda}[d^2+|\mathbf R+\boldsymbol\delta|^2]^{-p}.
\]

The per-atom full power sum is

\[
 S_p^{\rm full}={1\over2}\sum_{(l,m,n)\ne(0,0,0)}r_{mnl}^{-2p}
 ={1\over2}Z_\triangle(p)+\sum_{l=1}^{\infty}
 S_p^{2D}(la,\boldsymbol\Delta_l).
\]

The equality uses inversion symmetry: positive and negative layers are equal,
so their factor two cancels the pair-energy factor one half. The pair energy is

\[
 E_{\rm pair}=4\epsilon_{\rm LJ}
 [\sigma_{\rm LJ}^{12}S_6^{\rm full}
  -\sigma_{\rm LJ}^{6}S_3^{\rm full}].
\]

The same-plane triangular Epstein sum is evaluated as

\[
 Z_\triangle(p)=6\zeta(p)L_{-3}(p)b^{-2p},\quad
 L_{-3}(p)=3^{-p}[\zeta(p,1/3)-\zeta(p,2/3)].
\]

## 8. Two-dimensional Poisson/Bessel derivation

Use the Fourier convention

\[
 \widehat f(\mathbf G)=\int_{\mathbb R^2}f(\mathbf r)
 e^{-i\mathbf G\cdot\mathbf r}\,d^2r,
\quad
 \sum_{\mathbf R}f(\mathbf R+\boldsymbol\delta)
 ={1\over A_{\rm atomic\,cell}}\sum_{\mathbf G}
 \widehat f(\mathbf G)e^{i\mathbf G\cdot\boldsymbol\delta}.
\]

The Gamma representation

\[
 (d^2+r^2)^{-p}={1\over\Gamma(p)}\int_0^\infty
 t^{p-1}e^{-t(d^2+r^2)}dt
\]

and the 2D Gaussian transform produce, for `G>0`,

\[
 \widehat f_p(G)={2\pi\over\Gamma(p)}
 \left({G\over2d}\right)^{p-1}K_{p-1}(dG).
\]

The zero mode is evaluated separately:

\[
 \widehat f_p(0)={\pi\over p-1}d^{2-2p}.
\]

This is a two-dimensional transform and therefore uses `K_(p-1)`. The old
one-dimensional row construction correctly uses `K_(p-1/2)`; the formulas
must not be interchanged. Combining `+G/-G` produces real cosine corrugation:

\[
 S_p^{2D}={\pi d^{2-2p}\over A_{\rm atomic\,cell}(p-1)}
 +{1\over A_{\rm atomic\,cell}}
 \sum_{\mathbf G\ne0}\widehat f_p(G)
 \cos(\mathbf G\cdot\boldsymbol\delta).
\]

## 9. Exact mean layer sum and zeta powers

Summing the reciprocal-zero term over positive layers gives

\[
 \sum_{l=1}^{\infty}{\pi(la)^{2-2p}\over
 A_{\rm atomic\,cell}(p-1)}
 ={\pi a^{2-2p}\over A_{\rm atomic\,cell}(p-1)}
 \zeta(2p-2).
\]

Hence the attractive LJ term (`p=3`) contains `zeta(4)` and the repulsive
term (`p=6`) contains `zeta(10)`. A 1D row integral changes the radial measure
and yields powers `zeta(5)` and `zeta(11)`. Those are different dimensional
reductions, not alternative conventions for the same sum.

For nonzero reciprocal vectors, `K_nu(la|G|)` decays approximately as
`exp(-la|G|)/sqrt(la|G|)`. Complete reciprocal shells and layers are retained
until derivative-aware relative tolerances are met. Since
`Delta_(l+3)=Delta_l mod Lambda`, ABC phases repeat in three layer classes;
the implementation retains the explicit phase rather than imposing an
unhelpful closed form.

## 10. Full exponential environment density

Write the existing density kernel as

\[
 f_\rho(r)=\rho_0e^{-\beta_\rho(r/r_e-1)}
 =C_\rho e^{-\kappa_\rho r},\quad
 C_\rho=\rho_0e^{\beta_\rho},\quad\kappa_\rho=\beta_\rho/r_e.
\]

For a plane at distance `d`, differentiating the 2D Yukawa transform with
respect to `kappa` yields

\[
 \int_{\mathbb R^2}e^{-\kappa\sqrt{d^2+r^2}}
 e^{-i\mathbf G\cdot\mathbf r}d^2r
 ={2\pi\kappa e^{-dq}(1+dq)\over q^3},
 \qquad q=\sqrt{\kappa^2+G^2}.
\]

The full, self-excluded per-atom environment is

\[
 \rho_{\rm full}=\sum_{(l,m,n)\ne(0,0,0)}
 C_\rho e^{-\kappa_\rho r_{mnl}}
 =\rho_{\rm same\ plane}+2\sum_{l=1}^{\infty}
 \rho_{2D}(la,\boldsymbol\Delta_l).
\]

Critically, the embedding is applied once after this complete sum:

\[
 W_{\rm FCC}^{\rm full}=E_{\rm pair}+F(\rho_{\rm full}).
\]

It is **not** `sum_l F(rho_l)` and it is not a sum of old two-row hybrid
energies. This per-atom convention differs from the reduced two-row cell,
whose two symmetry-equivalent sites give `2F(rho)`.

## 11. Analytic derivatives

For either plane kernel, write a reciprocal coefficient `C_G(d)`:

\[
 S(d,\delta)=C_0(d)+\sum_{G\ne0}C_G(d)\cos(G\cdot\delta).
\]

Then

\[
 S_d=C_0'+\sum C_G'\cos,\quad
 \nabla_\delta S=-\sum C_G\sin(G\cdot\delta)G,
\]

\[
 S_{dd}=C_0''+\sum C_G''\cos,\quad
 \nabla_\delta S_d=-\sum C_G'\sin(G\cdot\delta)G,
\]

\[
 \nabla_\delta^2S=-\sum C_G\cos(G\cdot\delta)G\otimes G.
\]

Because `d_l=la` and `Delta_l=l(tau+s e_s)`, chain factors are
`d/da=l`, `d Delta/ds=l e_s`; second derivatives acquire `l^2`.
Power coefficients use analytic `K_nu` derivatives (`scipy.special.kvp`).
For the exponential coefficient

\[
 C_G={2\pi\kappa\over A q^3}e^{-dq}(1+dq),\quad
 C_G'=-{2\pi\kappa\over A}{d e^{-dq}\over q},\quad
 C_G''={2\pi\kappa\over A}{e^{-dq}(dq-1)\over q}.
\]

For scalar total density `rho`, embedding derivatives remain

\[
 W_i=E_{{pair},i}+F'(\rho)\rho_i,
\quad W_{ij}=E_{{pair},ij}+F''(\rho)\rho_i\rho_j+F'(\rho)\rho_{ij}.
\]

Thus `W_a,W_s,W_aa,W_as,W_ss` are analytic/semi-analytic reciprocal sums;
finite differences are validation only. The architecture retains vector
`grad_delta` and `hess_delta`, so a future two-component registry Hessian is
not blocked by the scalar path.

## 12. Equilibrium, Hessian, and stress mapping

At zero load the stable reference satisfies `W_a(a0,0)=0` and a positive
Hessian. With `c=(1,chi)^T`, the existing relaxed mapping is recomputed, not
copied:

\[
 \kappa_{\rm axial}={a_0\over c^T H_0^{-1}c},\qquad
 f^*=\kappa_{\rm axial}{\sigma\over E}.
\]

Neither `A_c` nor a mesh area enters this relation. `A_atomic_cell` appears
inside crystallographic lattice normalization only.

## 13. Static barriers and their limitations

The loaded surface is `G=W-f[(a-a0)+chi s]`. The registry barrier follows
the locally stable normal-relaxed homogeneous path. Its static spinodal solves

\[
 W_a=f,\quad W_s=\chi f,\quad
 W_{ss}-{W_{as}^2\over W_{aa}}=0.
\]

The normal spinodal at fixed registry is `max_a W_a(a,s)`. At zero force the
reported normal separation work is the energy difference between the bound
stack and infinitely separated individual planes. Under force, the saddle is
the larger root of `W_a=f`. This homogeneous dilation is not yet the work of a
single localized cleavage interface, so it must not be sold as a calibrated
cleavage curve.

For identical, unrefitted parameter values, the audit gives:

| model evaluated on full stack | a0/b | min Hessian eig. | kappa (chi=0.2) | registry spinodal | opening spinodal |
|---|---:|---:|---:|---:|---:|
| LJ reference | 0.6845871381 | 78.76610337 | 135.17737716 | 8.31726653 | 10.54407767 |
| hypothetical hybrid | 0.8983754061 | 0.929092275 | 4.002031175 | 0.198080934 | 0.295489421 |
| Al-target best feasible | no bound root | -- | -- | -- | -- |

At zero load, full-stack registry/opening barriers are respectively
`1.56084494/4.77709965` for LJ and `0.02831355/0.16632657` for the hypothetical
hybrid. Barrier ratios remain below one in the sampled stable range. These are
static landscape facts only; event ordering also requires mobility and time.

The important negative result is that the reduced Al-target best-feasible
parameters produce `W_a<0` over the searched full-stack range and approach the
isolated-plane limit without a stable normal root. Their reduced coordination
fit is therefore not transferable to full FCC coordination. No compensating
refit was performed.

## 14. Direct-sum validation and numerical status

`results/fcc111_full_stack` contains reproducible CSV tables. Representative
single-plane reciprocal/direct comparisons give maximum errors:

* exponential value/normal derivatives: `3.6e-15` absolute;
* `p=6` power quantities: `1.8e-13` absolute;
* `p=3` value: `2.24e-10`, equal to the direct-radius continuum tail estimate;
  its derivatives agree to about `1.5e-13`.

For the independent radius-60, 60-layer full sum, energy differences are
`2.13e-5` (LJ) and `7.52e-7` (hypothetical hybrid), both within the estimated
algebraic direct-sum tail. The reciprocal implementation uses typically 9
pair layers, 16 density layers, and at most 36 complete reciprocal shells for
the audited points. Analytic first/second derivatives agree with independent
finite differences at the recorded derivative-dependent tolerances.

Status labels are:

* exact analytic: crystallography, same-plane Epstein sum, Fourier transforms,
  zero-mode zeta layer sum;
* semi-analytic: tolerance-truncated reciprocal corrugation/layer series,
  branch roots and saddles;
* numerical: direct-sum validation and future PDE dynamics.

## 15. Relationship to probability dynamics and physical time

This task changes only the optional static energy reference. It does not alter
the Smoluchowski/Fokker--Planck equation, opening absorption, registry flux,
strain decomposition, statistical specimen area, or production defaults.
`M_a,phys`, `M_s,phys`, and `t0_seconds` remain unavailable. Consequently no
result here is a physical-Hz prediction.

The full stack is deliberately kept out of the desktop production selector:
the current best-feasible Al parameters are unbound under the fuller
coordination, and the homogeneous registry/opening coordinates require a new
mapping/refit before production use.

## 16. Limitations and next validation

The model is still a scalar homogeneous collective-coordinate reference. It
does not represent a localized dislocation core, a nonuniform plane opening,
all two-dimensional gamma-surface paths, free-surface relaxation, or physical
kinetics. The next scientifically defensible steps are (1) refit/identify the
analytic family using the full FCC coordination and consistent per-atom/
per-area targets, (2) validate both direct `<110>` and Shockley `<112>` paths
against DFT/EAM gamma surfaces, (3) introduce a localized interface coordinate
before claiming cleavage work, and (4) obtain independent kinetic mobility
data. None of those steps may be replaced by fatigue-life fitting.
