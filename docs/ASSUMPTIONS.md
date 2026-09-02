# Assumptions and approximations — active 1D theory
# 가정과 근사 — 활성 1D 이론

> **Normative status / 기준 상태:** This file follows the active 1D normal-only
> formulation. Historical `(a,s)` Smoluchowski/registry assumptions are retained only
> as a future extension and are **not** the governing assumptions of the present paper.
> 현재 기준 이론은 1D normal-only `P-u-Theta` 체계다. 과거 `(a,s)`
> Smoluchowski/registry 가정은 향후 확장 이력으로만 남으며 현재 논문의 활성
> 지배가정이 아니다.

Authoritative equation and notation sources:

- `README_EQUATION_INDEX.md`
- `docs/EQUATION_SUMMARY_1D_P_U_THETA.md`
- `docs/VARIABLE_INDEX_1D_P_U_THETA.md`
- `docs/MASTER_1D_P_U_THETA_FORMULATION.md`
- `docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`

## 0. Mandatory time-dependence notation / 시간의존성 표기 규칙

The active probability and moment fields are functions of spacing and time. Their
canonical nondimensional-time forms are

\[
\boxed{
P(\lambda,\tau),\qquad
u\text{ is not an active state symbol},\qquad
u\not\equiv u,
}
\]

and the actual active mean-rate/variance fields are

\[
\boxed{
 u(\lambda,\tau)=\mathbb E[c\mid\lambda,\tau],
\qquad
\Theta(\lambda,\tau)=\operatorname{Var}(c\mid\lambda,\tau).
}
\]

With physical time,

\[
\boxed{
P(\lambda,t)=P\!\left(\lambda,\tau=t/t_0\right),
}
\]

with analogous notation for `u` and `Theta`.

Bare `P`, `u`, or `Theta` may appear only as local shorthand after the full
functional dependence has already been declared. 정식 표기는 반드시
`P(lambda,t)` 또는 `P(lambda,tau)`, `u(lambda,t)` 또는 `u(lambda,tau)`,
`Theta(lambda,t)` 또는 `Theta(lambda,tau)`처럼 시간의존성을 포함한다.

## 1. Active physical assumptions / 활성 물리 가정

1. The target baseline is a pure single crystal represented by a one-dimensional
   normal chain under repeated uniaxial normal loading.
2. The active microscopic coordinates are node positions `x_j(tau)` or normalized
   nearest-neighbour spacings `lambda_i(tau)=x_i-x_{i-1}`. Registry/slip `s` is not
   required by the current normal-only paper mainline.
3. The active microscopic configurational energy is the nearest-neighbour
   generalized-LJ chain energy
   \[
   \boxed{V^*(\boldsymbol\lambda)=\sum_{i=1}^M\phi(\lambda_i).}
   \]
   The same energy must generate both the equations of motion and G2 if exact
   mechanical consistency is claimed.
4. The current finite-chain dynamics is deterministic and conservative apart from
   prescribed boundary work. No viscous damping, white noise, empirical damage,
   phonon bath, or stochastic diffusion is inserted into the active baseline.
5. The initial homogeneous ideal baseline may be
   \[
   \boxed{\lambda_i(0)=1,\qquad \dot\lambda_i(0)=0,}
   \]
   which gives a delta empirical state. Any broader initial realization measure
   `mu_0` must be physically declared rather than silently assumed.
6. The one-point probability state is mechanically generated. For one deterministic
   chain it is a spatial empirical counting measure over represented spacings; for
   an ensemble it is the push-forward of a declared full-state initial measure.
7. No named probability family (Gaussian, Weibull, Gibbs/Boltzmann, etc.) is imposed
   on `P(lambda,tau)`.
8. Smooth one-point fields are used only where the empirical measure admits a
   meaningful smooth/coarse representation and the required conditional moments
   exist.
9. Moment integration assumes sufficient decay in spacing-rate `c` so that the
   required velocity-space boundary terms vanish.
10. The divided density-shape formula is used only where
    \[
    \boxed{P(\lambda,\tau)>0,\qquad \Theta(\lambda,\tau)>0.}
    \]
    At `Theta=0`, the undivided transport/moment equations are used.
11. Neighbour independence is not assumed. The conditional acceleration and the
    `Theta` source retain neighbour joint statistics.
12. The exact general second-central-moment equation is
    \[
    \boxed{
    D_\tau\Theta
    +2\Theta\,\partial_\lambda u
    +\frac1P\partial_\lambda(PC_3)
    =2\Psi,
    }
    \]
    where
    \[
    \boxed{\Psi(\lambda,\tau)=\operatorname{Cov}(c,\ddot\lambda\mid\lambda,\tau).}
    \]
    Setting `Psi=0` or `C_3=0` is an additional closure assumption and is not active
    by default.
13. `Theta(lambda,tau)` is a conditional spacing-rate variance, not an empirical
    fatigue-damage scalar and not by itself the total kinetic-energy density.
14. Same-load loading/unloading non-retracing of `(P,u,Theta)` establishes dynamic
    history dependence, not irreversible dissipation.
15. The current conservative baseline has
    \[
    \boxed{\dot D_{\rm irr}=0,\qquad E_{\rm hyst}=0.}
    \]
    G3 requires a separately derived physical irreversible mechanism before a
    nonzero irreversible hysteresis energy can be claimed.
16. The operational local initiation threshold is the loss of positive tangent
    stiffness,
    \[
    \boxed{\phi''(\lambda_c)=0,}
    \]
    followed by first passage through `lambda_c`.
17. Local spatial first-passage fraction and specimen-to-specimen crack-initiation
    probability are distinct. Specimen probability requires a physically declared
    realization measure `mu_0` and represented correlation scale; independent-cell
    multiplication is not assumed.
18. The exact reduced differential equations are hierarchical rather than an
    autonomous three-field closure, but the closed finite LJ dynamics supplies
    exact push-forward/integral representations for `F`, `P`, `u`, `Theta`, `C_3`,
    `Psi`, G1, G2, and first-passage survival.
19. The microscopic time scale `t_0` must not be confused with laboratory fatigue
    cycling. A physical bridge from microscopic dynamics/history to laboratory
    Hz-scale fatigue accumulation remains open.
20. 2D/3D CAD/FEM may transport/visualize scalar normal mechanics but does not by
    itself activate multiaxial constitutive physics.

## 2. Assumptions that are NOT active / 현재 활성화하지 않는 가정

The active 1D paper does **not** assume:

- Boltzmann/Gibbs equilibrium as the fundamental spacing distribution;
- Gaussian/Weibull spacing or fatigue-life PDFs;
- Smoluchowski/Fokker--Planck mobility closure;
- Einstein fluctuation--dissipation relation as an already justified reduced law;
- arbitrary diffusion kernels or white noise;
- independent neighbouring spacings;
- independent statistical cells or independent FEM element failure probabilities;
- viscous damping or an empirical fatigue damage variable;
- FCC lattice reconstruction;
- registry coordinate `s` or unwrapped slip index `z` as a required current state;
- pure normal loading as an already proven low-frequency slip/dislocation mechanism.

## 3. Historical `(a,s)` registry theory / 과거 `(a,s)` registry 이론의 위치

The multiplicity-free multilayer energy `U_0(a,s)`, Bessel/polylog representations,
and unwrapped registry variable `z` are retained as mathematically defined extension
material. They are not deleted. However, current numerical/mechanical checks show
that the perfect-symmetry pure-normal baseline does not justify claiming an active
low-frequency `s` transition or plastic slip mechanism.

Thus, for the ideal symmetric pure-normal baseline, the conservative extension is

\[
\boxed{
P(a,s,t)=P(a,t)\,\delta(s-s_0),
}
\]

until a physical symmetry-breaking mechanism is introduced and justified.

Possible future mechanisms such as defects, dislocations, a physically derived
phonon bath, or another slow internal variable are OPEN and must be derived rather
than inserted ad hoc.

## 4. Remaining open physical inputs / 남아 있는 물리적 미해결 항목

- physical irreversible microscopic mechanism for G3;
- physically justified specimen/initial full-state measure `mu_0` and correlation scale;
- bridge from microscopic time/history to laboratory fatigue cycling;
- quantitative material calibration beyond the present reduced LJ bridge;
- experimental validation of first-passage initiation;
- any future symmetry-breaking/plasticity extension and its relation to dislocation mechanics.

These are physical open problems, not reasons to revert the active theory to an
assumed PDF family or an unjustified stochastic closure.
