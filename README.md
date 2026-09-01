# Al Fatigue Probability Theory

## Persistent research status / 영구 진행상황 기록

**Repository rule:** every commit that changes theory, physical scope,
numerical results, paper content, or validation status must update this section
before it is pushed. A new work session must read this section first. Completed,
in-progress, rejected, and physically unresolved work must never be inferred
only from chat history.

**저장소 유지 규칙:** 이론, 물리 범위, 수치 결과, 논문 또는 검증 상태가
바뀌는 모든 작업은 push 전에 반드시 이 절을 갱신한다. 새로운 작업 세션은
항상 이 절부터 읽는다. 채팅 기록에만 진행상황을 남기지 않는다.

Last updated: **2026-09-02 (multilayer derivation implemented and locally verified)**

- **Current integration task:** the normal spacing $a$ and one scalar
  crystallographic registry $s$ are being derived from one common intrinsic
  multilayer potential
  $U_0(a,s)=\sum_{k=1}^{\infty}W(ka,s)$.  The row--row $W$ is a kernel, not
  the final lattice energy.  There is no multiplicity prefactor $k$, and the
  same collective/unwrapped $s$ is used in every layer term (never $ks$).
- **Probability state:** the active fundamental state is the reduced
  two-coordinate density $P(a,s,t)$ for uniaxial tensile fatigue in a single
  crystal.  A declared loading axis and one declared slip system give
  $Q_a=A_0\sigma(t)$ and $Q_s=A_0M\sigma(t)$; there is no independent shear
  fatigue input or multiaxial fatigue criterion.
- **Plasticity criterion:** registry is unwrapped as $s=zb+\tilde s$.
  Intrawell phase lag alone is not plasticity; plastic deformation requires a
  residual change $\Delta\langle z\rangle\ne0$ after unloading/relaxation.
- **Verified current result:** a subcritical dimensionless resolved-shear
  pulse leaves $\langle z\rangle=0.4915266$ while the mean intrawell registry
  returns to $1.6\times10^{-10}$. Six symmetric zero-mean cycles leave only
  $\langle z\rangle=0.001004$. This is a mechanism demonstration, not an
  aluminum plastic-strain calibration.
- **Unified-energy rule:** normal opening and slip are the exact identity
  $\Delta U_0=\Delta U_n+V_{\rm slip}$ from the same $U_0$.  The prior
  collinear $U_\infty$ and single-row $W$ remain historical/reduced
  derivations but are not added as the fundamental total energy.
- **Still unresolved:** active slip-system selection, representative
  mechanical area, MD-derived mobility/memory,
  homogenization thickness, dislocation storage, hardening, and two-way
  normal--slip coupling. Therefore quantitative Al cyclic plasticity is not
  claimed.
- **FEM/UI:** geometry and scalar normal-stress visualization remain available.
  The new registry solver is not yet coupled to FEM elements; mesh dimension
  does not create multiaxial constitutive physics.
- **Verification:** direct $(k,p)$ sums agree with the exact
  $H_q$ Bessel--Lambert form for $q=6,12$ at tested points at relative errors
  of order $10^{-11}$ or smaller.  The independently derived 12--6 polylog
  closure agrees with the Bessel--Lambert evaluation at machine precision.
  The complete local suite passes (167 tests). The canonical paper was built
  locally with Tectonic 0.17.0; its log contains no LaTeX errors, missing
  glyphs, overflow warnings, or unresolved references. The tracked PDF is
  `output/pdf/slip_lattice_energy_derivation.pdf` (SHA-256
  `e433296316a04a1ae4239da4ce323476c81feb1e8559e304789fd00bdf9fe23a`).
  No GitHub Actions workflow or persistent automation is used or changed.

### Historical pre-multilayer branch description

The following branch-separation paragraphs are retained only as the immediate
pre-2026-09-02 research history. They are not the active governing model.

The project now has two deliberately separated reduced branches. The primary
fatigue branch remains the one-dimensional normal-tensile spacing model. The
optional plasticity branch is also one-dimensional in its lattice registry and
uses one scalar resolved-slip coordinate. Activating this branch does not
activate Rubin chains, an arbitrary shear-fatigue criterion, multiaxial
fatigue, or a full crystal-plasticity law. A joint future $(a,s)$ evolution
would have a two-coordinate probability state space, but it would still not be
a 2D continuum constitutive law.

현재 프로젝트는 두 개의 축약 branch를 분리해 사용한다. 주 균열개시 branch는
1D normal spacing 모델이고, 선택적 소성 branch는 하나의 slip system에 대한 1D
scalar registry 모델이다. Bessel 격자합은 이제 실제 활성 에너지로 사용된다.
다만 이는 전위 증식·경화까지 포함하는 정량적 단결정 소성모델은 아니다.

Reproduce the active registry demonstration and its focused verification with:

```powershell
py -3 -m simulations.run_registry_plasticity
py -3 -m simulations.verify_multilayer_lattice
py -3 -m pytest tests/test_multilayer_lattice.py tests/test_registry_plasticity.py -q
$env:MPLBACKEND='Agg'; py -3 -m pytest -q
```

The generated data and figures are under `results/data/registry_plasticity/`
and `results/figures/registry_plasticity/`.

## Active scope freeze / 현재 연구 범위

The active fundamental state is now the coupled reduced coordinate $(a,s)$
under one uniaxial tensile waveform. Its intrinsic energy is exclusively
$U_0(a,s)=\sum_{k\ge1}W(ka,s)$: no layer multiplicity $k$, no $ks$, and no
external-work term is included in $U_0$. Normal and slip excess energies are
an exact decomposition of this same potential. Plasticity is a residual
unwrapped-well transition, and crack initiation is normal-barrier first
passage. EAM/DFT remains future validation only. See
`docs/ACTIVE_IDEAL_REGISTRY_PLASTICITY.md` and `docs/ASSUMPTIONS.md`.

현재 활성 fundamental state는 하나의 단축 인장파형 아래 $(a,s)$를 쓰는
축약모델이다. intrinsic energy는 오직
$U_0(a,s)=\sum_{k\ge1}W(ka,s)$이며, layer multiplicity $k$, $ks$, 외력 일은
$U_0$에 넣지 않는다. normal/slip energy는 같은 $U_0$의 항등분해이고,
소성은 unwrapped well의 잔류 이동, 균열개시는 normal barrier first passage로
정의한다. EAM/DFT는 향후 검증용이다.

### Superseded pre-multilayer status (historical context only)

The paragraphs below in this subsection describe the previous separated
normal-chain/two-row baseline and are retained only to document the transition.
They are not the active governing model after 2026-09-02.

The primary crack-initiation theory is a one-dimensional normal-tensile spacing
model for pure single-crystal aluminum. The exact one-registry Bessel model is
now an optional active reduced-plasticity branch. Full FCC half-space,
Rubin-chain, multiaxial-fatigue, and conventional crystal-plasticity models are
not active. Even when a 2D/3D mesh is shown, the mesh does not create additional
constitutive physics.

주 균열개시 이론은 **순수 단결정 알루미늄의 1D normal tensile spacing
model**이다. 여기에 단일 slip system의 scalar registry를 쓰는 선택적 1D
Bessel 소성 branch가 활성화되어 있다. 완전한 FCC half-space, Rubin chain,
multiaxial fatigue 및 통상적인 crystal-plasticity 모델은 아직 활성화하지
않았다. 2D/3D mesh의 존재 자체가 구성방정식을 다축으로 만들지는 않는다.

The active coordinate is the local approximately homogeneous atomic spacing
before crack initiation, so the exact homogeneous zeta-lattice energy is the
active energy. The local crack-gap energy is retained only as a post-initiation
archive comparison in `theory/exact_lattice_energy.py`. The finite-volume kinetic model
with reflecting or absorbing escape is in `theory/smoluchowski_escape.py`.
The periodic absorbing solver and its one-cycle survival spectrum are in
`theory/smoluchowski_floquet.py`. Its principal multiplier is derived from the
governing PDE, not fitted as a fatigue-life coefficient. See
`docs/SMOLUCHOWSKI_FLOQUET_SURVIVAL.md`.
수식을 최소화한 한국어 연구진전 설명은
`docs/RESEARCH_PROGRESS_KO.md`에 정리되어 있다.
Uncalibrated runs are dimensionless demonstrations, not aluminum-life
predictions. In particular, the mechanical representative area $A_0$ is not a
correlation area $A_c$ and is not a FEM element area.

주기하중의 비가역 균열개시는 흡수형 Smoluchowski 방정식의 한 주기
생존연산자로 계산한다. 최대 고유값은 장시간 cycle 생존비이며 별도로 맞춘
피로계수가 아니다. 현재 수치는 무차원 예제일 뿐 단결정 알루미늄의 수명
예측값이 아니다.

The owner-supplied corrected two-row ideal-slip derivation and symbol index are
preserved under `research/source/` and audited under `libraries/shear/`. The
exact shifted Epstein--Hurwitz/Poisson--Bessel identity and scalar unwrapped
registry dynamics are active in `theory/registry_lattice.py` and
`theory/registry_plasticity.py`. Unsupported mixed patch energy, automatic
irreversibility, and quantitative-Al claims remain rejected. See
`docs/ACTIVE_IDEAL_REGISTRY_PLASTICITY.md`.

프로젝트 소유자가 제공한 23쪽 두 원자열 ideal-slip 유도자료는
`research/source/`에 원본으로 보존하고 `libraries/shear/`에서 검토한다.
이동된 Epstein--Hurwitz/Poisson--Bessel 항등식은 정확한 수학으로 보존하지만,
slip 변수, 혼합 patch energy 및 소성 주장은 활성화하지 않는다. 자세한 판정은
`libraries/shear/docs/SLIP_LATTICE_ENERGY_REVIEW.md`에 있다.

## `.ftgsim` project files

The desktop tensile app now reads and writes the open `.ftgsim` container.
It is a versioned ZIP/ZIP64 bundle containing checksummed JSON/CSV/PNG data,
never Python pickle or executable content. A project can be opened without OS
file association:

```powershell
py -3 -m simulations.fem_tension_app path/to/model.ftgsim
```

The format records the 1D mesh dimension, loading axis, material/loading
inputs, display state, and optional FEM results. Windows registry association
is intentionally not performed. See `docs/FTGSIM_FILE_FORMAT.md`.

## CAD-style mesh inspection

OBJ, binary/ASCII STL, ASCII PLY and legacy ASCII VTK geometry can be opened
before or after analysis in a common 1D/2D/3D viewport. It provides 3D orbit,
pan, wheel/right-drag zoom, reset, projection switching and the scalar normal
loading-axis marker:

```powershell
py -3 -m simulations.fem_tension_app examples/cube_3d.obj
```

STEP/IGES and automatic solid/volume meshing are not claimed because the
repository has no CAD-kernel/mesher backend. See `docs/MESH_VIEWPORT.md`.

When a project contains a separately solved `initiation_elements.csv` channel,
the FEM UI can color axial elements by cumulative initiation probability,
survival, or hazard. Missing probability results are shown as missing rather
than fabricated from FEM element count or display geometry. The active
definition is first passage through the zero-tangent-stiffness stretch; see
`docs/CRACK_INITIATION_DEFINITION.md`.

## Single-crystal loading direction

Every active run declares a nonzero cubic loading direction `[h k l]`. The GUI
accepts it as `Crystal axis [h k l]`. If all of $C_{11},C_{12},C_{44}$ are
provided programmatically or in `.ftgsim`, the application projects the scalar
directional Young modulus; otherwise `E_axis` is treated as a user-supplied
direction-specific value. The optional registry branch additionally requires a
declared slip-plane normal and in-plane slip direction and uses the signed
Schmid projection of the uniaxial load. The GUI does not yet solve this branch,
and no multiaxial fatigue criterion is activated. See
`docs/SINGLE_CRYSTAL_ORIENTATION.md`.
## Candidate kinetic probability and FEM coupling

The four working quantities are now connected by a separate, explicitly labeled kinetic post-processor:

$$
\int p_e(\lambda,t)d\lambda=1,
\qquad
\bar a_e=a_0\int\lambda p_e\,d\lambda,
$$

$$
u_{{\rm LJ},e}
=E\int[\phi(\lambda)-\phi(1)]p_e\,d\lambda,
\qquad
H_{e,k}=\oint_k\sigma_e\,d\bar\lambda_e.
$$

Each C FEM element supplies its local normal-stress history to a conditional-intact Smoluchowski solver. The implementation preserves normalization and the finite-volume Gibbs stationary state, exports mean spacing, variance, energy, loop work, and the tail above $\lambda_c$, and renders the actual 1D node/element mesh plus 2D/3D tensile-only views.

This is a **candidate kinetic extension**, not yet a calibrated aluminum
fatigue-life law. Reflecting runs use the critical tail only as an
instantaneous instability diagnostic. Separate absorbing runs define
irreversible initiation as first passage through $\lambda_c$ and report
survival, outgoing flux and hazard. These two observables are not mixed. See
`docs/MILESTONE16_PROBABILITY_ENERGY_HYSTERESIS_FEM.md` and
`docs/SMOLUCHOWSKI_FLOQUET_SURVIVAL.md`.

The periodic survival calculation is reproducible with

```powershell
py -3 -m simulations.run_smoluchowski_floquet
```

It writes verification plots and CSV/JSON data under
`results/figures/smoluchowski_floquet/` and
`results/data/smoluchowski_floquet/`. At the explicitly uncalibrated demo
point, $r\simeq0.9048$ explains the previously observed low-cycle probability
loss. The result also proves that conditioned survivor energy becomes
periodic; the active Markov model accumulates escaped probability, not an
arbitrarily retained fraction of hysteresis work.

With $E_0=E_{[hkl]}A_0a_0$, the reduced load is exactly
$Fa_0/E_0=\sigma/E_{[hkl]}$, but the inverse temperature and physical clock
remain $E_{[hkl]}A_0a_0/(k_BT)$ and
$t_r=\gamma a_0/(E_{[hkl]}A_0)$. Thus loading-axis stress alone cannot identify
life; the representative area and spacing mobility remain required physical
inputs and are never replaced by FEM mesh size.

The geometry layer now supports actual 1D line, 2D quad, and 3D hex connectivity, dependency-free STL/OBJ surface import, and optional STEP/IGES/BREP meshing through Gmsh. A lightweight NumPy/Matplotlib mesh UI exposes nodes, edges, opacity, and axial clipping without adding a VTK stack to the core installation. This does **not** change the active theory: every mesh cell receives only the declared tensile normal scalar $\sigma_{nn}$, or a scalar derived from its 1D $P(a,t)$. The current 1D-to-mesh mapping is explicitly a visualization/post-processing projection, not a 2D/3D elasticity solve. See `docs/MILESTONE17_2D_3D_MESH_CAD_NORMAL_ONLY.md`.

## Active 1D statistical-cell dependence scale

<!-- STATISTICAL_CELL_STATUS_EN -->

Probability aggregation is now explicitly separated into complete identical dependence, partial dependence, and true independence. For a second-order stationary 1D spacing process, the exact variance identity is

$$
\operatorname{Var}(\bar\lambda_M)
=\frac{\sigma_\lambda^2}{M}\tau_M,
\qquad
\tau_M=1+2\sum_{k=1}^{M-1}\left(1-\frac{k}{M}\right)\rho_k.
$$

This defines the variance-equivalent independent count and axial statistical length

$$
M_{\rm eff}=\frac{M}{\tau_M},
\qquad
\ell_{\rm stat}^{(2)}=a_0\tau_M.
$$

A single deterministic finite snapshot cannot use the all-lag population formula directly because sample-mean centering gives an exact weighted zero-sum identity. Such snapshots therefore use a separately labeled first-positive-lobe estimator. In the dynamically matched sweep $M=31,63,127,255$, the corrected estimate gives $M_{\rm eff}^{(+)}\approx2.93,2.99,3.03,3.05$ while $\ell_{\rm stat}^{(2,+)}/a_0\approx10.58,21.10,41.93,83.49$. Thus the tested protocol retains system-scale coherence rather than converging to a local material correlation length.

The mechanical calibration area $A_0$ is not identified with a transverse statistical independence area. The active scope remains strictly one-dimensional.

## Active correction — quasistatic limit of the deterministic correlation snapshot

<!-- QUASISTATIC_PROTOCOL_STATUS_EN -->

The earlier $M_{\rm eff}^{(+)}\approx3$ arithmetic remains a valid normalized-shape diagnostic for the selected deterministic snapshot, but its physical interpretation is now corrected. The committed Milestone 13 snapshot is taken at integer cycle 2 under zero-mean sinusoidal end loading, so the exact applied force at that phase is zero.

For the homogeneous force-controlled chain,

$$
\Pi=\sum_i[\phi(\lambda_i)-f\lambda_i],
\qquad
\phi'(\lambda_i)=f.
$$

On the stable branch $\phi''>0$, the root is unique and therefore every spacing is identical:

$$
\lambda_i=\lambda_s(f)\quad\forall i.
$$

At the sampled zero-force phase this gives $\lambda_i=1$ and $C_0=0$ exactly in the quasistatic state. A new $\alpha=\omega M$ sweep shows that the deterministic residual variance collapses rapidly toward zero as the drive is slowed, while the normalized correlation shape and positive-window $M_{\rm eff}^{(+)}$ remain near three. Therefore normalized residual correlation alone cannot define a material statistical-cell length.

The project now explicitly distinguishes a single-trajectory spatial empirical measure $P_{M,\mathrm{spatial}}^{\mathrm{traj}}$ from a future physically specified ensemble probability $P_{\rm ens}$. The next active target is a justified 1D initial ensemble under the same nonlinear cyclic mechanics, not a fitted named distribution.

## Active physical-statistical constraint on $P$

<!-- PHYSICAL_P_STATUS_EN -->

The current mechanics now gives a physically grounded hierarchy of possible spacing distributions without fitting a named family. With

$$
U(a)=E_0\phi(a/a_0)+C,
$$

the elastic calibration $E=(a_0/A_0)U''(a_0)$ and $\phi''(1)=1$ imply

$$
\boxed{E_0=EA_0a_0},
\qquad
\boxed{\chi=\frac{EA_0a_0}{k_BT}}.
$$

At zero-temperature homogeneous quasistatic equilibrium, the distribution is not broad:

$$
\boxed{P(\lambda\mid f)=\delta[\lambda-\lambda_s(f)]}.
$$

At fixed total normalized length, canonical equilibrium gives the exact finite-$M$ marginal

$$
\boxed{
P_M(\lambda\mid L,\chi)
=\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}.
}
$$

For constant tensile force $f>0$, the full-domain Gibbs integral diverges because $\phi(\lambda)-f\lambda\to-\infty$. For $0<f<f_c$, a Gibbs density conditioned on the intact basin $0<\lambda<\lambda_b(f)$ is therefore only a **controlled metastable/local-equilibrium approximation**, not a global equilibrium law or a fatigue-life law. See `docs/MILESTONE12_PHYSICAL_STATISTICAL_P.md`.

The representative layer area $A_0$ is now an explicit physical bottleneck: no numerical aluminum thermal distribution is claimed until $A_0$ is defined consistently with the coarse-grained layer interaction.

## Active-theory correction — exact nonlinear transport

**CURRENT ACTIVE STATUS:** the harmonic/single-mode and Taylor-expanded push-forward material retained later in this README is historical/local diagnostic work and is **not** the active global form of $P(\lambda,t)$.

The current derivation keeps the original nonlinear generalized-LJ force and starts from the exact empirical transport identities

$$
\partial_tP+\partial_\lambda J=0,
$$

with the spacing-velocity phase-space state $F_1(\lambda,v,t)$ and neighboring-spacing joint state $P_2$ entering the exact hierarchy. See `docs/MILESTONE10_EXACT_DISTRIBUTION_TRANSPORT.md`.

<!-- ACTIVE_TRANSPORT_CORRECTION_EN -->

Mechanics-first research framework for fatigue crack initiation under **one-dimensional normal cyclic loading** in high-purity / single-crystal aluminum.

## Legacy normal-chain scope (historical module)

The section below documents the earlier normal-only branch. The active common
energy and current scope are defined in the persistent status above and in
`docs/ACTIVE_IDEAL_REGISTRY_PLASTICITY.md`.

The active derivation is deliberately restricted to a one-dimensional stack of represented material layers. The microscopic reduced coordinate is the normal spacing

$$
a_i(t)>0,
$$

or the normalized spacing

$$
\lambda_i(t)=a_i(t)/a_0.
$$

The effective normal interaction between layers is represented by the calibrated generalized Lennard--Jones model. Three-dimensional FCC work remains archived under `libraries/fcc_normal/`. The corrected one-index registry/Bessel mechanism is active as a separate optional branch and is never added to the collinear normal energy.

Physical time $t$ is fundamental. Fatigue cycle count is not an independent state variable.

## Calibrated 1D layer-LJ mechanics

The normalized layer energy is

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
$$

with

$$
m=12.19,
\qquad
n=6.
$$

The calibration gives

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

The interior spacing equation is

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
}
$$

The current idealized tangent-instability stretch is

$$
\boxed{
\lambda_c
=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
\approx1.1077715386.
}
$$

## One-point probability state

The normalized one-point spacing density is

$$
\int_0^\infty p_\lambda(\lambda,t)\,d\lambda=1.
$$

Its mean and shifted configurational energy are

$$
\mu(t)=\int_0^\infty\lambda p_\lambda(\lambda,t)\,d\lambda,
$$

$$
\psi(\lambda)=\phi(\lambda)-\phi(1),
$$

and

$$
\mathcal E(t)=\int_0^\infty\psi(\lambda)p_\lambda(\lambda,t)\,d\lambda.
$$

The exact energy-feasibility work shows that normalization, mean, and energy alone cannot force a tensile tail because arbitrarily large energy can mathematically be hidden in the LJ compression branch as $\lambda\to0^+$. Any exact safe-energy ceiling therefore requires an independently justified compression-side constraint.

## Earlier two-moment distribution closure

A fixed-length/fixed-configurational-energy equal-base-measure assumption plus a large-$M$ saddle-point reduction gave the controlled approximation

$$
\boxed{
p_\lambda(\lambda,t)
=Z^{-1}\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)].
}
$$

The multipliers are determined by $\mu(t)$ and $\mathcal E(t)$ rather than fitted to histograms.

Direct deterministic tests showed that this closure can reproduce near-equilibrium variance reasonably well but does not reproduce the full driven distribution. In the tested 32-node states the Kolmogorov distance is about $0.15$, and the slower case shows a strong skewness mismatch. Therefore $\mu$ and $\mathcal E$ are not sufficient to determine the driven one-point distribution.

## Spatial-correlation result

The deterministic chain retains strong spatial ordering. Define

$$
C_k(t)
=
\frac{1}{M-k}
\sum_{i=1}^{M-k}
[\lambda_i-\mu][\lambda_{i+k}-\mu],
$$

and

$$
\rho_k=C_k/C_0.
$$

In the dynamically matched sweep with $\omega M=0.62$, the nearest-neighbor correlation rises from about $0.933$ at $M=31$ to about $0.991$ at $M=255$. The profile approximately collapses when plotted against $k/M$, with the first zero crossing near $0.35M$.

A one-point density is exactly invariant under permutation of the layer labels while $C_k$ is not. Therefore

$$
\boxed{p_\lambda(\lambda,t)\text{ cannot encode the complete spatial mechanical state}.}
$$

## New governing-equation clue for the form of $P$

The one-point density is also the exact spatial push-forward of the deterministic spacing field. With continuum layer label $\xi\in[0,1]$ and spacing field $\Lambda(\xi,t)$,

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1\delta[\lambda-\Lambda(\xi,t)]\,d\xi.
}
$$

For a piecewise monotone field,

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}.
}
$$

This provides a direct clue for the functional form of $P$: its shape is the kinematic image of the mechanically generated spatial waveform rather than an independently selected probability family.

Linearization around $\lambda=1$ gives

$$
\ddot u_i=u_{i+1}-2u_i+u_{i-1},
$$

with dispersion

$$
\boxed{\omega_q^2=4\sin^2(q/2).}
$$

A single coherent linear mode therefore generates the exact spatial push-forward

$$
\boxed{
p_{\rm 1mode}(\lambda)
=
\frac{\mathbf 1_{|\lambda-\mu|<|A|}}
{\pi\sqrt{A^2-(\lambda-\mu)^2}}.
}
$$

This arcsine form is not a fitted distribution; it follows from $\Lambda=\mu+A\cos\vartheta$ with uniform spatial phase.

The calibrated LJ force is strongly nonlinear:

$$
\boxed{
\phi'(1+u)
=u-10.595u^2+62.97935u^3+O(u^4).
}
$$

Thus the base mechanics itself generates harmonic distortion, skewness, and potentially multimodal one-point densities.

A first-plus-second-harmonic waveform

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta
$$

has

$$
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
$$

$$
\mu_3=\frac34A^2B,
$$

and the exact bound

$$
\boxed{|\gamma_1|\le\sqrt{2/3}\approx0.816497.}
$$

The slower deterministic snapshot has $\gamma_1\approx1.062$, so even two spatial harmonics are insufficient for that state. A single-mode arcsine reference also gives larger KS error than the earlier exponential closure. These are negative results that point toward richer mechanically generated mode content, not toward fitting another named probability family.

## Current research direction

The current hierarchy is

$$
\boxed{
\text{1D layer-LJ governing equation}
\rightarrow
\Lambda(\xi,t)\text{ / neighboring-spacing pair state}
\rightarrow
p_\lambda(\lambda,t)\text{ by push-forward}
\rightarrow
Q_c(t)
\rightarrow
\text{normal-opening first passage}.
}
$$

The next target is to derive the mechanically excited spatial-mode or pair-state evolution directly from the 1D governing equation and determine the smallest representation that reproduces both the one-point distribution $P_1$ and the observed $C_k$ without empirical relaxation laws or fitted probability families.

## Active files

- `theory/normal_lj_chain.py` — conservative 1D layer-LJ chain
- `theory/normal_lj_energy_feasibility.py` — exact energy-feasibility bounds under stated compression constraints
- `theory/normal_lj_distribution.py` — earlier two-moment large-$M$ closure
- `theory/normal_lj_spatial_correlation.py` — exact finite-chain correlation diagnostics
- `theory/normal_lj_pushforward.py` — spatial push-forward and mode-derived distribution identities
- `simulations/run_normal_lj_spatial_correlation.py` — dynamically matched correlation sweep
- `simulations/run_normal_lj_pushforward_clue.py` — single-mode falsification diagnostic
- `docs/MILESTONE8_SPATIAL_CORRELATION.md` — spatial-ordering result
- `docs/MILESTONE9_GOVERNING_EQUATION_PUSHFORWARD.md` — governing-equation clue for the form of $P$
- `results/data/result_manifest.json` — current machine-readable research state

## Reproduce active results

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

## Research labels

Important statements are classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**
- **NUMERICAL RESULT / DIAGNOSTIC**

No fitted Gaussian/Weibull distribution, cycle-dependent LJ parameter, empirical damage variable, or fitted correlation length is accepted as a mechanics derivation.

---

# 한국어 번역

## 후보 확률동역학 및 FEM 결합

현재 네 개의 작동식은 별도로 표시한 kinetic post-processor로 연결된다.

$$
\int p_e(\lambda,t)d\lambda=1,
\qquad
\bar a_e=a_0\int\lambda p_e\,d\lambda,
$$

$$
u_{{\rm LJ},e}
=E\int[\phi(\lambda)-\phi(1)]p_e\,d\lambda,
\qquad
H_{e,k}=\oint_k\sigma_e\,d\bar\lambda_e.
$$

각 C FEM element의 local normal-stress history를 conditional-intact Smoluchowski solver로 넘긴다. 구현은 normalization과 finite-volume Gibbs stationary state를 보존하고, 평균거리·분산·에너지·loop work·$\lambda_c$ 이상 tail을 출력하며 실제 1D node/element mesh 및 2D/3D tensile-only view를 그린다.

이는 **후보 kinetic extension**이며 아직 보정된 aluminum fatigue-life law가 아니다. 현재 upper boundary가 no-flux이므로 critical tail은 irreversible crack probability가 아니라 instability diagnostic이다. 자세한 내용은 `docs/MILESTONE16_PROBABILITY_ENERGY_HYSTERESIS_FEM.md`에 있다.

geometry 계층은 이제 실제 1D line, 2D quad, 3D hex connectivity, 외부 의존성이 없는 STL/OBJ surface import, 선택적 Gmsh 기반 STEP/IGES/BREP meshing을 지원한다. 경량 NumPy/Matplotlib mesh UI에서 node, edge, opacity 및 axial clipping을 조절할 수 있으며 core 설치에 무거운 VTK stack을 추가하지 않는다. 이것은 활성 이론을 바꾸지 않는다. 각 mesh cell에는 선언한 tensile normal scalar $\sigma_{nn}$ 또는 그 cell의 1D $P(a,t)$에서 유도한 scalar만 전달한다. 현재 1D-to-mesh mapping은 2D/3D elasticity solve가 아니라 명시적으로 구분한 visualization/post-processing projection이다. 자세한 내용은 `docs/MILESTONE17_2D_3D_MESH_CAD_NORMAL_ONLY.md`에 있다.

## 활성 1D 통계셀 종속성 척도

<!-- STATISTICAL_CELL_STATUS_KO -->

확률을 합칠 때 완전 동일 종속, 부분 종속, 진짜 독립을 명시적으로 구분한다. second-order stationary 1D spacing process에서는

$$
\operatorname{Var}(\bar\lambda_M)
=\frac{\sigma_\lambda^2}{M}\tau_M,
\qquad
\tau_M=1+2\sum_{k=1}^{M-1}\left(1-\frac{k}{M}\right)\rho_k
$$

가 정확하고, 이에 따라

$$
M_{\rm eff}=\frac{M}{\tau_M},
\qquad
\ell_{\rm stat}^{(2)}=a_0\tau_M
$$

를 정의한다.

하지만 하나의 deterministic finite snapshot을 자기 sample mean으로 center하면 모든 lag를 넣은 weighted correlation sum이 정확히 0이 되는 finite-sample identity가 있으므로 population 식을 그대로 plug-in하면 안 된다. 따라서 snapshot에는 별도로 표시한 first-positive-lobe estimator를 사용한다. Dynamically matched $M=31,63,127,255$ sweep에서는 corrected estimate가 $M_{\rm eff}^{(+)}\approx2.93,2.99,3.03,3.05$이고 $\ell_{\rm stat}^{(2,+)}/a_0\approx10.58,21.10,41.93,83.49$이다. 즉 tested protocol은 local material correlation length로 수렴하지 않고 system-scale coherence를 유지한다.

mechanical calibration area $A_0$를 transverse statistical independence area와 동일시하지 않는다. 활성 범위는 계속 엄격한 1D다.

## 활성 교정 — deterministic correlation snapshot의 준정적 극한

<!-- QUASISTATIC_PROTOCOL_STATUS_KO -->

기존 $M_{\rm eff}^{(+)}\approx3$ 계산 자체는 선택한 deterministic snapshot의 normalized-shape 진단값으로 유효하지만 물리적 해석은 교정되었다. Milestone 13 snapshot은 zero-mean sinusoidal end loading에서 정수 cycle 2에 저장되므로 그 정확한 위상에서 applied force는 0이다.

균질 force-controlled chain에서는

$$
\Pi=\sum_i[\phi(\lambda_i)-f\lambda_i],
\qquad
\phi'(\lambda_i)=f
$$

이고 안정 branch에서 $\phi''>0$이므로 안정 root가 유일하다. 따라서

$$
\lambda_i=\lambda_s(f)\quad\forall i
$$

이다.

현재 snapshot의 zero-force 위상에서는 준정적 상태가 정확히 $\lambda_i=1$, $C_0=0$이다. 새 $\alpha=\omega M$ sweep에서는 drive를 느리게 할수록 deterministic residual variance가 0으로 급격히 감소하지만 normalized correlation shape와 positive-window $M_{\rm eff}^{(+)}$는 약 3을 유지한다. 따라서 normalized residual correlation만으로 물질 고유 statistical-cell length를 정의할 수 없다.

이제 한 trajectory의 spatial empirical measure $P_{M,\mathrm{spatial}}^{\mathrm{traj}}$와 물리적으로 정의해야 할 ensemble probability $P_{\rm ens}$를 명시적으로 구분한다. 다음 활성 목표는 임의 named distribution fitting이 아니라 같은 nonlinear cyclic mechanics 위에서 물리적으로 정당한 1D initial ensemble을 만드는 것이다.

## $P$에 대한 활성 물리 통계역학 제약

<!-- PHYSICAL_P_STATUS_KO -->

현재 mechanics에서는 named distribution을 fitting하지 않고도 가능한 spacing distribution의 물리적 hierarchy를 얻는다.

$$
U(a)=E_0\phi(a/a_0)+C
$$

로 두면 elastic calibration $E=(a_0/A_0)U''(a_0)$와 $\phi''(1)=1$에서

$$
\boxed{E_0=EA_0a_0},
\qquad
\boxed{\chi=\frac{EA_0a_0}{k_BT}}
$$

가 나온다.

zero-temperature homogeneous quasistatic equilibrium에서는 broad distribution이 아니라

$$
\boxed{P(\lambda\mid f)=\delta[\lambda-\lambda_s(f)]}
$$

이다.

fixed total normalized length의 canonical equilibrium에서는 exact finite-$M$ marginal

$$
\boxed{
P_M(\lambda\mid L,\chi)
=\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}
}
$$

을 얻는다.

constant tensile force $f>0$에서는 $\phi(\lambda)-f\lambda\to-\infty$이므로 full-domain Gibbs integral이 발산한다. 따라서 $0<f<f_c$에서 intact basin $0<\lambda<\lambda_b(f)$에 조건부로 둔 Gibbs density는 **controlled metastable/local-equilibrium approximation**일 뿐 global equilibrium law나 fatigue-life law가 아니다. `docs/MILESTONE12_PHYSICAL_STATISTICAL_P.md`를 현재 물리 통계역학 기준으로 사용한다.

representative layer area $A_0$가 이제 명시적인 physical bottleneck이다. coarse-grained layer interaction과 일관된 $A_0$가 정해지기 전에는 numerical aluminum thermal distribution을 주장하지 않는다.

## 활성 이론 정정 — 정확한 비선형 transport

**현재 활성 상태:** 이 README 뒤쪽에 남아 있는 harmonic/single-mode 및 Taylor-expanded push-forward 내용은 historical/local diagnostic이며 $P(\lambda,t)$의 활성 전역 형식으로 사용하지 않는다.

현재 derivation은 원래 nonlinear generalized-LJ force를 그대로 유지하고 exact empirical transport identity

$$
\partial_tP+\partial_\lambda J=0
$$

에서 출발한다. 이후 spacing-velocity phase-space state $F_1(\lambda,v,t)$와 neighboring-spacing joint state $P_2$가 exact hierarchy에 들어간다. `docs/MILESTONE10_EXACT_DISTRIBUTION_TRANSPORT.md`를 현재 기준으로 사용한다.

<!-- ACTIVE_TRANSPORT_CORRECTION_KO -->

고순도 또는 단결정 알루미늄의 **1차원 수직 반복하중** 아래 피로 균열개시를 mechanics-first 방식으로 전개하는 연구 저장소다.

## 활성 범위

활성 derivation은 represented material layer의 1차원 stack으로 의도적으로 제한한다. microscopic reduced coordinate는 수직 layer spacing

$$
a_i(t)>0
$$

또는 normalized spacing

$$
\lambda_i(t)=a_i(t)/a_0
$$

이다.

layer 사이의 유효 normal interaction은 calibration된 generalized Lennard--Jones model로 표현한다. 3차원 FCC 연구는 `libraries/fcc_normal/` 아래 archive로 유지한다. 교정된 단일 registry/Bessel mechanism은 별도 선택 branch로 활성화하며 collinear normal energy에 더하지 않는다.

물리적 시간 $t$가 근본 evolution coordinate이며 fatigue cycle count는 독립 상태변수가 아니다.

## Calibration된 1D layer-LJ mechanics

normalized layer energy는

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
$$

이고

$$
m=12.19,
\qquad
n=6
$$

이다.

calibration은

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

을 준다.

interior spacing 지배방정식은

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1})
}
$$

이다.

현재 이상화된 tangent-instability stretch는

$$
\boxed{
\lambda_c
=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
\approx1.1077715386
}
$$

이다.

## One-point probability state

normalized one-point spacing density는

$$
\int_0^\infty p_\lambda(\lambda,t)\,d\lambda=1
$$

을 만족한다.

평균과 shifted configurational energy는

$$
\mu(t)=\int_0^\infty\lambda p_\lambda(\lambda,t)\,d\lambda,
$$

$$
\psi(\lambda)=\phi(\lambda)-\phi(1),
$$

그리고

$$
\mathcal E(t)=\int_0^\infty\psi(\lambda)p_\lambda(\lambda,t)\,d\lambda
$$

이다.

exact energy-feasibility 연구에서 normalization, mean, energy만으로는 tensile tail을 강제할 수 없다는 것이 나왔다. LJ compression branch가 $\lambda\to0^+$에서 수학적으로 임의로 큰 energy를 담을 수 있기 때문이다. 따라서 exact safe-energy ceiling을 쓰려면 독립적으로 정당화된 compression-side constraint가 필요하다.

## 이전 two-moment distribution closure

fixed-length/fixed-configurational-energy equal-base-measure assumption과 large-$M$ saddle-point reduction으로

$$
\boxed{
p_\lambda(\lambda,t)
=Z^{-1}\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
}
$$

이라는 controlled approximation을 얻었다.

$\alpha$, $\beta$는 histogram fitting이 아니라 $\mu(t)$와 $\mathcal E(t)$로 결정된다.

하지만 deterministic 1D layer-LJ 직접시험에서 이 closure는 near-equilibrium variance는 상당히 잘 맞춰도 driven distribution 전체 shape를 재현하지 못했다. tested 32-node state에서 Kolmogorov distance는 약 $0.15$이고, 느린 case에서는 skewness mismatch가 크게 남았다. 따라서 $\mu$와 $\mathcal E$만으로 driven one-point distribution을 정할 수는 없다.

## Spatial-correlation 결과

deterministic chain은 강한 spatial ordering을 유지한다. 다음을 정의한다.

$$
C_k(t)
=
\frac{1}{M-k}
\sum_{i=1}^{M-k}
[\lambda_i-\mu][\lambda_{i+k}-\mu],
$$

$$
\rho_k=C_k/C_0.
$$

$\omega M=0.62$로 dynamic similarity를 맞춘 sweep에서 nearest-neighbor correlation은 $M=31$의 약 $0.933$에서 $M=255$의 약 $0.991$까지 증가했다. $k/M$에 대해 그리면 profile이 거의 collapse하며 첫 zero crossing은 약 $0.35M$에 있다.

one-point density는 layer label을 permutation해도 정확히 불변이지만 $C_k$는 그렇지 않다. 따라서

$$
\boxed{p_\lambda(\lambda,t)\text{만으로 complete spatial mechanical state를 담을 수 없다}}
$$

는 구조적 결과를 얻었다.

## $P$의 형식에 대한 새로운 지배방정식 단서

one-point density는 deterministic spacing field의 정확한 spatial push-forward이기도 하다. continuum layer label $\xi\in[0,1]$과 spacing field $\Lambda(\xi,t)$를 쓰면

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1\delta[\lambda-\Lambda(\xi,t)]\,d\xi
}
$$

이다.

field가 piecewise monotone이면

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}
}
$$

이다.

즉 $P$의 함수형식은 독립적으로 고른 probability family가 아니라 mechanically generated spatial waveform의 kinematic image라는 직접적인 단서를 얻는다.

$\lambda=1$ 주변에서 선형화하면

$$
\ddot u_i=u_{i+1}-2u_i+u_{i-1}
$$

이고 dispersion은

$$
\boxed{\omega_q^2=4\sin^2(q/2)}
$$

이다.

하나의 coherent linear mode는 따라서 정확히

$$
\boxed{
p_{\rm 1mode}(\lambda)
=
\frac{\mathbf 1_{|\lambda-\mu|<|A|}}
{\pi\sqrt{A^2-(\lambda-\mu)^2}}
}
$$

이라는 spatial push-forward를 만든다.

이 arcsine form은 fitted distribution이 아니라 $\Lambda=\mu+A\cos\vartheta$와 uniform spatial phase에서 직접 나온다.

calibration된 LJ force는 강하게 nonlinear하다.

$$
\boxed{
\phi'(1+u)
=u-10.595u^2+62.97935u^3+O(u^4)
}
$$

이다.

따라서 base mechanics 자체가 harmonic distortion, skewness, 그리고 잠재적인 multimodal one-point density를 만든다.

first+second harmonic waveform

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta
$$

에서는

$$
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
$$

$$
\mu_3=\frac34A^2B,
$$

그리고 정확히

$$
\boxed{|\gamma_1|\le\sqrt{2/3}\approx0.816497}
$$

이다.

느린 deterministic snapshot의 $\gamma_1\approx1.062$는 이 bound도 넘기 때문에 two spatial harmonic만으로도 충분하지 않다. single-mode arcsine reference 역시 이전 exponential closure보다 KS error가 더 컸다. 이는 실패결과이며, 새로운 named probability family를 fitting해야 한다는 뜻이 아니라 mechanically generated mode content가 더 풍부하다는 뜻이다.

## 현재 연구방향

현재 hierarchy는

$$
\boxed{
\text{1D layer-LJ governing equation}
\rightarrow
\Lambda(\xi,t)\text{ / neighboring-spacing pair state}
\rightarrow
p_\lambda(\lambda,t)\text{ by push-forward}
\rightarrow
Q_c(t)
\rightarrow
\text{normal-opening first passage}
}
$$

이다.

다음 목표는 1D 지배방정식에서 mechanically excited spatial-mode 또는 pair-state evolution을 직접 유도하고 empirical relaxation law나 fitted probability family 없이 one-point distribution $P_1$과 관찰된 $C_k$를 동시에 재현하는 최소 representation을 찾는 것이다.

## 활성 파일

- `theory/normal_lj_chain.py` — conservative 1D layer-LJ chain
- `theory/normal_lj_energy_feasibility.py` — stated compression constraint 아래 exact energy-feasibility bound
- `theory/normal_lj_distribution.py` — 이전 two-moment large-$M$ closure
- `theory/normal_lj_spatial_correlation.py` — exact finite-chain correlation diagnostic
- `theory/normal_lj_pushforward.py` — spatial push-forward 및 mode-derived distribution identity
- `simulations/run_normal_lj_spatial_correlation.py` — dynamically matched correlation sweep
- `simulations/run_normal_lj_pushforward_clue.py` — single-mode falsification diagnostic
- `docs/MILESTONE8_SPATIAL_CORRELATION.md` — spatial-ordering result
- `docs/MILESTONE9_GOVERNING_EQUATION_PUSHFORWARD.md` — $P$ 형식에 대한 지배방정식 단서
- `results/data/result_manifest.json` — 현재 machine-readable research state

## 활성 결과 재현

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

## 연구 분류 라벨

중요한 statement는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**
- **NUMERICAL RESULT / DIAGNOSTIC**

fitted Gaussian/Weibull distribution, cycle-dependent LJ parameter, empirical damage variable, fitted correlation length은 mechanics derivation으로 인정하지 않는다.
