# Variable Definitions

This document defines the variables currently used in the theory, simulations, and research notes. The goal is to prevent the same symbol from changing meaning across derivations and to make every reduced coordinate traceable to a microscopic definition.

Important classifications used below:

- **EXACT / IDENTITY** — exact under the stated microscopic model.
- **DEFINITION** — a chosen mathematical definition.
- **ASSUMPTION** — a modeling assumption.
- **CONTROLLED APPROXIMATION** — a simplification that must be tested.
- **EMPIRICAL INPUT** — a measured or externally supplied material quantity.

## 1. Microscopic geometry and spacing variables

| Symbol | Definition | Physical meaning | Typical unit | Classification |
|---|---|---|---|---|
| $i,j,k$ | Integer indices | Atom, bond, neighbor, or lattice-site indices depending on context | dimensionless | DEFINITION |
| $N$ | Number of sampled local spacing states or atoms in a finite approximation | Finite-system size before the thermodynamic limit | dimensionless | DEFINITION |
| $N_{\rm RA}$ | Number of representative areas in an ensemble approximation | Finite approximation to a population of local structural regions | dimensionless | DEFINITION |
| $x_i(t)$ | Atomic or lattice-site coordinate | Position of site $i$ along a reduced one-dimensional chain | m | DEFINITION |
| $\mathbf r_i(t)$ | Atomic position vector | Full spatial position of atom $i$ | m | DEFINITION |
| $a_i(t)$ | Local spacing descriptor, often $a_i=x_{i+1}-x_i$ in a 1D chain | Instantaneous local interatomic spacing represented in the reduced theory | m | DEFINITION |
| $a$ | Continuous spacing-space coordinate | Argument of $P(a,t)$ | m | DEFINITION |
| $a_0$ | Reference/equilibrium spacing | Equilibrium spacing of the reference crystal state | m | EMPIRICAL INPUT or atomistic-calibration output |
| $a_c$ | Candidate critical spacing | A provisional spacing threshold associated with local instability | m | ASSUMPTION unless derived from a stated stability condition |
| $\bar a$ | $\langle a\rangle$ | Mean local spacing | m | DEFINITION |
| $R_i^{(k)}$ | $a_i+a_{i+1}+\cdots+a_{i+k-1}$ | Distance from site $i$ to its $k$-th neighbor in the reduced chain | m | EXACT / IDENTITY under the chain geometry |
| $r$ | Continuous pair-distance coordinate | Argument of pair-distance density $P_k(r,t)$ or pair potential $v(r)$ | m | DEFINITION |

## 2. Probability and state-density variables

### Thermodynamic-limit spacing density

For finite $N$,

$$
P_N(a,t)=\frac{1}{N}\sum_{i=1}^{N}\delta\!\left(a-a_i(t)\right).
$$

The central state density is

$$
\boxed{
P(a,t)=\lim_{N\to\infty}P_N(a,t)
}
$$

when the thermodynamic/statistical limit exists.

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $P_N(a,t)$ | Finite empirical spacing density | Finite-system distribution of local spacings | m$^{-1}$ |
| $P(a,t)$ | Thermodynamic-limit spacing density | Population density of local spacing states | m$^{-1}$ |
| $P_k(r,t)$ | Density of $k$-th-neighbor distance $R_i^{(k)}$ | Exact pair-distance hierarchy entry | m$^{-1}$ |
| $P^{*k}$ | $k$-fold convolution of $P$ | Approximation to $P_k$ if adjacent spacings are independent | m$^{-1}$ |
| $P_s(s,t)$ | Population density of non-affine slip coordinate | Ensemble density in slip/disregistry space | m$^{-1}$ if $s$ has units of length |
| $P(a,s,t)$ | Joint density of spacing and non-affine slip | Candidate minimal structural state retaining normal and shear information | m$^{-2}$ if both $a$ and $s$ are lengths |
| $F(a,c,t)$ | Joint density of spacing $a$ and spacing velocity $c$ | Phase-space lift required because microscopic mechanics is second order | s/m$^2$ |
| $Q_c(t)$ | $\int_{a_c}^{\infty}P(a,t)\,da$ | Instantaneous fraction beyond a spacing threshold | dimensionless |

Normalization conditions are

$$
\int P(a,t)\,da=1,
$$

$$
\iint P(a,s,t)\,da\,ds=1.
$$

## 3. Spacing-space kinematics

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $v(a,t)$ | $\langle\dot a_i\mid a_i=a\rangle$ | Conditional mean velocity in spacing space | m/s |
| $c$ | $\dot a$ | Local spacing velocity used in the phase-space density | m/s |
| $A(a,c,t)$ | $\langle\ddot a_i\mid a_i=a,\dot a_i=c\rangle$ | Conditional spacing acceleration | m/s$^2$ |
| $v_a$ | Conditional transport velocity in the $a$ direction of a joint state | Normal/spacing component of state-space transport | m/s |
| $v_s$ | Conditional transport velocity in the $s$ direction | Non-affine/shear component of state-space transport | m/s |

The exact spacing continuity equation is

$$
\boxed{
\partial_tP+\partial_a(Pv)=0.
}
$$

For a joint state,

$$
\boxed{
\partial_tP+\partial_a(Pv_a)+\partial_s(Pv_s)=0.
}
$$

Here the symbol $P$ in the second equation means the joint density $P(a,s,t)$.

## 4. Macroscopic loading variables

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $t$ | Time | Physical or nondimensional simulation time depending on the model | s unless explicitly nondimensional |
| $N$ in $P_N$ | Finite sample count | Do not confuse with fatigue cycle number | dimensionless |
| $N_{\rm cyc}$ or cycle index $N$ | Integer cycle count | Number of applied loading cycles | dimensionless |
| $T$ | $2\pi/\omega$ | Loading period | s |
| $f$ | Cyclic frequency | Cycles per second | Hz |
| $\omega$ | $2\pi f$ | Angular loading frequency | rad/s |
| $\sigma(t)$ | Applied macroscopic normal stress | Prescribed stress history | Pa |
| $\sigma_m$ | Mean stress | Mean value of a cyclic stress history | Pa |
| $\sigma_a$ | Stress amplitude | Half-range for a sinusoidal stress history | Pa |
| $\epsilon(t)$ | Macroscopic strain | Work-conjugate strain for $\sigma$ under the chosen measure | dimensionless |
| $\tau(t)$ | Resolved shear stress | Shear component acting on a slip system | Pa |
| $m_s$ | Schmid factor | Geometric factor relating uniaxial stress to resolved shear stress | dimensionless |
| $F(t)$ | Generalized force in reduced models | External force conjugate to a resolved coordinate such as $Q$ or $s$ | N or nondimensional, depending on model |
| $F_a$ | Generalized-force amplitude | Amplitude of $F(t)$ | same unit as $F$ |
| $F_m$ | Mean generalized force | Mean value of $F(t)$ | same unit as $F$ |

For sinusoidal stress,

$$
\sigma(t)=\sigma_m+\sigma_a\sin(\omega t).
$$

## 5. Affine and non-affine decomposition

A useful decomposition is

$$
a=\lambda(t)x,
$$

where $\lambda$ describes macroscopic affine stretch and $x$ describes the internal spacing coordinate after removing that affine stretch.

The corresponding density transformation is

$$
P(a,t)=\frac{1}{\lambda(t)}R\!\left(\frac{a}{\lambda(t)},t\right).
$$

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $\lambda(t)$ | Macroscopic stretch ratio | Affine part of deformation; for engineering strain one may use $\lambda=1+\epsilon$ | dimensionless |
| $x$ | $a/\lambda$ | Internal spacing coordinate with affine deformation removed | m |
| $R(x,t)$ | Internal spacing density | Shape/state density in the non-affine spacing coordinate | m$^{-1}$ |
| $w(x,t)$ | $\dot x$ | Non-affine internal spacing velocity | m/s |

Pure reversible affine elasticity corresponds to

$$
R(x,t)=R_0(x),
$$

whereas irreversible structural evolution requires the internal distribution itself to change.

## 6. Non-affine slip/disregistry coordinate $s$

### 6.1 Physical definition

The non-affine slip coordinate is **not** an empirical damage variable and is **not** introduced as a phenomenological plastic strain.

Choose a candidate FCC slip plane, typically a local $\{111\}$ plane, and define two atom groups immediately above and below that plane. Let their coarse displacement vectors relative to a reference lattice be

$$
\overline{\mathbf u}^{+}(t)
$$

and

$$
\overline{\mathbf u}^{-}(t).
$$

For slip system $\alpha$, let $\hat{\mathbf b}^{\alpha}$ be the unit vector in the corresponding crystallographic slip direction. Then define

$$
\boxed{
s^{\alpha}(t)
=
\left[
\overline{\mathbf u}^{+}(t)
-
\overline{\mathbf u}^{-}(t)
\right]
\cdot
\hat{\mathbf b}^{\alpha}.
}
$$

Thus $s^{\alpha}$ measures the relative tangential displacement across a selected slip plane after rigid-body translation has been removed by the use of relative displacement.

For a single reduced slip system the superscript may be dropped:

$$
s(t)\equiv s^{\alpha}(t).
$$

### 6.2 Why it is called non-affine

A homogeneous affine deformation can be represented by a deformation gradient $\mathbf F_{\rm mac}$. The displacement expected from purely affine deformation is subtracted from the local atomic displacement before constructing the disregistry. Therefore $s$ is intended to represent the local displacement that cannot be explained by the macroscopic affine deformation alone.

A more explicit atomistic implementation may use

$$
\mathbf u_i^{\rm na}
=
\mathbf r_i(t)
-
\mathbf F_{\rm mac}(t)\mathbf r_i^{0},
$$

followed by coarse averages of $\mathbf u_i^{\rm na}$ above and below the slip plane. Here $\mathbf r_i^{0}$ denotes a reference position.

This precise implementation remains to be finalized for the 3D FCC simulation. Until then, the scalar $s$ used in the current Hamiltonian proof-of-principle is a **CONTROLLED APPROXIMATION** to this atomistically defined disregistry coordinate.

### 6.3 Variables related to slip

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $s^{\alpha}$ | Relative non-affine displacement projected onto slip direction $\alpha$ | Local disregistry/slip coordinate | m |
| $s$ | Single-slip reduced version of $s^{\alpha}$ | Resolved non-affine coordinate used in current proof-of-principle simulation | m or nondimensional |
| $\mathbf s$ | Vector of one or more slip/disregistry coordinates | Multi-slip internal coordinate | m |
| $\hat{\mathbf b}^{\alpha}$ | Unit vector along slip direction | Crystallographic slip direction | dimensionless |
| $\mathbf b^{\alpha}$ | Burgers vector | Lattice translation associated with slip system $\alpha$ | m |
| $b=|\mathbf b|$ | Burgers-vector magnitude / periodicity scale | Period of an ideal slip energy landscape | m |
| $\overline{\mathbf u}^{+}$ | Mean non-affine displacement above a selected slip plane | Upper-side coarse displacement | m |
| $\overline{\mathbf u}^{-}$ | Mean non-affine displacement below a selected slip plane | Lower-side coarse displacement | m |
| $\alpha$ | Slip-system index | Labels crystallographic slip systems | dimensionless |

## 7. Gamma-surface and structural potential variables

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $\gamma_{\rm Al}(\mathbf s)$ | Generalized stacking-fault energy per area | Atomistic energy landscape for relative displacement across a slip plane | J/m$^2$ |
| $V_\gamma(s)$ | Reduced potential associated with slip | Energy landscape used for the resolved slip coordinate | J or nondimensional |
| $\Delta_\gamma$ | Barrier amplitude in the one-harmonic approximation | Artificial barrier scale in the current proof-of-principle model | J or nondimensional |
| $A_{\rm RA}$ | Representative slip area | Area over which a local gamma-surface energy is coarse-grained | m$^2$ |
| $\Phi(\mathbf s,t)$ | Driven structural energy landscape | Gamma-surface energy minus work of resolved shear loading | J |
| $\lambda_{\min}$ | Smallest eigenvalue of a Hessian | Local mechanical-stability indicator | depends on coordinate convention |

Current one-harmonic approximation:

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right].
$$

Desired Al-specific replacement:

$$
V_\gamma(\mathbf s)=A_{\rm RA}\gamma_{\rm Al}(\mathbf s).
$$

A candidate local stability-loss condition is

$$
\lambda_{\min}
\left[
\nabla_{\mathbf s}^{2}\Phi
\right]=0.
$$

## 8. Interatomic-potential and energy variables

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $v(r)$ | Pair potential | Microscopic interaction energy of a pair separated by $r$ | J |
| $U(a)$ | $\sum_{k=1}^{\infty}v(ka)$ for a uniform chain | Uniform-lattice energy per atom in the semi-infinite reduced model | J |
| $\mathcal U(t)$ | $\sum_k\int v(r)P_k(r,t)\,dr$ | Distribution-level configurational pair-potential energy per atom | J |
| $E_{\rm int}$ | Internal kinetic + configurational energy excluding external loading potential | Energy used in conservation tests | J or nondimensional |
| $W_{\rm ext}$ | $\int F\,dq$ or continuum equivalent | Mechanical work done by external loading | J or J/m$^3$ |
| $A_H$ | $\oint\sigma\,d\epsilon$ or reduced-coordinate equivalent | Hysteresis-loop work per cycle | J/m$^3$ in continuum stress-strain form |

For the reduced structural coordinate $s$,

$$
A_H=\oint F\,ds.
$$

For macroscopic stress and strain,

$$
A_H=\oint\sigma\,d\epsilon.
$$

## 9. Generalized Lennard-Jones variables

When the generalized Lennard-Jones baseline is used,

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $\varepsilon_{\rm LJ}$ | Energy scale | Pair-potential energy scale | J or eV |
| $\sigma_{\rm LJ}$ | Length scale | Characteristic LJ distance parameter | m |
| $m,n$ | Repulsive/attractive exponents with $m>n>1$ | Shape parameters of the generalized LJ interaction | dimensionless |
| $\zeta(\cdot)$ | Riemann zeta function | Appears when summing pair interactions in a semi-infinite uniform chain | dimensionless |
| $A_0$ | Reference effective cross-sectional area in the 1D-to-stress mapping | Converts generalized force to engineering stress | m$^2$ |
| $E$ | Young's modulus | Small-strain elastic modulus used in calibration | Pa |
| $\sigma_u$ | Ultimate/instability stress used in the old calibration | Calibration target in the previous 1D LJ model | Pa |
| $a_u$ | Spacing at the old instability/ultimate point | Location where $U''(a_u)=0$ in the old model | m |

## 10. Rubin-chain hysteresis variables

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $Q(t)$ | Resolved structural coordinate | Coarse coordinate coupled to a harmonic lattice bath | m or nondimensional |
| $P_Q$ | Momentum conjugate to $Q$ | Resolved-coordinate momentum | kg m/s |
| $M$ | Mass associated with the resolved coordinate | Effective inertia of $Q$ or $s$ | kg or nondimensional |
| $x_n$ | Harmonic-bath displacement of site $n$ | Unresolved lattice coordinate | m or nondimensional |
| $p_n$ | Momentum of bath site $n$ | Bath momentum | kg m/s |
| $m$ | Bath-site mass | Atomic/effective lattice mass | kg or nondimensional |
| $k$ | Harmonic nearest-neighbor spring constant | Bath coupling stiffness | N/m or nondimensional |
| $K_0$ | On-site stiffness of $Q$ | Local restoring stiffness of resolved coordinate | N/m or nondimensional |
| $q$ | Lattice wave number in the reduced chain | Phase increment per lattice site | rad/site |
| $\omega_D$ | $2\sqrt{k/m}$ | Upper frequency band edge of the monatomic harmonic chain | rad/s or nondimensional |
| $Z(\omega)$ | $\hat F/\hat Q$ | Complex dynamic stiffness seen by the resolved coordinate | N/m |
| $\hat Q$ | Complex harmonic amplitude of $Q$ | Frequency-domain response amplitude | m |
| $\hat F$ | Complex harmonic amplitude of $F$ | Frequency-domain forcing amplitude | N |
| $\phi$ | $\arg Z$ or the corresponding response phase lag | Phase lag between forcing and resolved displacement | rad |

Inside the propagating band, the imaginary part of $Z$ represents energy radiation into unresolved lattice modes; no viscous coefficient is inserted into the full microscopic model.

## 11. Hamiltonian slip-bath variables

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $u_j$ | Harmonic bath coordinate coupled to $s$ | Unresolved lattice displacement | m or nondimensional |
| $p_j$ | Momentum conjugate to $u_j$ | Bath momentum | kg m/s or nondimensional |
| $k_c$ | Coupling stiffness between $s$ and $u_1$ | Mechanical coupling of the resolved slip coordinate to the bath | N/m or nondimensional |
| $P_s$ | Momentum conjugate to $s$ | Resolved slip-coordinate momentum | kg m/s or nondimensional |
| $s_N$ | $s(NT)$ | Slip-coordinate state sampled at the end of loading cycle $N$ | same unit as $s$ |
| $\Delta s_{\rm cycle}$ | $s_{N+1}-s_N$ | Cycle-to-cycle structural drift | same unit as $s$ |

A periodic internal-friction state satisfies

$$
s_{N+1}=s_N,
$$

whereas a secular structural state satisfies

$$
s_{N+1}\neq s_N.
$$

## 12. Crack-initiation and first-passage variables

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $\tau_c$ | $\inf\{t>0:\text{state enters the chosen unstable set}\}$ | First-passage time to a crack-initiation criterion | s |
| $F_{\rm ci}(t)$ | $\Pr(\tau_c\le t)$ | Cumulative crack-initiation probability | dimensionless |
| $S(t)$ | $1-F_{\rm ci}(t)$ | Survival probability | dimensionless |
| $h(t)$ | $-\dot S/S$ | Hazard rate | s$^{-1}$ |
| $j_c(t)$ | Probability flux through an absorbing instability boundary | Instantaneous first-passage flux | s$^{-1}$ |
| $N_i$ | First cycle number satisfying the chosen crack-initiation condition | Predicted crack-initiation life in cycles | cycles |

These quantities are **definitions of a first-passage framework**, not yet a calibrated crack-initiation law for aluminum.

## 13. Moments and statistical descriptors

| Symbol | Definition | Physical meaning | Unit |
|---|---|---|---|
| $\langle g(a)\rangle$ | $\int g(a)P(a,t)\,da$ | Population average | depends on $g$ |
| $\mu_2$ | $\langle(a-\bar a)^2\rangle$ | Variance of local spacing | m$^2$ |
| $\mu_3$ | $\langle(a-\bar a)^3\rangle$ | Third central moment / asymmetry indicator | m$^3$ |
| $\operatorname{Var}(a)$ | $\mu_2$ | Spacing-distribution width | m$^2$ |
| $\operatorname{Cov}(a,v)$ | $\langle(a-\bar a)(v-\bar v)\rangle$ | Controls deterministic broadening or narrowing of $P$ | m$^2$/s |

An exact moment identity is

$$
\frac{d}{dt}\operatorname{Var}(a)
=2\operatorname{Cov}(a,v).
$$

## 14. Symbols that are intentionally not yet fixed

The following quantities require future microscopic definitions or calibration and must not be silently assigned arbitrary values:

- the Representative Area geometry and therefore $A_{\rm RA}$;
- the exact 3D atomistic definition of $s^{\alpha}$ for FCC Al;
- the Al-specific $\gamma_{\rm Al}(\mathbf s)$ surface used in the final model;
- the correct local instability set for crack initiation;
- the minimal correlation hierarchy required to close $P(a,s,t)$;
- the physically relevant mapping from reduced generalized force $F$ to macroscopic stress for a given RA;
- surface, defect, thermal, and multi-slip amplification mechanisms needed to bridge ordinary fatigue stresses to atomistic instability scales.

Any future document that introduces a new symbol should add it to this file in the same commit.

---

# 한국어 번역 — 변수 정의

이 문서는 현재 이론, 시뮬레이션, 연구노트에 사용되는 변수들을 정의한다. 목적은 서로 다른 유도에서 같은 기호의 의미가 바뀌는 것을 막고, 모든 축약 좌표가 가능한 한 미시적 정의까지 추적될 수 있도록 하는 것이다.

아래에서 사용하는 분류는 다음과 같다.

- **EXACT / IDENTITY** — 명시된 미시모델 아래에서 정확히 성립.
- **DEFINITION** — 선택한 수학적 정의.
- **ASSUMPTION** — 모델링 가정.
- **CONTROLLED APPROXIMATION** — 검증이 필요한 통제된 근사.
- **EMPIRICAL INPUT** — 실험 또는 외부 자료에서 주어지는 물성 입력.

## 1. 미시 기하와 원자간격 변수

| 기호 | 정의 | 물리적 의미 | 대표 단위 | 분류 |
|---|---|---|---|---|
| $i,j,k$ | 정수 index | 문맥에 따라 원자, bond, 이웃, lattice site index | 무차원 | DEFINITION |
| $N$ | 유한계에서의 원자 또는 spacing sample 수 | 열역학적 극한을 취하기 전의 finite-system size | 무차원 | DEFINITION |
| $N_{\rm RA}$ | 대표영역 population의 유한 sample 수 | 국부 구조영역 ensemble의 유한 근사 | 무차원 | DEFINITION |
| $x_i(t)$ | 원자 또는 lattice site의 1D 좌표 | 축약 사슬에서 site $i$의 위치 | m | DEFINITION |
| $\mathbf r_i(t)$ | 원자 위치 vector | 3차원에서 원자 $i$의 위치 | m | DEFINITION |
| $a_i(t)$ | 예: $a_i=x_{i+1}-x_i$ | 축약 이론에서 표현되는 순간 국부 원자간격 | m | DEFINITION |
| $a$ | 연속적인 spacing-space 좌표 | $P(a,t)$의 독립변수 | m | DEFINITION |
| $a_0$ | 기준/평형 spacing | 기준 결정의 평형 원자간격 | m | EMPIRICAL INPUT 또는 원자모델 calibration 결과 |
| $a_c$ | 후보 임계 spacing | 국부 기계적 불안정성과 연결되는 임시 threshold | m | 명시된 안정성 조건에서 유도되지 않으면 ASSUMPTION |
| $\bar a$ | $\langle a\rangle$ | 평균 국부 원자간격 | m | DEFINITION |
| $R_i^{(k)}$ | $a_i+\cdots+a_{i+k-1}$ | $i$ site에서 $k$번째 이웃까지 거리 | m | 사슬 기하에서 EXACT / IDENTITY |
| $r$ | 연속 pair-distance 좌표 | $P_k(r,t)$ 또는 $v(r)$의 독립변수 | m | DEFINITION |

## 2. 확률 및 상태밀도 변수

유한 $N$에서

$$
P_N(a,t)=\frac{1}{N}\sum_{i=1}^{N}\delta\!\left(a-a_i(t)\right)
$$

로 정의하고, 중심 상태밀도는

$$
\boxed{
P(a,t)=\lim_{N\to\infty}P_N(a,t)
}
$$

이다.

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $P_N(a,t)$ | finite empirical spacing density | 유한계의 국부 spacing 분포 | m$^{-1}$ |
| $P(a,t)$ | thermodynamic-limit spacing density | 국부 spacing 상태의 population density | m$^{-1}$ |
| $P_k(r,t)$ | $k$번째 이웃거리의 밀도 | 정확한 pair-distance hierarchy 항 | m$^{-1}$ |
| $P^{*k}$ | $P$의 $k$중 convolution | 인접 spacing이 독립일 때의 $P_k$ 근사 | m$^{-1}$ |
| $P_s(s,t)$ | non-affine slip coordinate의 population density | slip/disregistry 공간의 ensemble density | $s$가 길이면 m$^{-1}$ |
| $P(a,s,t)$ | spacing과 non-affine slip의 joint density | normal + shear 구조정보를 함께 보유하는 후보 최소상태 | 둘 다 길이면 m$^{-2}$ |
| $F(a,c,t)$ | spacing과 spacing velocity의 joint density | 원자역학의 2계 성질을 반영한 phase-space lift | s/m$^2$ |
| $Q_c(t)$ | $\int_{a_c}^{\infty}P(a,t)\,da$ | spacing threshold를 넘은 순간 population fraction | 무차원 |

정규화는

$$
\int P(a,t)\,da=1,
$$

$$
\iint P(a,s,t)\,da\,ds=1
$$

을 만족해야 한다.

## 3. spacing-space 운동학 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $v(a,t)$ | $\langle\dot a_i\mid a_i=a\rangle$ | spacing space에서의 conditional mean velocity | m/s |
| $c$ | $\dot a$ | phase-space density에 쓰이는 local spacing velocity | m/s |
| $A(a,c,t)$ | $\langle\ddot a_i\mid a_i=a,\dot a_i=c\rangle$ | conditional spacing acceleration | m/s$^2$ |
| $v_a$ | joint state의 $a$ 방향 transport velocity | normal/spacing 성분의 state-space transport | m/s |
| $v_s$ | joint state의 $s$ 방향 transport velocity | non-affine/shear 성분의 state-space transport | m/s |

정확한 spacing continuity equation은

$$
\boxed{
\partial_tP+\partial_a(Pv)=0
}
$$

이다.

joint state에서는

$$
\boxed{
\partial_tP+\partial_a(Pv_a)+\partial_s(Pv_s)=0
}
$$

가 된다. 두 번째 식의 $P$는 $P(a,s,t)$를 뜻한다.

## 4. 거시 하중 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $t$ | 시간 | 명시가 없으면 물리시간 | s |
| $N$ in $P_N$ | finite sample count | fatigue cycle number와 혼동 금지 | 무차원 |
| cycle index $N$ | 정수 cycle count | 반복하중 cycle 수 | cycle |
| $T$ | $2\pi/\omega$ | 하중주기 | s |
| $f$ | cyclic frequency | 초당 cycle 수 | Hz |
| $\omega$ | $2\pi f$ | angular frequency | rad/s |
| $\sigma(t)$ | applied normal stress | 외부 거시 응력 history | Pa |
| $\sigma_m$ | mean stress | 반복응력의 평균값 | Pa |
| $\sigma_a$ | stress amplitude | sinusoidal stress amplitude | Pa |
| $\epsilon(t)$ | macroscopic strain | 선택한 stress measure의 work-conjugate strain | 무차원 |
| $\tau(t)$ | resolved shear stress | slip system에 작용하는 전단응력 | Pa |
| $m_s$ | Schmid factor | uniaxial stress를 resolved shear로 변환하는 기하계수 | 무차원 |
| $F(t)$ | reduced model의 generalized force | $Q$ 또는 $s$와 work-conjugate한 외력 | N 또는 무차원 |
| $F_a$ | generalized-force amplitude | $F(t)$의 진폭 | $F$와 동일 |
| $F_m$ | mean generalized force | $F(t)$의 평균값 | $F$와 동일 |

sinusoidal stress는

$$
\sigma(t)=\sigma_m+\sigma_a\sin(\omega t)
$$

로 쓴다.

## 5. affine / non-affine 분해

$$
a=\lambda(t)x
$$

로 두어 $\lambda$는 거시 affine stretch, $x$는 affine 성분을 제거한 internal spacing coordinate로 해석한다.

분포는

$$
P(a,t)=\frac{1}{\lambda(t)}R\!\left(\frac{a}{\lambda(t)},t\right)
$$

로 변환된다.

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $\lambda(t)$ | macroscopic stretch ratio | 변형의 affine 부분 | 무차원 |
| $x$ | $a/\lambda$ | affine 변형을 제거한 internal spacing coordinate | m |
| $R(x,t)$ | internal spacing density | non-affine spacing coordinate에서의 상태밀도 | m$^{-1}$ |
| $w(x,t)$ | $\dot x$ | non-affine internal spacing velocity | m/s |

완전히 가역적인 affine elasticity에서는

$$
R(x,t)=R_0(x)
$$

이고, 비가역 구조진화가 있으려면 internal distribution 자체가 변해야 한다.

## 6. non-affine slip/disregistry coordinate $s$

### 6.1 물리적 정의

non-affine slip coordinate는 경험적 damage variable이 아니며 phenomenological plastic strain으로 정의한 변수도 아니다.

FCC의 후보 slip plane, 일반적으로 국부 $\{111\}$ 면을 선택하고, 그 면 바로 위와 아래에 있는 원자집단을 정의한다. 기준격자 대비 두 집단의 coarse displacement vector를

$$
\overline{\mathbf u}^{+}(t),
\qquad
\overline{\mathbf u}^{-}(t)
$$

라고 하자.

slip system $\alpha$의 slip direction unit vector를 $\hat{\mathbf b}^{\alpha}$라고 하면

$$
\boxed{
s^{\alpha}(t)
=
\left[
\overline{\mathbf u}^{+}(t)
-
\overline{\mathbf u}^{-}(t)
\right]
\cdot
\hat{\mathbf b}^{\alpha}
}
$$

로 정의한다.

즉 $s^{\alpha}$는 선택한 slip plane 양쪽 원자집단의 **상대적인 접선방향 displacement**, 즉 local disregistry를 나타낸다.

single-slip reduced model에서는

$$
s(t)\equiv s^{\alpha}(t)
$$

로 superscript를 생략한다.

### 6.2 왜 non-affine인가

거시적인 homogeneous affine deformation은 deformation gradient $\mathbf F_{\rm mac}$로 표현할 수 있다. local atomic displacement에서 이 affine prediction을 제거한 뒤 slip-plane 양쪽의 상대변위를 계산하면, 거시 affine deformation만으로 설명되지 않는 local displacement만 남는다.

보다 명시적인 3D 구현 후보는

$$
\mathbf u_i^{\rm na}
=
\mathbf r_i(t)
-
\mathbf F_{\rm mac}(t)\mathbf r_i^{0}
$$

를 먼저 정의한 뒤, slip plane 위/아래에서 $\mathbf u_i^{\rm na}$를 coarse average하는 방식이다.

여기서 $\mathbf r_i^{0}$는 기준 원자위치다.

이 3D FCC 구현은 아직 최종 확정 전이다. 따라서 현재 Hamiltonian proof-of-principle에서 사용하는 scalar $s$는 이 원자수준 disregistry coordinate의 **CONTROLLED APPROXIMATION**으로 분류한다.

### 6.3 slip 관련 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $s^{\alpha}$ | slip direction $\alpha$로 투영한 relative non-affine displacement | local disregistry/slip coordinate | m |
| $s$ | single-slip reduced $s^{\alpha}$ | 현재 proof-of-principle simulation의 resolved non-affine coordinate | m 또는 무차원 |
| $\mathbf s$ | 하나 이상의 slip/disregistry coordinate vector | multi-slip internal coordinate | m |
| $\hat{\mathbf b}^{\alpha}$ | slip direction unit vector | 결정학적 slip 방향 | 무차원 |
| $\mathbf b^{\alpha}$ | Burgers vector | slip system의 lattice translation | m |
| $b$ | $|\mathbf b|$ | Burgers-vector magnitude 또는 ideal slip periodicity | m |
| $\overline{\mathbf u}^{+}$ | slip plane 위쪽의 평균 non-affine displacement | upper-side coarse displacement | m |
| $\overline{\mathbf u}^{-}$ | slip plane 아래쪽의 평균 non-affine displacement | lower-side coarse displacement | m |
| $\alpha$ | slip-system index | 결정학적 slip system label | 무차원 |

## 7. gamma-surface 및 구조 potential 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $\gamma_{\rm Al}(\mathbf s)$ | generalized stacking-fault energy per area | slip plane 상대변위에 대한 Al의 원자수준 energy landscape | J/m$^2$ |
| $V_\gamma(s)$ | slip 관련 reduced potential | resolved slip coordinate의 energy landscape | J 또는 무차원 |
| $\Delta_\gamma$ | one-harmonic approximation의 barrier amplitude | 현재 원리증명 모델의 인공 barrier scale | J 또는 무차원 |
| $A_{\rm RA}$ | representative slip area | gamma-surface energy를 coarse-grain하는 국부 면적 | m$^2$ |
| $\Phi(\mathbf s,t)$ | driven structural energy landscape | gamma-surface energy에서 resolved shear work를 뺀 energy | J |
| $\lambda_{\min}$ | Hessian의 최소 eigenvalue | 국부 기계적 안정성 indicator | 좌표 정의에 따라 달라짐 |

현재 one-harmonic approximation은

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right]
$$

이고, 최종 목표는

$$
V_\gamma(\mathbf s)=A_{\rm RA}\gamma_{\rm Al}(\mathbf s)
$$

로 교체하는 것이다.

후보 안정성 상실 조건은

$$
\lambda_{\min}\left[\nabla_{\mathbf s}^{2}\Phi\right]=0
$$

이다.

## 8. 원자간 potential 및 energy 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $v(r)$ | pair potential | 거리 $r$인 원자쌍의 microscopic interaction energy | J |
| $U(a)$ | $\sum_kv(ka)$ | 균일 semi-infinite reduced lattice의 원자당 energy | J |
| $\mathcal U(t)$ | $\sum_k\int v(r)P_k(r,t)\,dr$ | distribution-level configurational pair-potential energy per atom | J |
| $E_{\rm int}$ | external loading potential을 제외한 kinetic + configurational energy | 보존검증에 사용되는 내부에너지 | J 또는 무차원 |
| $W_{\rm ext}$ | $\int F\,dq$ 또는 continuum equivalent | 외부하중이 한 mechanical work | J 또는 J/m$^3$ |
| $A_H$ | $\oint\sigma\,d\epsilon$ 또는 reduced equivalent | cycle당 hysteresis-loop work | continuum에서는 J/m$^3$ |

resolved coordinate $s$에서는

$$
A_H=\oint F\,ds
$$

이고 macroscopic stress-strain에서는

$$
A_H=\oint\sigma\,d\epsilon
$$

이다.

## 9. generalized Lennard-Jones 변수

baseline에서

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
$$

를 사용하면 다음 변수들이 등장한다.

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $\varepsilon_{\rm LJ}$ | energy scale | pair-potential energy scale | J 또는 eV |
| $\sigma_{\rm LJ}$ | length scale | LJ characteristic distance parameter | m |
| $m,n$ | $m>n>1$ | generalized LJ interaction shape parameter | 무차원 |
| $\zeta$ | Riemann zeta function | semi-infinite uniform chain의 pair sum에서 등장 | 무차원 |
| $A_0$ | reference effective cross-sectional area | 1D generalized force를 engineering stress로 변환 | m$^2$ |
| $E$ | Young's modulus | small-strain elastic modulus | Pa |
| $\sigma_u$ | old calibration의 ultimate/instability stress | 과거 1D LJ calibration target | Pa |
| $a_u$ | $U''(a_u)=0$이 되는 spacing | 과거 모델의 instability location | m |

## 10. Rubin-chain 히스테리시스 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $Q(t)$ | resolved structural coordinate | harmonic lattice bath와 결합된 관심 구조좌표 | m 또는 무차원 |
| $P_Q$ | $Q$의 conjugate momentum | resolved-coordinate momentum | kg m/s |
| $M$ | resolved coordinate mass | $Q$ 또는 $s$의 effective inertia | kg 또는 무차원 |
| $x_n$ | harmonic-bath site displacement | unresolved lattice coordinate | m 또는 무차원 |
| $p_n$ | bath-site momentum | bath momentum | kg m/s |
| $m$ | bath-site mass | atomic/effective lattice mass | kg 또는 무차원 |
| $k$ | nearest-neighbor spring constant | harmonic bath coupling stiffness | N/m 또는 무차원 |
| $K_0$ | $Q$의 on-site stiffness | resolved coordinate restoring stiffness | N/m 또는 무차원 |
| $q$ | lattice wave number | site당 phase increment | rad/site |
| $\omega_D$ | $2\sqrt{k/m}$ | monatomic harmonic-chain upper band edge | rad/s 또는 무차원 |
| $Z(\omega)$ | $\hat F/\hat Q$ | resolved coordinate가 보는 complex dynamic stiffness | N/m |
| $\hat Q$ | $Q$의 complex harmonic amplitude | frequency-domain response amplitude | m |
| $\hat F$ | $F$의 complex harmonic amplitude | frequency-domain force amplitude | N |
| $\phi$ | phase lag | forcing과 resolved displacement 사이 위상차 | rad |

propagating band 내부에서 $\operatorname{Im}Z>0$은 전체 미시모델에 점성항을 넣지 않아도 unresolved lattice mode로 energy가 방사될 수 있음을 나타낸다.

## 11. Hamiltonian slip-bath 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $u_j$ | $s$와 연결된 harmonic bath coordinate | unresolved lattice displacement | m 또는 무차원 |
| $p_j$ | $u_j$의 conjugate momentum | bath momentum | kg m/s 또는 무차원 |
| $k_c$ | $s$와 $u_1$ 사이 coupling stiffness | resolved slip-bath mechanical coupling | N/m 또는 무차원 |
| $P_s$ | $s$의 conjugate momentum | resolved slip-coordinate momentum | kg m/s 또는 무차원 |
| $s_N$ | $s(NT)$ | cycle $N$ 종료 시의 slip coordinate state | $s$와 동일 |
| $\Delta s_{\rm cycle}$ | $s_{N+1}-s_N$ | cycle-to-cycle structural drift | $s$와 동일 |

완전한 periodic internal-friction 상태는

$$
s_{N+1}=s_N
$$

이고 secular structural evolution은

$$
s_{N+1}\neq s_N
$$

이다.

## 12. crack-initiation 및 first-passage 변수

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $\tau_c$ | $\inf\{t>0:\text{state enters unstable set}\}$ | crack-initiation criterion으로의 first-passage time | s |
| $F_{\rm ci}(t)$ | $\Pr(\tau_c\le t)$ | cumulative crack-initiation probability | 무차원 |
| $S(t)$ | $1-F_{\rm ci}(t)$ | survival probability | 무차원 |
| $h(t)$ | $-\dot S/S$ | hazard rate | s$^{-1}$ |
| $j_c(t)$ | absorbing instability boundary를 통과하는 probability flux | instantaneous first-passage flux | s$^{-1}$ |
| $N_i$ | crack-initiation criterion을 처음 만족하는 cycle | 예측 crack-initiation life | cycle |

이 변수들은 first-passage framework의 **정의**이며, 아직 Al에 calibration된 crack-initiation law가 아니다.

## 13. moment와 통계 descriptor

| 기호 | 정의 | 물리적 의미 | 단위 |
|---|---|---|---|
| $\langle g(a)\rangle$ | $\int g(a)P(a,t)\,da$ | population average | $g$에 따라 결정 |
| $\mu_2$ | $\langle(a-\bar a)^2\rangle$ | spacing variance | m$^2$ |
| $\mu_3$ | $\langle(a-\bar a)^3\rangle$ | third central moment / asymmetry indicator | m$^3$ |
| $\operatorname{Var}(a)$ | $\mu_2$ | spacing-distribution width | m$^2$ |
| $\operatorname{Cov}(a,v)$ | $\langle(a-\bar a)(v-\bar v)\rangle$ | $P$가 deterministic하게 넓어지거나 좁아지는 정도를 결정 | m$^2$/s |

정확한 identity는

$$
\frac{d}{dt}\operatorname{Var}(a)
=2\operatorname{Cov}(a,v)
$$

이다.

## 14. 아직 의도적으로 확정하지 않은 변수/정의

다음 항목은 이후 미시역학에서 정의 또는 calibration되어야 하며 임의값을 넣으면 안 된다.

- Representative Area의 실제 기하와 $A_{\rm RA}$;
- FCC Al에서 $s^{\alpha}$의 최종 3D 원자수준 구현;
- 최종 모델에 사용할 Al-specific $\gamma_{\rm Al}(\mathbf s)$;
- crack initiation을 위한 정확한 local instability set;
- $P(a,s,t)$를 닫기 위해 필요한 최소 correlation hierarchy;
- 특정 RA에서 reduced generalized force $F$와 거시응력을 연결하는 정확한 mapping;
- 일반적인 피로응력을 atomistic instability scale까지 연결하는 surface, defect, thermal, multi-slip amplification mechanism.

앞으로 새로운 기호를 도입하는 문서는 같은 commit에서 이 파일에도 해당 기호를 추가하는 것을 원칙으로 한다.
