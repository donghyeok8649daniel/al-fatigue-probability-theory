# Assumptions and Approximations

## Statistical-mechanical ensemble rule

<!-- PHYSICAL_P_ASSUMPTIONS_EN -->

A smooth thermal $P$ is not assumed merely because the coordinate is treated probabilistically. The following physical distinctions are mandatory:

- the zero-temperature homogeneous quasistatic state is a delta distribution on the stable mechanical branch;
- a microcanonical distribution requires an isolated conservative equilibrium interpretation and includes the kinetic density-of-states factor after momenta are integrated out;
- a canonical fixed-length distribution requires a thermal reservoir and fixed total length;
- a tensile intact-basin Gibbs density is only a metastable/local-equilibrium approximation and requires intrabasin equilibration to be fast relative to loading/escape;
- no Kramers or Arrhenius escape prefactor is introduced without an independently derived bath/friction/phonon time scale;
- $A_0$ is a physical coarse-graining input and may not be tuned to obtain a desired tail probability or fatigue life.

This file records the assumptions of the active **one-dimensional normal-LJ / continuous-time** mainline.

## Active assumptions

1. The research target is high-purity / single-crystal aluminum under primarily normal cyclic loading.
2. The active theory is intentionally one-dimensional. Three-dimensional FCC calculations are archived and are not part of the default mainline.
3. The primary microscopic coordinate is the local normal interatomic spacing $a$ or its normalized form $\lambda=a/a_0$.
4. The central state density is $P(a,t)$ or equivalently $P(\lambda,t)$.
5. Physical time $t$ is the fundamental evolution coordinate. Fatigue cycle count is not an independent state variable.
6. The microscopic energy baseline is a fixed generalized Lennard-Jones interaction when explicitly stated.
7. LJ parameters do not evolve with loading history merely to imitate damage.
8. No named probability family is assumed for $P$.
9. Crack initiation is pursued through normal opening / normal stability loss, not through an inserted empirical damage law.
10. A 2D/3D specimen mesh is permitted as a geometry, storage, and visualization structure, but the active microscopic state and probability dynamics remain 1D normal-only.
11. Until a multidimensional continuum solver is separately validated, copying a 1D axial field to 2D/3D cells is labeled as a visualization/post-processing projection rather than a multidimensional mechanics result.

## Active Smoluchowski/Floquet assumptions

1. Eliminated atomic and phonon coordinates act as an isothermal bath on the
   resolved timescale, and spacing-velocity relaxation is fast enough for the
   overdamped reduction.
2. Mobility is constant. No coordinate-dependent mobility, stochastic
   convention correction, or independently fitted diffusivity is active.
3. The tangent-instability point $\lambda_c$ is an operational absorbing
   initiation boundary. Its mathematical use is exact after declaration; its
   equivalence to observed crack initiation remains a physical validation
   problem.
4. The imposed load is periodic when the one-cycle spectrum is used. Cycle
   count is only a stroboscopic observation of the continuous-time PDE, not a
   newly introduced evolution variable.
5. The generator itself does not degrade between cycles. Therefore the
   long-cycle survivor-conditioned density and energy are periodic, while
   irreversible accumulation occurs through escaped probability.

## Exact result versus physical constraint

The continuous-time energy-feasibility theorem is exact only after the admissible support is stated.

A central additional condition is

$$
\boxed{
\lambda\ge\lambda_L(t)>0.
}
$$

The mathematical theorem treats $\lambda_L(t)$ as a given hard lower bound. Its physical value is **not yet derived** and must not be fitted to obtain a desired fatigue life.

Without a compression-side constraint, normalization, mean, and LJ energy alone cannot force a tensile tail because LJ repulsion diverges as $\lambda\to0^+$.

## Controlled approximations that may be tested

- finite 1D chain as a numerical representation of the active 1D theory;
- nearest-neighbor truncation when explicitly used;
- finite empirical density as an approximation to $P(a,t)$;
- any later Markov, Fokker-Planck, or moment closure only after a derivation justifies it.

## Forbidden shortcuts unless explicitly justified

- Fitting a Weibull, Gaussian, or another named family to $P(a,t)$ merely for convenience.
- Using fatigue cycle count as the fundamental evolution coordinate.
- Changing LJ parameters with time or loading history to create degradation.
- Introducing an empirical damage variable and calling it mechanics-derived.
- Introducing damping, relaxation times, transition rates, kernels, or thresholds solely to obtain a desired fatigue curve.
- Assuming that all external work is retained as configurational energy.
- Using $\mathcal E_{\rm safe}^{\max}(t)$ without stating the compression-side constraint.
- Treating illustrative $\lambda_L$ values as aluminum material constants.
- Calling an atomic-frequency dynamic instability a 20 Hz fatigue prediction without a derived bridge.
- Tuning away the reversible 100 MPa null result merely because fatigue accumulation was expected.

---

# 한국어 번역 — 가정과 근사

## 통계역학 ensemble 규칙

<!-- PHYSICAL_P_ASSUMPTIONS_KO -->

coordinate를 probabilistic하게 다룬다는 이유만으로 smooth thermal $P$를 가정하지 않는다. 다음 물리적 구분을 반드시 유지한다.

- zero-temperature homogeneous quasistatic state는 stable mechanical branch 위의 delta distribution이다.
- microcanonical distribution은 isolated conservative equilibrium 해석을 요구하며 momentum 적분 뒤 kinetic density-of-states factor를 포함한다.
- canonical fixed-length distribution은 thermal reservoir와 fixed total length를 요구한다.
- tensile intact-basin Gibbs density는 metastable/local-equilibrium approximation일 뿐이며 intrabasin equilibration이 loading/escape보다 빨라야 한다.
- independently derived bath/friction/phonon time scale 없이 Kramers 또는 Arrhenius escape prefactor를 넣지 않는다.
- $A_0$는 physical coarse-graining input이며 원하는 tail probability나 fatigue life를 얻도록 tuning하면 안 된다.

이 문서는 활성 **1차원 normal-LJ / 연속시간** mainline의 가정을 기록한다.

## 활성 가정

1. 연구대상은 주로 normal cyclic loading을 받는 고순도 또는 단결정 알루미늄이다.
2. 활성 이론은 의도적으로 1차원이다. 3차원 FCC 계산은 archive에 보존하며 기본 mainline에는 포함하지 않는다.
3. 주 microscopic coordinate는 국부 수직 원자간격 $a$ 또는 normalized form $\lambda=a/a_0$이다.
4. 중심 상태밀도는 $P(a,t)$ 또는 동등한 $P(\lambda,t)$이다.
5. 물리적 시간 $t$가 근본 evolution coordinate다. fatigue cycle count는 독립적인 state variable이 아니다.
6. 명시적으로 사용할 때 microscopic energy baseline은 고정 generalized Lennard-Jones interaction이다.
7. damage를 흉내내기 위해 loading history에 따라 LJ parameter를 바꾸지 않는다.
8. $P$에 특정 named probability family를 가정하지 않는다.
9. crack initiation은 empirical damage law가 아니라 normal opening / normal stability loss로 다룬다.
10. 2D/3D 시편 mesh는 geometry, storage 및 visualization 구조로 허용하지만 활성 microscopic state와 probability dynamics는 1D normal-only로 유지한다.
11. multidimensional continuum solver를 별도로 검증하기 전까지 1D axial field를 2D/3D cell에 복사한 결과는 multidimensional mechanics result가 아니라 visualization/post-processing projection으로 표시한다.

## 활성 Smoluchowski/Floquet 가정

1. 생략한 atomic/phonon 좌표는 해석 시간척도에서 isothermal bath로
   작용하고, spacing velocity relaxation은 overdamped 축약이 가능할 만큼
   빠르다.
2. mobility는 상수다. 위치의존 mobility, stochastic convention 보정,
   독립 fitting diffusion은 활성화하지 않는다.
3. 변곡점 $\lambda_c$를 operational absorbing initiation boundary로 둔다.
   선언 뒤의 수학은 정확하지만 실제 균열개시와 같은지는 별도 물리 검증
   대상이다.
4. 한 주기 spectrum을 쓸 때 하중은 주기적이다. cycle count는 연속시간
   PDE를 stroboscopic하게 관찰한 것이며 새로운 근본 evolution variable이
   아니다.
5. generator 자체는 cycle에 따라 열화되지 않는다. 따라서 장기 생존조건부
   분포와 에너지는 주기적이고, 비가역 누적은 유출확률로 발생한다.

## Exact result와 physical constraint의 구분

continuous-time energy-feasibility theorem은 admissible support가 명시된 뒤에만 exact하다.

핵심 추가조건은

$$
\boxed{
\lambda\ge\lambda_L(t)>0
}
$$

이다.

수학 theorem은 $\lambda_L(t)$를 주어진 hard lower bound로 취급한다. 그 물리적 값은 **아직 유도되지 않았으며** 원하는 fatigue life를 얻기 위한 fitting parameter로 사용하면 안 된다.

compression-side constraint가 없으면 LJ repulsion이 $\lambda\to0^+$에서 발산하므로 정규화, 평균, LJ energy만으로 tensile tail을 강제할 수 없다.

## 검증 가능한 controlled approximation

- 활성 1D theory의 numerical representation으로 사용하는 finite 1D chain;
- 명시적으로 선언한 nearest-neighbor truncation;
- $P(a,t)$의 approximation으로 사용하는 finite empirical density;
- 향후 Markov, Fokker-Planck, moment closure는 derivation이 정당화한 뒤에만 사용.

## 명시적 정당화 없이 금지되는 shortcut

- 편의를 위해 $P(a,t)$에 Weibull, Gaussian 또는 다른 named family를 fitting하는 것.
- fatigue cycle count를 근본 evolution coordinate로 사용하는 것.
- degradation을 만들기 위해 시간이나 loading history에 따라 LJ parameter를 바꾸는 것.
- empirical damage variable을 넣고 mechanics-derived라고 부르는 것.
- 원하는 fatigue curve를 얻기 위해 damping, relaxation time, transition rate, kernel, threshold를 임의로 넣는 것.
- 외부 work 전체가 configurational energy로 저장된다고 가정하는 것.
- compression-side constraint를 밝히지 않고 $\mathcal E_{\rm safe}^{\max}(t)$를 사용하는 것.
- illustrative $\lambda_L$ 값을 aluminum material constant로 취급하는 것.
- derived bridge 없이 atomic-frequency dynamic instability를 20 Hz fatigue prediction이라고 부르는 것.
- 피로누적을 기대했다는 이유로 reversible 100 MPa null result를 tuning으로 없애는 것.
