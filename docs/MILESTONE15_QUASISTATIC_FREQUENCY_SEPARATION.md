# Milestone 15 — Quasistatic Frequency Separation and a No-Go Result for Pure Conservative Fatigue Accumulation

## Status

**ACTIVE 1D THEORY CONSTRAINT / FALSIFICATION RESULT**

This milestone clarifies the role of loading frequency in the current one-dimensional normal layer-LJ research. It does not add an empirical fatigue law.

The key distinction is between:

1. **laboratory fatigue frequency** such as 10–100 Hz;
2. the microscopic layer/atomic mechanical time scale;
3. artificial or controlled finite-chain frequencies used in numerical diagnostics.

These must not be conflated.

## Dimensionless loading frequency

The microscopic reduced angular frequency is

$$
\omega^* = 2\pi f t_0,
$$

where $f$ is the laboratory loading frequency and $t_0$ is the microscopic normal-mechanics time scale.

When

$$
\boxed{\omega^*\ll 1},
$$

the microscopic mechanical coordinates can be treated as quasistatic relative to the laboratory load variation, provided no other slow internal state is introduced.

This does **not** mean that all fatigue mechanisms are frequency-independent. It means only that ordinary laboratory frequency should not be interpreted as inertial resonance of the microscopic conservative LJ coordinate unless a separate time-scale calculation supports that interpretation.

## Exact stable homogeneous branch

For the normalized generalized-LJ layer energy $\phi(\lambda)$, define the dimensionless tensile force $f^*$. Below the ideal normal instability,

$$
0\le f^*<f_c,
$$

the stable homogeneous branch satisfies

$$
\boxed{\phi'(\lambda_s)=f^*.}
$$

On the stable tensile branch $1\le\lambda_s<\lambda_c$, this relation is single-valued.

For a tension-only cyclic load written in phase form,

$$
f^*(\theta)=f_m^*+f_a^*\sin\theta,
\qquad
0\le\theta\le2\pi,
$$

the quasistatic response is therefore

$$
\lambda_s(\theta)=\lambda_s\!\left(f^*(\theta)\right).
$$

Physical frequency does not appear in this equilibrium path. Frequency only changes the map between phase and laboratory time,

$$
\theta=2\pi f t.
$$

## Exact no-hysteresis result under stated assumptions

Assume:

1. the mechanics is purely one-dimensional normal generalized-LJ;
2. the state remains on the unique stable homogeneous branch;
3. the load remains subcritical;
4. the response is quasistatic;
5. no additional irreversible internal variable, stochastic escape process, defect state, or absorbing crack state is present.

Then the response is a single-valued function of the instantaneous force.

Therefore

$$
f^*(2\pi)=f^*(0)
\quad\Rightarrow\quad
\lambda_s(2\pi)=\lambda_s(0),
$$

and loading and unloading retrace the same equilibrium curve.

Hence, under these assumptions,

$$
\boxed{\text{hysteresis area}=0,}
$$

$$
\boxed{\text{cycle-to-cycle state accumulation}=0.}
$$

At zero temperature the one-point distribution is

$$
P(\lambda\mid f^*)=\delta\!\left(\lambda-\lambda_s(f^*)\right),
$$

so the entire distribution also returns exactly after every closed subcritical cycle.

Likewise, if a finite-temperature distribution is assumed to instantaneously equilibrate to a unique $P_{\rm eq}(\lambda\mid f^*,T)$, it is again single-valued in the instantaneous load and retraces the same family of distributions. Instantaneous equilibrium alone cannot produce fatigue hysteresis.

## Consequence for the current deterministic chain simulations

Earlier dynamically matched sweeps used

$$
\omega^* M=\text{constant}
$$

to keep the loading period comparable to the finite represented-chain acoustic transit time while changing $M$.

That construction remains a useful **CONTROLLED NUMERICAL FINITE-SIZE DIAGNOSTIC**.

It must not be interpreted as evidence that real fatigue strength is controlled by microscopic inertia at the same dimensionless frequency. Finite $\omega^*$ in the conservative chain can generate wave structure, phase lag, and system-scale coherence that vanish or change in the quasistatic limit.

Therefore the previously observed approximately constant positive-window effective count

$$
\widehat M_{\rm eff}^{(+)}\approx 3
$$

must be tested for survival as $\omega^*\to0$ before it can be considered a material statistical structure.

## Relation to experimental fatigue-frequency observations

For many metallic high-cycle fatigue tests at ordinary room-temperature frequencies, cycles-to-failure can be much less sensitive to frequency than a purely time-controlled activated-damage model would predict. Aluminum studies exist in which no statistically significant fatigue-life change was observed over substantial frequency ranges, although frequency effects can appear when temperature rise, environment, creep, diffusion, strain-rate sensitivity, or other time-dependent mechanisms become important.

Therefore weak laboratory-frequency sensitivity is consistent with treating the immediate elastic layer response as quasistatic, but it does **not** by itself identify the missing irreversible fatigue mechanism.

## Research implication: the missing mechanism is not ordinary microscopic inertia

The no-go result leaves a sharp question:

> What physically justified state or process makes $P(\lambda,t)$ fail to retrace under a slow closed load cycle?

Candidate classes must be kept distinct:

### A. Stochastic first passage / metastable escape

A thermally or otherwise stochastically driven layer state may escape an intact basin under a periodically modulated barrier. This can introduce irreversibility through an absorbing or first-passage condition.

However, a conventional continuous-time escape rate generally introduces an explicit time-per-cycle dependence. Therefore its predicted frequency dependence must be checked against experimental frequency insensitivity rather than assumed acceptable.

### B. Additional irreversible microscopic state

A hidden internal coordinate can create path dependence even when the normal elastic coordinate is quasistatic. Examples in real metals include defect/dislocation structure, but introducing such a state would go beyond the present pure spacing-only model and requires a first-principles justification rather than an arbitrary damage variable.

### C. Spatial heterogeneity or boundary-induced localization

Nonuniform initial states, defects, or boundary conditions can create a nontrivial distribution and local first-passage behavior. These may be investigated in 1D without moving to 2D/3D, but they must be physically defined rather than inserted as fitted noise.

### D. Finite-frequency conservative dynamics

This route can produce phase lag and nonuniformity numerically, but it is now classified as an unlikely primary explanation for ordinary low-frequency metal fatigue unless a separate physical time-scale argument demonstrates otherwise.

## Next numerical test

The immediate 1D test should be a fixed-$M$, fixed-force-path sweep

$$
\omega^*\downarrow0
$$

with snapshots compared at identical cycle phase rather than identical laboratory time.

Track:

- spacing variance;
- skewness;
- $\rho_k$;
- $\widehat M_{\rm eff}^{(+)}$;
- deviation from the homogeneous stable branch;
- cycle-to-cycle residual state.

If all nonuniformity and apparent statistical-cell structure collapse toward the homogeneous branch as $\omega^*\to0$, the earlier coherence is a finite-rate/finite-chain artifact.

If a nontrivial distribution survives the quasistatic limit, its origin must be identified from the initial/boundary/ensemble assumptions before assigning it material meaning.

## Classification summary

- $\phi'(\lambda_s)=f^*$ on the stable branch: **EXACT within the stated 1D potential**.
- Single-valued closed-cycle retracing on that branch: **EXACT**.
- Zero hysteresis and zero accumulation for a pure quasistatic conservative state with no hidden variables: **EXACT under the stated assumptions**.
- $\omega^*\ll1$ for a specific experiment: **requires a physical time-scale calibration**.
- Weak frequency sensitivity of real fatigue: **EMPIRICAL and material/environment dependent**.
- The missing irreversible mechanism: **OPEN RESEARCH PROBLEM**.

---

# 마일스톤 15 — 준정적 주파수 분리와 순수 보존계 피로누적의 No-Go 결과

## 상태

**활성 1D 이론 제약 / 반증 결과**

이 마일스톤은 현재 1차원 normal layer-LJ 연구에서 하중 주파수의 역할을 명확히 한다. 경험적 피로 법칙을 새로 넣지 않는다.

가장 중요한 것은 다음 세 주파수/시간척도를 구분하는 것이다.

1. 10–100 Hz 같은 **실제 피로시험 주파수**;
2. 원자/층 normal mechanics의 microscopic time scale;
3. 유한 chain 수치진단에서 인위적 또는 통제적으로 사용하는 주파수.

이 셋을 같은 의미로 해석하면 안 된다.

## 무차원 하중 주파수

microscopic reduced angular frequency는

$$
\omega^* = 2\pi f t_0
$$

이다. 여기서 $f$는 실험실 하중주파수, $t_0$는 microscopic normal-mechanics 시간척도다.

만약

$$
\boxed{\omega^*\ll1}
$$

이면, 별도의 느린 내부상태가 없다는 조건 아래 microscopic mechanical coordinate는 실험 하중 변화에 비해 준정적으로 볼 수 있다.

이것은 모든 실제 피로기구가 주파수와 무관하다는 뜻이 아니다. 별도의 시간척도 계산 없이 실제 시험주파수를 microscopic conservative LJ coordinate의 관성 공진으로 해석하면 안 된다는 뜻이다.

## 정확한 homogeneous stable branch

normalized generalized-LJ layer energy를 $\phi(\lambda)$라고 하고 dimensionless tensile force를 $f^*$라 하자. 이상적인 normal instability보다 아래에서는

$$
0\le f^*<f_c
$$

이고 stable homogeneous branch는

$$
\boxed{\phi'(\lambda_s)=f^*}
$$

를 만족한다.

stable tensile branch $1\le\lambda_s<\lambda_c$에서는 이 관계가 단일값이다.

순수 인장 cyclic load를 phase로

$$
f^*(\theta)=f_m^*+f_a^*\sin\theta,
\qquad
0\le\theta\le2\pi
$$

라고 쓰면 준정적 응답은

$$
\lambda_s(\theta)=\lambda_s\!\left(f^*(\theta)\right)
$$

이다.

이 equilibrium path에는 실제 주파수 자체가 나타나지 않는다. 주파수는 phase와 실험실 시간 사이의 관계

$$
\theta=2\pi f t
$$

만 바꾼다.

## 명시한 가정 아래의 exact no-hysteresis 결과

다음을 가정한다.

1. mechanics는 순수 1D normal generalized-LJ다;
2. 상태는 unique stable homogeneous branch에 머문다;
3. 하중은 subcritical이다;
4. 응답은 준정적이다;
5. 별도의 irreversible internal variable, stochastic escape, defect state, absorbing crack state가 없다.

그러면 응답은 현재 force만의 단일값 함수다.

따라서

$$
f^*(2\pi)=f^*(0)
\quad\Rightarrow\quad
\lambda_s(2\pi)=\lambda_s(0)
$$

이고 loading과 unloading은 동일 equilibrium curve를 정확히 되짚는다.

그러므로 위 가정 아래에서는

$$
\boxed{\text{hysteresis area}=0}
$$

이고

$$
\boxed{\text{cycle-to-cycle state accumulation}=0}
$$

이다.

$T=0$에서 one-point distribution은

$$
P(\lambda\mid f^*)=\delta\!\left(\lambda-\lambda_s(f^*)\right)
$$

이므로 닫힌 subcritical cycle 뒤 전체 분포도 정확히 원래 상태로 돌아온다.

finite-temperature distribution을 매 순간 unique $P_{\rm eq}(\lambda\mid f^*,T)$로 즉시 평형화한다고 가정해도 마찬가지다. 현재 하중의 단일값 함수이므로 동일한 distribution family를 되짚는다. instantaneous equilibrium만으로는 피로 hysteresis가 생기지 않는다.

## 현재 deterministic chain simulation에 대한 의미

과거 dynamically matched sweep에서는 $M$을 바꿀 때 유한 represented chain의 acoustic transit time과 loading period의 비를 대략 유지하기 위해

$$
\omega^* M=\text{constant}
$$

을 사용했다.

이 방법은 여전히 **통제된 수치적 finite-size diagnostic**으로 유효하다.

그러나 이를 실제 fatigue strength가 동일한 dimensionless microscopic inertia에 의해 결정된다는 근거로 해석하면 안 된다. 보존적 finite-$\omega^*$ chain에서는 wave structure, phase lag, system-scale coherence가 생길 수 있고, 이는 준정적 limit에서 사라지거나 변할 수 있다.

따라서 과거에 관측한

$$
\widehat M_{\rm eff}^{(+)}\approx3
$$

이라는 결과도 $\omega^*\to0$에서 살아남는지 확인하기 전에는 물질의 statistical structure라고 볼 수 없다.

## 실제 피로 주파수 관측과의 관계

많은 금속의 상온 high-cycle fatigue에서는 일반적인 시험주파수 범위에서 cycles-to-failure가 순수한 시간제어 activated-damage 모델이 예측할 만큼 강한 주파수 의존성을 보이지 않는 경우가 있다. 알루미늄에서도 상당한 주파수 범위에 걸쳐 통계적으로 유의한 fatigue-life 변화가 관측되지 않은 연구들이 있다. 반면 온도상승, 환경, creep, diffusion, strain-rate sensitivity 같은 시간의존 메커니즘이 중요해지면 주파수 효과가 나타날 수 있다.

따라서 약한 실험 주파수 민감도는 즉각적인 elastic layer response를 준정적으로 보는 것과 양립한다. 그러나 이것만으로 빠져 있는 irreversible fatigue mechanism이 무엇인지 정해지는 것은 아니다.

## 연구적 의미: 빠져 있는 메커니즘은 일반적인 microscopic inertia가 아니다

no-go 결과가 남기는 질문은 명확하다.

> 느린 닫힌 하중 cycle에서 무엇이 $P(\lambda,t)$를 원래 경로로 돌아가지 못하게 만드는가?

가능한 후보들은 반드시 구분해야 한다.

### A. stochastic first passage / metastable escape

thermal 또는 다른 stochastic fluctuation으로 intact basin에서 빠져나가는 과정은 barrier가 주기적으로 바뀌는 상황에서 absorbing/first-passage irreversibility를 만들 수 있다.

그러나 일반적인 continuous-time escape rate는 cycle당 실제 소요시간에 직접 의존한다. 따라서 이런 모델이 예측하는 주파수 의존성은 실제 실험의 약한 주파수 민감도와 반드시 비교해야 한다.

### B. 추가적인 irreversible microscopic state

normal elastic coordinate가 준정적이어도 hidden internal coordinate가 있으면 path dependence가 생길 수 있다. 실제 금속에서는 defect/dislocation structure 같은 것이 예시지만, 이를 넣으면 현재 pure spacing-only model을 넘어가므로 임의 damage variable이 아니라 first-principles 근거가 필요하다.

### C. spatial heterogeneity 또는 boundary-induced localization

비균일 초기상태, defect, boundary condition은 nontrivial distribution과 local first-passage behavior를 만들 수 있다. 이는 2D/3D로 가지 않고 1D에서 연구할 수 있지만 fitted noise로 넣는 것이 아니라 물리적으로 정의해야 한다.

### D. finite-frequency conservative dynamics

수치적으로 phase lag와 nonuniformity를 만들 수는 있지만, 별도의 물리 시간척도 근거가 없는 한 일반적인 저주파 금속 피로의 주된 원인으로 보는 것은 현재 우선순위에서 내린다.

## 다음 수치시험

즉시 할 1D 시험은 fixed-$M$, fixed-force-path 조건에서

$$
\omega^*\downarrow0
$$

sweep을 하는 것이다.

동일한 실험실 시간이 아니라 동일 cycle phase에서 snapshot을 비교해야 한다.

추적할 값은 다음과 같다.

- spacing variance;
- skewness;
- $\rho_k$;
- $\widehat M_{\rm eff}^{(+)}$;
- homogeneous stable branch로부터의 deviation;
- cycle-to-cycle residual state.

$\omega^*\to0$에서 모든 nonuniformity와 apparent statistical-cell structure가 homogeneous branch로 무너지면 과거 coherence는 finite-rate/finite-chain artifact다.

반대로 quasistatic limit에서도 nontrivial distribution이 남는다면, 물질적 의미를 부여하기 전에 그것이 어떤 initial/boundary/ensemble assumption에서 나온 것인지 밝혀야 한다.

## 분류 요약

- stable branch에서 $\phi'(\lambda_s)=f^*$: **명시한 1D potential 안에서 EXACT**.
- 그 branch에서 닫힌 cycle의 single-valued retracing: **EXACT**.
- hidden variable 없는 순수 준정적 보존상태의 zero hysteresis / zero accumulation: **명시한 가정 아래 EXACT**.
- 특정 실험에서 $\omega^*\ll1$: **physical time-scale calibration 필요**.
- 실제 피로의 약한 frequency sensitivity: **재료/환경 의존 EMPIRICAL**.
- 빠져 있는 irreversible mechanism: **OPEN RESEARCH PROBLEM**.
