# Milestone 2 — Normal-deformation-driven structural evolution

## Status

This document defines the **main research direction** of the project.

The project is focused primarily on cyclic **normal stress, normal bond stretching, and normal opening instability**, not on shear-slip evolution as the principal fatigue mechanism.

The reason for selecting high-purity / single-crystal Al in this project is to study a material/system where the target failure route can be isolated and analyzed through normal deformation as cleanly as possible. Whether a given Al orientation or specimen is quantitatively weaker in normal opening than in shear must be checked independently; it is not assumed as a universal material law.

## 1. Primary loading variable

The main imposed loading is a cyclic normal stress

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t).
$$

The corresponding macroscopic normal strain is

$$
\epsilon_n(t).
$$

The primary hysteresis quantity is therefore

$$
\boxed{
A_H=\oint \sigma_n\,d\epsilon_n.
}
$$

Shear stress $\tau$ is not the principal loading variable of the mainline theory.

## 2. Primary microscopic coordinate

The central microscopic coordinate remains the local normal interatomic-spacing descriptor

$$
a_i(t).
$$

The thermodynamic-limit state density is

$$
\boxed{
P(a,t)=\lim_{N\to\infty}\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right).
}
$$

The main theoretical target is

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{normal hysteresis}
\rightarrow
\text{cycle-to-cycle distribution evolution}
\rightarrow
\text{normal-opening instability}.
}
$$

## 3. Energy model hierarchy

The preferred microscopic starting point is a fixed interatomic potential, with the generalized Lennard-Jones pair law used as the principal analytic baseline:

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

For an atomistic configuration $\{\mathbf r_i\}$, the pair-energy baseline is

$$
\boxed{
U_{\rm pair}=\frac12\sum_{i\ne j}v(|\mathbf r_i-\mathbf r_j|).
}
$$

For the reduced pair-distance hierarchy,

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
}
$$

The potential parameters are not allowed to evolve with fatigue damage. Structural evolution must appear through the microscopic configuration and the distributions $P$, $P_k$, or their necessary correlation hierarchy.

## 4. Exact spacing dynamics and closure problem

For a deterministic spacing ensemble,

$$
\boxed{
\partial_tP+\partial_a(Pv_a)=0,
}
$$

with

$$
v_a(a,t)=\langle\dot a_i\mid a_i=a\rangle.
$$

This continuity equation is exact as a kinematic identity, but it is not generally closed.

For a nearest-neighbor chain,

$$
m\ddot a_i
=v'(a_{i+1})-2v'(a_i)+v'(a_{i-1}),
$$

so neighboring-spacing correlations and/or a phase-space lift may be required.

The main closure question is therefore not "which damage law should be fitted?" but

$$
\boxed{
\text{what is the minimum microscopic state required to determine }v_a?
}
$$

## 5. Main Milestone-2 target

A purely affine reversible solution has

$$
P_{N+1}(a)=P_N(a)
$$

at identical cycle phase.

Fatigue accumulation requires a mechanically generated failure of complete recovery:

$$
\boxed{
P_{N+1}(a)\neq P_N(a).
}
$$

The next model must obtain this under cyclic normal loading without inserting an empirical damage variable, arbitrary transition kernel, fitted damping law, or prescribed probability family.

Candidate mechanisms are to be tested only if they can be defined from microscopic mechanics, for example:

- anharmonic energy transfer among normal lattice modes;
- spatially nonuniform normal spacing and its correlations;
- surface-normal relaxation and free-surface opening;
- finite-temperature phase-space distributions;
- localized normal instability produced by geometry or existing microscopic heterogeneity.

## 6. Crack initiation target

A local idealized normal-instability criterion under a one-coordinate lattice baseline can be written as

$$
U''(a_c)=0.
$$

This is only an idealized reduced-lattice stability condition.

The more general objective is to formulate crack initiation as either

1. a first-passage event into a mechanically unstable normal-opening state, or
2. a loss of stability of the relevant distribution/atomistic energy landscape.

The instantaneous tail

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

must not automatically be equated with cumulative crack-initiation probability.

## 7. Role of the shear/slip model

The existing Hamiltonian slip-bath calculation is retained as an **auxiliary proof-of-principle only**.

It established that conservative microscopic dynamics plus a nonlinear structural coordinate can produce

$$
s_{N+1}\neq s_N
$$

without an empirical fatigue law.

It does **not** define the main physical mechanism of this project.

The mainline model must now reproduce the analogous result directly in the normal-spacing sector:

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

under cyclic normal stress.

---

# 한국어 번역 — 수직변형 기반 구조진화

## 상태

이 문서는 프로젝트의 **주 연구방향**을 정의한다.

이 프로젝트는 전단 slip 진화를 주 피로메커니즘으로 두는 것이 아니라, 반복적인 **수직응력, 원자결합의 수직 신장, 수직 opening instability**에 집중한다.

고순도/단결정 Al을 선택한 연구 의도는 가능한 한 수직변형 경로를 깨끗하게 분리하여 분석하는 것이다. 다만 특정 Al 결정방향이나 시편이 정량적으로 전단보다 수직 opening에 더 약하다는 사실은 별도로 검증해야 하며, 이를 보편적 재료법칙으로 가정하지 않는다.

## 1. 주 하중변수

주 외력은 반복 수직응력이다.

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t).
$$

이에 대응하는 거시적 수직변형률은

$$
\epsilon_n(t)
$$

이다.

따라서 주 히스테리시스 양은

$$
\boxed{
A_H=\oint \sigma_n\,d\epsilon_n
}
$$

이다.

전단응력 $\tau$는 메인 이론의 주 하중변수가 아니다.

## 2. 주 미시좌표

중심 미시좌표는 계속 국부적인 수직 원자간격 descriptor

$$
a_i(t)
$$

이다.

열역학적 극한의 상태밀도는

$$
\boxed{
P(a,t)=\lim_{N\to\infty}\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right)
}
$$

이다.

메인 이론의 목표 흐름은

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{수직 히스테리시스}
\rightarrow
\text{cycle별 분포진화}
\rightarrow
\text{수직 opening instability}
}
$$

이다.

## 3. 에너지모델 계층

미시역학의 출발점은 고정된 원자간 potential로 두며, generalized Lennard-Jones pair law를 주 해석 baseline으로 사용한다.

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

원자배열 $\{\mathbf r_i\}$에 대한 pair-energy baseline은

$$
\boxed{
U_{\rm pair}=\frac12\sum_{i\ne j}v(|\mathbf r_i-\mathbf r_j|)
}
$$

이다.

축약된 pair-distance hierarchy에서는

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr
}
$$

이다.

피로손상에 따라 potential parameter 자체를 변화시키지 않는다. 구조진화는 원자배열과 $P$, $P_k$, 또는 필요한 correlation hierarchy의 변화로 나타나야 한다.

## 4. 정확한 spacing dynamics와 closure 문제

결정론적 spacing ensemble에 대해

$$
\boxed{
\partial_tP+\partial_a(Pv_a)=0
}
$$

이며,

$$
v_a(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

이 continuity equation은 운동학적 항등식으로 정확하지만 일반적으로 닫혀 있지 않다.

최근접 사슬에서는

$$
m\ddot a_i
=v'(a_{i+1})-2v'(a_i)+v'(a_{i-1})
$$

이므로 인접 spacing correlation 또는 phase-space 확장이 필요할 수 있다.

따라서 핵심 closure 문제는 "어떤 damage law를 fitting할 것인가"가 아니라

$$
\boxed{
\text{$v_a$를 결정하기 위해 필요한 최소 미시상태는 무엇인가?}
}
$$

이다.

## 5. Milestone 2의 주 목표

순수 affine 가역응답에서는 동일한 cycle phase에서

$$
P_{N+1}(a)=P_N(a)
$$

이다.

피로누적을 설명하려면 역학적으로 완전복원이 깨져

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

가 되어야 한다.

다음 모델은 반복 수직하중 아래에서 empirical damage variable, arbitrary transition kernel, fitted damping law, prescribed probability family 없이 이 결과를 만들어야 한다.

검토할 수 있는 후보 메커니즘도 미시역학에서 정의 가능한 경우에만 허용한다. 예시는 다음과 같다.

- 수직 lattice mode 사이의 anharmonic energy transfer;
- 공간적으로 불균일한 normal spacing과 그 correlation;
- 자유표면의 normal relaxation 및 opening;
- 유한온도 phase-space distribution;
- 형상 또는 기존 미시적 이질성에서 발생하는 국부 수직 instability.

## 6. 균열개시 목표

단일좌표 격자 baseline에서 이상적인 국부 수직 instability는

$$
U''(a_c)=0
$$

과 같이 둘 수 있다.

다만 이것은 이상화된 축약격자 안정성 기준이다.

더 일반적인 목표는 균열개시를

1. 기계적으로 불안정한 수직 opening 상태로의 first-passage 사건 또는
2. 관련 분포/원자론적 energy landscape의 안정성 상실

로 정식화하는 것이다.

순간 tail

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

를 누적 crack-initiation probability와 자동으로 동일시하면 안 된다.

## 7. 전단/slip 모델의 역할

기존 Hamiltonian slip-bath 계산은 **보조적인 원리증명**으로만 유지한다.

그 모델은 경험적 피로법칙 없이 보존적인 미시역학과 비선형 구조좌표만으로

$$
s_{N+1}\neq s_N
$$

이 가능함을 보였다.

하지만 이것은 이 프로젝트의 주 물리메커니즘을 정의하지 않는다.

이제 메인 이론은 동일한 종류의 cycle-to-cycle 비복원성을 수직 spacing sector에서 직접 만들어야 한다.

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

즉 반복 수직응력에서 바로 이 결과를 유도하는 것이 다음 단계다.
