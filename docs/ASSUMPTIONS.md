# Assumptions and Approximations

This file must be updated whenever a new modeling assumption is introduced.

## Research scope

The main physical focus is cyclic **normal stress and normal interatomic deformation**. The primary microscopic structural variable is the local spacing $a_i(t)$ and the primary state density is $P(a,t)$.

Shear-slip variables are auxiliary unless future microscopic derivation shows that they are indispensable for closure of the normal-deformation problem.

The project motivation is to study a high-purity / single-crystal Al system in which normal-opening mechanics can be isolated and analyzed as cleanly as possible. The statement that a particular Al orientation/specimen is quantitatively weaker in normal opening than in shear is **not assumed as a universal material law** and must be checked independently.

## Current working assumptions

1. The material of interest is high-purity / single-crystal aluminum.
2. The primary reduced structural coordinate is an interatomic-spacing descriptor $a$ associated with normal separation/stretching.
3. The primary imposed loading variable is cyclic normal stress $\sigma_n(t)$.
4. The state density $P(a,t)$ is interpreted as a thermodynamic-limit population density, not merely a finite histogram.
5. A generalized Lennard-Jones-type pair interaction is the principal analytic microscopic baseline when explicitly stated.
6. The microscopic potential parameters are fixed; fatigue evolution must arise through configuration/distribution evolution rather than by changing LJ parameters with damage.
7. Macroscopic affine stretch may be separated from internal structural evolution by writing $a=\lambda x$.
8. Existing Rubin-chain and Hamiltonian slip-bath models are proof-of-principle auxiliary models; they are not the main physical model of Al fatigue in this project.

## Controlled approximations that may be tested

- Independent adjacent spacings: $P_k\approx P^{*k}$. This is not exact and must be validated against correlated simulations.
- Finite-chain approximation to a semi-infinite lattice before reflected waves return.
- Markov closure in spacing space. Allowed only if memory is shown to be negligible on the scale of interest.
- Fokker–Planck truncation. Allowed only after a small-jump/Kramers–Moyal argument.
- Moment closure. Last-resort reduced model, not a starting axiom.
- One-dimensional normal-chain models. Allowed as derivation/null-test baselines, not automatically as quantitative 3D Al predictions.

## Forbidden shortcuts unless explicitly justified

- Fitting a Weibull distribution to $P(a,t)$ merely because fatigue data often look Weibull-like.
- Inserting an empirical hysteresis loop law and then claiming it was derived from mechanics.
- Introducing damping, damage variables, transition rates, barriers, or kernels only to obtain the desired fatigue curve.
- Treating a single reversible LJ coordinate as a complete fatigue model.
- Lowering an atomistic barrier solely so that a desired macroscopic fatigue stress causes instability.
- Making shear/dislocation evolution the central mechanism when the model has not shown that it is required for the normal-opening problem.
- Confusing instantaneous unstable-tail occupancy with first-passage crack initiation.

---

# 한국어 번역 — 가정과 근사

새로운 모델링 가정이 도입될 때마다 이 파일을 반드시 갱신한다.

## 연구 범위

주 물리적 관심은 반복적인 **수직응력과 수직 원자간 변형**이다. 중심 미시 구조변수는 국부 spacing $a_i(t)$이고, 주 상태밀도는 $P(a,t)$이다.

전단/slip 변수는 향후 미시역학적 유도에서 수직변형 문제의 closure에 반드시 필요하다고 확인되지 않는 한 보조변수로만 취급한다.

이 프로젝트에서 고순도/단결정 Al을 선택한 의도는 수직 opening mechanics를 가능한 한 깨끗하게 분리하여 분석하기 위해서다. 특정 Al 결정방향이나 시편이 정량적으로 전단보다 수직 opening에 더 약하다는 명제는 **보편적인 재료법칙으로 가정하지 않으며** 별도로 검증해야 한다.

## 현재의 작업 가정

1. 연구 대상 재료는 고순도 또는 단결정 알루미늄이다.
2. 가장 기본적인 축약 구조좌표는 수직 분리/신장을 나타내는 원자간격 descriptor $a$이다.
3. 주 외부하중 변수는 반복 수직응력 $\sigma_n(t)$이다.
4. 상태밀도 $P(a,t)$는 단순한 유한 표본 히스토그램이 아니라 열역학적 극한의 population density로 해석한다.
5. generalized Lennard-Jones pair interaction을 명시적인 주 해석 미시 baseline으로 사용한다.
6. microscopic potential parameter는 고정한다. 피로진화는 LJ parameter의 damage-dependent 변화가 아니라 원자배열/분포의 변화로 나타나야 한다.
7. 거시적인 affine stretch와 내부 구조진화를 $a=\lambda x$로 분리할 수 있다.
8. 기존 Rubin-chain 및 Hamiltonian slip-bath 모델은 원리증명용 보조모델이며, 프로젝트의 실제 Al 피로 주모델이 아니다.

## 검증 가능한 controlled approximation

- 인접 원자간격의 독립성: $P_k\approx P^{*k}$. 이는 정확식이 아니며 상관관계 simulation과 비교해 검증해야 한다.
- 반사파가 돌아오기 전까지 유한 사슬로 준무한 격자를 근사하는 것.
- spacing space에서의 Markov closure. 관심 시간척도에서 memory가 무시 가능하다는 것이 확인된 경우에만 허용한다.
- Fokker–Planck truncation. 작은 jump에 대한 Kramers–Moyal 논증 뒤에만 허용한다.
- moment closure. 출발 공리가 아니라 마지막 단계의 축약모델로만 사용한다.
- 1차원 수직사슬 모델. 유도와 null test baseline으로는 허용하지만 자동으로 정량적인 3D Al 예측으로 취급하지 않는다.

## 명시적인 정당화 없이는 금지되는 지름길

- 피로 데이터가 Weibull 형태로 보인다는 이유만으로 $P(a,t)$에 Weibull 분포를 fitting하는 것.
- 경험적인 hysteresis loop 법칙을 넣은 뒤 역학에서 유도했다고 주장하는 것.
- 원하는 피로곡선을 얻기 위해 damping, damage variable, transition rate, barrier, kernel 등을 임의로 도입하는 것.
- 하나의 가역적인 LJ 좌표를 완전한 피로모델로 취급하는 것.
- 원하는 macroscopic fatigue stress에서 instability가 발생하도록 atomistic barrier를 임의로 낮추는 것.
- 수직 opening 문제에 필수적이라는 역학적 근거 없이 shear/dislocation evolution을 중심 메커니즘으로 승격하는 것.
- 순간적인 unstable-tail 점유율을 누적 crack-initiation probability와 혼동하는 것.
