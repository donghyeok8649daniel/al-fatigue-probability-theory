# Assumptions and Approximations

This file must be updated whenever a new modeling assumption is introduced.

## Current working assumptions

1. The material of interest is high-purity / single-crystal aluminum.
2. The primary reduced structural coordinate is an interatomic-spacing descriptor $a$.
3. The state density $P(a,t)$ is interpreted as a thermodynamic-limit population density, not merely a finite histogram.
4. A generalized Lennard-Jones-type interaction may be used as a baseline microscopic pair potential when explicitly stated.
5. Macroscopic affine stretch may be separated from internal/non-affine structural evolution by writing $a=\lambda x$.

## Controlled approximations that may be tested

- Independent adjacent spacings: $P_k\approx P^{*k}$. This is not exact and must be validated against correlated simulations.
- Markov closure in spacing space. Allowed only if memory is shown to be negligible on the scale of interest.
- Fokker–Planck truncation. Allowed only after a small-jump/Kramers–Moyal argument.
- Moment closure. Last-resort reduced model, not a starting axiom.

## Forbidden shortcuts unless explicitly justified

- Fitting a Weibull distribution to $P(a,t)$ merely because fatigue data often look Weibull-like.
- Inserting an empirical hysteresis loop law and then claiming it was derived from mechanics.
- Introducing damping, damage variables, transition rates, barriers, or kernels only to obtain the desired fatigue curve.
- Treating a single reversible LJ coordinate as a complete fatigue model.
- Confusing instantaneous unstable-tail occupancy with first-passage crack initiation.

---

# 한국어 번역 — 가정과 근사

새로운 모델링 가정이 도입될 때마다 이 파일을 반드시 갱신한다.

## 현재의 작업 가정

1. 연구 대상 재료는 고순도 또는 단결정 알루미늄이다.
2. 가장 기본적인 축약 구조좌표는 원자간격을 나타내는 변수 $a$이다.
3. 상태밀도 $P(a,t)$는 단순한 유한 표본 히스토그램이 아니라 열역학적 극한에서의 population density로 해석한다.
4. 명시적으로 선언하는 경우 generalized Lennard-Jones 형태의 상호작용을 기준 미시 pair potential로 사용할 수 있다.
5. 거시적인 affine stretch와 내부 비아핀 구조진화를 $a=\lambda x$와 같이 분리하여 표현할 수 있다.

## 검증 가능한 controlled approximation

- 인접 원자간격의 독립성: $P_k\approx P^{*k}$. 이는 정확식이 아니며 상관관계를 포함한 simulation과 비교해 검증해야 한다.
- spacing space에서의 Markov closure. 관심 시간척도에서 memory가 무시 가능하다는 것이 확인된 경우에만 허용한다.
- Fokker–Planck truncation. 작은 jump에 대한 Kramers–Moyal 논증을 거친 뒤에만 허용한다.
- moment closure. 출발 공리가 아니라 마지막 단계의 축약모델로만 사용한다.

## 명시적인 정당화 없이는 금지되는 지름길

- 피로 데이터가 Weibull 형태로 보인다는 이유만으로 $P(a,t)$에 Weibull 분포를 fitting하는 것.
- 경험적인 hysteresis loop 법칙을 넣은 뒤 역학에서 유도했다고 주장하는 것.
- 원하는 피로곡선을 얻기 위해서만 damping, damage variable, transition rate, barrier, kernel 등을 임의로 도입하는 것.
- 하나의 가역적인 LJ 좌표를 완전한 피로모델로 취급하는 것.
- 순간적인 unstable-tail 점유율을 누적 crack-initiation probability와 혼동하는 것.
