# Assumptions and Approximations

This file must be updated whenever a new modeling assumption is introduced.

## Active mainline assumptions

1. The research target is high-purity / single-crystal aluminum under primarily **normal cyclic loading**.
2. The primary reduced structural coordinate is the local normal interatomic-spacing descriptor $a$.
3. The central state density is $P(a,t)$, interpreted as a thermodynamic-limit population density rather than a finite histogram.
4. The active microscopic energy baseline is a fixed generalized Lennard-Jones pair interaction when explicitly stated.
5. Fatigue evolution must arise from atomic configuration, spacing distributions, correlations, memory, or stability loss. The LJ parameters are not allowed to evolve with cycle count merely to create damage.
6. Shear/slip coordinates are not part of the active mainline. Earlier shear-oriented work is preserved under `libraries/shear/` as an auxiliary research library.

## Controlled approximations that may be tested

- One-dimensional normal chain as a reduction of the full 3D crystal.
- Independent adjacent spacings: $P_k\approx P^{*k}$. This is not exact and must be checked against correlated simulations.
- Markov closure in spacing space, only if projected memory is negligible on the target time scale.
- Fokker–Planck truncation, only after a small-jump/Kramers–Moyal argument.
- Moment closure only as a late-stage reduction.

## Forbidden shortcuts unless explicitly justified

- Fitting a Weibull, Gaussian, or other named distribution to $P(a,t)$ merely for convenience.
- Inserting an empirical hysteresis or fatigue-damage law and claiming it was derived from mechanics.
- Changing LJ parameters with cycle count to imitate damage.
- Introducing damping, barriers, transition rates, kernels, or thresholds solely to obtain a desired fatigue curve.
- Calling a high-frequency lattice resonance result a 20 Hz fatigue prediction without an explicit time-scale bridge.
- Confusing instantaneous unstable-tail occupancy with first-passage crack initiation.

---

# 한국어 번역 — 가정과 근사

새로운 모델링 가정이 도입될 때마다 이 파일을 반드시 갱신한다.

## 활성 mainline 가정

1. 연구대상은 주로 **수직 반복하중**을 받는 고순도 또는 단결정 알루미늄이다.
2. 가장 기본적인 축약 구조좌표는 국부 수직 원자간격 descriptor $a$이다.
3. 중심 상태밀도는 $P(a,t)$이며, 유한 histogram이 아니라 열역학적 극한의 population density로 해석한다.
4. 명시적으로 선언하는 경우 활성 미시에너지 baseline은 고정된 generalized Lennard-Jones pair interaction이다.
5. 피로진화는 원자배열, spacing distribution, correlation, memory 또는 stability loss에서 나와야 한다. damage를 만들기 위해 cycle에 따라 LJ parameter를 변화시키는 것은 허용하지 않는다.
6. shear/slip 좌표는 활성 mainline에 포함하지 않는다. 기존 전단 지향 연구는 `libraries/shear/` 아래 보조 연구 라이브러리로 보존한다.

## 검증 가능한 controlled approximation

- full 3D crystal을 축약한 1차원 normal chain.
- 인접 원자간격 독립성 $P_k\approx P^{*k}$. 이는 정확하지 않으며 correlated simulation과 비교해야 한다.
- 목표 시간척도에서 projected memory가 무시 가능할 때만 사용하는 spacing-space Markov closure.
- 작은 jump에 대한 Kramers–Moyal 논증 이후에만 사용하는 Fokker–Planck truncation.
- 마지막 단계의 축약으로만 사용하는 moment closure.

## 명시적 정당화 없이 금지되는 지름길

- 편의를 위해 $P(a,t)$에 Weibull, Gaussian 등 특정 분포를 fitting하는 것.
- 경험적 hysteresis 또는 fatigue-damage law를 넣고 역학에서 유도했다고 주장하는 것.
- damage를 흉내내기 위해 cycle 수에 따라 LJ parameter를 변경하는 것.
- 원하는 피로곡선을 얻기 위해 damping, barrier, transition rate, kernel, threshold를 임의로 넣는 것.
- 명시적인 시간척도 연결 없이 고주파 lattice resonance 결과를 20 Hz 피로예측이라고 부르는 것.
- 순간적인 unstable-tail 점유와 first-passage crack initiation을 혼동하는 것.
