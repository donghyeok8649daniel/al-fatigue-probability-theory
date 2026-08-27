# Assumptions and Approximations

This file records the active assumptions of the **normal-deformation / normal-opening** mainline only.

## Active assumptions

1. The research target is high-purity / single-crystal aluminum under primarily **normal cyclic loading**.
2. The primary reduced structural coordinate is the local normal interatomic spacing $a$.
3. The central state density is $P(a,t)$, interpreted as a thermodynamic-limit population density rather than a finite histogram.
4. The principal analytic microscopic baseline is a fixed generalized Lennard-Jones pair interaction when explicitly stated.
5. Fatigue evolution must arise from atomic configuration, spacing distributions, correlations, projected memory, or mechanical stability loss.
6. Lennard-Jones parameters are not allowed to evolve with cycle count merely to imitate damage.
7. The current 1D normal chain is a reduced model and is not claimed to be an exact 3D description of FCC aluminum.

## Controlled approximations that may be tested

- One-dimensional normal chain as a reduction of the full 3D crystal.
- Nearest-neighbor truncation when explicitly used.
- Independent adjacent spacings, $P_k\approx P^{*k}$, only after comparison with correlated simulations.
- Markov closure in spacing space only if projected memory is negligible on the target time scale.
- Fokker–Planck truncation only after a small-jump/Kramers–Moyal argument.
- Moment closure only as a late-stage reduction.

## Forbidden shortcuts unless explicitly justified

- Fitting a Weibull, Gaussian, or other named family to $P(a,t)$ merely for convenience.
- Inserting an empirical hysteresis law or fatigue-damage evolution equation and then calling it mechanics-derived.
- Changing LJ parameters with cycle count to create degradation.
- Introducing damping, relaxation times, transition rates, kernels, or thresholds solely to obtain a desired fatigue curve.
- Calling an atomic-frequency dynamic instability a 20 Hz fatigue prediction without a physically derived time-scale bridge.
- Confusing instantaneous unstable-tail occupancy with cumulative first-passage crack initiation.
- Tuning away a reversible 100 MPa null result merely because fatigue accumulation was expected.

---

# 한국어 번역 — 가정과 근사

이 문서는 활성 **수직변형 / normal-opening** mainline의 가정만 기록한다.

## 활성 가정

1. 연구대상은 주로 **수직 반복하중**을 받는 고순도 또는 단결정 알루미늄이다.
2. 가장 기본적인 축약 구조좌표는 국부 수직 원자간격 $a$이다.
3. 중심 상태밀도는 $P(a,t)$이며, 유한 histogram이 아니라 열역학적 극한의 population density로 해석한다.
4. 명시적으로 선언하는 경우 주된 해석적 미시 baseline은 고정 generalized Lennard-Jones pair interaction이다.
5. 피로진화는 원자배열, spacing distribution, correlation, projected memory 또는 mechanical stability loss에서 나와야 한다.
6. damage를 흉내내기 위해 cycle에 따라 Lennard-Jones parameter를 변화시키지 않는다.
7. 현재 1D normal chain은 축약모델이며 FCC aluminum의 정확한 3D 표현이라고 주장하지 않는다.

## 검증 가능한 controlled approximation

- full 3D crystal을 축약한 1차원 normal chain.
- 명시적으로 사용하는 nearest-neighbor truncation.
- correlated simulation과 비교한 뒤에만 허용하는 adjacent-spacing independence $P_k\approx P^{*k}$.
- 목표 시간척도에서 projected memory가 무시 가능할 때만 사용하는 spacing-space Markov closure.
- 작은 jump에 대한 Kramers–Moyal 논증 이후에만 사용하는 Fokker–Planck truncation.
- 마지막 단계의 축약으로만 사용하는 moment closure.

## 명시적인 정당화 없이 금지되는 지름길

- 편의를 위해 $P(a,t)$에 Weibull, Gaussian 등 특정 분포 family를 fitting하는 것.
- 경험적 hysteresis law 또는 fatigue-damage evolution equation을 넣고 mechanics-derived라고 부르는 것.
- degradation을 만들기 위해 cycle 수에 따라 LJ parameter를 변경하는 것.
- 원하는 피로곡선을 만들 목적으로 damping, relaxation time, transition rate, kernel, threshold를 임의로 넣는 것.
- 물리적으로 유도된 time-scale bridge 없이 atomic-frequency instability를 20 Hz 피로예측이라고 부르는 것.
- 순간 unstable-tail occupancy와 누적 first-passage crack initiation을 혼동하는 것.
- 피로누적을 기대했다는 이유만으로 가역적인 100 MPa null result를 tuning으로 없애는 것.
