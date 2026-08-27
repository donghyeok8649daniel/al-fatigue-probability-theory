# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation in high-purity / single-crystal aluminum.

## Research goal

Derive cyclic hysteresis, irreversible fatigue accumulation, and crack initiation from microscopic mechanics and distribution-valued state variables, while minimizing empirical fatigue laws and uncontrolled fitting.

## Core state

The primary reduced state is the thermodynamic-limit spacing density

$$
P(a,t)=\lim_{N\to\infty}\frac{1}{N}\sum_{i=1}^{N}\delta(a-a_i(t)).
$$

where $a_i(t)$ is a local interatomic-spacing descriptor. The long-term goal is to derive the evolution of $P(a,t)$ from mechanics rather than prescribe a probability family.

## Current backbone

- Microscopic mechanics: Newton/Hamiltonian dynamics
- Deterministic lattice baseline: generalized Lennard-Jones interaction
- Distribution energy: pair-distance distributions and correlation hierarchy
- Exact density conservation: $\partial_tP+\partial_a(Pv)=0$
- Hysteresis target: obtain $\oint \sigma\,d\epsilon>0$ without inserting an empirical hysteresis law
- Fatigue target: obtain secular cycle-to-cycle evolution $P_{N+1}\neq P_N$
- Crack initiation: formulate as loss of mechanical stability / first-passage into an unstable state

## Current progress

Milestone 1 now has a mechanics proof-of-principle using a resolved structural coordinate coupled to a semi-infinite harmonic Rubin chain. The full microscopic model is conservative, but outgoing unresolved lattice modes generate an exact phase lag and a nonzero reduced hysteresis loop.

Reference nondimensional result:

$$
A_H^{\mathrm{analytic}}=0.0152091700,
$$

$$
A_H^{\mathrm{numeric}}=0.0152088400.
$$

The relative loop-area error is about $2.17\times10^{-5}$ and the relative external-work / internal-energy balance error is about $1.25\times10^{-5}$.

See `docs/MILESTONE1_RUBIN_CHAIN.md` and `results/REFERENCE_RUN.md`.

This is **not yet fatigue accumulation**. The highest-priority theoretical problem remains deriving a mechanically justified slow structural state that gives

$$
P_{N+1}\neq P_N.
$$

## Run the current simulation

From the repository root:

```bash
python -m pip install -r requirements.txt
python -m simulations.run_rubin_hysteresis
python -m unittest tests.test_rubin_chain
```

## Repository structure

- `docs/` — theory notes, assumptions, derivations, failed approaches
- `theory/` — analytic and numerical model code
- `simulations/` — numerical experiments
- `tests/` — conservation, limiting-case, and falsification tests
- `firmware/` — eventual C/C++ reduced model and fatigue-tester integration
- `results/` — generated summaries and figures; raw experimental data should remain outside normal Git history

## Research rule

Every important result should be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that only reproduces a known fatigue curve by fitting is not considered a successful derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 피로 균열 개시를 미시역학에서부터 설명하기 위한 mechanics-first 연구 framework이다.

## 연구 목표

경험적인 피로법칙과 통제되지 않은 fitting을 가능한 한 배제하면서, 미시역학과 분포값 상태변수로부터 반복 히스테리시스, 비가역적 피로누적, 균열개시를 유도하는 것이 목표다.

## 핵심 상태변수

가장 기본적인 축약 상태는 열역학적 극한의 원자간격 밀도다.

$$
P(a,t)=\lim_{N\to\infty}\frac{1}{N}\sum_{i=1}^{N}\delta(a-a_i(t)).
$$

여기서 $a_i(t)$는 국부 원자간격을 나타내는 descriptor다. 장기 목표는 특정 확률분포 family를 미리 가정하는 대신 $P(a,t)$의 시간진화를 역학으로부터 유도하는 것이다.

## 현재 이론의 기본 골격

- 미시역학: Newton/Hamiltonian dynamics
- 결정론적 lattice baseline: generalized Lennard-Jones interaction
- 분포 에너지: pair-distance distribution 및 correlation hierarchy
- 정확한 밀도 보존식: $\partial_tP+\partial_a(Pv)=0$
- 히스테리시스 목표: 경험적 hysteresis law 없이 $\oint \sigma\,d\epsilon>0$를 얻는 것
- 피로누적 목표: cycle에 따라 $P_{N+1}\neq P_N$인 secular evolution을 얻는 것
- 균열개시: mechanical stability loss 또는 unstable state로의 first-passage 문제로 정식화하는 것

## 현재 진행상황

Milestone 1에서는 관심 구조좌표를 준무한 harmonic Rubin chain에 결합한 역학적 원리증명 모델을 구성했다. 전체 미시모델은 보존계지만, 외부로 전파되는 unresolved lattice mode 때문에 정확한 phase lag와 0이 아닌 축약 히스테리시스 loop가 발생한다.

기준 무차원 결과는

$$
A_H^{\mathrm{analytic}}=0.0152091700,
$$

$$
A_H^{\mathrm{numeric}}=0.0152088400
$$

이다.

loop area의 해석값-수치값 상대오차는 약 $2.17\times10^{-5}$이고, 외부 일과 내부에너지 수지의 상대오차는 약 $1.25\times10^{-5}$이다.

자세한 내용은 `docs/MILESTONE1_RUBIN_CHAIN.md`와 `results/REFERENCE_RUN.md`에 기록한다.

그러나 이것은 **아직 피로누적 자체가 아니다.** 현재 가장 우선적인 이론 문제는 다음 조건을 만드는 역학적으로 정당화된 느린 구조상태를 유도하는 것이다.

$$
P_{N+1}\neq P_N.
$$

## 현재 simulation 실행 방법

repository root에서 다음을 실행한다.

```bash
python -m pip install -r requirements.txt
python -m simulations.run_rubin_hysteresis
python -m unittest tests.test_rubin_chain
```

## Repository 구조

- `docs/` — 이론 노트, 가정, 유도, 실패한 접근
- `theory/` — 해석 및 수치 모델 코드
- `simulations/` — 수치실험
- `tests/` — 보존법칙, limiting case, 반증 test
- `firmware/` — 향후 C/C++ 축약모델 및 피로시험기 연동
- `results/` — 생성된 요약 및 figure; 원본 실험데이터는 일반 Git history 밖에서 관리하는 것을 원칙으로 함

## 연구 규칙

모든 중요한 결과는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

기존 피로곡선을 fitting으로 재현했을 뿐인 모델은 성공적인 이론 유도로 간주하지 않는다.
