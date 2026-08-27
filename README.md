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

# 한국어 요약

이 repository의 목표는 순수/단결정 Al의 반복 하중에서 히스테리시스, 피로 누적, 균열 개시를 가능한 한 미시역학에서 직접 유도하는 것이다.

현재 Milestone 1에서는 **경험적 damping을 넣지 않은 전체 보존 Newton 사슬**을 이용해, 숨은 격자 모드로의 에너지 전달만으로 관심 구조좌표에서 $A_H>0$가 생길 수 있음을 해석식과 수치적분으로 확인했다.

그러나 이것은 아직 실제 Al 피로수명 모델이 아니다. 다음 핵심 문제는 한 cycle 뒤 구조가 완전히 복원되지 않는

$$
P_{N+1}\neq P_N
$$

의 미시역학적 원인을 유도하는 것이다.
