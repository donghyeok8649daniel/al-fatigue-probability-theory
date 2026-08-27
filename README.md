# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation under **normal cyclic loading** in high-purity / single-crystal aluminum.

## Main research direction

The active theory is strictly centered on normal deformation and normal-opening instability:

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{normal hysteresis / memory}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{normal-opening instability}
}
$$

Earlier non-normal proof-of-principle work has been removed from the active repository. The current code, tests, results, and documentation are organized around normal interatomic spacing only.

## Core state

For local normal spacings $a_i(t)$,

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right),
$$

and

$$
\boxed{
P(a,t)=\lim_{N\to\infty}P_N(a,t).
}
$$

For deterministic microscopic trajectories,

$$
\boxed{
\partial_tP+\partial_a(Pv)=0,
}
$$

where

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle.
$$

The central mathematical problem is closure: determine the minimum microscopic state needed to derive $v(a,t)$ without inserting empirical fatigue evolution laws.

## Microscopic energy baseline

The main analytic baseline is a fixed generalized Lennard-Jones pair interaction

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

The LJ parameters do **not** change with cycle count. Fatigue evolution must arise from the atomic configuration, spacing distribution, correlations, memory, or mechanical instability.

The exact pair-distance energy hierarchy is

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
}
$$

## Active normal-chain model

The present reduced normal chain uses

$$
V=\sum_i\phi(\lambda_i),
\qquad
\lambda_i=\frac{a_i}{a_0},
$$

with generalized-LJ exponents

$$
m=12.19,
\qquad
n=6.
$$

The normalized potential is chosen so that

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

The idealized local normal stability-loss condition

$$
\phi''(\lambda_c)=0
$$

gives

$$
\boxed{\lambda_c\approx1.10777154.}
$$

Using the existing $E=69\,\mathrm{GPa}$ mapping gives an idealized 1D normal instability scale of about

$$
2.555\,\mathrm{GPa}.
$$

## Current numerical result

A 32-atom perfect-chain null test at

$$
\sigma_a=100\,\mathrm{MPa}
$$

does not produce false fatigue accumulation. No spacing crosses $\lambda_c$, and the work-energy relative error is approximately

$$
1.24\times10^{-10}.
$$

This is an intended falsification result and must not be tuned away.

A larger but statically subcritical dimensionless forcing can produce dynamic local opening at atomic-scale frequencies, showing that internal mode structure and history matter. However, this is **not** a 20 Hz fatigue prediction. The major unresolved problem is the physically derived time-scale bridge from fast atomistic LJ dynamics to slow cycle-to-cycle evolution of $P(a,t)$.

See `docs/MILESTONE2_NORMAL_DEFORMATION.md` and `results/reports/NORMAL_LJ_RESULTS.md`.

## Variable definitions

All active theory and simulation variables are defined in:

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`
- `firmware/VARIABLE_DEFINITIONS.md`

Any new variable must be added to the appropriate dictionary in the same commit that introduces it.

## Run the active simulation

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest tests.test_normal_lj_chain
```

## Repository structure

- `docs/` — active normal-theory derivations, assumptions, open problems, variable definitions, and firmware architecture
- `theory/` — active normal-deformation theory code
- `simulations/` — active normal numerical experiments
- `tests/` — normal-model conservation and falsification tests
- `results/data/` — machine-readable normal simulation results
- `results/figures/` — normal-deformation figures
- `results/reports/` — bilingual interpretation of normal results
- `firmware/` — hardware-independent fatigue-tester controller core
- `tools/` — PC-side telemetry / analysis helpers

## Research rule

Every important result must be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that reproduces a fatigue curve only by fitting is not considered a successful derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 **수직 반복하중** 아래 피로 균열개시를 미시역학으로부터 유도하기 위한 mechanics-first 연구 저장소다.

## 메인 연구방향

활성 이론은 수직변형과 normal-opening instability만을 중심으로 한다.

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{수직 히스테리시스 / memory}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{normal-opening instability}
}
$$

기존의 non-normal 원리증명 연구는 활성 저장소에서 제거했다. 현재 code, test, result, documentation은 수직 원자간격을 중심으로 정리한다.

## 핵심 상태변수

국부 수직 원자간격 $a_i(t)$에 대해

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right)
$$

를 정의하고,

$$
\boxed{
P(a,t)=\lim_{N\to\infty}P_N(a,t)
}
$$

를 중심 상태밀도로 사용한다.

결정론적 미시 trajectory에 대해서는

$$
\boxed{
\partial_tP+\partial_a(Pv)=0
}
$$

이 정확하게 성립하며,

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

핵심 수학 문제는 경험적인 fatigue evolution law를 넣지 않고 $v(a,t)$를 유도하기 위해 필요한 최소 미시상태를 결정하는 것이다.

## 미시에너지 baseline

주된 해석적 baseline은 고정된 generalized Lennard-Jones pair interaction이다.

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

LJ parameter는 cycle 수에 따라 변하지 않는다. 피로진화는 원자배열, spacing distribution, correlation, memory 또는 mechanical instability에서 나와야 한다.

정확한 pair-distance energy hierarchy는

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr
}
$$

이다.

## 활성 normal-chain 모델

현재 축약 normal chain은

$$
V=\sum_i\phi(\lambda_i),
\qquad
\lambda_i=\frac{a_i}{a_0}
$$

를 사용하고 generalized-LJ exponent는

$$
m=12.19,
\qquad
n=6
$$

이다.

정규화된 potential은

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

을 만족하도록 잡는다.

이상화된 국부 수직 안정성 상실조건

$$
\phi''(\lambda_c)=0
$$

으로부터

$$
\boxed{\lambda_c\approx1.10777154}
$$

를 얻는다.

기존 $E=69\,\mathrm{GPa}$ mapping을 사용하면 이상화된 1D normal instability scale은 약

$$
2.555\,\mathrm{GPa}
$$

이다.

## 현재 수치결과

32-atom 완전사슬에서

$$
\sigma_a=100\,\mathrm{MPa}
$$

수직 반복하중 null test를 수행하면 인공적인 피로누적이 생기지 않는다. 어떤 spacing도 $\lambda_c$를 넘지 않으며 work-energy 상대오차는 약

$$
1.24\times10^{-10}
$$

이다.

이것은 tuning으로 없애야 하는 실패가 아니라 의도된 반증결과다.

더 큰 statically subcritical dimensionless forcing에서는 atomic-scale frequency에서 dynamic local opening이 생길 수 있어 내부 mode 구조와 history의 중요성을 보여준다. 하지만 이것은 **20 Hz 피로예측이 아니다.** 현재 가장 큰 미해결 문제는 빠른 atomistic LJ dynamics에서 느린 cycle-to-cycle $P(a,t)$ 진화까지 이어지는 물리적으로 유도된 time-scale bridge다.

자세한 내용은 `docs/MILESTONE2_NORMAL_DEFORMATION.md`와 `results/reports/NORMAL_LJ_RESULTS.md`에 있다.

## 변수정의

활성 theory와 simulation의 모든 변수는 다음 문서에 정의한다.

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`
- `firmware/VARIABLE_DEFINITIONS.md`

새로운 변수를 도입하면 같은 commit에서 해당 변수사전을 반드시 갱신한다.

## 활성 simulation 실행

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest tests.test_normal_lj_chain
```

## Repository 구조

- `docs/` — 활성 normal 이론유도, 가정, open problem, 변수정의, firmware architecture
- `theory/` — 활성 normal-deformation theory code
- `simulations/` — 활성 normal 수치실험
- `tests/` — normal model의 conservation/falsification test
- `results/data/` — machine-readable normal simulation 결과
- `results/figures/` — normal-deformation figure
- `results/reports/` — normal result의 영문/한국어 해석
- `firmware/` — hardware-independent 피로시험기 controller core
- `tools/` — PC telemetry / analysis helper

## 연구 규칙

모든 중요한 결과는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

기존 피로곡선을 단순 fitting으로 재현했을 뿐인 모델은 성공적인 이론 유도로 간주하지 않는다.
