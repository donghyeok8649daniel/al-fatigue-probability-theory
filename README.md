# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation under **normal cyclic loading** in high-purity / single-crystal aluminum.

## Main research direction

The active theory is centered on normal deformation and normal-opening instability:

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

The repository root is the active normal-deformation mainline. Earlier Rubin-chain, non-affine slip, gamma-surface, and shear-oriented proof-of-principle work is **preserved**, not deleted, under `libraries/shear/`.

## Core state

For local normal spacings $a_i(t)$,

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right),
$$

and

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t).}
$$

For deterministic microscopic trajectories,

$$
\boxed{\partial_tP+\partial_a(Pv)=0,}
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

The LJ parameters do **not** change with cycle count. Fatigue evolution must arise from atomic configuration, spacing distribution, correlations, projected memory, or mechanical instability.

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

with

$$
m=12.19,
\qquad
n=6.
$$

The normalized potential satisfies

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

Using the current $E=69\,\mathrm{GPa}$ mapping gives an idealized 1D normal instability scale of about $2.555\,\mathrm{GPa}$.

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

A larger but statically subcritical dimensionless forcing can produce dynamic local opening at atomic-scale frequencies. That result shows that mode structure and history matter, but it is **not** a 20 Hz fatigue prediction. The current highest-priority problem is the physically derived time-scale bridge from fast LJ atomistic dynamics to slow cycle-to-cycle evolution of $P(a,t)$.

See `docs/MILESTONE2_NORMAL_DEFORMATION.md` and `results/reports/NORMAL_LJ_RESULTS.md`.

## Preserved shear / auxiliary library

All earlier non-normal work is grouped under

```text
libraries/shear/
├── README.md
├── docs/
├── theory/
├── simulations/
├── tests/
└── results/
```

This includes the Rubin-chain hysteresis proof, Hamiltonian slip-bath model, gamma-surface notes, pure-Al shear constraints, historical tests, data, and figures. These files are retained for comparison and future coupled-mode studies but are not imported by the default normal workflow.

## Variable definitions

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — active normal-LJ variables
- `firmware/VARIABLE_DEFINITIONS.md` — fatigue-tester firmware variables
- `libraries/shear/docs/VARIABLE_DEFINITIONS.md` — preserved historical/shear variables

Any new active variable must be added to the appropriate active dictionary in the same commit.

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
- `results/` — active normal data, figures, and reports
- `libraries/shear/` — preserved Rubin/shear/non-affine auxiliary research library
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

활성 이론은 수직변형과 normal-opening instability를 중심으로 한다.

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

repository root는 활성 normal-deformation mainline이다. 기존 Rubin-chain, non-affine slip, gamma-surface, shear-oriented proof-of-principle 연구는 삭제하지 않고 `libraries/shear/` 아래에 보존한다.

## 핵심 상태변수

국부 수직 원자간격 $a_i(t)$에 대해

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right)
$$

를 정의하고,

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t)}
$$

를 중심 상태밀도로 사용한다.

결정론적 미시 trajectory에 대해서는

$$
\boxed{\partial_tP+\partial_a(Pv)=0}
$$

이 정확하게 성립하며,

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

핵심 수학 문제는 경험적 fatigue evolution law를 넣지 않고 $v(a,t)$를 유도하기 위해 필요한 최소 미시상태를 결정하는 것이다.

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

LJ parameter는 cycle 수에 따라 변하지 않는다. 피로진화는 원자배열, spacing distribution, correlation, projected memory 또는 mechanical instability에서 나와야 한다.

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

를 사용하고

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

을 만족한다.

이상화된 국부 수직 안정성 상실조건

$$
\phi''(\lambda_c)=0
$$

으로부터

$$
\boxed{\lambda_c\approx1.10777154}
$$

를 얻는다.

현재 $E=69\,\mathrm{GPa}$ mapping을 사용하면 이상화된 1D normal instability scale은 약 $2.555\,\mathrm{GPa}$이다.

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

더 큰 statically subcritical dimensionless forcing에서는 atomic-scale frequency에서 dynamic local opening이 생길 수 있다. 이 결과는 mode 구조와 history의 중요성을 보여주지만 **20 Hz 피로예측은 아니다.** 현재 최우선 문제는 빠른 LJ atomistic dynamics에서 느린 cycle-to-cycle $P(a,t)$ 진화로 가는 물리적으로 유도된 time-scale bridge다.

자세한 내용은 `docs/MILESTONE2_NORMAL_DEFORMATION.md`와 `results/reports/NORMAL_LJ_RESULTS.md`에 있다.

## 보존된 전단 / 보조 library

기존 non-normal 연구는 전부

```text
libraries/shear/
├── README.md
├── docs/
├── theory/
├── simulations/
├── tests/
└── results/
```

아래에 묶어 보존한다.

여기에는 Rubin-chain hysteresis proof, Hamiltonian slip-bath model, gamma-surface note, pure-Al shear constraint, 과거 test/data/figure가 포함된다. 비교연구와 향후 coupled-mode 연구를 위해 남겨두지만 기본 normal workflow에서는 import하지 않는다.

## 변수정의

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — 활성 normal-LJ 변수
- `firmware/VARIABLE_DEFINITIONS.md` — 피로시험기 firmware 변수
- `libraries/shear/docs/VARIABLE_DEFINITIONS.md` — 보존된 과거/shear 변수

새로운 활성변수는 같은 commit에서 해당 활성 변수사전에 추가한다.

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
- `tests/` — normal model conservation/falsification test
- `results/` — 활성 normal data, figure, report
- `libraries/shear/` — 보존된 Rubin/shear/non-affine 보조 연구 library
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
