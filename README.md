# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation under **normal cyclic loading** in high-purity / single-crystal aluminum.

## Main research direction

The active mainline of this repository is normal deformation and normal-opening instability:

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

The primary microscopic energy baseline is a fixed generalized Lennard-Jones pair potential. Fatigue evolution must arise from atomic configuration and distribution evolution, not by changing the potential parameters with cycle count.

## Active normal model

The present normal-chain model uses

$$
V=\sum_i \phi(\lambda_i),
\qquad
\lambda_i=\frac{a_i}{a_0},
$$

with generalized LJ exponents $m=12.19$ and $n=6$. The normalized potential is chosen so that the reference state satisfies $\phi'(1)=0$ and $\phi''(1)=1$.

The local normal stability-loss condition is

$$
\phi''(\lambda_c)=0,
$$

which gives

$$
\lambda_c\approx1.10777154.
$$

Using the existing $E=69\,\mathrm{GPa}$ stress mapping, the corresponding idealized normal stress scale is about $2.555\,\mathrm{GPa}$.

A 100 MPa cyclic normal-loading null test does not produce false fatigue accumulation in the current perfect-chain model. This is an intended falsification result, not a failure to be tuned away.

See `docs/MILESTONE2_NORMAL_DEFORMATION.md`, `theory/normal_lj_chain.py`, and `results/reports/NORMAL_LJ_RESULTS.md`.

## Shear / auxiliary research library

Earlier Rubin-chain hysteresis, non-affine slip, gamma-surface, shear-oriented simulations, tests, and results are preserved under

`libraries/shear/`

They are retained as an auxiliary research library and historical proof-of-principle work. They are **not part of the active normal-deformation mainline** and are not run by the default normal test workflow.

## Active variable dictionaries

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — active normal-LJ variables
- `firmware/VARIABLE_DEFINITIONS.md` — fatigue-tester firmware fields
- `libraries/shear/docs/VARIABLE_DEFINITIONS.md` — archived broad/shear-era variables

## Run the active normal simulation

From the repository root:

```bash
python -m pip install -r requirements.txt
python -m simulations.run_normal_lj_chain
python -m unittest tests.test_normal_lj_chain
```

## Repository structure

- `docs/` — active normal-theory notes, assumptions, derivations, open problems, architecture, and style rules
- `theory/` — active normal-theory code
- `simulations/` — active normal numerical experiments
- `tests/` — active normal falsification/conservation tests
- `results/` — active normal numerical data, figures, and reports
- `libraries/shear/` — preserved auxiliary shear/Rubin/non-affine research library
- `firmware/` — hardware-independent fatigue-tester controller core
- `tools/` — PC-side telemetry / analysis helpers

## Research rule

Every important result should be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that reproduces a fatigue curve only by fitting is not considered a successful derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 **수직 반복하중** 아래 피로 균열개시를 미시역학에서 유도하기 위한 mechanics-first 연구 저장소다.

## 메인 연구방향

현재 저장소의 활성 mainline은 수직변형과 normal-opening instability다.

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

주된 미시 에너지 baseline은 고정된 generalized Lennard-Jones pair potential이다. 피로진화는 cycle 수에 따라 potential parameter를 바꾸는 방식이 아니라 원자배열과 분포의 진화에서 나와야 한다.

## 활성 normal 모델

현재 normal-chain 모델은

$$
V=\sum_i \phi(\lambda_i),
\qquad
\lambda_i=\frac{a_i}{a_0}
$$

를 사용하며 generalized LJ 지수는 $m=12.19$, $n=6$이다. 정규화된 potential은 기준상태에서 $\phi'(1)=0$, $\phi''(1)=1$을 만족하도록 잡는다.

국부 수직 안정성 상실조건은

$$
\phi''(\lambda_c)=0
$$

이고,

$$
\lambda_c\approx1.10777154
$$

를 준다.

기존 $E=69\,\mathrm{GPa}$ stress mapping을 사용하면 대응하는 이상화 normal stress scale은 약 $2.555\,\mathrm{GPa}$이다.

100 MPa 수직 반복하중 null test에서는 현재 완전사슬 모델이 가짜 피로누적을 만들지 않는다. 이것은 tuning으로 없애야 할 실패가 아니라 의도된 반증결과다.

자세한 내용은 `docs/MILESTONE2_NORMAL_DEFORMATION.md`, `theory/normal_lj_chain.py`, `results/reports/NORMAL_LJ_RESULTS.md`에 있다.

## 전단 / 보조 연구 라이브러리

기존 Rubin-chain 히스테리시스, non-affine slip, gamma-surface, 전단 지향 simulation, test, result는

`libraries/shear/`

아래에 보존한다.

이는 보조 연구 라이브러리와 역사적 proof-of-principle 작업으로 유지하지만, **활성 normal-deformation mainline에는 포함하지 않는다.** 기본 normal test workflow에서도 실행하지 않는다.

## 활성 변수사전

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — 활성 normal-LJ 변수
- `firmware/VARIABLE_DEFINITIONS.md` — 피로시험기 firmware 변수
- `libraries/shear/docs/VARIABLE_DEFINITIONS.md` — 보존된 과거 broad/shear 변수사전

## 활성 normal simulation 실행

repository root에서 다음을 실행한다.

```bash
python -m pip install -r requirements.txt
python -m simulations.run_normal_lj_chain
python -m unittest tests.test_normal_lj_chain
```

## Repository 구조

- `docs/` — 활성 normal 이론노트, 가정, 유도, open problem, architecture, style rule
- `theory/` — 활성 normal 이론코드
- `simulations/` — 활성 normal 수치실험
- `tests/` — 활성 normal 반증/보존 test
- `results/` — 활성 normal 수치데이터, figure, report
- `libraries/shear/` — 보존된 shear/Rubin/non-affine 보조 연구 라이브러리
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
