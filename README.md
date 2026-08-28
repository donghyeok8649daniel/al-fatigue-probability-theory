# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation under **normal cyclic loading** in high-purity / single-crystal aluminum.

## Active research direction

The repository root is the active normal-deformation mainline:

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{normal memory / hysteresis}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{normal-opening instability}
}
$$

Earlier Rubin-chain, slip, gamma-surface, and shear-oriented work is preserved under `libraries/shear/` and is not part of the default active workflow.

## Core state

For local normal spacings $a_i(t)$,

$$
P_N(a,t)=\frac1N\sum_{i=1}^{N}\delta\!\left(a-a_i(t)\right),
$$

and the central state density is

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

The main closure problem is to derive $v(a,t)$, or the minimum enlarged state required to determine it, from microscopic mechanics rather than from an empirical fatigue evolution law.

## Microscopic interaction baseline

The active analytic baseline is a fixed generalized Lennard-Jones pair interaction

$$
\boxed{
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
}
$$

The LJ parameters do **not** evolve with cycle count. Structural evolution must come from atomic configuration, spacing distributions, correlations, memory, or stability loss.

The exact pair-distance energy hierarchy is

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
}
$$

## Current active models

### 1. Reduced 1D normal chain

`theory/normal_lj_chain.py` is the current direct cyclic-dynamics null model.

It uses

$$
V=\sum_i\phi(\lambda_i),
\qquad
\lambda_i=\frac{a_i}{a_0},
$$

with

$$
m=12.19,\qquad n=6.
$$

A 32-atom perfect-chain calculation at $100$ MPa normal stress amplitude does **not** produce artificial fatigue accumulation. No local spacing crosses the idealized LJ tangent-instability stretch and the work-energy balance error is approximately

$$
1.24\times10^{-10}.
$$

This is an intended null/falsification result.

See `results/reports/NORMAL_LJ_RESULTS.md`.

### 2. Three-dimensional FCC normal lattice sum

`theory/fcc_normal_lj.py` removes the nearest-neighbor 1D geometry and evaluates

$$
\boxed{
U(\mathbf F)
=\frac12\sum_{\mathbf R\neq0}
v(|\mathbf F\mathbf R|)
}
$$

directly on an FCC lattice.

For [001] normal loading,

$$
\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n),
$$

and $\lambda_t$ is relaxed by

$$
\frac{\partial U}{\partial\lambda_t}=0.
$$

Using external Al reference values

$$
C_{11}=107\ \text{GPa},\qquad C_{12}=61\ \text{GPa},
$$

gives

$$
E_{[001]}\approx62.7024\ \text{GPa},
\qquad
\nu_{[001]}\approx0.363095.
$$

With $(m,n)=(12.19,6)$ and one LJ energy-scale calibration to $E_{[001]}$, the FCC pair model predicts

$$
C_{11}^{\rm LJ}\approx107.169\ \text{GPa},
$$

$$
C_{12}^{\rm LJ}\approx61.180\ \text{GPa},
$$

so the **normal elastic sector is reproduced very closely**.

The same model predicts an unfitted relaxed [001] ideal engineering strength of approximately

$$
\boxed{9.045\ \text{GPa}.}
$$

A first-principles reference reports about $10.63$ GPa for pure Al [001] tension, so the strength scale is reasonable without fitting the peak.

However, two exact/quantitative limitations are now explicit:

1. a cubic central pair potential obeys the Cauchy relation
   $$
   \boxed{C_{12}=C_{44}},
   $$
   while real Al has approximately $C_{12}=61$ GPa and $C_{44}=29$ GPa;
2. fitting the normal modulus gives a cohesive energy of only about
   $$
   0.976\ \text{eV/atom},
   $$
   whereas the experimental scale is about $3.43$ eV/atom. Fitting cohesion instead makes $E_{[001]}\approx220.5$ GPa.

Therefore the generalized LJ pair law is currently a strong **effective normal-mechanics baseline**, but it is not yet a quantitatively valid thermodynamic cohesive-energy model.

See `docs/MILESTONE3_FCC_NORMAL_LJ.md`.

## Consequence for the next theory step

The current result prevents a premature shortcut.

A thermal first-passage model containing

$$
\exp\!\left(-\frac{\Delta U}{k_BT}\right)
$$

cannot be quantitatively trusted if the absolute separation-energy scale is wrong, even when the normal elastic response is good.

So the next active problem is

$$
\boxed{
\text{preserve successful LJ normal mechanics}
\; + \;
\text{derive the minimum physically required cohesive/many-body correction}
}
$$

before introducing thermal escape rates or fatigue-life predictions.

## Variable definitions

Active variables are defined in:

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — spacing, probability, 1D normal dynamics, first passage;
- `docs/VARIABLE_DEFINITIONS_FCC_NORMAL_LJ.md` — FCC geometry, $\mathbf F$, $\lambda_n$, $\lambda_t$, lattice sums, cubic elastic constants and calibration variables;
- `firmware/VARIABLE_DEFINITIONS.md` — tester firmware fields and fault flags.

Any new variable must be documented in the appropriate dictionary in the same change that introduces it.

## Reproduce active results

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

The active result generator runs both the reduced normal chain and the FCC normal-LJ validation.

## Repository structure

- `docs/` — active normal-theory derivations, assumptions, variable definitions and open problems
- `theory/` — active normal mechanics code
- `simulations/` — active normal numerical experiments
- `tests/` — conservation, calibration and falsification tests
- `results/data/` — machine-readable results
- `results/figures/` — generated normal-deformation figures
- `results/reports/` — bilingual result interpretation
- `firmware/` — hardware-independent fatigue-tester controller core
- `tools/` — PC telemetry / analysis helpers
- `libraries/shear/` — preserved auxiliary shear/Rubin/slip research library

## Research rule

Every important result must be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that only reproduces a fatigue curve through fitted damage parameters is not considered a successful mechanics derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 **수직 반복하중** 아래 피로 균열개시를 미시역학에서부터 유도하기 위한 mechanics-first 연구 저장소다.

## 활성 연구방향

repository root의 active mainline은 수직변형이다.

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{수직 memory / hysteresis}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{normal-opening instability}
}
$$

기존 Rubin-chain, slip, gamma-surface, shear-oriented 연구는 삭제하지 않고 `libraries/shear/`에 보존하며 기본 active workflow에는 포함하지 않는다.

## 핵심 상태변수

국부 수직 원자간격 $a_i(t)$에 대해

$$
P_N(a,t)=\frac1N\sum_{i=1}^{N}\delta\!\left(a-a_i(t)\right)
$$

를 정의하고 중심 상태밀도는

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t)}
$$

이다.

결정론적 microscopic trajectory에서는

$$
\boxed{\partial_tP+\partial_a(Pv)=0}
$$

이 정확하게 성립하며

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

핵심 closure 문제는 경험적 fatigue evolution law 없이 $v(a,t)$ 또는 이를 결정하기 위해 필요한 최소 enlarged state를 microscopic mechanics에서 유도하는 것이다.

## 미시 상호작용 baseline

활성 해석 baseline은 고정 generalized Lennard-Jones pair interaction이다.

$$
\boxed{
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
}
$$

LJ parameter는 cycle에 따라 변하지 않는다. 구조진화는 atomic configuration, spacing distribution, correlation, memory 또는 stability loss에서 나와야 한다.

정확한 pair-distance energy hierarchy는

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr
}
$$

이다.

## 현재 활성 모델

### 1. 축약 1D normal chain

`theory/normal_lj_chain.py`는 현재 직접 cyclic dynamics를 확인하는 null model이다.

$$
V=\sum_i\phi(\lambda_i),
\qquad
\lambda_i=\frac{a_i}{a_0}
$$

을 사용하고

$$
m=12.19,\qquad n=6
$$

이다.

32-atom 완전사슬에 $100$ MPa normal stress amplitude를 가한 계산에서는 인공적인 fatigue accumulation이 생기지 않는다. 어떤 local spacing도 이상화된 LJ tangent-instability stretch를 넘지 않았고 work-energy balance error는 약

$$
1.24\times10^{-10}
$$

이다.

이것은 의도된 null/falsification result다.

자세한 결과는 `results/reports/NORMAL_LJ_RESULTS.md`에 있다.

### 2. 3차원 FCC normal lattice sum

`theory/fcc_normal_lj.py`는 1D nearest-neighbor geometry를 제거하고

$$
\boxed{
U(\mathbf F)
=\frac12\sum_{\mathbf R\neq0}v(|\mathbf F\mathbf R|)
}
$$

를 FCC lattice에 직접 계산한다.

[001] normal loading에서는

$$
\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n)
$$

이고 $\lambda_t$는

$$
\frac{\partial U}{\partial\lambda_t}=0
$$

으로 relaxation한다.

외부 Al reference

$$
C_{11}=107\ \text{GPa},\qquad C_{12}=61\ \text{GPa}
$$

로부터

$$
E_{[001]}\approx62.7024\ \text{GPa},
\qquad
\nu_{[001]}\approx0.363095
$$

를 얻는다.

$(m,n)=(12.19,6)$을 유지하고 LJ energy scale 하나만 $E_{[001]}$에 맞추면

$$
C_{11}^{\rm LJ}\approx107.169\ \text{GPa},
$$

$$
C_{12}^{\rm LJ}\approx61.180\ \text{GPa}
$$

가 되어 **normal elastic sector를 매우 잘 재현**한다.

그리고 peak strength를 fitting하지 않았는데 relaxed [001] ideal engineering strength는 약

$$
\boxed{9.045\ \text{GPa}}
$$

가 나온다.

first-principles reference의 pure Al [001] ideal strength 약 $10.63$ GPa와 비교하면 strength scale도 상당히 근접하다.

하지만 두 가지 한계도 정확히 드러난다.

1. cubic central pair potential은
   $$
   \boxed{C_{12}=C_{44}}
   $$
   라는 Cauchy relation을 만족해야 하지만 실제 Al은 대략 $C_{12}=61$ GPa, $C_{44}=29$ GPa다.
2. normal modulus를 맞추면 cohesive energy는 약
   $$
   0.976\ \text{eV/atom}
   $$
   에 불과하지만 experimental scale은 약 $3.43$ eV/atom이다. 반대로 cohesion을 맞추면 $E_{[001]}\approx220.5$ GPa가 된다.

따라서 generalized LJ pair law는 현재 **effective normal-mechanics baseline**으로는 매우 강하지만, 아직 quantitative thermodynamic cohesive-energy model은 아니다.

자세한 내용은 `docs/MILESTONE3_FCC_NORMAL_LJ.md`에 있다.

## 다음 이론단계에 주는 의미

현재 결과는 성급한 shortcut 하나를 막아준다.

$$
\exp\!\left(-\frac{\Delta U}{k_BT}\right)
$$

같은 thermal first-passage factor는 absolute energy scale에 지수적으로 민감하다. normal elasticity가 잘 맞더라도 separation energy가 틀리면 quantitative thermal escape rate에 그대로 사용할 수 없다.

그래서 다음 active problem은

$$
\boxed{
\text{성공적인 LJ normal mechanics 유지}
\; + \;
\text{최소한의 물리적으로 필요한 cohesive/many-body correction 유도}
}
$$

이다.

thermal escape rate나 fatigue-life prediction은 그 다음에 넣는다.

## 변수정의

활성 변수는 다음 문서에 정의한다.

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — spacing, probability, 1D normal dynamics, first passage;
- `docs/VARIABLE_DEFINITIONS_FCC_NORMAL_LJ.md` — FCC geometry, $\mathbf F$, $\lambda_n$, $\lambda_t$, lattice sum, cubic elastic constants, calibration variables;
- `firmware/VARIABLE_DEFINITIONS.md` — tester firmware field와 fault flag.

새 변수를 도입하면 같은 변경에서 적절한 변수사전을 갱신한다.

## 활성 결과 재생성

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

active result generator는 축약 normal chain과 FCC normal-LJ validation을 모두 실행한다.

## Repository 구조

- `docs/` — active normal theory derivation, assumption, variable definition, open problem
- `theory/` — active normal mechanics code
- `simulations/` — active normal numerical experiment
- `tests/` — conservation, calibration, falsification test
- `results/data/` — machine-readable result
- `results/figures/` — normal-deformation figure
- `results/reports/` — bilingual result interpretation
- `firmware/` — hardware-independent fatigue-tester controller core
- `tools/` — PC telemetry / analysis helper
- `libraries/shear/` — 보존된 auxiliary shear/Rubin/slip 연구 library

## 연구 규칙

모든 중요한 결과는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

fitted damage parameter만으로 fatigue curve를 재현하는 모델은 성공적인 mechanics derivation으로 간주하지 않는다.
