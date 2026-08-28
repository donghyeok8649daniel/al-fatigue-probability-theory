# Al Fatigue Probability Theory

Mechanics-first framework for fatigue crack initiation under **one-dimensional normal cyclic loading** in high-purity / single-crystal aluminum.

## Active research direction

The active repository mainline is deliberately restricted to

$$
\boxed{
\sigma_n(t)
\rightarrow
P(a,t)
\rightarrow
\mu(t),\mathcal E(t)
\rightarrow
\text{crack-free energy feasibility}
\rightarrow
\text{normal-opening tail}
\rightarrow
\tau_c.
}
$$

The fundamental evolution variable is physical time $t$, not fatigue cycle count.

If the loading frequency is constant, an equivalent cycle count may later be reported as

$$
N=f t,
$$

but $N$ is not used as the state-evolution coordinate.

All active theory is **1D, normal-only, generalized-LJ based**. Earlier shear/Rubin work is preserved under `libraries/shear/`. Earlier 3D FCC normal-LJ work is preserved under `libraries/fcc_normal/`. Neither archive is part of the default active workflow.

## Core probability state

For finite local spacings $a_i(t)$,

$$
P_M(a,t)
=
\frac1M\sum_{i=1}^M\delta(a-a_i(t)),
$$

where $M$ denotes a finite number of represented spacings, not fatigue cycles.

The continuum/thermodynamic state is

$$
\boxed{
P(a,t).
}
$$

The density is not assumed to be Gaussian, Weibull, or any other named family.

Its normalization and mean are

$$
\boxed{
\int P(a,t)\,da=1,
}
$$

$$
\boxed{
\bar a(t)=\int aP(a,t)\,da.
}
$$

For deterministic microscopic trajectories,

$$
\boxed{
\partial_tP+\partial_a(Pv)=0,
}
$$

with

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle.
$$

## Fixed generalized Lennard-Jones baseline

The active microscopic interaction is

$$
\boxed{
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}r\right)^m
-
\left(\frac{\sigma_{\rm LJ}}r\right)^n
\right].
}
$$

The potential parameters do not evolve with time or fatigue history.

For the reduced normalized 1D chain,

$$
\lambda=\frac{a}{a_0},
$$

and

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
$$

with

$$
m=12.19,
\qquad
n=6.
$$

The normalization gives

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

The idealized normal tangent-stability limit is

$$
\boxed{
\phi''(\lambda_c)=0,
}
$$

with

$$
\boxed{
\lambda_c\approx1.1077715386.
}
$$

## Current null result

The existing perfect 32-atom 1D reference simulation at a 100 MPa normal stress amplitude remains an important null test. It does not generate artificial fatigue accumulation and no local spacing crosses $\lambda_c$ in the reference run.

The work-energy balance error is approximately

$$
1.24\times10^{-10}.
$$

This result must not be tuned away.

## New continuous-time energy formulation

Define the shifted LJ energy

$$
\boxed{
\psi(\lambda)=\phi(\lambda)-\phi(1),
}
$$

and

$$
\boxed{
\mu(t)=\int\lambda P(\lambda,t)\,d\lambda,
}
$$

$$
\boxed{
\mathcal E(t)
=
\int\psi(\lambda)P(\lambda,t)\,d\lambda.
}
$$

At nearly conserved mean, energy above the homogeneous value is exactly

$$
\boxed{
\mathcal E(t)-\psi(\mu(t))
=
\int
\left[
\psi(\lambda)-\psi(\mu)-\psi'(\mu)(\lambda-\mu)
\right]
P(\lambda,t)\,d\lambda.
}
$$

Inside the convex LJ region this quantity is non-negative and measures distributional spreading away from the mean.

## Exact impossibility result

Normalization, mean, and energy alone do **not** force a tensile tail.

Because

$$
\psi(\lambda)\rightarrow\infty
\qquad
(\lambda\rightarrow0^+),
$$

arbitrarily large energy can mathematically be stored by reverse compression while all probability remains below $\lambda_c$.

Therefore one more physical condition is necessary.

The current minimal condition is a mechanically justified lower support bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0.
}
$$

This is not a fitting parameter. Deriving or independently constraining $\lambda_L(t)$ is now the main physical problem.

## Exact crack-free energy ceiling

Assume

$$
\operatorname{supp}P(\lambda,t)
\subset
[\lambda_L(t),\lambda_c].
$$

Since the generalized LJ energy is convex throughout this interval, the exact minimum and maximum possible energies at a fixed mean are

$$
\boxed{
\mathcal E_{\rm safe}^{\min}(t)=\psi(\mu(t)),
}
$$

and

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_c).
}
$$

The maximum is attained by an endpoint two-point measure, so this is an exact extremum over all admissible probability distributions, not a chosen distribution ansatz.

Define

$$
\boxed{
M_E(t)=\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t).
}
$$

If a valid hard compression bound has been established and

$$
M_E(t)<0,
$$

then no crack-free probability distribution can satisfy normalization, mean, energy, and support simultaneously. Therefore probability mass beyond $\lambda_c$ becomes unavoidable.

The corresponding continuous-time first passage is

$$
\boxed{
\tau_E
=
\inf\{t\ge0:M_E(t)<0\}.
}
$$

See `docs/MILESTONE4_TIME_ENERGY_FEASIBILITY.md`.

## Current next step

The highest-priority problem is now

$$
\boxed{
\text{derive }\lambda_L(t)\text{ from 1D normal-LJ mechanics.}
}
$$

No 3D model, shear coordinate, fitted damage law, fitted relaxation time, or assumed probability family is required for this step.

## Active code

- `theory/normal_lj_chain.py` — conservative 1D LJ normal-chain dynamics
- `theory/normal_lj_energy_feasibility.py` — exact probability-measure energy bounds
- `theory/normal_lj_timescale.py` — retained 1D time-scale falsification diagnostic
- `simulations/run_normal_lj_chain.py` — direct 1D null/dynamic simulation
- `simulations/run_normal_lj_energy_feasibility.py` — exact energy-ceiling parameter sweep
- `tests/test_normal_lj_chain.py`
- `tests/test_normal_lj_energy_feasibility.py`

## Variable definitions

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`
- `docs/VARIABLE_DEFINITIONS_ENERGY_FEASIBILITY.md`
- `docs/VARIABLE_DEFINITIONS_NORMAL_TIMESCALE.md`
- `firmware/VARIABLE_DEFINITIONS.md`

## Reproduce active results

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

The default workflow runs only active 1D normal-LJ calculations.

## Repository structure

- `docs/` — active 1D normal theory and variable definitions
- `theory/` — active 1D normal-LJ mathematical code
- `simulations/` — active 1D numerical/analytical runners
- `tests/` — active 1D falsification and theorem tests
- `results/` — active 1D result data, figures, and reports
- `firmware/` — fatigue-tester controller core
- `libraries/shear/` — preserved shear/Rubin/slip archive
- `libraries/fcc_normal/` — preserved 3D FCC normal-LJ archive

## Research rule

Every important statement must be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**

A result obtained only by fitting a fatigue curve is not considered a successful mechanics derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 **1차원 수직 반복하중** 아래 피로 균열개시를 mechanics-first 방식으로 유도하기 위한 연구 framework다.

## 활성 연구방향

현재 repository mainline은 의도적으로 다음 경로로 제한한다.

$$
\boxed{
\sigma_n(t)
\rightarrow
P(a,t)
\rightarrow
\mu(t),\mathcal E(t)
\rightarrow
\text{crack-free energy feasibility}
\rightarrow
\text{normal-opening tail}
\rightarrow
\tau_c
}
$$

근본 evolution variable은 fatigue cycle count가 아니라 물리적 시간 $t$다.

loading frequency가 일정할 때 필요하면 나중에

$$
N=ft
$$

로 equivalent cycle count를 표시할 수 있지만 $N$을 state-evolution coordinate로 사용하지 않는다.

모든 활성 이론은 **1D, normal-only, generalized-LJ 기반**이다. 기존 shear/Rubin 연구는 `libraries/shear/`에 보존하고, 기존 3D FCC normal-LJ 연구는 `libraries/fcc_normal/`에 보존한다. 두 archive 모두 기본 active workflow에는 들어가지 않는다.

## 핵심 확률상태

유한한 local spacing $a_i(t)$에 대해

$$
P_M(a,t)
=
\frac1M\sum_{i=1}^M\delta(a-a_i(t))
$$

를 정의한다. 여기서 $M$은 fatigue cycle 수가 아니라 유한한 represented spacing 수다.

continuum/thermodynamic state는

$$
\boxed{P(a,t)}
$$

이다.

이 density를 Gaussian, Weibull 또는 다른 named family로 가정하지 않는다.

정규화와 평균은

$$
\boxed{
\int P(a,t)\,da=1
}
$$

및

$$
\boxed{
\bar a(t)=\int aP(a,t)\,da
}
$$

이다.

deterministic microscopic trajectory에서는

$$
\boxed{
\partial_tP+\partial_a(Pv)=0
}
$$

이고

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

## 고정 generalized Lennard-Jones baseline

활성 microscopic interaction은

$$
\boxed{
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}r\right)^m
-
\left(\frac{\sigma_{\rm LJ}}r\right)^n
\right]
}
$$

이다.

potential parameter는 시간이나 fatigue history에 따라 변하지 않는다.

축약 normalized 1D chain에서는

$$
\lambda=\frac{a}{a_0}
$$

이고

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
$$

이며

$$
m=12.19,
\qquad
n=6
$$

이다.

normalization에 의해

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

이다.

이상화된 normal tangent-stability limit는

$$
\boxed{
\phi''(\lambda_c)=0
}
$$

이고

$$
\boxed{
\lambda_c\approx1.1077715386
}
$$

이다.

## 현재 null result

기존 perfect 32-atom 1D reference simulation에서 normal stress amplitude 100 MPa를 가한 결과는 중요한 null test로 유지한다. reference run에서는 인공적인 fatigue accumulation이 생기지 않았고 어떤 local spacing도 $\lambda_c$를 넘지 않았다.

work-energy balance error는 약

$$
1.24\times10^{-10}
$$

이다.

이 결과를 tuning으로 없애면 안 된다.

## 새로운 연속시간 에너지 formulation

shifted LJ energy를

$$
\boxed{
\psi(\lambda)=\phi(\lambda)-\phi(1)
}
$$

로 정의하고

$$
\boxed{
\mu(t)=\int\lambda P(\lambda,t)\,d\lambda
}
$$

및

$$
\boxed{
\mathcal E(t)
=
\int\psi(\lambda)P(\lambda,t)\,d\lambda
}
$$

를 사용한다.

평균이 거의 보존될 때 homogeneous value를 초과하는 에너지는 정확히

$$
\boxed{
\mathcal E(t)-\psi(\mu(t))
=
\int
\left[
\psi(\lambda)-\psi(\mu)-\psi'(\mu)(\lambda-\mu)
\right]
P(\lambda,t)\,d\lambda
}
$$

이다.

convex LJ 영역에서는 이 값이 음수가 아니며 평균에서 분포가 퍼지는 정도와 연결된다.

## 정확한 impossibility result

정규화, 평균, 에너지만으로는 **tensile tail을 강제할 수 없다.**

$$
\psi(\lambda)\rightarrow\infty
\qquad
(\lambda\rightarrow0^+)
$$

이므로 모든 확률질량을 $\lambda_c$ 아래에 유지하면서 reverse compression으로 임의로 큰 에너지를 수학적으로 저장할 수 있다.

따라서 하나의 추가 물리조건이 필요하다.

현재 최소조건은 역학적으로 정당화된 lower support bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0
}
$$

이다.

이 값은 fitting parameter가 아니다. $\lambda_L(t)$를 1D mechanics로부터 유도하거나 독립적으로 제약하는 것이 현재 중심 물리문제다.

## 정확한 crack-free energy ceiling

$$
\operatorname{supp}P(\lambda,t)
\subset
[\lambda_L(t),\lambda_c]
$$

를 가정한다.

generalized LJ energy는 이 interval 전체에서 convex하므로 고정 평균에서 가능한 최소 및 최대 energy는 정확히

$$
\boxed{
\mathcal E_{\rm safe}^{\min}(t)=\psi(\mu(t))
}
$$

및

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_c)
}
$$

이다.

최댓값은 endpoint two-point measure에서 실제로 달성되므로 특정 distribution ansatz가 아니라 모든 admissible probability distribution에 대한 정확한 extremum이다.

$$
\boxed{
M_E(t)=\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t)
}
$$

를 정의한다.

유효한 hard compression bound가 확보된 상태에서

$$
M_E(t)<0
$$

이면 정규화, 평균, 에너지 및 support를 동시에 만족하는 crack-free probability distribution은 존재할 수 없다. 따라서 $\lambda_c$를 넘는 probability mass가 필수가 된다.

이에 대응하는 연속시간 first passage는

$$
\boxed{
\tau_E
=
\inf\{t\ge0:M_E(t)<0\}
}
$$

이다.

자세한 내용은 `docs/MILESTONE4_TIME_ENERGY_FEASIBILITY.md`에 있다.

## 현재 다음 단계

가장 우선적인 문제는

$$
\boxed{
\text{1D normal-LJ mechanics로부터 }\lambda_L(t)\text{를 유도하는 것}
}
$$

이다.

이 단계에는 3D model, shear coordinate, fitted damage law, fitted relaxation time 또는 assumed probability family가 필요하지 않다.

## 활성 code

- `theory/normal_lj_chain.py` — conservative 1D LJ normal-chain dynamics
- `theory/normal_lj_energy_feasibility.py` — 정확한 probability-measure energy bound
- `theory/normal_lj_timescale.py` — 유지되는 1D time-scale falsification diagnostic
- `simulations/run_normal_lj_chain.py` — 직접 1D null/dynamic simulation
- `simulations/run_normal_lj_energy_feasibility.py` — exact energy-ceiling parameter sweep
- `tests/test_normal_lj_chain.py`
- `tests/test_normal_lj_energy_feasibility.py`

## 변수정의

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`
- `docs/VARIABLE_DEFINITIONS_ENERGY_FEASIBILITY.md`
- `docs/VARIABLE_DEFINITIONS_NORMAL_TIMESCALE.md`
- `firmware/VARIABLE_DEFINITIONS.md`

## 활성 결과 재현

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

기본 workflow는 active 1D normal-LJ calculation만 실행한다.

## Repository 구조

- `docs/` — 활성 1D normal theory 및 변수정의
- `theory/` — 활성 1D normal-LJ mathematical code
- `simulations/` — 활성 1D numerical/analytical runner
- `tests/` — 활성 1D falsification 및 theorem test
- `results/` — 활성 1D result data, figure, report
- `firmware/` — fatigue-tester controller core
- `libraries/shear/` — 보존된 shear/Rubin/slip archive
- `libraries/fcc_normal/` — 보존된 3D FCC normal-LJ archive

## 연구 규칙

모든 중요한 statement는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**

fatigue curve fitting만으로 얻은 결과는 성공적인 mechanics derivation으로 보지 않는다.
