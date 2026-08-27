# Al Fatigue Probability Theory

Mechanics-first framework for fatigue crack initiation in high-purity / single-crystal aluminum, with the main physical focus on **cyclic normal stress, normal interatomic stretching, and normal-opening instability**.

## Research goal

The main target is

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{normal hysteresis}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{normal-opening instability / crack initiation}.
}
$$

The project tries to derive this chain from microscopic mechanics with as few phenomenological assumptions as possible. Empirical damage variables, fitted hysteresis laws, arbitrary transition kernels, and prescribed probability families are not accepted as starting axioms.

Shear/slip calculations are retained only as auxiliary mechanism tests. They are not the main physical hypothesis of the project.

## Why aluminum is used

High-purity / single-crystal Al is used as a target system for studying a normal-deformation-driven route as cleanly as possible.

This is a research-design choice, not a universal claim that every Al orientation is always weaker in normal opening than in shear. Any quantitative normal-versus-shear comparison must be checked separately for the chosen orientation, temperature, and microscopic model.

## Core state variable

For local normal interatomic spacings $a_i(t)$,

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta\!\left(a-a_i(t)\right),
$$

and the main distribution-valued state is

$$
\boxed{
P(a,t)=\lim_{N\to\infty}P_N(a,t).
}
$$

For deterministic trajectories,

$$
\boxed{
\partial_tP+\partial_a(Pv_a)=0,
}
$$

where

$$
v_a(a,t)=\langle\dot a_i\mid a_i=a\rangle.
$$

This is an exact kinematic identity. The central closure problem is to derive $v_a$ from the minimum necessary microscopic state.

## Primary energy baseline

The preferred analytic microscopic baseline is a fixed generalized Lennard-Jones pair interaction

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

The potential parameters do **not** evolve with fatigue. Structural evolution must come from the microscopic configuration and the induced distributions/correlations.

For exact $k$-th-neighbor distance densities $P_k(r,t)$,

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
}
$$

## Current progress

### Milestone 1 — mechanics-derived reduced hysteresis

A resolved coordinate coupled to a semi-infinite harmonic Rubin chain gives a nonzero reduced hysteresis loop even though the full microscopic model is conservative and contains no fitted viscous damping.

Reference nondimensional result:

$$
A_H^{\rm analytic}=0.0152091700,
$$

$$
A_H^{\rm numeric}=0.0152088400.
$$

Relative loop-area error:

$$
2.17\times10^{-5}.
$$

Relative work-energy error:

$$
1.25\times10^{-5}.
$$

This is an auxiliary existence proof for hidden-mode phase lag, not the final Al fatigue model.

### Milestone 2 — mainline normal generalized-LJ chain

The first direct normal-deformation simulation now exists in `theory/normal_lj_chain.py`.

The finite model is

$$
V=\sum_i\phi(\lambda_i),
\qquad
\lambda_i=x_{i+1}-x_i,
$$

with normalized generalized LJ energy

$$
\boxed{
\phi(\lambda)=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

and

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

For $m=12.19$, $n=6$,

$$
\boxed{
\lambda_c=1.1077715386
}
$$

from

$$
\phi''(\lambda_c)=0.
$$

The corresponding dimensionless static critical force is

$$
\boxed{
f_c=0.03703426967.
}
$$

Using the earlier mapping $f=\sigma/E$ with $E=69$ GPa gives the idealized 1D normal instability scale

$$
\sigma_c\approx2.5554\ \mathrm{GPa}.
$$

#### 100 MPa null test

For

$$
\sigma_a=100\ \mathrm{MPa},
$$

$$
f_a=\sigma_a/E=1.44927536\times10^{-3}.
$$

In the 32-atom 12-cycle reference run, no local spacing crossed $\lambda_c$.

The final recorded spacing variance was approximately

$$
7.55\times10^{-12},
$$

and the global work-energy relative error was

$$
\boxed{1.24\times10^{-10}}.
$$

This is an important null result: the perfect normal LJ chain does not invent fatigue merely because cyclic loading is present.

#### Dynamic sub-static-critical crossing

For a larger but still statically subcritical amplitude

$$
f_a=0.03<f_c,
$$

a direct conservative simulation at

$$
\omega^*=0.02
$$

first reached

$$
\lambda_{\max}\ge\lambda_c
$$

at approximately

$$
\boxed{N=2.25074\ \text{cycles}}.
$$

This crossing is strongly frequency dependent, showing that internal mode structure, phase, and history matter. It is not a single scalar stress-threshold problem.

However, this atomic-scale dynamical result is **not a 20 Hz fatigue prediction**.

Using the previous $a_0$, $A_0$, Al atomic mass, and $E=69$ GPa gives an atomic time scale of about

$$
t_0\approx5.55\times10^{-14}\ \mathrm{s}.
$$

Thus

$$
\omega^*=0.02
$$

corresponds to roughly

$$
5.73\times10^{10}\ \mathrm{Hz},
$$

while 20 Hz corresponds to only

$$
\boxed{
\omega^*_{20\mathrm{Hz}}\approx6.97\times10^{-12}.
}
$$

The resulting time-scale separation is now one of the central theoretical problems.

See:

- `docs/MILESTONE2_NORMAL_DEFORMATION.md`
- `results/reports/NORMAL_LJ_RESULTS.md`
- `results/data/normal_lj_summary.json`
- `results/data/normal_lj_cycle_history.csv`
- `results/figures/normal_lj_traction_stretch.svg`
- `results/figures/normal_lj_cycle_max_spacing.svg`
- `results/figures/normal_lj_spacing_distribution.svg`
- `results/figures/normal_lj_frequency_sweep.svg`

## Current theoretical bottleneck

The main unresolved problem is no longer whether an LJ normal model can be written. It can.

The main problem is to derive a slow, experimentally relevant evolution law from the exact microscopic dynamics across the enormous separation

$$
\text{atomic time scale}\ll\text{20 Hz cycle time}.
$$

The target remains

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

at the same cycle phase **without** inserting a fitted damage law.

Priority directions:

1. exact projected normal-spacing dynamics and memory;
2. fast/slow separation between phonon motion and slow structural state;
3. spacing-correlation hierarchy;
4. free-surface and geometry-defined normal opening;
5. finite-temperature phase-space ensembles;
6. first-passage formulation for normal instability.

## Reproduce the current simulations

```bash
python -m pip install -r requirements.txt
python -m simulations.run_normal_lj_chain
python -m unittest tests.test_normal_lj_chain
python -m simulations.run_rubin_hysteresis
python -m unittest tests.test_rubin_chain
```

The normal-LJ runner writes numerical data to `results/data/` and plots to `results/figures/`.

## Fatigue-tester firmware

A hardware-independent C99 axial-fatigue controller core exists under `firmware/`.

Implemented at the core level:

- sine / triangle normal-stress reference generation;
- stress-to-force conversion;
- cycle counting;
- PI load-cell force-loop structure and anti-windup;
- force, displacement, sensor-validity and E-stop fault handling;
- zero actuator command on fault;
- target-cycle stop;
- MCU hardware-abstraction boundary.

It is not yet a board-complete flash image because the final MCU, actuator drive, sensor electronics, and validated controller gains are not fixed.

## Variable definitions

Theory and simulation symbols are tracked in `docs/VARIABLE_DEFINITIONS.md`. Firmware fields and flags are tracked in `firmware/VARIABLE_DEFINITIONS.md`.

## Research rule

Every important statement should be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that reproduces a fatigue curve only by fitting is not considered a successful derivation.

---

# 한국어 번역

고순도/단결정 알루미늄의 피로 균열개시를 미시역학에서 설명하기 위한 mechanics-first framework이며, 주 물리적 관심은 **반복 수직응력, 수직 원자간 신장, 수직 opening instability**다.

## 연구 목표

주 목표는

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{수직 히스테리시스}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{수직 opening instability / crack initiation}
}
$$

를 가능한 한 적은 phenomenological assumption으로 microscopic mechanics에서 직접 유도하는 것이다.

경험적 damage variable, fitted hysteresis law, arbitrary transition kernel, prescribed probability family는 출발 공리로 사용하지 않는다.

전단/slip 계산은 보조 mechanism test로만 유지한다.

## 왜 Al을 쓰는가

고순도/단결정 Al을 수직변형 기반 failure route를 가능한 한 깨끗하게 연구하기 위한 대상계로 사용한다.

다만 모든 Al orientation이 언제나 전단보다 수직 opening에 약하다고 보편적으로 가정하지 않는다. 정량비교는 결정방향, 온도, 하중상태, microscopic model에 대해 별도로 검증해야 한다.

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

를 주 distribution-valued state로 둔다.

결정론적 trajectory에서는

$$
\boxed{
\partial_tP+\partial_a(Pv_a)=0
}
$$

이며

$$
v_a(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

이 식은 정확한 운동학적 항등식이고, 핵심 closure 문제는 필요한 최소 microscopic state에서 $v_a$를 유도하는 것이다.

## 주 에너지 baseline

주 해석 microscopic baseline은 고정된 generalized Lennard-Jones pair interaction이다.

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

피로가 진행된다고 potential parameter를 바꾸지 않는다. 구조진화는 microscopic configuration과 distribution/correlation의 변화로 나타나야 한다.

정확한 $k$-th-neighbor distance density $P_k(r,t)$에 대해서는

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr
}
$$

이다.

## 현재 진행상황

### Milestone 1 — mechanics-derived reduced hysteresis

관심 좌표를 준무한 harmonic Rubin chain에 결합하면 전체 미시계가 보존계이고 fitted viscous damping이 없어도 0이 아닌 reduced hysteresis가 생긴다.

기준 무차원 결과는

$$
A_H^{\rm analytic}=0.0152091700,
$$

$$
A_H^{\rm numeric}=0.0152088400
$$

이고, loop-area 상대오차는

$$
2.17\times10^{-5},
$$

work-energy 상대오차는

$$
1.25\times10^{-5}
$$

이다.

이것은 hidden-mode phase lag의 보조 existence proof이며 최종 Al 피로모델은 아니다.

### Milestone 2 — 메인 수직 generalized-LJ chain

첫 직접 수직변형 simulation을 `theory/normal_lj_chain.py`에 추가했다.

finite model은

$$
V=\sum_i\phi(\lambda_i),
\qquad
\lambda_i=x_{i+1}-x_i
$$

이고 normalized generalized LJ energy는

$$
\boxed{
\phi(\lambda)=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

이다.

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

이 되도록 normalization했다.

$m=12.19$, $n=6$에서

$$
\boxed{
\lambda_c=1.1077715386
}
$$

이고,

$$
\phi''(\lambda_c)=0
$$

이다.

이에 대응하는 무차원 static critical force는

$$
\boxed{
f_c=0.03703426967
}
$$

이다.

기존 $f=\sigma/E$, $E=69$ GPa mapping을 쓰면 idealized 1D normal instability scale은

$$
\sigma_c\approx2.5554\ \mathrm{GPa}
$$

이다.

#### 100 MPa null test

$$
\sigma_a=100\ \mathrm{MPa}
$$

이면

$$
f_a=\sigma_a/E=1.44927536\times10^{-3}
$$

이다.

32-atom, 12-cycle 기준 계산에서 어떤 local spacing도 $\lambda_c$를 넘지 않았다.

마지막 기록 spacing variance는 약

$$
7.55\times10^{-12}
$$

이고 global work-energy 상대오차는

$$
\boxed{1.24\times10^{-10}}
$$

이었다.

즉 perfect normal LJ chain은 cyclic loading이 있다는 이유만으로 가짜 fatigue를 만들지 않는다.

#### 정적 임계보다 낮은 dynamic crossing

더 큰 하중이지만 여전히

$$
f_a=0.03<f_c
$$

인 경우, 직접 보존동역학 simulation에서

$$
\omega^*=0.02
$$

일 때 처음으로

$$
\lambda_{\max}\ge\lambda_c
$$

가 되는 시점이 약

$$
\boxed{N=2.25074\ \text{cycle}}
$$

이었다.

이 결과는 frequency dependence가 매우 크다. 즉 internal mode structure, phase, history가 중요하고 단일 scalar stress threshold만으로 finite dynamics를 설명할 수 없다.

하지만 이것은 **20 Hz fatigue prediction이 아니다.**

기존 $a_0$, $A_0$, Al atomic mass, $E=69$ GPa를 쓰면 atomic time scale은

$$
t_0\approx5.55\times10^{-14}\ \mathrm{s}
$$

이고,

$$
\omega^*=0.02
$$

는 약

$$
5.73\times10^{10}\ \mathrm{Hz}
$$

에 해당한다.

반대로 20 Hz는

$$
\boxed{
\omega^*_{20\mathrm{Hz}}\approx6.97\times10^{-12}
}
$$

정도다.

따라서 atomic time scale과 실제 fatigue-test time scale 사이의 엄청난 분리가 이제 핵심 이론문제다.

관련 파일:

- `docs/MILESTONE2_NORMAL_DEFORMATION.md`
- `results/reports/NORMAL_LJ_RESULTS.md`
- `results/data/normal_lj_summary.json`
- `results/data/normal_lj_cycle_history.csv`
- `results/figures/normal_lj_traction_stretch.svg`
- `results/figures/normal_lj_cycle_max_spacing.svg`
- `results/figures/normal_lj_spacing_distribution.svg`
- `results/figures/normal_lj_frequency_sweep.svg`

## 현재 이론적 병목

이제 문제는 LJ normal model을 만들 수 있느냐가 아니다. 이미 만들었다.

남은 핵심은

$$
\text{atomic time scale}\ll\text{20 Hz cycle time}
$$

이라는 거대한 time-scale separation을 가로질러 실험적으로 의미 있는 slow evolution을 유도하는 것이다.

목표는 여전히 동일한 cycle phase에서

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

를 empirical damage law 없이 얻는 것이다.

우선순위는 다음과 같다.

1. exact projected normal-spacing dynamics와 memory;
2. phonon motion과 slow structural state의 fast/slow separation;
3. spacing-correlation hierarchy;
4. free-surface / geometry-defined normal opening;
5. finite-temperature phase-space ensemble;
6. normal instability의 first-passage formulation.

## 현재 simulation 재현

```bash
python -m pip install -r requirements.txt
python -m simulations.run_normal_lj_chain
python -m unittest tests.test_normal_lj_chain
python -m simulations.run_rubin_hysteresis
python -m unittest tests.test_rubin_chain
```

normal-LJ runner는 수치데이터를 `results/data/`, 그래프를 `results/figures/`에 저장한다.

## 피로시험기 firmware

`firmware/`에는 수직 축방향 fatigue test용 hardware-independent C99 controller core가 있다.

현재 구현된 core 기능은 normal-stress reference 생성, stress-to-force 변환, cycle count, load-cell PI force loop, anti-windup, force/displacement/sensor/E-stop fault handling, fault 시 zero command, target-cycle stop, MCU HAL boundary다.

최종 MCU, actuator drive, sensor electronics, 검증된 controller gain이 확정되지 않아 아직 board-complete flash image는 아니다.

## 변수정의

이론/simulation 기호는 `docs/VARIABLE_DEFINITIONS.md`, firmware field/flag는 `firmware/VARIABLE_DEFINITIONS.md`에서 관리한다.

## 연구 규칙

중요한 문장은 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

기존 fatigue curve를 fitting으로만 재현한 모델은 성공적인 유도로 간주하지 않는다.
