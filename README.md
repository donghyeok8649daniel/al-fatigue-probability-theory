# Al Fatigue Probability Theory

Mechanics-first framework for fatigue crack initiation in high-purity / single-crystal aluminum, with the main physical focus on **cyclic normal stress, normal interatomic stretching, and normal-opening instability**.

## Research goal

Derive normal hysteresis, irreversible cycle-to-cycle structural evolution, and crack initiation from microscopic mechanics and distribution-valued state variables while minimizing empirical fatigue laws and uncontrolled fitting.

The principal mechanics chain is

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

Shear/slip models are retained only as auxiliary proof-of-principle calculations unless future derivation shows that they are mathematically required for closure of the normal-deformation problem.

## Why aluminum is used in this project

The project uses high-purity / single-crystal Al as a target system for isolating and studying normal-deformation-driven failure as cleanly as possible.

This is a research design choice. The repository does **not** assume as a universal material law that every Al orientation or specimen is always weaker in normal opening than in shear. Any such quantitative comparison must be checked separately for the chosen crystal orientation, temperature, loading state, and microscopic model.

## Core state

The primary reduced state is the thermodynamic-limit spacing density

$$
P(a,t)=\lim_{N\to\infty}\frac{1}{N}\sum_{i=1}^{N}\delta(a-a_i(t)),
$$

where $a_i(t)$ is a local normal interatomic-spacing descriptor.

The long-term objective is to derive the evolution of $P(a,t)$ from microscopic mechanics rather than prescribe a probability family.

## Primary microscopic energy baseline

The principal analytic baseline is a fixed generalized Lennard-Jones pair potential

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

The potential itself does not change with fatigue damage. Structural evolution must appear through the microscopic configuration and the induced state distributions.

For exact pair-distance densities $P_k(r,t)$,

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
}
$$

For deterministic spacing trajectories, the exact kinematic conservation law is

$$
\boxed{
\partial_tP+\partial_a(Pv_a)=0,
}
$$

where

$$
v_a(a,t)=\langle\dot a_i\mid a_i=a\rangle.
$$

The central closure problem is to derive $v_a$ from the minimum necessary microscopic state without inserting a fitted constitutive law.

## Current progress

### Milestone 1 — reduced hysteresis mechanism

A resolved coordinate coupled to a semi-infinite harmonic Rubin chain gives a mechanics-derived phase lag and nonzero reduced hysteresis while the full microscopic system remains conservative.

Reference nondimensional result:

$$
A_H^{\rm analytic}=0.0152091700,
$$

$$
A_H^{\rm numeric}=0.0152088400.
$$

The relative loop-area error is approximately $2.17\times10^{-5}$ and the energy-balance relative error is approximately $1.25\times10^{-5}$.

This is an auxiliary existence proof. The mainline task is now to derive the analogous hysteresis directly in the normal-spacing sector.

### Milestone 2 — normal spacing accumulation

The current main target is

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

at identical cycle phase under cyclic normal stress.

A previous nonlinear Hamiltonian slip-bath simulation demonstrated $s_{N+1}\neq s_N$ without an empirical damage law. That result is retained as an auxiliary proof that conservative nonlinear microscopic dynamics can generate a nontrivial cycle map, but **shear slip is not the principal physical mechanism of this project**.

The next mainline simulation must obtain secular evolution directly from normal spacing, normal-mode coupling, correlations, free-surface opening, or another microscopically defined normal-deformation mechanism.

See:

- `docs/MILESTONE2_NORMAL_DEFORMATION.md` — current mainline theory
- `docs/OPEN_PROBLEMS.md` — current unresolved mechanics
- `docs/MILESTONE2_HAMILTONIAN_SLIP.md` — auxiliary shear/slip proof-of-principle

## Crack initiation

The principal crack-initiation route is a normal-opening stability / first-passage problem.

An idealized reduced-lattice baseline is

$$
U''(a_c)=0,
$$

but the final initiation theory must distinguish instantaneous tail occupancy

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

from cumulative first-passage probability.

## Variable dictionaries

Theory and simulation symbols are defined in `docs/VARIABLE_DEFINITIONS.md`. Firmware fields and fault flags are defined in `firmware/VARIABLE_DEFINITIONS.md`.

Any source or document that introduces a new research symbol must update the relevant variable dictionary in the same commit.

## Reproduce current proof-of-principle results

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest tests.test_rubin_chain
python -m unittest tests.test_hamiltonian_slip_bath
```

Generated numerical data are written under `results/data/` and figures under `results/figures/`.

These simulations are currently mechanism tests, not quantitative Al fatigue-life predictions.

## Fatigue-tester firmware status

A target-independent C99 real-time controller core exists under `firmware/` for normal axial fatigue loading.

Implemented at the core level:

- sine / triangle normal-stress reference generation;
- stress-to-force conversion;
- cycle counting;
- PI load-cell force-loop structure and anti-windup;
- force, displacement, sensor-validity and E-stop fault handling;
- zero command on fault;
- target-cycle stop;
- MCU hardware-abstraction boundary.

The firmware is not yet a board-complete flash image because the final MCU, actuator drive, sensor electronics, and validated control gains have not been fixed.

## Repository structure

- `docs/` — theory notes, assumptions, exact derivations, failed approaches, normal-deformation milestone, architecture, and variable definitions
- `theory/` — analytic and numerical model code
- `simulations/` — numerical experiments and reproducible result generation
- `tests/` — conservation, limiting-case, and falsification tests
- `firmware/` — hardware-independent axial-fatigue real-time controller core
- `tools/` — PC-side telemetry / analysis helpers
- `results/data/` — machine-readable simulation results
- `results/figures/` — generated plots
- `results/reports/` — bilingual interpretation of simulation results

## Research rule

Every important result must be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that reproduces a known fatigue curve only by fitting is not considered a successful derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 피로 균열개시를 미시역학에서 설명하기 위한 mechanics-first framework이며, 주 물리적 관심은 **반복 수직응력, 수직 원자간 신장, 수직 opening instability**이다.

## 연구 목표

경험적인 피로법칙과 통제되지 않은 fitting을 가능한 한 배제하면서 microscopic mechanics와 distribution-valued state variable로부터 수직 히스테리시스, cycle-to-cycle 비가역 구조진화, 균열개시를 유도한다.

주 역학 chain은

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

이다.

전단/slip 모델은 향후 수직변형 문제의 closure에 수학적으로 필요하다고 유도되지 않는 한 보조 원리증명으로만 유지한다.

## 이 프로젝트에서 Al을 사용하는 이유

고순도/단결정 Al을 수직변형 기반 파괴를 가능한 한 깨끗하게 분리하여 연구하기 위한 대상계로 사용한다.

이것은 연구설계상의 선택이다. 저장소에서는 모든 Al 결정방향이나 시편이 항상 전단보다 수직 opening에 더 약하다는 것을 보편적인 재료법칙으로 가정하지 않는다. 그러한 정량비교는 선택된 결정방향, 온도, 하중상태, microscopic model에 대해 별도로 검증해야 한다.

## 핵심 상태

가장 기본적인 축약상태는 열역학적 극한의 spacing density

$$
P(a,t)=\lim_{N\to\infty}\frac{1}{N}\sum_{i=1}^{N}\delta(a-a_i(t))
$$

이다.

여기서 $a_i(t)$는 국부 수직 원자간격 descriptor다. 특정 probability family를 미리 가정하지 않고 $P(a,t)$의 진화를 microscopic mechanics에서 유도하는 것이 장기 목표다.

## 주 미시에너지 baseline

주 해석 baseline은 고정된 generalized Lennard-Jones pair potential이다.

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

피로손상에 따라 potential 자체를 변화시키지 않는다. 구조진화는 microscopic configuration과 그로부터 유도되는 state distribution의 변화로 나타나야 한다.

정확한 pair-distance density $P_k(r,t)$에 대해

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr
}
$$

이다.

결정론적 spacing trajectory에 대한 정확한 운동학적 보존식은

$$
\boxed{
\partial_tP+\partial_a(Pv_a)=0
}
$$

이고,

$$
v_a(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

핵심 closure 문제는 fitted constitutive law를 넣지 않고 필요한 최소 microscopic state에서 $v_a$를 유도하는 것이다.

## 현재 진행상황

### Milestone 1 — 축약 히스테리시스 메커니즘

관심 좌표를 준무한 harmonic Rubin chain에 결합하면 전체 미시계가 보존계인 상태에서 mechanics-derived phase lag와 0이 아닌 축약 히스테리시스가 발생한다.

기준 무차원 결과는

$$
A_H^{\rm analytic}=0.0152091700,
$$

$$
A_H^{\rm numeric}=0.0152088400
$$

이다.

loop-area 상대오차는 약 $2.17\times10^{-5}$이고 energy-balance 상대오차는 약 $1.25\times10^{-5}$이다.

이 결과는 보조적인 existence proof다. 이제 메인 과제는 동일한 히스테리시스를 수직 spacing sector에서 직접 유도하는 것이다.

### Milestone 2 — normal spacing 누적진화

현재 메인 목표는 반복 수직응력에서 동일한 cycle phase에

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

를 얻는 것이다.

이전 nonlinear Hamiltonian slip-bath simulation은 empirical damage law 없이 $s_{N+1}\neq s_N$이 가능함을 보였다. 이 결과는 보존적인 비선형 microscopic dynamics가 nontrivial cycle map을 만들 수 있다는 보조 증거로 남기지만, **전단 slip은 이 프로젝트의 주 물리메커니즘이 아니다.**

다음 mainline simulation은 normal spacing, normal-mode coupling, spacing correlation, free-surface opening 또는 다른 microscopically defined normal-deformation mechanism에서 secular evolution을 직접 만들어야 한다.

관련 문서:

- `docs/MILESTONE2_NORMAL_DEFORMATION.md` — 현재 메인 이론
- `docs/OPEN_PROBLEMS.md` — 현재 미해결 역학문제
- `docs/MILESTONE2_HAMILTONIAN_SLIP.md` — 보조 shear/slip 원리증명

## 균열개시

주 crack-initiation route는 수직 opening의 안정성 / first-passage 문제다.

이상화된 reduced-lattice baseline은

$$
U''(a_c)=0
$$

이지만, 최종 initiation theory는 순간 tail occupancy

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

와 cumulative first-passage probability를 구분해야 한다.

## 변수사전

이론 및 simulation 기호는 `docs/VARIABLE_DEFINITIONS.md`, firmware field와 fault flag는 `firmware/VARIABLE_DEFINITIONS.md`에 정의한다.

새로운 연구기호를 도입하는 source/document는 같은 commit에서 해당 변수정의 파일을 갱신한다.

## 현재 원리증명 simulation 재현

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest tests.test_rubin_chain
python -m unittest tests.test_hamiltonian_slip_bath
```

수치데이터는 `results/data/`, 그래프는 `results/figures/`에 생성된다.

현재 simulation은 mechanism test이며 정량적인 Al fatigue-life prediction은 아니다.

## 피로시험기 firmware 상태

`firmware/` 아래에는 수직 축방향 피로하중을 위한 target-independent C99 real-time controller core가 있다.

현재 core 수준에서 구현된 기능은 다음과 같다.

- sine / triangle normal-stress reference 생성;
- stress-to-force 변환;
- cycle counting;
- PI load-cell force-loop 구조 및 anti-windup;
- force, displacement, sensor-validity, E-stop fault 처리;
- fault 시 command 0;
- target-cycle stop;
- MCU hardware-abstraction boundary.

최종 MCU, actuator drive, sensor electronics, 검증된 control gain이 확정되지 않았기 때문에 아직 board-complete flash image는 아니다.

## Repository 구조

- `docs/` — 이론노트, 가정, 정확한 유도, 실패한 접근, 수직변형 milestone, architecture, 변수정의
- `theory/` — 해석 및 수치모델 코드
- `simulations/` — 수치실험과 재현 가능한 결과 생성
- `tests/` — 보존법칙, limiting-case, 반증 test
- `firmware/` — hardware-independent axial-fatigue real-time controller core
- `tools/` — PC telemetry / analysis helper
- `results/data/` — machine-readable simulation 결과
- `results/figures/` — 생성된 그래프
- `results/reports/` — simulation 결과의 영문/한국어 해석

## 연구 규칙

모든 중요한 결과는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

기존 피로곡선을 fitting으로만 재현하는 모델은 성공적인 이론유도로 보지 않는다.
