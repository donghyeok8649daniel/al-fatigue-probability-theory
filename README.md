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

The current minimal extension for non-affine structure is the joint state

$$
P(a,s,t),
$$

where $s$ is an atomistically traceable slip/disregistry coordinate and

$$
P(a,t)=\int P(a,s,t)\,ds.
$$

## Current backbone

- Microscopic mechanics: Newton/Hamiltonian dynamics
- Deterministic lattice baseline: generalized Lennard-Jones interaction
- Distribution energy: pair-distance distributions and correlation hierarchy
- Exact density conservation: $\partial_tP+\partial_a(Pv)=0$
- Hysteresis target: obtain $\oint \sigma\,d\epsilon>0$ without inserting an empirical hysteresis law
- Fatigue target: obtain secular cycle-to-cycle evolution $P_{N+1}\neq P_N$
- Crack initiation: formulate as loss of mechanical stability / first-passage into an unstable state

## Variable dictionaries

Theory and simulation symbols are defined in `docs/VARIABLE_DEFINITIONS.md`. Firmware fields and fault flags are defined in `firmware/VARIABLE_DEFINITIONS.md`.

Any future document or source change that introduces a new research symbol should update the relevant variable dictionary in the same commit.

## Current progress

### Milestone 1 — mechanics-derived reduced hysteresis

A resolved structural coordinate coupled to a semi-infinite harmonic Rubin chain produces an exact phase lag and nonzero reduced hysteresis while the full microscopic model remains conservative.

Reference nondimensional result:

$$
A_H^{\mathrm{analytic}}=0.0152091700,
$$

$$
A_H^{\mathrm{numeric}}=0.0152088400.
$$

The relative loop-area error is about $2.17\times10^{-5}$ and the relative external-work / internal-energy balance error is about $1.25\times10^{-5}$.

### Milestone 2 — proof of cycle-to-cycle structural evolution

A nonlinear periodic non-affine coordinate coupled to the conservative lattice bath produces three numerical regimes: bounded intra-basin response, finite relocation followed by periodic response, and a running inter-basin state.

For the current nondimensional $F_a=0.50$ proof-of-principle case,

$$
s_{N+1}-s_N\approx-1
$$

in the late-cycle regime, while the global energy-balance relative error is about

$$
1.77\times10^{-7}.
$$

This is **not yet a calibrated Al fatigue-life prediction**. Real Al still requires an atomistic $\gamma(\mathbf s)$ input and a mechanically derived explanation for the gap between ordinary fatigue stresses and perfect-crystal ideal shear strength.

See `results/reports/SIMULATION_RESULTS.md` for figures and interpretation.

## Reproduce the current simulation results

From the repository root:

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest tests.test_rubin_chain
python -m unittest tests.test_hamiltonian_slip_bath
```

Generated numerical data are written under `results/data/` and figures under `results/figures/`.

## Fatigue-tester firmware status

A target-independent C99 real-time controller core now exists under `firmware/`.

Implemented at the core level:

- sine / triangle stress-reference generation;
- stress-to-force conversion;
- cycle counting;
- PI load-cell force-loop structure and anti-windup;
- force, displacement, sensor-validity and E-stop fault handling;
- zero command on fault;
- target-cycle stop;
- MCU hardware-abstraction boundary.

The core compiled locally with strict warnings and its host logic/safety test passed. It is **not yet a board-complete flash image** because the actual MCU, actuator drive, ADC/load-cell chain, displacement sensor, DCPD hardware, and real controller gains still need to be fixed and validated.

See `docs/FIRMWARE_ARCHITECTURE.md` and `firmware/README.md`.

## Repository structure

- `docs/` — theory notes, assumptions, derivations, failed approaches, architecture, and variable definitions
- `theory/` — analytic and numerical model code
- `simulations/` — numerical experiments and reproducible result generation
- `tests/` — conservation, limiting-case, and falsification tests
- `firmware/` — hardware-independent C real-time controller core and MCU porting boundary
- `tools/` — PC-side telemetry / analysis helpers
- `results/data/` — machine-readable simulation results
- `results/figures/` — generated plots
- `results/reports/` — bilingual interpretation of simulation results

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

가장 기본적인 축약상태는 열역학적 극한의 원자간격 밀도다.

$$
P(a,t)=\lim_{N\to\infty}\frac{1}{N}\sum_{i=1}^{N}\delta(a-a_i(t)).
$$

여기서 $a_i(t)$는 국부 원자간격 descriptor다. 특정 확률분포 family를 미리 가정하지 않고 $P(a,t)$의 시간진화를 역학에서 유도하는 것이 장기 목표다.

현재 비아핀 구조를 포함하기 위한 최소 확장 후보는

$$
P(a,s,t)
$$

이며, $s$는 원자좌표로 추적 가능한 slip/disregistry coordinate다. 기존 spacing state는

$$
P(a,t)=\int P(a,s,t)\,ds
$$

라는 정확한 marginal로 남는다.

## 현재 이론의 기본 골격

- 미시역학: Newton/Hamiltonian dynamics
- 결정론적 lattice baseline: generalized Lennard-Jones interaction
- 분포에너지: pair-distance distribution 및 correlation hierarchy
- 정확한 밀도보존식: $\partial_tP+\partial_a(Pv)=0$
- 히스테리시스 목표: 경험적 hysteresis law 없이 $\oint \sigma\,d\epsilon>0$
- 피로누적 목표: $P_{N+1}\neq P_N$인 secular evolution
- 균열개시: mechanical stability loss 또는 first-passage 문제

## 변수 사전

이론 및 simulation 기호는 `docs/VARIABLE_DEFINITIONS.md`에 정의하고, firmware field와 fault flag는 `firmware/VARIABLE_DEFINITIONS.md`에 정의한다.

앞으로 새로운 연구기호를 도입하는 문서나 source change는 같은 commit에서 해당 변수정의 파일을 함께 갱신한다.

## 현재 진행상황

### Milestone 1 — 역학으로부터 유도된 축약 히스테리시스

관심 구조좌표를 준무한 harmonic Rubin chain에 결합하면 전체 미시계는 보존계인 상태에서 정확한 phase lag와 0이 아닌 축약 히스테리시스가 발생한다.

기준 무차원 결과는

$$
A_H^{\mathrm{analytic}}=0.0152091700,
$$

$$
A_H^{\mathrm{numeric}}=0.0152088400
$$

이다.

loop area 해석값-수치값 상대오차는 약 $2.17\times10^{-5}$이고, 외부 일 / 내부에너지 수지 상대오차는 약 $1.25\times10^{-5}$이다.

### Milestone 2 — cycle-to-cycle 구조변화 원리증명

비선형 주기적 비아핀 좌표를 보존적 lattice bath에 결합한 모델에서 세 영역이 나타났다. 하나의 basin 안에서 제한된 응답, 유한한 relocation 후 주기응답, 그리고 basin 사이를 계속 이동하는 running state다.

현재 무차원 $F_a=0.50$ 원리증명 계산의 후반부에서는

$$
s_{N+1}-s_N\approx-1
$$

이고, 전체 energy-balance 상대오차는 약

$$
1.77\times10^{-7}
$$

이다.

하지만 이것은 **아직 실제 Al 피로수명을 보정한 예측모델이 아니다.** 실제 Al의 atomistic $\gamma(\mathbf s)$ 입력과, 일반적인 피로응력과 완전결정 ideal shear strength 사이의 간극을 줄이는 미시역학적 메커니즘이 아직 필요하다.

그래프와 해석은 `results/reports/SIMULATION_RESULTS.md`에 정리한다.

## 현재 simulation 결과 재생성

repository root에서 다음을 실행한다.

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest tests.test_rubin_chain
python -m unittest tests.test_hamiltonian_slip_bath
```

수치데이터는 `results/data/`, 그래프는 `results/figures/`에 생성된다.

## 피로시험기 firmware 상태

`firmware/` 아래에 target-independent C99 실시간 controller core를 추가했다.

현재 core 수준에서 구현된 기능은 다음과 같다.

- sine / triangle stress reference 생성;
- stress-to-force 변환;
- cycle counting;
- PI load-cell force-loop 구조와 anti-windup;
- force, displacement, sensor-validity, E-stop fault 처리;
- fault 시 command 0;
- target-cycle stop;
- MCU hardware-abstraction boundary.

이 core는 strict warning 옵션으로 로컬 compile했고 host logic/safety test를 통과했다. 다만 실제 MCU, actuator drive, ADC/load-cell chain, displacement sensor, DCPD hardware, 실제 controller gain이 확정되지 않았으므로 **아직 특정 보드에 바로 flash하는 완성 image는 아니다.**

자세한 구조는 `docs/FIRMWARE_ARCHITECTURE.md`와 `firmware/README.md`에 정리한다.

## Repository 구조

- `docs/` — 이론노트, 가정, 유도, 실패한 접근, architecture, 변수정의
- `theory/` — 해석 및 수치모델 코드
- `simulations/` — 수치실험 및 재현 가능한 결과 생성
- `tests/` — 보존법칙, limiting case, 반증 test
- `firmware/` — hardware-independent C 실시간 controller core와 MCU 포팅 경계
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

기존 피로곡선을 fitting으로 재현했을 뿐인 모델은 성공적인 이론 유도로 간주하지 않는다.
