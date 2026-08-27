# Fatigue Tester Firmware Architecture

## Design principle

The real-time machine controller and the research theory solver are separated.

The MCU is responsible for deterministic normal-load control and safety:

$$
\sigma_{\rm ref}(t)
\rightarrow
F_{\rm ref}(t)
\rightarrow
\text{closed-loop force control}.
$$

The PC is responsible for logging, visualization, and computationally expensive analysis of the active normal-spacing theory $P(a,t)$.

The full atomistic/generalized-LJ solver is **not placed inside the hard real-time control loop**. The MCU controls measured force; the PC later compares the measured normal stress/strain history with the theory.

## Real-time reference generation

For a sine test,

$$
\sigma_{\rm ref}(t)
=
\sigma_m+\sigma_a\sin(2\pi f t),
$$

and

$$
\boxed{F_{\rm ref}(t)=A_{\rm specimen}\sigma_{\rm ref}(t).}
$$

The firmware tracks $F_{\rm ref}$ using the measured load-cell force.

## Recorded quantities

The intended telemetry includes at least:

- reference normal stress and force;
- measured normal force;
- displacement;
- normal strain when available;
- temperature;
- DCPD voltage when enabled;
- cycle count;
- actuator command;
- fault flags.

These measurements are experimental inputs/observables. They are not substitutes for the microscopic state $P(a,t)$.

## Safety model

The firmware core commands zero actuator output if any of the following occurs:

- emergency stop or travel-limit input;
- invalid sensor sample;
- measured force exceeds the configured machine limit;
- displacement exceeds the configured travel limit;
- the requested force itself exceeds the configured machine limit;
- requested target cycle count is reached.

Software checks do not replace independent hardwired safety relays, drive limits, mechanical stops, or emergency-stop circuitry.

## Controller status

The included PI controller is a framework only. Real actuator gains are not supplied because $K_p$ and $K_i$ depend on the actuator, power amplifier, fixture stiffness, load cell, sampling rate, and specimen.

Arbitrary controller gains must not be copied into a real machine. They have to be identified and validated on the actual hardware under conservative force/travel limits.

## Files

- `firmware/include/fatigue_controller.h` — target-independent controller API
- `firmware/src/fatigue_controller.c` — waveform, cycle counting, PI force control, anti-windup, safety state
- `firmware/include/fatigue_hal.h` — MCU-specific hardware abstraction boundary
- `firmware/src/main_loop_example.c` — fixed-period control-task integration example
- `firmware/src/host_test.c` — host-side logic and safety tests
- `firmware/CMakeLists.txt` — host compilation test
- `tools/fatigue_pc_bridge.py` — initial PC-side telemetry/log helper

## Porting sequence

1. choose the MCU and actuator drive;
2. implement `fatigue_hal.h` for ADC/load cell, displacement, DCPD, temperature, actuator output, watchdog, and E-stop;
3. calibrate sensor units independently;
4. identify the closed-loop plant at low force;
5. tune and validate $K_p,K_i$ with conservative machine limits;
6. verify zero-output behavior for every fault;
7. run a dummy specimen before high-purity Al;
8. only then enable full test amplitudes.

---

# 한국어 번역 — 피로시험기 펌웨어 구조

## 설계 원칙

실시간 장비 제어기와 연구 이론 solver를 분리한다.

MCU는 deterministic normal-load control과 안전을 담당한다.

$$
\sigma_{\rm ref}(t)
\rightarrow
F_{\rm ref}(t)
\rightarrow
\text{폐루프 force control}
$$

PC는 데이터 저장, 시각화, 그리고 활성 normal-spacing theory $P(a,t)$의 계산량이 큰 분석을 담당한다.

full atomistic/generalized-LJ solver를 hard real-time 제어 loop 안에 직접 넣지 않는다. MCU는 측정 force를 제어하고, PC가 이후 측정된 normal stress/strain history를 이론과 비교한다.

## 실시간 기준값 생성

사인 피로시험에서는

$$
\sigma_{\rm ref}(t)
=
\sigma_m+\sigma_a\sin(2\pi f t)
$$

이고,

$$
\boxed{F_{\rm ref}(t)=A_{\rm specimen}\sigma_{\rm ref}(t)}
$$

이다.

펌웨어는 load cell에서 측정한 힘을 이용해 $F_{\rm ref}$를 추종한다.

## 기록량

telemetry에는 최소한 다음을 포함하는 것을 목표로 한다.

- reference normal stress와 force;
- measured normal force;
- displacement;
- 가능하면 normal strain;
- temperature;
- 사용하는 경우 DCPD voltage;
- cycle count;
- actuator command;
- fault flag.

이 측정량들은 experimental input/observable이며 microscopic state $P(a,t)$ 자체를 대신하지 않는다.

## 안전 구조

다음 조건에서는 firmware core가 actuator command를 0으로 만든다.

- emergency stop 또는 travel-limit 입력;
- 잘못된 sensor sample;
- 측정힘이 설정된 장비 한계를 초과;
- 변위가 설정된 travel 한계를 초과;
- 명령한 force reference 자체가 장비 한계를 초과;
- 목표 cycle 수 도달.

software check는 독립적인 hardwired safety relay, drive limit, mechanical stop, emergency-stop 회로를 대체하지 않는다.

## 제어기 상태

포함된 PI force controller는 framework다. 실제 $K_p,K_i$는 actuator, power amplifier, fixture stiffness, load cell, sampling rate, specimen에 따라 달라지므로 임의값을 제공하지 않는다.

실제 장비에 arbitrary controller gain을 복사해 넣으면 안 된다. 실제 hardware에서 보수적인 force/travel limit 아래 plant를 식별하고 검증해야 한다.

## 파일

- `firmware/include/fatigue_controller.h` — MCU 독립 controller API
- `firmware/src/fatigue_controller.c` — waveform, cycle count, PI force control, anti-windup, safety state
- `firmware/include/fatigue_hal.h` — MCU-specific hardware abstraction 경계
- `firmware/src/main_loop_example.c` — 고정주기 control task 연결 예
- `firmware/src/host_test.c` — PC에서 실행하는 logic/safety test
- `firmware/CMakeLists.txt` — host compile 검증
- `tools/fatigue_pc_bridge.py` — 초기 PC telemetry/log helper

## 실제 장비 포팅 순서

1. MCU와 actuator drive를 확정한다.
2. ADC/load cell, displacement, DCPD, temperature, actuator output, watchdog, E-stop에 대해 `fatigue_hal.h`를 구현한다.
3. sensor 단위를 독립적으로 calibration한다.
4. 낮은 force에서 closed-loop plant를 식별한다.
5. 보수적인 machine limit 아래에서 $K_p,K_i$를 tuning 및 검증한다.
6. 모든 fault에서 output이 0이 되는지 확인한다.
7. 실제 고순도 Al 전에 dummy specimen으로 시험한다.
8. 그 다음에만 full test amplitude를 허용한다.
