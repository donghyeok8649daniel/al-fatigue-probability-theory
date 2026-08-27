# Fatigue Tester Firmware Architecture

## Design principle

The real-time machine controller and the research theory solver are separated.

The MCU is responsible for deterministic control and safety:

$$
\sigma_{\rm ref}(t)
\rightarrow
F_{\rm ref}(t)
\rightarrow
\text{closed-loop force control}.
$$

The PC is responsible for logging, visualization, and computationally expensive theory analysis.

The current proof-of-principle $P(a,s,t)$ / Hamiltonian solvers are **not placed inside the hard real-time control loop**.

## Real-time reference generation

For a sine test,

$$
\sigma_{\rm ref}(t)
=
\sigma_m+\sigma_a\sin(2\pi f t),
$$

and

$$
F_{\rm ref}(t)=A_{\rm specimen}\sigma_{\rm ref}(t).
$$

The firmware tracks $F_{\rm ref}$ using measured load-cell force.

## Safety model

The firmware core immediately commands zero actuator output if any of the following occurs:

- emergency stop / travel-limit input;
- invalid sensor sample;
- measured force exceeds the configured machine limit;
- displacement exceeds the configured travel limit;
- the requested force itself exceeds the configured machine limit;
- requested target cycle count is reached.

These software checks do not replace independent hardwired safety relays, drive limits, mechanical stops, or emergency-stop circuitry.

## Controller status

The included PI controller is a framework only. No real actuator gains are supplied because $K_p$ and $K_i$ depend on the actual actuator, power amplifier, fixture stiffness, load cell, sampling rate, and specimen.

Setting gains by copying arbitrary values would violate the mechanics-first / no-hidden-fitting philosophy and can be unsafe.

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

MCU는 deterministic control과 안전을 담당한다.

$$
\sigma_{\rm ref}(t)
\rightarrow
F_{\rm ref}(t)
\rightarrow
\text{폐루프 force control}
$$

PC는 데이터 저장, 시각화, 계산량이 큰 이론분석을 담당한다.

현재 원리증명 단계의 $P(a,s,t)$ / Hamiltonian solver를 hard real-time 제어 loop 안에 직접 넣지 않는다.

## 실시간 기준값 생성

사인 피로시험이면

$$
\sigma_{\rm ref}(t)
=
\sigma_m+\sigma_a\sin(2\pi f t)
$$

이고,

$$
F_{\rm ref}(t)
=
A_{\rm specimen}\sigma_{\rm ref}(t)
$$

이다.

펌웨어는 load cell에서 측정한 힘을 이용해 $F_{\rm ref}$를 추종한다.

## 안전 구조

다음 조건에서는 펌웨어 core가 actuator command를 즉시 0으로 만든다.

- emergency stop 또는 travel-limit 입력;
- 잘못된 sensor sample;
- 측정힘이 설정된 장비 한계를 초과;
- 변위가 설정된 travel 한계를 초과;
- 명령한 force reference 자체가 장비 한계를 초과;
- 목표 cycle 수 도달.

이 software check는 독립적인 hardwired safety relay, drive limit, mechanical stop, emergency-stop 회로를 대체하지 않는다.

## 제어기 상태

포함된 PI force controller는 framework다. 실제 $K_p,K_i$는 넣지 않았다. 실제 gain은 actuator, power amplifier, fixture stiffness, load cell, sampling rate, specimen에 따라 달라진다.

임의의 gain을 복사해 넣는 것은 이 연구의 no-hidden-fitting 원칙에도 맞지 않고 장비 안전 측면에서도 부적절하다.

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
