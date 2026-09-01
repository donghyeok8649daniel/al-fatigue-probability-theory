# Fatigue Tester Firmware Architecture

## Design principle

The real-time machine controller and the research theory solver are separated.

The MCU is responsible for deterministic normal-load control, synchronized sensing, safety, and telemetry:

$$
\sigma_{\rm ref}(t)
\rightarrow
F_{\rm ref}(t)
\rightarrow
\text{closed-loop force control}.
$$

The PC is responsible for logging, visualization, and computationally expensive analysis of the active probability theory $P(a,s,t)$.

The full generalized-LJ lattice-sum solver, Smoluchowski evolution, plastic unwrapped-slip bookkeeping, and the four governing-equation outputs are **not placed inside the hard real-time control loop**. The MCU controls measured force and streams synchronized measurements; the PC may solve the probability model online or replay the logged history offline.

The detailed hardware stack and sizing rules are defined in `docs/FATIGUE_TESTER_HARDWARE.md`.

## Hardware boundary

The real tester requires the following minimum hardware blocks:

- deterministic 32-bit MCU with hardware timer and watchdog;
- dynamic axial actuator and drive sized from force/stroke/velocity/frequency requirements;
- axial load cell;
- bridge excitation, instrumentation AFE, and sufficiently fast high-resolution ADC;
- displacement/travel sensing and preferably an independent specimen-strain/extensometer channel;
- DCPD constant-current source, Kelvin voltage taps, low-noise differential AFE, and ADC;
- specimen-near temperature sensing;
- hardwired E-stop, travel limits, drive-enable interruption, and power isolation;
- timestamped PC communication.

Machine force sizing starts from

$$
F_{\rm req}
=
A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|),
$$

and sinusoidal actuator velocity is estimated from

$$
x_a\approx L_g\epsilon_a,
\qquad
v_{\rm pk}\approx 2\pi f x_a.
$$

A component is not acceptable merely because its static force rating is large enough; dynamic duty cycle, stroke, velocity, bandwidth, thermal limits, and fixture stiffness must also be checked.

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

- MCU timestamp and packet sequence number;
- cycle count and phase;
- reference normal stress and force;
- measured normal force;
- displacement;
- normal strain when available;
- temperature;
- DCPD voltage when enabled;
- actuator command;
- fault flags.

A useful minimum record is therefore

$$
\boxed{
\{t,N,\phi,\sigma_{\rm ref},F_{\rm ref},F,u,\epsilon,T,V_{\rm DCPD},u_{\rm act},\text{faults}\}
}
$$

These measurements are experimental inputs/observables. They are not substitutes for the microscopic state $P(a,s,t)$.

The stable telemetry contract allows the probability solver to be imported or replaced without modifying the MCU real-time controller.

## Probability-theory boundary

The intended PC-side chain is

$$
\boxed{
\text{MCU telemetry}
\rightarrow
P(a,s,t)
\rightarrow
\{\bar a,\bar U,E_{\rm hyst},S\}
\rightarrow
\varepsilon_p,P_{\rm crack}
}
$$

where the four governing-equation outputs are:

1. mean interlayer spacing $\bar a$;
2. mean intrinsic lattice energy $\bar U$;
3. cumulative irreversible/hysteretic energy $E_{\rm hyst}$;
4. normalization/survival $S$.

The real-time control loop must remain operational and safe even if the PC solver is disconnected or crashes.

## Safety model

The firmware core commands zero actuator output if any of the following occurs:

- emergency stop or travel-limit input;
- invalid sensor sample;
- measured force exceeds the configured machine limit;
- displacement exceeds the configured travel limit;
- the requested force itself exceeds the configured machine limit;
- requested target cycle count is reached.

Software checks do not replace independent hardwired safety relays, drive limits, mechanical stops, or emergency-stop circuitry.

The E-stop/drive-disable path must remove actuator authority without waiting for a PC response.

## Controller status

The included PI controller is a framework only. Real actuator gains are not supplied because $K_p$ and $K_i$ depend on the actuator, power amplifier, fixture stiffness, load cell, ADC/filter delay, sampling rate, and specimen.

Arbitrary controller gains must not be copied into a real machine. They have to be identified and validated on the actual hardware under conservative force/travel limits.

The 10 kHz period shown in the example code is an integration example, not a frozen hardware requirement. The final control/sample rate must be chosen after identifying the mechanical/control bandwidth and ADC/filter delay.

## HAL mapping

- `ft_hal_read_sample()` — load-cell ADC, displacement/strain, temperature, DCPD, E-stop, and limit inputs;
- `ft_hal_write_actuator()` — selected drive command interface;
- `ft_hal_send_telemetry()` — timestamped PC transport;
- `ft_hal_watchdog_kick()` — independent MCU watchdog;
- `ft_hal_safe_shutdown()` — drive-enable removal / safe-state output.

## Files

- `docs/FATIGUE_TESTER_HARDWARE.md` — hardware sizing, acquisition chains, safety, telemetry, and procurement order
- `firmware/include/fatigue_controller.h` — target-independent controller API
- `firmware/src/fatigue_controller.c` — waveform, cycle counting, PI force control, anti-windup, safety state
- `firmware/include/fatigue_hal.h` — MCU-specific hardware abstraction boundary
- `firmware/src/main_loop_example.c` — fixed-period control-task integration example
- `firmware/src/host_test.c` — host-side logic and safety tests
- `firmware/CMakeLists.txt` — host compilation test
- `tools/fatigue_pc_bridge.py` — initial PC-side telemetry/log helper

## Porting sequence

1. freeze specimen geometry and maximum stress envelope;
2. calculate required force, displacement, velocity, frequency, and duty cycle;
3. choose the actuator and drive;
4. choose the load cell and force AFE/ADC;
5. choose displacement/extensometer sensing;
6. choose the MCU and communication interface;
7. design/select DCPD and temperature acquisition;
8. implement `fatigue_hal.h`;
9. calibrate sensor units independently;
10. identify the closed-loop plant at low force;
11. tune and validate $K_p,K_i$ with conservative machine limits;
12. verify zero-output behavior for every fault;
13. run a dummy specimen before high-purity Al;
14. only then enable full test amplitudes.

---

# 한국어 번역 — 피로시험기 펌웨어 구조

## 설계 원칙

실시간 장비 제어기와 연구 이론 solver를 분리한다.

MCU는 deterministic normal-load control, 동기화된 센싱, 안전, telemetry를 담당한다.

$$
\sigma_{\rm ref}(t)
\rightarrow
F_{\rm ref}(t)
\rightarrow
\text{폐루프 force control}
$$

PC는 데이터 저장, 시각화, 그리고 $P(a,s,t)$ 확률이론 계산을 담당한다.

full generalized-LJ lattice sum, Smoluchowski evolution, unwrapped slip에 의한 소성 bookkeeping, 4방정식 계산을 hard real-time MCU loop 안에 넣지 않는다. MCU는 측정 force를 제어하고 동기화된 실험 데이터를 PC로 보내며, PC는 이를 실시간 또는 사후 replay 방식으로 확률 solver에 넣는다.

구체적인 하드웨어 stack과 sizing 규칙은 `docs/FATIGUE_TESTER_HARDWARE.md`를 기준으로 한다.

## 필요한 하드웨어 경계

실제 시험기에는 최소한 다음이 필요하다.

- deterministic hardware timer/watchdog을 갖는 32-bit MCU;
- force/stroke/velocity/frequency 요구조건으로 선정한 동적 축방향 actuator + drive;
- axial load cell;
- bridge excitation + instrumentation AFE + 충분히 빠른 고해상도 ADC;
- displacement/travel sensor와 가능하면 독립 specimen strain/extensometer channel;
- DCPD constant-current source + Kelvin voltage tap + 저잡음 differential AFE/ADC;
- specimen-near temperature sensor;
- hardwired E-stop, travel limit, drive-enable 차단 및 power isolation;
- timestamp가 포함된 PC 통신.

machine force 요구량은

$$
F_{\rm req}
=
A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|)
$$

에서 시작해 정하고, 사인 하중의 actuator 속도는

$$
x_a\approx L_g\epsilon_a,
\qquad
v_{\rm pk}\approx2\pi f x_a
$$

로 추정한다.

정적 force rating만 맞는다고 actuator를 선택하면 안 되며 stroke, velocity, bandwidth, thermal duty, fixture stiffness까지 확인해야 한다.

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

telemetry에는 최소한 다음을 포함한다.

- MCU timestamp / packet sequence;
- cycle count / phase;
- reference stress / force;
- measured force;
- displacement;
- normal strain;
- temperature;
- DCPD voltage;
- actuator command;
- fault flags.

즉 최소 telemetry record는

$$
\boxed{
\{t,N,\phi,\sigma_{\rm ref},F_{\rm ref},F,u,\epsilon,T,V_{\rm DCPD},u_{\rm act},\text{faults}\}
}
$$

로 잡는다.

이 측정량들은 experimental input/observable이며 microscopic state $P(a,s,t)$ 자체를 대신하지 않는다.

telemetry contract를 안정적으로 유지하면 PC 확률 solver는 MCU firmware를 바꾸지 않고 가져오거나 교체할 수 있다.

## 확률이론 경계

PC에서는

$$
\boxed{
\text{MCU telemetry}
\rightarrow
P(a,s,t)
\rightarrow
\{\bar a,\bar U,E_{\rm hyst},S\}
\rightarrow
\varepsilon_p,P_{\rm crack}
}
$$

를 계산한다.

여기서 4방정식 출력은

1. 평균 층간거리 $\bar a$;
2. 평균 intrinsic lattice energy $\bar U$;
3. 누적 비가역/hysteresis 에너지 $E_{\rm hyst}$;
4. 정규화/생존확률 $S$.

PC solver가 끊기거나 오류가 나도 MCU의 실시간 제어와 안전은 독립적으로 유지되어야 한다.

## 안전 구조

다음 조건에서는 firmware core가 actuator command를 0으로 만든다.

- emergency stop 또는 travel-limit 입력;
- 잘못된 sensor sample;
- 측정힘이 설정된 장비 한계를 초과;
- 변위가 설정된 travel 한계를 초과;
- 명령한 force reference 자체가 장비 한계를 초과;
- 목표 cycle 수 도달.

software check는 독립적인 hardwired safety relay, drive limit, mechanical stop, emergency-stop 회로를 대체하지 않는다. E-stop/drive-disable은 PC 응답을 기다리지 않고 actuator authority를 제거할 수 있어야 한다.

## 제어기 상태

포함된 PI force controller는 framework다. 실제 $K_p,K_i$는 actuator, power amplifier, fixture stiffness, load cell, ADC/filter delay, sampling rate, specimen에 따라 달라진다.

예제 코드의 10 kHz는 integration example이지 확정된 hardware requirement가 아니다. 최종 control/sample rate는 실제 plant bandwidth와 ADC/filter delay를 식별한 뒤 정한다.

## HAL mapping

- `ft_hal_read_sample()` — load cell, displacement/strain, temperature, DCPD, E-stop/limit 입력
- `ft_hal_write_actuator()` — 선택한 actuator drive 명령
- `ft_hal_send_telemetry()` — timestamped PC transport
- `ft_hal_watchdog_kick()` — MCU independent watchdog
- `ft_hal_safe_shutdown()` — drive-enable 제거 / safe state

## 실제 장비 포팅 순서

1. 시편 geometry와 최대 stress envelope를 확정한다.
2. force/displacement/velocity/frequency/duty requirement를 계산한다.
3. actuator + drive를 선택한다.
4. load cell + force AFE/ADC를 선택한다.
5. displacement/extensometer sensor를 선택한다.
6. MCU와 communication interface를 선택한다.
7. DCPD/temperature acquisition을 설계한다.
8. `fatigue_hal.h`를 구현한다.
9. sensor calibration을 수행한다.
10. 낮은 force에서 plant를 식별한다.
11. 보수적인 limit 아래 $K_p,K_i$를 tuning한다.
12. 모든 fault의 zero-output을 검증한다.
13. dummy specimen으로 시험한다.
14. 그 다음 full amplitude를 허용한다.
