# Firmware

This directory contains the hardware-independent real-time core intended for later upload to the fatigue tester MCU.

It is **not a board-complete firmware image yet**, because the MCU, actuator drive, ADC/load-cell interface, displacement sensor, and DCPD electronics have not been fixed.

The concrete hardware requirements, sizing equations, sensor chain, DCPD chain, safety hardware, and HAL mapping are documented in:

- `docs/FATIGUE_TESTER_HARDWARE.md`

The architecture is intentionally split as follows:

- **MCU:** deterministic cyclic loading, force feedback, sensing, hard real-time safety, timestamped telemetry;
- **PC:** logging/visualization and the imported probability-fatigue solver `P(a,s,t)` with the four governing-equation outputs and crack-initiation inference.

The generalized-LJ/probability solver is therefore not required inside the hard real-time MCU loop.

## Host verification

```bash
cd firmware
cmake -S . -B build
cmake --build build
./build/firmware_host_test
```

Expected output:

```text
firmware core host tests: PASS
```

## What is already implemented

- sine and triangle stress-reference generation;
- stress-to-force conversion;
- cycle counting;
- PI load-cell force-loop core;
- conditional-integration anti-windup;
- force, displacement, sensor-validity and E-stop fault handling;
- zero actuator command on fault;
- target-cycle stop;
- hardware abstraction interface.

## Required hardware blocks

The board-specific implementation shall provide at least:

- a deterministic 32-bit MCU with watchdog and suitable timer/SPI/communication peripherals;
- a dynamic axial actuator and drive sized from the specimen force/stroke/frequency requirements;
- an axial load cell;
- bridge excitation + instrumentation AFE + sufficiently fast high-resolution ADC for force feedback;
- displacement/travel sensing and preferably an independent specimen-strain/extensometer channel;
- DCPD constant-current source, Kelvin voltage taps, low-noise differential AFE and ADC;
- specimen-near temperature sensing;
- hardwired E-stop, travel limits, drive-enable interruption and power isolation;
- a timestamped PC telemetry link.

Actuator/load-cell sizing starts from

$$
F_{\rm req}=A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|),
$$

not from an arbitrary fixed machine force rating.

See `docs/FATIGUE_TESTER_HARDWARE.md` before selecting specific parts.

## What must be added for the real tester

- select the actual MCU and actuator/drive after machine bandwidth sizing;
- MCU-specific timer/RTOS scheduling;
- load-cell AFE/ADC implementation and calibration;
- actuator DAC/PWM/CAN/EtherCAT output adapter as required by the chosen drive;
- independent hardware E-stop path;
- displacement/strain acquisition;
- temperature acquisition;
- DCPD current source and acquisition electronics;
- persistent logging / serial, USB, CAN or Ethernet telemetry;
- identified and validated controller gains.

Never enable a real actuator with placeholder gains or unverified sensor polarity.

---

# 한국어 번역 — 펌웨어

이 폴더에는 향후 피로시험기 MCU에 올릴 hardware-independent 실시간 core가 들어 있다.

아직 MCU, actuator drive, ADC-load-cell interface, displacement sensor, DCPD electronics가 확정되지 않았기 때문에 **특정 보드에 바로 flash하는 완성 binary 단계는 아니다.**

구체적인 하드웨어 요구조건, sizing 식, sensor chain, DCPD chain, 안전회로 및 HAL mapping은 다음 문서에 정리했다.

- `docs/FATIGUE_TESTER_HARDWARE.md`

역할은 다음처럼 분리한다.

- **MCU:** 반복하중 생성, force feedback, 센싱, hard real-time 안전, timestamp telemetry;
- **PC:** 데이터 저장/시각화와 가져온 `P(a,s,t)` 확률 피로 solver, 4방정식 출력 및 균열개시 추론.

따라서 generalized-LJ/확률 solver를 MCU hard real-time loop 안에 넣을 필요는 없다.

## PC에서 logic 검증

```bash
cd firmware
cmake -S . -B build
cmake --build build
./build/firmware_host_test
```

정상 결과:

```text
firmware core host tests: PASS
```

## 현재 구현된 기능

- sine/triangle stress reference 생성;
- stress-to-force 변환;
- cycle counter;
- load-cell 기반 PI force control core;
- conditional-integration anti-windup;
- force/displacement/sensor/E-stop fault 처리;
- fault 발생 시 actuator command 0;
- 목표 cycle 도달 시 정지;
- hardware abstraction interface.

## 필요한 하드웨어 블록

실제 장비에는 최소한 다음이 필요하다.

- deterministic timer와 watchdog을 갖는 32-bit MCU;
- 시편의 힘/변위/주파수 요구조건으로 선정한 동적 축방향 actuator + drive;
- axial load cell;
- load-cell bridge excitation + instrumentation AFE + 충분히 빠른 고해상도 ADC;
- displacement/travel sensor와 가능하면 독립 specimen strain/extensometer channel;
- DCPD constant-current source + Kelvin voltage tap + 저잡음 differential AFE/ADC;
- 시편 근처 temperature sensor;
- hardwired E-stop, travel limit, drive-enable 차단, power isolation;
- timestamp가 포함되는 PC telemetry interface.

actuator/load cell 용량은 임의로 고정하지 않고

$$
F_{\rm req}=A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|)
$$

에서 시작해 정한다.

특정 부품을 고르기 전에 `docs/FATIGUE_TESTER_HARDWARE.md`를 기준으로 삼는다.

## 실제 시험기에 추가해야 하는 부분

- machine bandwidth 계산 후 MCU와 actuator/drive 확정;
- MCU-specific timer 또는 RTOS scheduling;
- load cell AFE/ADC 및 calibration;
- 선택한 drive에 맞는 DAC/PWM/CAN/EtherCAT 출력;
- 독립적인 hardware E-stop;
- displacement/strain acquisition;
- temperature acquisition;
- DCPD current source 및 acquisition electronics;
- persistent logging 및 serial/USB/CAN/Ethernet telemetry;
- 실제 plant에서 식별하고 검증한 controller gain.

placeholder gain이나 sensor polarity가 검증되지 않은 상태에서는 실제 actuator를 enable하면 안 된다.
