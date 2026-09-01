# Firmware

This directory contains the hardware-independent real-time core intended for later upload to the fatigue tester MCU.

It is **not a board-complete firmware image yet**. The real machine still needs the selected servo drive, load-cell AFE/ADC, displacement/strain sensor, DCPD electronics, board-specific HAL and validated controller gains.

Hardware and procurement information now lives under:

- `../hardware/README.md`
- `../hardware/bom/fatigue_tester_bom.csv`

The architecture is intentionally split as follows:

- **MCU:** deterministic cyclic loading, force feedback, sensing, hard real-time safety, timestamped telemetry;
- **PC:** logging/visualization and the imported probability-fatigue solver `P(a,s,t)` with the four governing-equation outputs and crack-initiation inference.

The generalized-LJ/probability solver is therefore not required inside the hard real-time MCU loop.

## Host verification

```bash
cd fatigue_tester/firmware
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

- STM32-class deterministic MCU;
- dynamic axial servo actuator + drive;
- axial load cell;
- load-cell instrumentation amplifier/AFE;
- displacement/travel and specimen-strain channels;
- DCPD current source + Kelvin voltage acquisition;
- temperature channel;
- hardwired E-stop, travel limits, breaker/contactor and drive disable;
- timestamped PC telemetry link.

Actuator/load-cell sizing starts from

\[
F_{\rm req}=A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|),
\]

not from an arbitrary fixed machine force rating.

## What must be added for the real tester

- MCU-specific timer/RTOS scheduling;
- board-specific HAL;
- sensor calibration;
- actuator command adapter for the selected drive;
- independent hardware E-stop path;
- persistent telemetry/logging;
- identified and validated controller gains.

Never enable a real actuator with placeholder gains or unverified sensor polarity.

---

# 한국어

이 폴더에는 피로시험기에 올릴 hardware-independent 실시간 제어 core가 들어 있다.

확률 시뮬레이션 `P(a,s,t)`와 4방정식 계산은 MCU에 넣지 않고 PC에서 처리한다. MCU는 반복 인장하중 생성, force feedback, 센싱, 안전, timestamp telemetry에 집중한다.

부품/BOM은 `../hardware/bom/fatigue_tester_bom.csv`를 기준으로 관리한다.
