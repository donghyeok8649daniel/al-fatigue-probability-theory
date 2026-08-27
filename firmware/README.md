# Firmware

This directory contains the hardware-independent real-time core intended for later upload to the fatigue tester MCU.

It is **not a board-complete firmware image yet**, because the MCU, actuator drive, ADC/load-cell interface, displacement sensor, and DCPD electronics have not been fixed.

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

## What must be added for the real tester

- MCU-specific timer/RTOS scheduling;
- ADC/load-cell conversion and calibration;
- actuator DAC/PWM/CAN/EtherCAT output;
- independent hardware E-stop path;
- displacement/strain acquisition;
- temperature acquisition;
- DCPD acquisition;
- persistent logging / serial or USB telemetry;
- identified and validated controller gains.

Never enable a real actuator with placeholder gains or unverified sensor polarity.

---

# 한국어 번역 — 펌웨어

이 폴더에는 향후 피로시험기 MCU에 올릴 hardware-independent 실시간 core가 들어 있다.

아직 MCU, actuator drive, ADC-load-cell interface, displacement sensor, DCPD electronics가 확정되지 않았기 때문에 **특정 보드에 바로 flash하는 완성 binary 단계는 아니다.**

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

## 실제 시험기에 추가해야 하는 부분

- MCU-specific timer 또는 RTOS scheduling;
- ADC/load cell 변환 및 calibration;
- actuator DAC/PWM/CAN/EtherCAT 출력;
- 독립적인 hardware E-stop;
- displacement/strain acquisition;
- temperature acquisition;
- DCPD acquisition;
- persistent logging 및 serial/USB telemetry;
- 실제 plant에서 식별하고 검증한 controller gain.

placeholder gain이나 sensor polarity가 검증되지 않은 상태에서는 실제 actuator를 enable하면 안 된다.
