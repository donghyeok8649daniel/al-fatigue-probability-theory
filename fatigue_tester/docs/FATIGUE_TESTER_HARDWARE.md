# Fatigue Tester Hardware Stack

## Purpose

This document fixes the **minimum hardware architecture** for the uniaxial tensile-fatigue tester. The real-time machine controller and the probability/fatigue theory solver remain separated:

\[
\boxed{\text{MCU} = \text{loading + sensing + safety + telemetry}}
\]

\[
\boxed{\text{PC} = P(a,s,t)\text{ solver + four governing equations + crack-initiation inference}}
\]

The MCU does **not** need to solve the generalized-LJ lattice sums, Smoluchowski PDE, or the four governing equations in hard real time. It must instead provide deterministic loading and synchronized experimental data of sufficient quality for the PC-side probability solver.

---

## 1. Machine sizing equations

The actuator and force sensor shall be selected from the requested specimen and loading condition, not from a fixed arbitrary force rating.

For specimen cross-sectional area \(A\), mean stress \(\sigma_m\), and stress amplitude \(\sigma_a\),

\[
\sigma_{\max}=\sigma_m+\sigma_a,
\qquad
\sigma_{\min}=\sigma_m-\sigma_a.
\]

The required absolute force range is

\[
\boxed{F_{\rm req}=A\max(|\sigma_{\max}|,|\sigma_{\min}|)}.
\]

The load cell, actuator, fixture, and mechanical frame shall all have a rated capacity above this value with an engineering safety margin determined during the detailed mechanical design.

For a gauge length \(L_g\) and approximate cyclic axial strain amplitude \(\epsilon_a\), the corresponding displacement amplitude is

\[
\boxed{x_a\approx L_g\epsilon_a}
\]

and for sinusoidal loading the approximate peak actuator velocity requirement is

\[
\boxed{v_{\rm pk}\approx 2\pi f x_a}.
\]

A candidate actuator is acceptable only if it can simultaneously satisfy the required force, stroke, velocity, acceleration, and duty-cycle/thermal limits at the requested fatigue frequency.

---

## 2. Required hardware blocks

### 2.1 Real-time controller

Required class:

- 32-bit MCU with deterministic hardware timers;
- hardware FPU strongly preferred;
- multiple SPI buses or sufficiently fast shared SPI;
- ADC/DAC/PWM/CAN/RS-485 capability as required by the selected drive;
- hardware watchdog;
- USB or isolated serial/CAN path to the PC.

A practical baseline is an **STM32G4/H7-class controller** or an equivalent MCU. The exact part is intentionally not frozen until the actuator drive and sensor interfaces are selected.

The controller runs the hard real-time force loop and must never depend on the PC for emergency shutdown.

### 2.2 Axial actuator and drive

Required functions:

- bidirectional or tension-biased cyclic axial loading as required by the test protocol;
- force/stroke/velocity capability satisfying Section 1;
- command interface compatible with the MCU/HAL;
- independent drive enable and fault input;
- preferably an independent drive-side current/torque/force limit.

Possible implementation classes include an electromechanical servo actuator/servo motor + screw mechanism, a suitable linear servo actuator, or another dynamic actuator whose continuous cyclic rating is verified for the target frequency. The final selection must be made from force-stroke-frequency calculations rather than static load rating alone.

### 2.3 Axial load cell

Required:

- tensile/compressive force transducer or tension-only transducer compatible with the fixture;
- rated force above \(F_{\rm req}\);
- low creep and adequate fatigue life;
- bandwidth comfortably above the mechanical test frequency and force-control bandwidth;
- overload capacity suitable for machine protection.

For the present prototype scale, the intended class is approximately **1--2 kN** only when the specimen/load calculation confirms that this covers \(F_{\rm req}\) with margin. This is not a universal fixed rating.

### 2.4 Load-cell analog front end and ADC

Do **not** use a slow scale-oriented ADC as the primary closed-loop force sensor if its output rate cannot support the force-control bandwidth.

Required chain:

\[
\boxed{\text{load cell}\rightarrow\text{instrumentation amplifier/AFE}\rightarrow\text{high-resolution ADC}\rightarrow\text{MCU}}
\]

Requirements:

- differential low-noise input;
- bridge excitation appropriate to the selected load cell;
- 20--24 bit class conversion preferred for metrology;
- synchronized sample rate sufficient for the chosen control loop;
- anti-aliasing/filtering characterized rather than hidden;
- known digital-filter group delay;
- calibration coefficients stored separately from controller gains.

Initial design target: synchronized force samples in the **multi-kS/s to tens-of-kS/s class**, with the final rate determined after plant bandwidth identification.

### 2.5 Displacement / strain measurement

At least one independent displacement or strain channel is required even when force is the controlled variable.

Preferred hierarchy:

1. extensometer or gauge-length displacement sensor for specimen strain;
2. actuator/fixture displacement sensor for machine travel and safety;
3. both channels when practical.

Possible sensor classes:

- LVDT;
- optical or magnetic linear encoder;
- high-resolution linear potentiometer for early low-cost prototypes only;
- strain-gauge extensometer with its own AFE.

The safety travel channel must remain valid even if the research strain channel fails.

### 2.6 DCPD crack-initiation channel

Because crack initiation is intended to be experimentally detected by resistance/potential-drop change, reserve a dedicated DCPD chain:

\[
\boxed{\text{stable current source}\rightarrow\text{specimen}\rightarrow\text{Kelvin voltage taps}\rightarrow\text{low-noise differential AFE}\rightarrow\text{ADC}}
\]

Required hardware:

- stable constant-current source;
- four-wire/Kelvin connection to the specimen;
- low-noise differential instrumentation amplifier;
- ADC channel isolated or filtered as required by the power electronics environment;
- reference measurement for current-source drift when practical;
- shielded/twisted wiring and careful grounding.

DCPD does not need to run inside the force-control feedback path. It shall be timestamped and streamed to the PC for crack-initiation inference and experimental validation.

### 2.7 Temperature measurement

Required because both resistance-based crack sensing and material response are temperature sensitive.

Minimum:

- specimen-near temperature sensor;
- preferably a second sensor near the DCPD/current-source electronics or fixture;
- synchronized timestamp in telemetry.

Sensor class may be RTD, thermocouple, or a suitable precision semiconductor sensor depending on the final temperature range.

### 2.8 Safety hardware

Software fault flags are not sufficient. Required independent safety hardware:

- physical emergency-stop switch;
- safety relay or hardwired drive-enable interruption;
- upper/lower travel limit switches;
- mechanical travel stops where appropriate;
- drive over-current/over-torque protection;
- fuse/breaker and appropriately rated power isolation;
- watchdog that defaults to actuator disable;
- force and displacement software limits as a secondary layer.

The E-stop path must be able to remove actuator authority without waiting for the PC or probability solver.

### 2.9 Communications and logging

Required PC link:

- USB CDC, isolated UART/RS-485, CAN/CAN-FD, or Ethernet depending on the selected MCU and noise environment;
- deterministic packet sequence number;
- MCU timestamp;
- cycle count and phase;
- fault flags;
- explicit units/version identifier.

At minimum each telemetry record should carry

\[
\boxed{
\{t,\,N,\,\phi,\,\sigma_{\rm ref},\,F_{\rm ref},\,F,\,u,\,\epsilon,\,T,\,V_{\rm DCPD},\,u_{\rm act},\,\text{faults}\}
}
\]

so the PC can reconstruct the loading history and feed the probability solver without guessing timing.

---

## 3. Recommended signal-rate separation

The channels do not all require the same bandwidth.

- **force feedback:** highest-priority synchronized real-time channel;
- **actuator displacement / travel:** synchronized with the control loop or sufficiently fast for safety;
- **specimen strain/extensometer:** high-rate measurement for hysteresis loops;
- **DCPD:** lower rate is acceptable than force feedback, but it must be timestamped and sufficiently fast to localize crack-initiation changes;
- **temperature:** much slower rate is acceptable;
- **PC probability solver:** may process streamed data online or replay logged data offline; it is not part of the hard real-time loop.

Do not force every sensor into the 10 kHz example rate. The final sample schedule shall be derived from sensor bandwidth, control bandwidth, anti-aliasing requirements, and storage/telemetry capacity.

---

## 4. Hardware-to-HAL mapping

The existing target-independent firmware interface maps to hardware as follows.

| HAL / firmware quantity | Required hardware source |
|---|---|
| `measured_force_n` | load cell + bridge excitation + AFE + ADC |
| `displacement_m` | LVDT / linear encoder / actuator travel sensor |
| `strain` | extensometer or calibrated gauge-length displacement |
| `temperature_c` | RTD / thermocouple / temperature sensor |
| `dcpd_v` | DCPD constant-current + Kelvin voltage AFE + ADC |
| `estop_active` | hardwired E-stop / safety relay status |
| `travel_limit_active` | physical limit switches |
| `ft_hal_write_actuator()` | servo drive analog/PWM/CAN/EtherCAT command adapter |
| `ft_hal_send_telemetry()` | USB/UART/RS-485/CAN/Ethernet adapter |
| `ft_hal_watchdog_kick()` | MCU independent watchdog |
| `ft_hal_safe_shutdown()` | drive-enable removal / safe-state output |

---

## 5. Probability-theory integration boundary

The experimental pipeline shall be

\[
\boxed{
\text{tester hardware}
\rightarrow
\text{MCU telemetry}
\rightarrow
\text{PC logger}
\rightarrow
P(a,s,t)
\rightarrow
\{\bar a,\bar U,E_{\rm hyst},S\}
\rightarrow
P_{\rm crack}
}
\]

The probability solver can therefore be imported/reused independently as long as the telemetry contract is stable. This is preferable to embedding the theory solver into the MCU.

---

## 6. Procurement/design freeze order

Freeze hardware in this order:

1. specimen geometry and maximum stress envelope;
2. required force/stroke/velocity/frequency;
3. actuator + drive;
4. load cell;
5. load-cell AFE/ADC;
6. displacement/extensometer sensing;
7. MCU and communication interface;
8. DCPD current source + AFE/ADC;
9. temperature sensing;
10. E-stop, travel limits, power isolation, enclosure and grounding;
11. only then finalize the board-specific HAL.

This order prevents selecting an MCU or ADC before the actual machine bandwidth and drive interface are known.

---

# 한국어 요약

피로시험기에서 MCU가 해야 할 일은 **확률 시뮬레이션이 아니라 실제 반복 인장하중을 안정적으로 만들고 측정하는 것**이다.

핵심 하드웨어는 다음과 같다.

- 실시간 MCU;
- 동적 인장 actuator + drive;
- load cell;
- load-cell instrumentation amplifier + 고속 고해상도 ADC;
- displacement/strain sensor;
- DCPD constant-current source + Kelvin 측정 AFE/ADC;
- temperature sensor;
- E-stop / limit switch / safety relay / power isolation;
- PC telemetry interface.

actuator와 load cell 용량은 임의로 정하지 않고

\[
F_{\rm req}=A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|)
\]

에서 정한다.

확률 이론은 PC에서

\[
P(a,s,t)\rightarrow\bar a,\bar U,E_{\rm hyst},S\rightarrow P_{\rm crack}
\]

를 계산하고, MCU는 정확한 timestamp가 붙은 force/strain/DCPD/temperature 데이터를 제공하는 구조로 고정한다.
