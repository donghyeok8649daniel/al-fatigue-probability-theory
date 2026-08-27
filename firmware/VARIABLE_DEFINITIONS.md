# Firmware Variable Definitions

This file defines the variables introduced by the hardware-independent fatigue-tester firmware core. Theory variables remain in `docs/VARIABLE_DEFINITIONS.md`.

## Configuration variables

| Symbol / field | Definition | Unit | Meaning |
|---|---|---:|---|
| `stress_mean_pa` | $\sigma_m$ | Pa | Mean commanded engineering stress |
| `stress_amplitude_pa` | $\sigma_a$ | Pa | Cyclic stress amplitude |
| `frequency_hz` | $f$ | Hz | Command waveform frequency |
| `specimen_area_m2` | $A_{\rm specimen}$ | m$^2$ | Reference specimen cross-sectional area used for stress-to-force conversion |
| `target_cycles` | $N_{\rm target}$ | cycle | Requested stopping cycle count; zero means no cycle-count stop |
| `kp` | $K_p$ | actuator-command/N | Proportional force-loop gain; real value must be identified on the actual machine |
| `ki` | $K_i$ | actuator-command/(N s) | Integral force-loop gain; real value must be identified on the actual machine |
| `max_abs_force_n` | $F_{\max}$ | N | Software force safety limit |
| `max_abs_displacement_m` | $x_{\max}$ | m | Software travel/displacement safety limit |
| `max_abs_actuator_command` | $u_{\max}$ | normalized or drive-specific | Absolute actuator command saturation |
| `control_dt_s` | $\Delta t_c$ | s | Fixed real-time controller update period |

The reference generator is

$$
\sigma_{\rm ref}(t)=\sigma_m+\sigma_a w(\phi),
$$

where $w$ is the selected normalized waveform and

$$
F_{\rm ref}(t)=A_{\rm specimen}\sigma_{\rm ref}(t).
$$

## Sample variables

| Field | Unit | Meaning |
|---|---:|---|
| `measured_force_n` | N | Load-cell force measurement |
| `displacement_m` | m | Measured actuator/specimen displacement channel |
| `strain` | dimensionless | Measured or derived specimen strain |
| `temperature_c` | °C | Specimen/environment temperature channel |
| `dcpd_v` | V | Direct-current-potential-drop crack-sensing channel |
| `estop_active` | boolean | Emergency-stop input state |
| `travel_limit_active` | boolean | Hardware travel-limit state |
| `sensors_valid` | boolean | Aggregate validity flag for required sensor samples |

## Controller state and output variables

| Field | Definition | Meaning |
|---|---|---|
| `phase_rad` | $\phi$ | Current waveform phase in radians |
| `integrator` | $I_c$ | PI-controller integral state |
| `cycle_count` | $N$ | Completed waveform cycles |
| `stress_reference_pa` | $\sigma_{\rm ref}$ | Current stress command |
| `force_reference_n` | $F_{\rm ref}$ | Current force command |
| `actuator_command` | $u$ | Saturated normalized or drive-specific actuator request |
| `fault_flags` | bit mask | Latched controller fault/status bits |
| `running` | boolean | Whether the real-time controller is actively generating commands |

## Fault flags

- `FT_FAULT_ESTOP`: E-stop or travel-limit path active.
- `FT_FAULT_FORCE_LIMIT`: measured or requested force exceeds configured limit.
- `FT_FAULT_DISPLACEMENT_LIMIT`: displacement exceeds configured limit.
- `FT_FAULT_SENSOR_INVALID`: required sensor sample is missing/non-finite/invalid.
- `FT_FAULT_CONFIG_INVALID`: controller configuration failed validation.
- `FT_FAULT_CYCLE_COMPLETE`: requested cycle count has been completed.

These software variables and checks do not replace hardwired machine safety.

---

# 한국어 번역 — 펌웨어 변수 정의

이 파일은 hardware-independent 피로시험기 firmware core에서 새로 도입한 변수를 정의한다. 이론 변수는 `docs/VARIABLE_DEFINITIONS.md`에 유지한다.

## 설정 변수

| 기호 / field | 정의 | 단위 | 의미 |
|---|---|---:|---|
| `stress_mean_pa` | $\sigma_m$ | Pa | 명령 engineering stress의 평균값 |
| `stress_amplitude_pa` | $\sigma_a$ | Pa | 반복응력 진폭 |
| `frequency_hz` | $f$ | Hz | 명령 waveform 주파수 |
| `specimen_area_m2` | $A_{\rm specimen}$ | m$^2$ | stress-to-force 변환에 사용하는 시편 기준 단면적 |
| `target_cycles` | $N_{\rm target}$ | cycle | 목표 정지 cycle 수; 0이면 cycle-count 자동정지 없음 |
| `kp` | $K_p$ | actuator-command/N | force loop 비례 gain; 실제 장비에서 식별해야 함 |
| `ki` | $K_i$ | actuator-command/(N s) | force loop 적분 gain; 실제 장비에서 식별해야 함 |
| `max_abs_force_n` | $F_{\max}$ | N | software force 안전한계 |
| `max_abs_displacement_m` | $x_{\max}$ | m | software travel/displacement 안전한계 |
| `max_abs_actuator_command` | $u_{\max}$ | normalized 또는 drive-specific | actuator command 절대 포화값 |
| `control_dt_s` | $\Delta t_c$ | s | 고정 real-time controller update 주기 |

기준응력은

$$
\sigma_{\rm ref}(t)=\sigma_m+\sigma_a w(\phi)
$$

이고, $w$는 선택한 정규화 waveform이다. 기준힘은

$$
F_{\rm ref}(t)=A_{\rm specimen}\sigma_{\rm ref}(t)
$$

이다.

## 측정 sample 변수

| Field | 단위 | 의미 |
|---|---:|---|
| `measured_force_n` | N | load-cell 힘 측정값 |
| `displacement_m` | m | actuator/specimen 변위 채널 |
| `strain` | 무차원 | 측정 또는 계산된 시편 변형률 |
| `temperature_c` | °C | 시편/환경 온도 채널 |
| `dcpd_v` | V | DCPD 균열감지 채널 |
| `estop_active` | boolean | emergency-stop 입력상태 |
| `travel_limit_active` | boolean | hardware travel-limit 상태 |
| `sensors_valid` | boolean | 필수 센서 sample의 종합 유효성 flag |

## controller 상태 및 출력 변수

| Field | 정의 | 의미 |
|---|---|---|
| `phase_rad` | $\phi$ | 현재 waveform 위상 |
| `integrator` | $I_c$ | PI controller 적분상태 |
| `cycle_count` | $N$ | 완료된 waveform cycle 수 |
| `stress_reference_pa` | $\sigma_{\rm ref}$ | 현재 stress command |
| `force_reference_n` | $F_{\rm ref}$ | 현재 force command |
| `actuator_command` | $u$ | 포화가 적용된 actuator 요청값 |
| `fault_flags` | bit mask | latched fault/status bit |
| `running` | boolean | 실시간 controller가 명령을 생성 중인지 여부 |

## fault flag

- `FT_FAULT_ESTOP`: E-stop 또는 travel-limit 경로 활성.
- `FT_FAULT_FORCE_LIMIT`: 측정힘 또는 요청힘이 설정한계를 초과.
- `FT_FAULT_DISPLACEMENT_LIMIT`: 변위가 설정한계를 초과.
- `FT_FAULT_SENSOR_INVALID`: 필수 sensor sample이 없거나 non-finite/invalid.
- `FT_FAULT_CONFIG_INVALID`: controller configuration 검증 실패.
- `FT_FAULT_CYCLE_COMPLETE`: 요청 cycle 수 완료.

이 software 변수와 check는 hardwired machine safety를 대체하지 않는다.
