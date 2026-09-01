# Firmware Variable Definitions

This file defines the variables introduced by the hardware-independent fatigue-tester firmware core. Theory variables remain in `../../docs/VARIABLE_DEFINITIONS.md`.

## Configuration variables

| Symbol / field | Definition | Unit | Meaning |
|---|---|---:|---|
| `stress_mean_pa` | $\sigma_m$ | Pa | Mean commanded engineering stress |
| `stress_amplitude_pa` | $\sigma_a$ | Pa | Cyclic stress amplitude |
| `frequency_hz` | $f$ | Hz | Command waveform frequency |
| `specimen_area_m2` | $A_{\rm specimen}$ | m$^2$ | Reference specimen cross-sectional area used for stress-to-force conversion |
| `target_cycles` | $N_{\rm target}$ | cycle | Requested stopping cycle count; zero means no cycle-count stop |
| `kp` | $K_p$ | actuator-command/N | Proportional force-loop gain; identify on the actual machine |
| `ki` | $K_i$ | actuator-command/(N s) | Integral force-loop gain; identify on the actual machine |
| `max_abs_force_n` | $F_{\max}$ | N | Software force safety limit |
| `max_abs_displacement_m` | $x_{\max}$ | m | Software travel/displacement safety limit |
| `max_abs_actuator_command` | $u_{\max}$ | normalized or drive-specific | Absolute actuator command saturation |
| `control_dt_s` | $\Delta t_c$ | s | Fixed real-time controller update period |

The reference generator is

\[
\sigma_{\rm ref}(t)=\sigma_m+\sigma_a w(\phi),
\qquad
F_{\rm ref}(t)=A_{\rm specimen}\sigma_{\rm ref}(t).
\]

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

## Controller state and outputs

| Field | Definition | Meaning |
|---|---|---|
| `phase_rad` | $\phi$ | Current waveform phase in radians |
| `integrator` | $I_c$ | PI-controller integral state |
| `cycle_count` | $N$ | Completed waveform cycles |
| `stress_reference_pa` | $\sigma_{\rm ref}$ | Current stress command |
| `force_reference_n` | $F_{\rm ref}$ | Current force command |
| `actuator_command` | $u$ | Saturated actuator request |
| `fault_flags` | bit mask | Latched controller fault/status bits |
| `running` | boolean | Whether the controller is actively generating commands |

Software checks do not replace hardwired machine safety.
