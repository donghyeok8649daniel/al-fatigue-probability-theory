#ifndef FATIGUE_CONTROLLER_H
#define FATIGUE_CONTROLLER_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    FT_WAVE_SINE = 0,
    FT_WAVE_TRIANGLE = 1
} ft_waveform_t;

enum {
    FT_FAULT_NONE = 0u,
    FT_FAULT_ESTOP = 1u << 0,
    FT_FAULT_FORCE_LIMIT = 1u << 1,
    FT_FAULT_DISPLACEMENT_LIMIT = 1u << 2,
    FT_FAULT_SENSOR_INVALID = 1u << 3,
    FT_FAULT_CONFIG_INVALID = 1u << 4,
    FT_FAULT_CYCLE_COMPLETE = 1u << 5
};

typedef struct {
    double stress_mean_pa;
    double stress_amplitude_pa;
    double frequency_hz;
    double specimen_area_m2;
    uint64_t target_cycles;
    ft_waveform_t waveform;

    /* PI force-loop gains. Must be calibrated for the real actuator/load cell. */
    double kp;
    double ki;

    /* Safety limits. Set from machine and specimen limits, not theory output. */
    double max_abs_force_n;
    double max_abs_displacement_m;
    double max_abs_actuator_command;

    /* Fixed controller period. */
    double control_dt_s;
} ft_config_t;

typedef struct {
    double measured_force_n;
    double displacement_m;
    double strain;
    double temperature_c;
    double dcpd_v;
    bool estop_active;
    bool travel_limit_active;
    bool sensors_valid;
} ft_sample_t;

typedef struct {
    double stress_reference_pa;
    double force_reference_n;
    double actuator_command;
    double phase_rad;
    uint64_t cycle_count;
    uint32_t fault_flags;
    bool running;
} ft_output_t;

typedef struct {
    ft_config_t cfg;
    double phase_rad;
    double integrator;
    uint64_t cycle_count;
    uint32_t fault_flags;
    bool running;
} ft_controller_t;

bool ft_config_valid(const ft_config_t *cfg);
void ft_init(ft_controller_t *ctl, const ft_config_t *cfg);
void ft_start(ft_controller_t *ctl);
void ft_stop(ft_controller_t *ctl);
ft_output_t ft_step(ft_controller_t *ctl, const ft_sample_t *sample);

#ifdef __cplusplus
}
#endif

#endif
