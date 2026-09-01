/* === 한국어 파일 안내 시작 ===
파일 역할: 파형 생성, force reference, 제어 상태 갱신 및 안전 로직을 구현한다.
주요 함수: clampd, waveform_value, ft_config_valid, ft_init, ft_start, ft_stop, ft_step
주의: 실제 하드웨어 구동에는 별도의 HAL 및 limit/interlock 검증이 필요하다.
=== 한국어 파일 안내 끝 === */
#include "fatigue_controller.h"

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double clampd(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double waveform_value(ft_waveform_t waveform, double phase) {
    if (waveform == FT_WAVE_TRIANGLE) {
        const double u = phase / (2.0 * M_PI);
        return 4.0 * fabs(u - floor(u + 0.5)) - 1.0;
    }
    return sin(phase);
}

bool ft_config_valid(const ft_config_t *cfg) {
    if (cfg == NULL) return false;
    if (!isfinite(cfg->stress_mean_pa) ||
        !isfinite(cfg->stress_amplitude_pa) ||
        !isfinite(cfg->frequency_hz) ||
        !isfinite(cfg->specimen_area_m2) ||
        !isfinite(cfg->kp) ||
        !isfinite(cfg->ki) ||
        !isfinite(cfg->control_dt_s)) {
        return false;
    }
    if (cfg->frequency_hz <= 0.0 ||
        cfg->specimen_area_m2 <= 0.0 ||
        cfg->control_dt_s <= 0.0 ||
        cfg->max_abs_force_n <= 0.0 ||
        cfg->max_abs_displacement_m <= 0.0 ||
        cfg->max_abs_actuator_command <= 0.0) {
        return false;
    }
    if (cfg->stress_amplitude_pa < 0.0) return false;
    return true;
}

void ft_init(ft_controller_t *ctl, const ft_config_t *cfg) {
    if (ctl == NULL) return;
    ctl->phase_rad = 0.0;
    ctl->integrator = 0.0;
    ctl->cycle_count = 0u;
    ctl->fault_flags = FT_FAULT_NONE;
    ctl->running = false;

    if (cfg == NULL || !ft_config_valid(cfg)) {
        ctl->fault_flags = FT_FAULT_CONFIG_INVALID;
        return;
    }
    ctl->cfg = *cfg;
}

void ft_start(ft_controller_t *ctl) {
    if (ctl == NULL) return;
    if (ctl->fault_flags != FT_FAULT_NONE) return;
    ctl->phase_rad = 0.0;
    ctl->integrator = 0.0;
    ctl->cycle_count = 0u;
    ctl->running = true;
}

void ft_stop(ft_controller_t *ctl) {
    if (ctl == NULL) return;
    ctl->running = false;
    ctl->integrator = 0.0;
}

ft_output_t ft_step(ft_controller_t *ctl, const ft_sample_t *sample) {
    ft_output_t out = {0};

    if (ctl == NULL || sample == NULL) {
        out.fault_flags = FT_FAULT_SENSOR_INVALID;
        return out;
    }

    if (!ctl->running) {
        out.cycle_count = ctl->cycle_count;
        out.fault_flags = ctl->fault_flags;
        out.phase_rad = ctl->phase_rad;
        out.running = false;
        return out;
    }

    if (sample->estop_active || sample->travel_limit_active) {
        ctl->fault_flags |= FT_FAULT_ESTOP;
    }
    if (!sample->sensors_valid ||
        !isfinite(sample->measured_force_n) ||
        !isfinite(sample->displacement_m)) {
        ctl->fault_flags |= FT_FAULT_SENSOR_INVALID;
    }
    if (fabs(sample->measured_force_n) > ctl->cfg.max_abs_force_n) {
        ctl->fault_flags |= FT_FAULT_FORCE_LIMIT;
    }
    if (fabs(sample->displacement_m) > ctl->cfg.max_abs_displacement_m) {
        ctl->fault_flags |= FT_FAULT_DISPLACEMENT_LIMIT;
    }

    if (ctl->fault_flags != FT_FAULT_NONE) {
        ctl->running = false;
        ctl->integrator = 0.0;
        out.cycle_count = ctl->cycle_count;
        out.fault_flags = ctl->fault_flags;
        out.phase_rad = ctl->phase_rad;
        out.running = false;
        out.actuator_command = 0.0;
        return out;
    }

    const double wave = waveform_value(ctl->cfg.waveform, ctl->phase_rad);
    out.stress_reference_pa =
        ctl->cfg.stress_mean_pa + ctl->cfg.stress_amplitude_pa * wave;
    out.force_reference_n =
        out.stress_reference_pa * ctl->cfg.specimen_area_m2;

    if (fabs(out.force_reference_n) > ctl->cfg.max_abs_force_n) {
        ctl->fault_flags |= FT_FAULT_FORCE_LIMIT;
        ctl->running = false;
        out.fault_flags = ctl->fault_flags;
        out.running = false;
        out.actuator_command = 0.0;
        return out;
    }

    const double error = out.force_reference_n - sample->measured_force_n;
    const double proposed_integrator =
        ctl->integrator + ctl->cfg.ki * error * ctl->cfg.control_dt_s;
    const double unsat = ctl->cfg.kp * error + proposed_integrator;
    const double sat = clampd(
        unsat,
        -ctl->cfg.max_abs_actuator_command,
        ctl->cfg.max_abs_actuator_command
    );

    if (sat == unsat ||
        (sat > 0.0 && error < 0.0) ||
        (sat < 0.0 && error > 0.0)) {
        ctl->integrator = proposed_integrator;
    }
    out.actuator_command = sat;

    const double dphase =
        2.0 * M_PI * ctl->cfg.frequency_hz * ctl->cfg.control_dt_s;
    ctl->phase_rad += dphase;
    if (ctl->phase_rad >= 2.0 * M_PI) {
        ctl->phase_rad = fmod(ctl->phase_rad, 2.0 * M_PI);
        ctl->cycle_count += 1u;
    }

    if (ctl->cfg.target_cycles > 0u &&
        ctl->cycle_count >= ctl->cfg.target_cycles) {
        ctl->fault_flags |= FT_FAULT_CYCLE_COMPLETE;
        ctl->running = false;
        out.actuator_command = 0.0;
    }

    out.phase_rad = ctl->phase_rad;
    out.cycle_count = ctl->cycle_count;
    out.fault_flags = ctl->fault_flags;
    out.running = ctl->running;
    return out;
}
