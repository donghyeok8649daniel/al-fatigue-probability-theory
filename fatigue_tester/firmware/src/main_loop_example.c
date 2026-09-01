/* === 한국어 파일 안내 시작 ===
파일 역할: MCU main loop에서 HAL과 fatigue controller를 연결하는 사용 예시다.
=== 한국어 파일 안내 끝 === */
#include "fatigue_controller.h"
#include "fatigue_hal.h"

static ft_controller_t g_controller;

void fatigue_control_init(void) {
    const ft_config_t cfg = {
        .stress_mean_pa = 0.0,
        .stress_amplitude_pa = 100.0e6,
        .frequency_hz = 20.0,
        .specimen_area_m2 = 7.068583470577034e-6,
        .target_cycles = 1000000u,
        .waveform = FT_WAVE_SINE,
        .kp = 0.0,
        .ki = 0.0,
        .max_abs_force_n = 1000.0,
        .max_abs_displacement_m = 0.005,
        .max_abs_actuator_command = 1.0,
        .control_dt_s = 0.0001
    };

    ft_init(&g_controller, &cfg);
}

void fatigue_control_start(void) {
    ft_start(&g_controller);
}

void fatigue_control_tick(void) {
    ft_sample_t sample = {0};

    if (!ft_hal_read_sample(&sample)) {
        sample.sensors_valid = false;
    }

    const ft_output_t out = ft_step(&g_controller, &sample);

    if (!out.running && out.fault_flags != FT_FAULT_NONE) {
        ft_hal_write_actuator(0.0);
        ft_hal_safe_shutdown();
    } else {
        ft_hal_write_actuator(out.actuator_command);
    }

    ft_hal_send_telemetry(&out, &sample);
    ft_hal_watchdog_kick();
}
