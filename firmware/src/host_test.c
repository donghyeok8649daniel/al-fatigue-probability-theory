/* === 한국어 파일 안내 시작 ===
파일 역할: 실제 MCU 없이 host 환경에서 firmware 제어 로직의 기본 동작을 확인한다.
주요 자료형: 없음 또는 다른 헤더에서 정의
주요 함수: main
주의: 이 안내는 코드 탐색용이며 실제 하드웨어 구동에는 별도의 HAL 및 limit/interlock 검증이 필요하다.
=== 한국어 파일 안내 끝 === */
#include "fatigue_controller.h"
#include <assert.h>
#include <stdio.h>

int main(void) {
    const ft_config_t cfg = {
        .stress_mean_pa = 0.0,
        .stress_amplitude_pa = 100.0e6,
        .frequency_hz = 20.0,
        .specimen_area_m2 = 1.0e-6,
        .target_cycles = 2u,
        .waveform = FT_WAVE_SINE,
        .kp = 1.0e-3,
        .ki = 1.0e-2,
        .max_abs_force_n = 1000.0,
        .max_abs_displacement_m = 0.01,
        .max_abs_actuator_command = 1.0,
        .control_dt_s = 1.0e-4
    };
    ft_controller_t ctl;
    ft_init(&ctl, &cfg);
    assert(ctl.fault_flags == FT_FAULT_NONE);
    ft_start(&ctl);

    ft_sample_t sample = {
        .measured_force_n = 0.0,
        .displacement_m = 0.0,
        .strain = 0.0,
        .temperature_c = 25.0,
        .dcpd_v = 0.0,
        .estop_active = false,
        .travel_limit_active = false,
        .sensors_valid = true
    };

    unsigned long steps = 0;
    while (ctl.running && steps < 200000UL) {
        const ft_output_t out = ft_step(&ctl, &sample);
        /* Simple idealized plant for a host logic test only. */
        sample.measured_force_n = out.force_reference_n;
        steps++;
    }

    assert(ctl.cycle_count >= 2u);
    assert((ctl.fault_flags & FT_FAULT_CYCLE_COMPLETE) != 0u);

    ft_init(&ctl, &cfg);
    ft_start(&ctl);
    sample.measured_force_n = cfg.max_abs_force_n * 1.1;
    const ft_output_t fault = ft_step(&ctl, &sample);
    assert(!fault.running);
    assert((fault.fault_flags & FT_FAULT_FORCE_LIMIT) != 0u);
    assert(fault.actuator_command == 0.0);

    puts("firmware core host tests: PASS");
    return 0;
}
