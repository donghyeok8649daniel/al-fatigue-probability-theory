/* === 한국어 파일 안내 시작 ===
파일 역할: load cell, actuator, 시간, 안전정지 등 하드웨어 의존 기능을 제어 로직에서 분리하는 HAL 인터페이스다.
주요 자료형: 없음 또는 다른 헤더에서 정의
주요 함수: ft_hal_read_sample, ft_hal_write_actuator, ft_hal_send_telemetry, ft_hal_watchdog_kick, ft_hal_safe_shutdown
주의: 이 안내는 코드 탐색용이며 실제 하드웨어 구동에는 별도의 HAL 및 limit/interlock 검증이 필요하다.
=== 한국어 파일 안내 끝 === */
#ifndef FATIGUE_HAL_H
#define FATIGUE_HAL_H

#include "fatigue_controller.h"

/*
 * Board-specific adapter.
 * Implement these functions for the selected MCU/ADC/DAC/PWM/load-cell stack.
 * The core controller does not know about STM32, ESP32, Teensy, etc.
 */
bool ft_hal_read_sample(ft_sample_t *sample);
void ft_hal_write_actuator(double normalized_command);
void ft_hal_send_telemetry(const ft_output_t *out, const ft_sample_t *sample);
void ft_hal_watchdog_kick(void);
void ft_hal_safe_shutdown(void);

#endif
