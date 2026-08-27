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
