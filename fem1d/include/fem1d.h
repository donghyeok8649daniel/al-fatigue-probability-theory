/*
 * === 한국어 파일 안내 시작 ===
 * - 파일 역할: 1차원 선형 bar FEM의 모델 구조체와 해석/CSV 출력 API를 선언한다.
 * - 주요 자료형: Fem1DModel
 * - 주요 함수: fem1d_validate_model, fem1d_solve_quasistatic_step, fem1d_write_history, fem1d_self_test
 * - 현재 범위: small-strain, 선형탄성, 2절점 1D bar element, 좌단 고정/우단 축응력, 준정적 cyclic loading.
 * - 주의: P(lambda,t), 균열확률, damage law는 이 파일에 포함하지 않는다. FEM과 확률이론의 인터페이스만 준비한다.
 * === 한국어 파일 안내 끝 ===
 */
#ifndef FEM1D_H
#define FEM1D_H

#include <stddef.h>

typedef struct {
    int elements;
    double length_m;
    double area_m2;
    double young_pa;
    double stress_mean_pa;
    double stress_amplitude_pa;
    double frequency_hz;
    int cycles;
    int steps_per_cycle;
} Fem1DModel;

int fem1d_validate_model(const Fem1DModel *model, char *message, size_t message_size);

int fem1d_solve_quasistatic_step(
    const Fem1DModel *model,
    double applied_stress_pa,
    double *nodal_displacement_m,
    double *element_strain,
    double *element_stress_pa
);

int fem1d_write_history(const Fem1DModel *model, const char *output_dir);

int fem1d_self_test(void);

#endif
