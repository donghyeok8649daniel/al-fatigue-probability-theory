/*
 * === 한국어 파일 안내 시작 ===
 * - 파일 역할: 1D bar FEM 실행파일의 CLI를 제공하고 모델 입력을 받아 cyclic 준정적 해석을 실행한다.
 * - 주요 기능: 기본 Al-like demo parameter, 명령행 옵션 파싱, self-test, CSV 출력 디렉터리 지정.
 * - 출력: nodes.csv, elements.csv, metadata.csv.
 * - 주의: 확률밀도함수/피로손상은 아직 결합하지 않는다. 이 실행기는 순수 1D FEM scaffold다.
 * === 한국어 파일 안내 끝 ===
 */
#include "fem1d.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_help(const char *program) {
    printf("Usage: %s [options]\n", program);
    printf("\n1D quasistatic linear bar FEM scaffold.\n\n");
    printf("Options:\n");
    printf("  --elements N                 number of linear bar elements (default 40)\n");
    printf("  --length-m X                 bar length in m (default 0.05)\n");
    printf("  --area-m2 X                  cross-sectional area in m^2 (default 1e-6)\n");
    printf("  --young-pa X                 Young's modulus in Pa (default 69e9)\n");
    printf("  --stress-mean-mpa X          mean applied axial stress in MPa (default 50)\n");
    printf("  --stress-amplitude-mpa X     stress amplitude in MPa (default 100)\n");
    printf("  --frequency-hz X             cyclic frequency in Hz (default 20)\n");
    printf("  --cycles N                   number of cycles (default 2)\n");
    printf("  --steps-per-cycle N          temporal samples per cycle (default 80)\n");
    printf("  --outdir PATH                output directory (default fem1d_output)\n");
    printf("  --self-test                  verify uniform-bar analytical solution\n");
    printf("  --help                       show this help\n");
}

static int require_value(int argc, char **argv, int *index, const char **value) {
    if (*index + 1 >= argc) {
        fprintf(stderr, "Missing value after %s\n", argv[*index]);
        return -1;
    }
    *index += 1;
    *value = argv[*index];
    return 0;
}

int main(int argc, char **argv) {
    Fem1DModel model;
    const char *outdir = "fem1d_output";
    int run_self_test = 0;
    int i;
    char validation_message[128];

    model.elements = 40;
    model.length_m = 0.05;
    model.area_m2 = 1.0e-6;
    model.young_pa = 69.0e9;
    model.stress_mean_pa = 50.0e6;
    model.stress_amplitude_pa = 100.0e6;
    model.frequency_hz = 20.0;
    model.cycles = 2;
    model.steps_per_cycle = 80;

    for (i = 1; i < argc; ++i) {
        const char *value = NULL;

        if (strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        }
        if (strcmp(argv[i], "--self-test") == 0) {
            run_self_test = 1;
            continue;
        }
        if (strcmp(argv[i], "--elements") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.elements = (int)strtol(value, NULL, 10);
        } else if (strcmp(argv[i], "--length-m") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.length_m = strtod(value, NULL);
        } else if (strcmp(argv[i], "--area-m2") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.area_m2 = strtod(value, NULL);
        } else if (strcmp(argv[i], "--young-pa") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.young_pa = strtod(value, NULL);
        } else if (strcmp(argv[i], "--stress-mean-mpa") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.stress_mean_pa = 1.0e6 * strtod(value, NULL);
        } else if (strcmp(argv[i], "--stress-amplitude-mpa") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.stress_amplitude_pa = 1.0e6 * strtod(value, NULL);
        } else if (strcmp(argv[i], "--frequency-hz") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.frequency_hz = strtod(value, NULL);
        } else if (strcmp(argv[i], "--cycles") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.cycles = (int)strtol(value, NULL, 10);
        } else if (strcmp(argv[i], "--steps-per-cycle") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            model.steps_per_cycle = (int)strtol(value, NULL, 10);
        } else if (strcmp(argv[i], "--outdir") == 0) {
            if (require_value(argc, argv, &i, &value) != 0) return 2;
            outdir = value;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_help(argv[0]);
            return 2;
        }
    }

    if (run_self_test) {
        return fem1d_self_test() == 0 ? 0 : 1;
    }

    if (fem1d_validate_model(&model, validation_message, sizeof(validation_message)) != 0) {
        fprintf(stderr, "Invalid model: %s\n", validation_message);
        return 2;
    }

    if (fem1d_write_history(&model, outdir) != 0) {
        fprintf(stderr, "FEM solve/output failed.\n");
        return 1;
    }

    printf("1D FEM history written to %s\n", outdir);
    printf("Scope: linear elastic, quasistatic, no fatigue/probability coupling yet.\n");
    return 0;
}
