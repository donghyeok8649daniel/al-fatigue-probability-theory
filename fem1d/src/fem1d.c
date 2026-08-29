/*
 * === 한국어 파일 안내 시작 ===
 * - 파일 역할: 2절점 1D bar element를 조립해 준정적 축인장 문제를 풀고 시간별 node/element CSV를 생성한다.
 * - 주요 기능: tridiagonal stiffness solve, element strain/stress recovery, cyclic stress history export, analytical uniform-bar self-test.
 * - 현재 물리: epsilon=du/dx, sigma=E*epsilon, K_e=(EA/le)[[1,-1],[-1,1]].
 * - 주의: 이 코드는 피로 damage나 P(lambda,t)를 임의로 넣지 않는다. 현재는 reversible linear-elastic null/reference solver다.
 * === 한국어 파일 안내 끝 ===
 */
#include "fem1d.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define FEM1D_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define FEM1D_MKDIR(path) mkdir(path, 0777)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int ensure_directory(const char *path) {
    if (path == NULL || path[0] == '\0') {
        return -1;
    }
    if (FEM1D_MKDIR(path) == 0 || errno == EEXIST) {
        return 0;
    }
    return -1;
}

static int solve_tridiagonal(
    int n,
    const double *lower,
    const double *diag,
    const double *upper,
    const double *rhs,
    double *solution
) {
    double *cprime;
    double *dprime;
    int i;

    if (n <= 0 || diag == NULL || rhs == NULL || solution == NULL) {
        return -1;
    }

    cprime = (double *)calloc((size_t)n, sizeof(double));
    dprime = (double *)calloc((size_t)n, sizeof(double));
    if (cprime == NULL || dprime == NULL) {
        free(cprime);
        free(dprime);
        return -1;
    }

    if (fabs(diag[0]) < 1.0e-300) {
        free(cprime);
        free(dprime);
        return -1;
    }

    cprime[0] = (n > 1) ? upper[0] / diag[0] : 0.0;
    dprime[0] = rhs[0] / diag[0];

    for (i = 1; i < n; ++i) {
        const double denom = diag[i] - lower[i - 1] * cprime[i - 1];
        if (fabs(denom) < 1.0e-300) {
            free(cprime);
            free(dprime);
            return -1;
        }
        cprime[i] = (i < n - 1) ? upper[i] / denom : 0.0;
        dprime[i] = (rhs[i] - lower[i - 1] * dprime[i - 1]) / denom;
    }

    solution[n - 1] = dprime[n - 1];
    for (i = n - 2; i >= 0; --i) {
        solution[i] = dprime[i] - cprime[i] * solution[i + 1];
    }

    free(cprime);
    free(dprime);
    return 0;
}

int fem1d_validate_model(const Fem1DModel *model, char *message, size_t message_size) {
    const char *error = NULL;

    if (model == NULL) {
        error = "model is null";
    } else if (model->elements < 1) {
        error = "elements must be >= 1";
    } else if (!(model->length_m > 0.0)) {
        error = "length_m must be positive";
    } else if (!(model->area_m2 > 0.0)) {
        error = "area_m2 must be positive";
    } else if (!(model->young_pa > 0.0)) {
        error = "young_pa must be positive";
    } else if (!(model->frequency_hz > 0.0)) {
        error = "frequency_hz must be positive";
    } else if (model->cycles < 1) {
        error = "cycles must be >= 1";
    } else if (model->steps_per_cycle < 4) {
        error = "steps_per_cycle must be >= 4";
    }

    if (message != NULL && message_size > 0) {
        if (error == NULL) {
            snprintf(message, message_size, "ok");
        } else {
            snprintf(message, message_size, "%s", error);
        }
    }
    return error == NULL ? 0 : -1;
}

int fem1d_solve_quasistatic_step(
    const Fem1DModel *model,
    double applied_stress_pa,
    double *nodal_displacement_m,
    double *element_strain,
    double *element_stress_pa
) {
    const int nfree = model != NULL ? model->elements : 0;
    const int nodes = nfree + 1;
    double *lower = NULL;
    double *diag = NULL;
    double *upper = NULL;
    double *rhs = NULL;
    double *ufree = NULL;
    double element_length;
    double element_stiffness;
    int i;
    char message[128];

    if (fem1d_validate_model(model, message, sizeof(message)) != 0 ||
        nodal_displacement_m == NULL || element_strain == NULL ||
        element_stress_pa == NULL) {
        return -1;
    }

    lower = (double *)calloc((size_t)(nfree > 1 ? nfree - 1 : 1), sizeof(double));
    diag = (double *)calloc((size_t)nfree, sizeof(double));
    upper = (double *)calloc((size_t)(nfree > 1 ? nfree - 1 : 1), sizeof(double));
    rhs = (double *)calloc((size_t)nfree, sizeof(double));
    ufree = (double *)calloc((size_t)nfree, sizeof(double));
    if (lower == NULL || diag == NULL || upper == NULL || rhs == NULL || ufree == NULL) {
        free(lower);
        free(diag);
        free(upper);
        free(rhs);
        free(ufree);
        return -1;
    }

    element_length = model->length_m / (double)model->elements;
    element_stiffness = model->young_pa * model->area_m2 / element_length;

    for (i = 0; i < nfree; ++i) {
        diag[i] = (i == nfree - 1) ? element_stiffness : 2.0 * element_stiffness;
    }
    for (i = 0; i < nfree - 1; ++i) {
        lower[i] = -element_stiffness;
        upper[i] = -element_stiffness;
    }
    rhs[nfree - 1] = applied_stress_pa * model->area_m2;

    if (solve_tridiagonal(nfree, lower, diag, upper, rhs, ufree) != 0) {
        free(lower);
        free(diag);
        free(upper);
        free(rhs);
        free(ufree);
        return -1;
    }

    nodal_displacement_m[0] = 0.0;
    for (i = 1; i < nodes; ++i) {
        nodal_displacement_m[i] = ufree[i - 1];
    }
    for (i = 0; i < model->elements; ++i) {
        element_strain[i] =
            (nodal_displacement_m[i + 1] - nodal_displacement_m[i]) / element_length;
        element_stress_pa[i] = model->young_pa * element_strain[i];
    }

    free(lower);
    free(diag);
    free(upper);
    free(rhs);
    free(ufree);
    return 0;
}

int fem1d_write_history(const Fem1DModel *model, const char *output_dir) {
    const int nodes = model != NULL ? model->elements + 1 : 0;
    const int total_steps = model != NULL ? model->cycles * model->steps_per_cycle + 1 : 0;
    double *u = NULL;
    double *strain = NULL;
    double *stress = NULL;
    FILE *node_file = NULL;
    FILE *element_file = NULL;
    FILE *metadata_file = NULL;
    char node_path[1024];
    char element_path[1024];
    char metadata_path[1024];
    char message[128];
    int step;

    if (fem1d_validate_model(model, message, sizeof(message)) != 0 || output_dir == NULL) {
        return -1;
    }
    if (ensure_directory(output_dir) != 0) {
        return -1;
    }

    snprintf(node_path, sizeof(node_path), "%s/nodes.csv", output_dir);
    snprintf(element_path, sizeof(element_path), "%s/elements.csv", output_dir);
    snprintf(metadata_path, sizeof(metadata_path), "%s/metadata.csv", output_dir);

    node_file = fopen(node_path, "w");
    element_file = fopen(element_path, "w");
    metadata_file = fopen(metadata_path, "w");
    if (node_file == NULL || element_file == NULL || metadata_file == NULL) {
        if (node_file != NULL) fclose(node_file);
        if (element_file != NULL) fclose(element_file);
        if (metadata_file != NULL) fclose(metadata_file);
        return -1;
    }

    u = (double *)calloc((size_t)nodes, sizeof(double));
    strain = (double *)calloc((size_t)model->elements, sizeof(double));
    stress = (double *)calloc((size_t)model->elements, sizeof(double));
    if (u == NULL || strain == NULL || stress == NULL) {
        fclose(node_file);
        fclose(element_file);
        fclose(metadata_file);
        free(u);
        free(strain);
        free(stress);
        return -1;
    }

    fprintf(node_file, "time_s,step,node,x_m,displacement_m,applied_stress_pa\n");
    fprintf(element_file, "time_s,step,element,x_mid_m,strain,stress_pa,applied_stress_pa\n");
    fprintf(metadata_file, "key,value\n");
    fprintf(metadata_file, "elements,%d\n", model->elements);
    fprintf(metadata_file, "length_m,%.17g\n", model->length_m);
    fprintf(metadata_file, "area_m2,%.17g\n", model->area_m2);
    fprintf(metadata_file, "young_pa,%.17g\n", model->young_pa);
    fprintf(metadata_file, "stress_mean_pa,%.17g\n", model->stress_mean_pa);
    fprintf(metadata_file, "stress_amplitude_pa,%.17g\n", model->stress_amplitude_pa);
    fprintf(metadata_file, "frequency_hz,%.17g\n", model->frequency_hz);
    fprintf(metadata_file, "cycles,%d\n", model->cycles);
    fprintf(metadata_file, "steps_per_cycle,%d\n", model->steps_per_cycle);
    fprintf(metadata_file, "solver,quasistatic_linear_bar_fem\n");
    fprintf(metadata_file, "probability_coupling,none_scaffold_only\n");

    for (step = 0; step < total_steps; ++step) {
        const double t = ((double)step / (double)model->steps_per_cycle) / model->frequency_hz;
        const double applied = model->stress_mean_pa + model->stress_amplitude_pa *
            sin(2.0 * M_PI * model->frequency_hz * t);
        const double dx = model->length_m / (double)model->elements;
        int i;

        if (fem1d_solve_quasistatic_step(model, applied, u, strain, stress) != 0) {
            fclose(node_file);
            fclose(element_file);
            fclose(metadata_file);
            free(u);
            free(strain);
            free(stress);
            return -1;
        }

        for (i = 0; i < nodes; ++i) {
            fprintf(
                node_file,
                "%.17g,%d,%d,%.17g,%.17g,%.17g\n",
                t,
                step,
                i,
                dx * (double)i,
                u[i],
                applied
            );
        }
        for (i = 0; i < model->elements; ++i) {
            fprintf(
                element_file,
                "%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g\n",
                t,
                step,
                i,
                dx * ((double)i + 0.5),
                strain[i],
                stress[i],
                applied
            );
        }
    }

    fclose(node_file);
    fclose(element_file);
    fclose(metadata_file);
    free(u);
    free(strain);
    free(stress);
    return 0;
}

int fem1d_self_test(void) {
    Fem1DModel model;
    const double applied = 123.4e6;
    const double stress_scale = fabs(applied) > 1.0 ? fabs(applied) : 1.0;
    double *u;
    double *strain;
    double *stress;
    double max_u_error = 0.0;
    double max_stress_error = 0.0;
    int i;

    model.elements = 16;
    model.length_m = 0.20;
    model.area_m2 = 2.0e-6;
    model.young_pa = 69.0e9;
    model.stress_mean_pa = 0.0;
    model.stress_amplitude_pa = applied;
    model.frequency_hz = 20.0;
    model.cycles = 1;
    model.steps_per_cycle = 16;

    u = (double *)calloc((size_t)(model.elements + 1), sizeof(double));
    strain = (double *)calloc((size_t)model.elements, sizeof(double));
    stress = (double *)calloc((size_t)model.elements, sizeof(double));
    if (u == NULL || strain == NULL || stress == NULL) {
        free(u);
        free(strain);
        free(stress);
        return -1;
    }

    if (fem1d_solve_quasistatic_step(&model, applied, u, strain, stress) != 0) {
        free(u);
        free(strain);
        free(stress);
        return -1;
    }

    for (i = 0; i <= model.elements; ++i) {
        const double x = model.length_m * (double)i / (double)model.elements;
        const double exact_u = (applied / model.young_pa) * x;
        const double error = fabs(u[i] - exact_u);
        if (error > max_u_error) max_u_error = error;
    }
    for (i = 0; i < model.elements; ++i) {
        const double error = fabs(stress[i] - applied);
        if (error > max_stress_error) max_stress_error = error;
    }

    free(u);
    free(strain);
    free(stress);

    if (max_u_error > 1.0e-11 * model.length_m + 1.0e-15) {
        fprintf(stderr, "self-test displacement error too large: %.17g\n", max_u_error);
        return -1;
    }
    if (max_stress_error / stress_scale > 1.0e-10) {
        fprintf(stderr, "self-test stress relative error too large: %.17g\n", max_stress_error / stress_scale);
        return -1;
    }

    printf("FEM1D self-test passed: max_u_error=%.6e m, max_stress_rel_error=%.6e\n",
           max_u_error, max_stress_error / stress_scale);
    return 0;
}
