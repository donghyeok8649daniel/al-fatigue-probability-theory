# Crack-initiation research modules

Research modules are allocated separately so that work on geometry, actual material response, spatial coupling and initiation probability can be reviewed independently.

| Module | Assigned branch | Scope |
|---|---|---|
| cell_work_units | `research/20261018-cell-work-units` | 원래 a/s의 기하·가상일·길이/에너지·셀 면적/부피 정규화 |
| real_stress | `research/20261018-real-stress` | 이상강도·결함/소스 작동·실제 균열 개시의 분리와 검증 경로 |
| fvm_coupling | `research/20261018-fvm-coupling` | 단위와 외부 일을 보존하는 one-axis adapter, 미정값 guard, trial/commit |
| initiation_probability | `research/20261018-initiation-probability` | 셀 공동 생존과 독립 곱의 차이·오차와 평균 개시시간 조건 |

The common starting commit is `c8ccff7553cbb2e66b0c37a76ee8a24808478824` on `probability-pde-solver-v1`. This records the GitHub base only; unpublished local research is not implicitly included or overwritten.

Each module owns only its new directory. Units, virtual work, flux conservation, conditional probability assumptions and missing physical inputs must be reviewed before integration. The original production model and both main branches are preserved. No calibrated aluminium stress, specimen initiation probability or physical fatigue life is claimed by these allocation files.
