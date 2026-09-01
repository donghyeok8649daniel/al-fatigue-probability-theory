# Normal-LJ Physical Statistical Distribution Diagnostic

## Classification

**DIMENSIONLESS MODEL DIAGNOSTIC / CONTROLLED METASTABLE APPROXIMATION.**

This report evaluates the full nonlinear reduced layer-LJ potential. It does not assign a physical aluminum temperature because the representative coarse-grained area $A_0$ has not yet been fixed.

The critical reduced force and tangent-instability spacing are

$$
f_c=0.0370342696708,
\qquad
\lambda_c=1.10777153855.
$$

For each $f/f_c$, the stable point $\lambda_s$, unstable barrier point $\lambda_b$, and

$$
\Delta w=w_f(\lambda_b)-w_f(\lambda_s)
$$

are computed without a Taylor expansion. Metastable densities are conditioned on $0<\lambda<\lambda_b$ and evaluated for illustrative dimensionless $\chi$ values.

The CSV and JSON files contain the complete numerical table. The important qualitative checks are:

- $\Delta w>0$ for every tested $0<f<f_c$;
- the barrier decreases as $f	o f_c^-$;
- increasing $\chi$ concentrates the intact-basin density near $\lambda_s$;
- the reported $Q_c$ is an instantaneous basin population, not an escape rate or fatigue life.

---

# 한국어 번역 — Normal-LJ 물리 통계분포 진단

## 분류

**DIMENSIONLESS MODEL DIAGNOSTIC / CONTROLLED METASTABLE APPROXIMATION.**

이 보고서는 full nonlinear reduced layer-LJ potential을 계산한다. representative coarse-grained area $A_0$가 아직 물리적으로 정해지지 않았으므로 실제 aluminum temperature를 부여하지 않는다.

critical reduced force와 tangent-instability spacing은

$$
f_c=0.0370342696708,
\qquad
\lambda_c=1.10777153855
$$

이다.

각 $f/f_c$에 대해 stable point $\lambda_s$, unstable barrier point $\lambda_b$ 및

$$
\Delta w=w_f(\lambda_b)-w_f(\lambda_s)
$$

를 Taylor expansion 없이 계산한다. metastable density는 $0<\lambda<\lambda_b$에 조건부로 두고 illustrative dimensionless $\chi$ 값에서 계산한다.

전체 수치표는 CSV와 JSON에 저장한다. 핵심 qualitative check는 다음과 같다.

- 모든 tested $0<f<f_c$에서 $\Delta w>0$;
- $f	o f_c^-$이면 barrier가 감소;
- $\chi$가 증가하면 intact-basin density가 $\lambda_s$ 근처로 집중;
- $Q_c$는 instantaneous basin population이며 escape rate나 fatigue life가 아님.
