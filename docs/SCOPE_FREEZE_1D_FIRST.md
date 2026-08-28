# Active Scope Freeze — One-Dimensional Tension First

## Decision

The active research program is frozen to **one-dimensional normal tension** until the probability-density mechanics is internally complete and numerically validated.

The current active chain is

$$
\boxed{
\text{1D layer-LJ mechanics}
\rightarrow
P(\lambda,t)
\rightarrow
\text{tail / first-passage criterion}
}
$$

or, when spatial position along a one-dimensional specimen must be retained,

$$
\boxed{
\sigma(x,t)
\rightarrow
P(\lambda;x,t)
\rightarrow
Q_c(x,t)\text{ or a 1D first-passage quantity}.
}
$$

## Explicitly deferred

The following are **not active tasks** at this stage:

- 2D or 3D finite-element models;
- FEM Gauss-point coupling;
- specimen-scale hazard integration;
- multiaxial stress criteria;
- principal-stress projection rules;
- probabilistic element deletion, XFEM, cohesive-zone coupling, or crack propagation;
- a three-dimensional statistical mini-mesh;
- orientation distributions or FCC/shear reactivation.

These topics may be revisited only after the one-dimensional theory is closed enough to define and validate $P(\lambda,t)$ under one-dimensional tensile loading.

## Immediate priority

The immediate priority is the **functional form and evolution of the one-dimensional spacing density**. In particular:

1. retain the full nonlinear calibrated layer-LJ potential;
2. avoid global Taylor and finite-harmonic ansatzes for the full distribution;
3. distinguish exact mechanics, statistical-mechanical ensemble assumptions, controlled approximations, and empirical inputs;
4. determine which physically derived form of $P$ applies to the actual 1D driven chain;
5. test that form directly against deterministic 1D layer-LJ dynamics;
6. only after this succeeds, define a physical coarse-graining length/cell for a 1D specimen.

A future discretized one-dimensional bar is allowed only as a **1D extension of the same theory**, not as a shortcut into multidimensional FEM.

---

# 활성 범위 고정 — 우선 1차원 인장만 수행

## 결정

확률밀도 mechanics가 내부적으로 완성되고 수치적으로 검증될 때까지 활성 연구 범위를 **1차원 수직 인장**으로 고정한다.

현재 활성 연결은

$$
\boxed{
\text{1D layer-LJ mechanics}
\rightarrow
P(\lambda,t)
\rightarrow
\text{tail / first-passage criterion}
}
$$

이다.

1차원 시편 내부 위치 $x$까지 필요한 경우에만

$$
\boxed{
\sigma(x,t)
\rightarrow
P(\lambda;x,t)
\rightarrow
Q_c(x,t)\text{ 또는 1D first-passage quantity}
}
$$

로 확장한다.

## 명시적으로 보류하는 항목

현 단계에서는 다음을 활성 과제로 다루지 않는다.

- 2D/3D 유한요소 모델;
- FEM Gauss-point coupling;
- 시편 전체 hazard 적분;
- 다축 응력 기준;
- 주응력 투영 규칙;
- 확률적 element deletion, XFEM, cohesive-zone coupling, crack propagation;
- 3차원 statistical mini-mesh;
- orientation distribution 또는 FCC/shear 재활성화.

이 항목들은 1차원 인장하중 아래 $P(\lambda,t)$를 정의하고 검증할 수 있을 정도로 1D 이론이 닫힌 뒤에만 다시 검토한다.

## 즉시 우선순위

지금의 최우선 목표는 **1차원 layer-spacing density의 함수형과 진화**다.

1. full nonlinear calibrated layer-LJ potential을 그대로 유지한다.
2. 전체 분포에 대한 global Taylor 전개와 finite-harmonic ansatz를 사용하지 않는다.
3. exact mechanics, statistical-mechanical ensemble assumption, controlled approximation, empirical input을 엄격히 구분한다.
4. 실제 1D driven chain에 어떤 물리적으로 유도된 $P$ 형식이 적용되는지 결정한다.
5. 해당 형식을 deterministic 1D layer-LJ dynamics와 직접 비교한다.
6. 이것이 성공한 뒤에만 1D 시편용 physical coarse-graining length/cell을 정의한다.

향후 discretized 1D bar는 허용하지만, 그것도 동일한 이론의 **1D 확장**이어야 하며 multidimensional FEM으로 넘어가기 위한 우회로로 사용하지 않는다.
