# Milestone 14 — Standalone 1D FEM Scaffold Before Probability Coupling

## Status

This milestone deliberately prepares the numerical continuum side **before** the layer-spacing probability law is finalized. The active scope remains strictly one-dimensional.

The FEM code must not force a form of $P(\lambda,t)$, a damage variable, an empirical fatigue law, or a crack-probability rule. It only produces verified local mechanical histories that can later become inputs to the probability model.

## 1. 1D continuum problem

For an axial bar with displacement $u(x,t)$, the general 1D balance law is

$$
\rho A\ddot u
=\frac{\partial}{\partial x}(A\sigma)+b.
$$

The first scaffold intentionally takes the quasistatic, zero-body-force limit

$$
\boxed{
\frac{\partial}{\partial x}(A\sigma)=0
}
$$

with small-strain linear elasticity

$$
\epsilon=\frac{\partial u}{\partial x},
\qquad
\boxed{\sigma=E\epsilon}.
$$

Classification: **ASSUMPTION / REFERENCE MODEL** for the constitutive choice; **EXACT** continuum balance within the stated 1D setting.

## 2. Two-node finite element

For a constant-area element of length $\ell_e$,

$$
\boxed{
K_e
=\frac{EA}{\ell_e}
\begin{bmatrix}
1&-1\\
-1&1
\end{bmatrix}
}
$$

and the element strain/stress recovery is

$$
\boxed{
\epsilon_e=\frac{u_{e+1}-u_e}{\ell_e},
\qquad
\sigma_e=E\epsilon_e.
}
$$

The C implementation assembles the resulting tridiagonal system after fixing the left node.

## 3. Cyclic loading used only as a mechanical history

The right-end traction is

$$
\boxed{
\sigma_{\rm app}(t)
=\sigma_m+\sigma_a\sin(2\pi f t).
}
$$

Every sampled time is solved quasistatically. No cycle-dependent state is accumulated by the FEM scaffold.

Therefore a uniform linear-elastic bar is a required reversible null case.

## 4. Analytical null solution

For a uniform bar with the left end fixed and right-end axial traction,

$$
\boxed{
\sigma(x,t)=\sigma_{\rm app}(t)
}
$$

and

$$
\boxed{
u(x,t)=\frac{\sigma_{\rm app}(t)}{E}x.}
$$

The C executable contains a `--self-test` that compares the finite-element result with this exact solution.

## 5. Data interface reserved for the probability theory

The solver exports:

- `nodes.csv`: time, node, axial coordinate, displacement, applied stress;
- `elements.csv`: time, element, midpoint, strain, stress, applied stress;
- `metadata.csv`: geometry, material, loading, discretization, solver classification.

The future coupling must be explicit, for example

$$
\{\sigma_e(t),\epsilon_e(t)\}
\longrightarrow
\text{physically justified local inputs to }P(\lambda,t),
$$

but **this arrow is not implemented yet**.

The current FEM element size is a numerical discretization length. It is not identified with $\ell_{\rm stat}^{(2)}$, a tail-clustering length, an atomic spacing, or any future statistical mini-cell length.

## 6. Python visualizer

`simulations/visualize_fem1d.py` reads the C output and generates four independent diagnostics:

1. displacement along the bar at a selected snapshot;
2. strain along the bar;
3. stress along the bar;
4. stress history at the center element.

No crack probability is drawn until the 1D probability law is physically fixed.

## 7. Deferred extensions

The scaffold is intentionally prepared for, but does not yet activate:

- nonlinear constitutive response;
- transient inertia and a consistent/lumped mass matrix;
- spatially varying area/material properties;
- $P(\lambda,t)$ state coupling;
- first-passage or crack-initiation probability;
- statistical-cell aggregation.

All of those require separate derivation/validation.

---

# 한국어 번역 — 확률결합 전 독립형 1D FEM 스캐폴드

## 상태

이번 단계에서는 layer-spacing 확률법칙이 완성되기 전에 **연속체 수치해석 쪽 기반만 미리 만든다.** 활성 범위는 계속 엄격한 1차원이다.

FEM 코드 안에 $P(\lambda,t)$의 함수형, damage variable, 경험적 피로법칙, 균열확률식을 미리 넣지 않는다. 이후 확률모델의 입력이 될 검증된 기계 이력만 만든다.

## 1. 1D 연속체 문제

축방향 bar의 변위를 $u(x,t)$라 하면 일반적인 1D 운동방정식은

$$
\rho A\ddot u
=\frac{\partial}{\partial x}(A\sigma)+b
$$

이다.

첫 scaffold에서는 의도적으로 준정적, body-force zero로 제한해

$$
\boxed{
\frac{\partial}{\partial x}(A\sigma)=0
}
$$

을 쓰고 small-strain 선형탄성

$$
\epsilon=\frac{\partial u}{\partial x},
\qquad
\boxed{\sigma=E\epsilon}
$$

을 사용한다.

분류: constitutive choice는 **ASSUMPTION / REFERENCE MODEL**, 명시된 1D 범위의 balance equation은 **EXACT**다.

## 2. 2절점 유한요소

길이 $\ell_e$인 일정 단면 element에 대해

$$
\boxed{
K_e
=\frac{EA}{\ell_e}
\begin{bmatrix}
1&-1\\
-1&1
\end{bmatrix}
}
$$

이고

$$
\boxed{
\epsilon_e=\frac{u_{e+1}-u_e}{\ell_e},
\qquad
\sigma_e=E\epsilon_e
}
$$

로 strain/stress를 복원한다.

C 구현은 좌단 고정 후 생기는 tridiagonal system을 실제로 푼다.

## 3. cyclic loading은 기계 이력으로만 사용

우단 traction은

$$
\boxed{
\sigma_{\rm app}(t)
=\sigma_m+\sigma_a\sin(2\pi f t)
}
$$

이다.

각 sampled time을 준정적으로 풀며 FEM 자체에는 cycle-dependent state가 없다. 따라서 uniform linear-elastic bar는 반드시 reversible null case여야 한다.

## 4. 해석적 null solution

균일 bar의 좌단을 고정하고 우단에 axial traction을 가하면

$$
\boxed{
\sigma(x,t)=\sigma_{\rm app}(t)
}
$$

이고

$$
\boxed{
u(x,t)=\frac{\sigma_{\rm app}(t)}{E}x}
$$

이다.

C 실행파일의 `--self-test`는 FEM 결과를 이 exact solution과 비교한다.

## 5. 확률이론 결합을 위한 데이터 인터페이스

solver는 다음을 출력한다.

- `nodes.csv`: 시간, node, 축좌표, displacement, applied stress;
- `elements.csv`: 시간, element, midpoint, strain, stress, applied stress;
- `metadata.csv`: 형상, 재료, 하중, discretization, solver 분류.

향후

$$
\{\sigma_e(t),\epsilon_e(t)\}
\longrightarrow
\text{$P(\lambda,t)$의 물리적으로 정당화된 local input}
$$

으로 연결할 수 있지만 **현재는 이 화살표를 구현하지 않는다.**

그리고 현재 FEM element 크기는 순수 numerical discretization length다. 이를 $\ell_{\rm stat}^{(2)}$, tail clustering length, 원자간격, 미래의 statistical mini-cell length와 동일시하지 않는다.

## 6. Python visualizer

`simulations/visualize_fem1d.py`는 C 출력 CSV를 읽어 다음 네 가지를 각각 그린다.

1. 선택 snapshot의 bar displacement;
2. strain 분포;
3. stress 분포;
4. 중앙 element의 stress history.

1D 확률법칙이 물리적으로 확정되기 전에는 crack probability를 그리지 않는다.

## 7. 보류된 확장

현재 scaffold는 다음 확장을 염두에 두지만 활성화하지 않는다.

- nonlinear constitutive response;
- transient inertia와 consistent/lumped mass matrix;
- 위치에 따른 area/material;
- $P(\lambda,t)$ state coupling;
- first-passage/crack-initiation probability;
- statistical-cell probability aggregation.

각 항목은 별도의 유도와 검증 뒤에만 추가한다.
