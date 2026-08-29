# 1D FEM Scaffold

This directory contains a deliberately minimal finite-element scaffold for the active **one-dimensional normal-tension** research program.

## Purpose

The FEM solver is kept independent of the unfinished probability theory. It provides a verified continuum-scale stress/strain history that can later be passed to the layer-spacing probability model without hiding assumptions inside the FEM code.

Current chain:

$$
\sigma_{\rm app}(t)
\rightarrow
u(x,t),\epsilon(x,t),\sigma(x,t)
\rightarrow
\text{CSV interface}
\rightarrow
\text{future }P(\lambda,t)\text{ coupling}.
$$

The probability coupling is intentionally **not implemented yet**.

## Current mechanical model

A uniform small-strain linear-elastic bar is discretized with two-node elements. For an element of length $\ell_e$,

$$
K_e
=\frac{EA}{\ell_e}
\begin{bmatrix}
1&-1\\
-1&1
\end{bmatrix}.
$$

The left end is fixed and the right end carries the prescribed axial stress history

$$
\sigma_{\rm app}(t)
=\sigma_m+\sigma_a\sin(2\pi f t).
$$

Each time sample is solved quasistatically. The recovered element fields are

$$
\epsilon_e=\frac{u_{e+1}-u_e}{\ell_e},
\qquad
\sigma_e=E\epsilon_e.
$$

Classification: **EXACT finite-element discretization of the stated linear-elastic quasistatic boundary-value problem**. It is not yet a fatigue constitutive model.

## Build

```bash
make -C fem1d
make -C fem1d self-test
```

## Demo

```bash
./fem1d/bin/fem1d_solver \
  --elements 40 \
  --length-m 0.05 \
  --area-m2 1e-6 \
  --young-pa 69e9 \
  --stress-mean-mpa 50 \
  --stress-amplitude-mpa 100 \
  --frequency-hz 20 \
  --cycles 2 \
  --steps-per-cycle 80 \
  --outdir results/data/fem1d_demo

python simulations/visualize_fem1d.py \
  --input-dir results/data/fem1d_demo \
  --output-dir results/figures/fem1d_demo
```

The solver writes `nodes.csv`, `elements.csv`, and `metadata.csv`.

## Required null test

For a uniform bar under end traction, the exact continuum solution is

$$
u(x,t)=\frac{\sigma_{\rm app}(t)}{E}x,
\qquad
\sigma(x,t)=\sigma_{\rm app}(t).
$$

`--self-test` checks the finite-element result against this solution. This reference case must remain reversible and must never generate fatigue or probability accumulation by itself.

---

# 한국어 번역 — 1D FEM 스캐폴드

이 디렉터리는 현재 연구의 **1차원 수직 인장 문제**만을 위한 최소 유한요소해석 스캐폴드다.

## 목적

아직 확정되지 않은 확률이론을 FEM 코드 안에 미리 집어넣지 않는다. 먼저 검증된 연속체 응력/변형률 이력을 계산하고, 나중에 layer-spacing 확률모델에 명시적인 인터페이스로 넘긴다.

현재 연결은

$$
\sigma_{\rm app}(t)
\rightarrow
u(x,t),\epsilon(x,t),\sigma(x,t)
\rightarrow
\text{CSV 인터페이스}
\rightarrow
\text{향후 }P(\lambda,t)\text{ 결합}
$$

이다.

확률 결합은 의도적으로 **아직 구현하지 않았다**.

## 현재 기계 모델

균일한 small-strain 선형탄성 bar를 2절점 element로 나눈다. 길이가 $\ell_e$인 element의 stiffness는

$$
K_e
=\frac{EA}{\ell_e}
\begin{bmatrix}
1&-1\\
-1&1
\end{bmatrix}
$$

이다.

왼쪽 끝은 고정하고 오른쪽 끝에는

$$
\sigma_{\rm app}(t)
=\sigma_m+\sigma_a\sin(2\pi f t)
$$

의 축응력을 가한다.

각 시간점은 준정적으로 풀며 element field는

$$
\epsilon_e=\frac{u_{e+1}-u_e}{\ell_e},
\qquad
\sigma_e=E\epsilon_e
$$

로 복원한다.

분류: **명시된 선형탄성 준정적 경계값 문제에 대한 정확한 FEM discretization**이다. 아직 fatigue constitutive model은 아니다.

## 빌드

```bash
make -C fem1d
make -C fem1d self-test
```

## 데모

```bash
./fem1d/bin/fem1d_solver \
  --elements 40 \
  --length-m 0.05 \
  --area-m2 1e-6 \
  --young-pa 69e9 \
  --stress-mean-mpa 50 \
  --stress-amplitude-mpa 100 \
  --frequency-hz 20 \
  --cycles 2 \
  --steps-per-cycle 80 \
  --outdir results/data/fem1d_demo

python simulations/visualize_fem1d.py \
  --input-dir results/data/fem1d_demo \
  --output-dir results/figures/fem1d_demo
```

solver는 `nodes.csv`, `elements.csv`, `metadata.csv`를 출력한다.

## 필수 null test

균일 bar의 끝에 traction을 가하면 정확해는

$$
u(x,t)=\frac{\sigma_{\rm app}(t)}{E}x,
\qquad
\sigma(x,t)=\sigma_{\rm app}(t)
$$

이다.

`--self-test`는 FEM 결과를 이 해석해와 비교한다. 이 reference case는 계속 reversible해야 하며, 그 자체로 피로나 확률 누적을 만들어서는 안 된다.
