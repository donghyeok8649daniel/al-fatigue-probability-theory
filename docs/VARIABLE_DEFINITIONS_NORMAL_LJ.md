# Variable Definitions — Active Normal-LJ Mainline

This document is the active variable dictionary for the normal-deformation theory. Earlier shear/non-affine variables are preserved separately under `libraries/shear/docs/VARIABLE_DEFINITIONS.md`.

## Classifications

- **EXACT / IDENTITY** — exact under the stated microscopic model.
- **DEFINITION** — chosen mathematical definition.
- **ASSUMPTION** — modeling assumption.
- **CONTROLLED APPROXIMATION** — simplification that must be tested.
- **EMPIRICAL INPUT** — measured or externally supplied material quantity.

## 1. Normal microscopic geometry

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $x_i(t)$ | atomic coordinate along the loading axis | normal atomic position in the reduced chain | m | DEFINITION |
| $a_i(t)$ | $x_{i+1}-x_i$ | local normal interatomic spacing | m | DEFINITION |
| $a_0$ | reference spacing | equilibrium/reference normal spacing | m | EMPIRICAL INPUT or calibration output |
| $a$ | continuous spacing coordinate | argument of $P(a,t)$ | m | DEFINITION |
| $\lambda_i$ | $a_i/a_0$ | local normalized normal spacing | dimensionless | DEFINITION |
| $\lambda$ | $a/a_0$ | continuous normalized spacing coordinate | dimensionless | DEFINITION |
| $\lambda_c$ | solution of $\phi''(\lambda_c)=0$ | local idealized normal stability-loss stretch | dimensionless | EXACT under stated LJ model |
| $a_c$ | $a_0\lambda_c$ | corresponding critical normal spacing | m | DEFINITION once $\lambda_c$ is fixed |

## 2. Normal spacing distribution

For finite $N$,

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta(a-a_i(t)).
$$

The thermodynamic-limit state is

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t)}.
$$

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $P_N(a,t)$ | finite empirical spacing density | finite-system normal-spacing distribution | m$^{-1}$ |
| $P(a,t)$ | thermodynamic-limit spacing density | central structural state | m$^{-1}$ |
| $\bar a$ | $\int aP(a,t)\,da$ | mean normal spacing | m |
| $\operatorname{Var}(a)$ | $\int(a-\bar a)^2P\,da$ | normal-spacing variance | m$^2$ |
| $Q_c(t)$ | $\int_{a_c}^{\infty}P(a,t)\,da$ | instantaneous unstable-tail occupancy | dimensionless |

Normalization:

$$
\int P(a,t)\,da=1.
$$

## 3. Spacing-space dynamics

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $v(a,t)$ | $\langle\dot a_i\mid a_i=a\rangle$ | conditional transport velocity in spacing space | m/s |
| $c$ | $\dot a$ | spacing velocity | m/s |
| $F(a,c,t)$ | joint density in spacing and spacing velocity | phase-space lift | s/m$^2$ |
| $A(a,c,t)$ | $\langle\ddot a_i\mid a_i=a,\dot a_i=c\rangle$ | conditional spacing acceleration | m/s$^2$ |

Exact kinematic identity:

$$
\boxed{\partial_tP+\partial_a(Pv)=0.}
$$

Phase-space form:

$$
\partial_tF+\partial_a(cF)+\partial_c(AF)=0.
$$

## 4. Normal loading variables

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $\sigma_n(t)$ | applied normal stress | cyclic stress along active loading axis | Pa |
| $\sigma_m$ | mean normal stress | cycle mean | Pa |
| $\sigma_a$ | normal stress amplitude | cyclic amplitude | Pa |
| $\epsilon_n(t)$ | normal strain | work-conjugate normal strain | dimensionless |
| $f$ | cyclic frequency | laboratory cycles per second | Hz |
| $\omega$ | $2\pi f$ | angular frequency | rad/s |
| $T$ | $2\pi/\omega$ | loading period | s |
| $N_{\rm cyc}$ | cycle index/count | number of completed cycles | dimensionless |

For a sinusoidal normal test,

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t).
$$

## 5. Generalized Lennard-Jones baseline

The active pair potential is

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $v(r)$ | pair interaction energy | microscopic normal interaction baseline | J |
| $\varepsilon_{\rm LJ}$ | energy scale | LJ energy parameter | J or eV |
| $\sigma_{\rm LJ}$ | length scale | LJ distance parameter | m |
| $m,n$ | generalized LJ exponents | repulsive/attractive shape exponents | dimensionless |
| $\phi(\lambda)$ | normalized LJ potential used in the active chain | dimensionless energy function | dimensionless |
| $\phi'(\lambda)$ | derivative of $\phi$ | normalized tensile force/stress coordinate | dimensionless |
| $\phi''(\lambda)$ | second derivative | local tangent stiffness | dimensionless |

The current normalized form is

$$
\boxed{
\phi(\lambda)=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

so that

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

The current calibrated exponents are

$$
m=12.19,
\qquad
n=6.
$$

## 6. Chain energy and force

For the nearest-neighbor active normal chain,

$$
V=\sum_i\phi(\lambda_i).
$$

The force follows by differentiation; no fatigue-force law is inserted independently.

The local idealized stability condition is

$$
\boxed{\phi''(\lambda_c)=0.}
$$

For the current exponents,

$$
\boxed{\lambda_c\approx1.1077715386.}
$$

## 7. Energy accounting

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $E_{\rm int}$ | kinetic + LJ configurational energy | internal energy in direct chain simulation | J or nondimensional |
| $W_{\rm ext}$ | time-integrated external power | external mechanical work | J or nondimensional |
| $A_H$ | $\oint\sigma_n\,d\epsilon_n$ | normal stress-strain hysteresis work per cycle | J/m$^3$ |

A numerical model must satisfy the appropriate work-energy balance to within time-discretization error.

## 8. Pair-distance hierarchy

For more than nearest-neighbor structure,

$$
R_i^{(k)}=a_i+a_{i+1}+\cdots+a_{i+k-1}.
$$

If $P_k(r,t)$ is the density of $R_i^{(k)}$,

$$
\mathcal U(t)=\sum_{k=1}^{\infty}\int v(r)P_k(r,t)\,dr.
$$

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $R_i^{(k)}$ | $k$-th neighbor distance | correlated normal pair distance | m |
| $P_k(r,t)$ | density of $R_i^{(k)}$ | exact pair-distance hierarchy entry | m$^{-1}$ |
| $P^{*k}$ | $k$-fold convolution of $P$ | independence approximation to $P_k$ | m$^{-1}$ |
| $\mathcal U$ | distribution-level configurational energy | pair-energy functional | J per atom under chosen convention |

The replacement $P_k\approx P^{*k}$ is a **CONTROLLED APPROXIMATION**, not an identity.

## 9. Time-scale variables

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $t_0$ | $\sqrt{m_{\rm Al}a_0/(EA_0)}$ in the current 1D scaling | atomic mechanical time scale | s |
| $\omega^*$ | $\omega t_0$ | dimensionless angular frequency | dimensionless |
| $m_{\rm Al}$ | atomic mass of Al | inertial scale | kg |
| $E$ | reference Young's modulus used in current mapping | stress scale | Pa |
| $A_0$ | effective 1D-to-stress area | force/stress conversion scale | m$^2$ |

A laboratory-frequency claim must explicitly map between $\omega$ and $\omega^*$.

## 10. Crack-initiation variables

| Symbol | Definition | Meaning | Unit |
|---|---|---|---|
| $\tau_c$ | first time a specified normal mechanical-instability condition occurs | crack-initiation first-passage time | s |
| $S(t)$ | $P(\tau_c>t)$ | survival probability | dimensionless |
| $F_{\rm ci}(t)$ | $1-S(t)$ | cumulative initiation probability | dimensionless |
| $h(t)$ | $-\dot S/S$ | initiation hazard if defined | s$^{-1}$ |

The active theory does not identify $Q_c(t)$ with $F_{\rm ci}(t)$ without a first-passage argument.

---

# 한국어 번역 — 활성 Normal-LJ 변수정의

이 문서는 수직변형 이론의 활성 변수사전이다. 과거 shear/non-affine 변수는 `libraries/shear/docs/VARIABLE_DEFINITIONS.md`에 별도로 보존한다.

## 분류

- **EXACT / IDENTITY** — 명시된 미시모델 아래 정확함.
- **DEFINITION** — 선택한 수학적 정의.
- **ASSUMPTION** — 모델링 가정.
- **CONTROLLED APPROXIMATION** — 검증해야 하는 축약.
- **EMPIRICAL INPUT** — 측정 또는 외부에서 공급한 재료량.

## 1. 수직 미시기하

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $x_i(t)$ | loading axis 방향 atomic coordinate | 축약사슬의 수직 원자위치 | m | DEFINITION |
| $a_i(t)$ | $x_{i+1}-x_i$ | 국부 수직 원자간격 | m | DEFINITION |
| $a_0$ | 기준 spacing | 평형/기준 수직 원자간격 | m | EMPIRICAL INPUT 또는 calibration output |
| $a$ | 연속 spacing coordinate | $P(a,t)$의 argument | m | DEFINITION |
| $\lambda_i$ | $a_i/a_0$ | 국부 normalized normal spacing | dimensionless | DEFINITION |
| $\lambda$ | $a/a_0$ | 연속 normalized spacing coordinate | dimensionless | DEFINITION |
| $\lambda_c$ | $\phi''(\lambda_c)=0$의 해 | 이상화된 국부 normal stability-loss stretch | dimensionless | stated LJ model 아래 EXACT |
| $a_c$ | $a_0\lambda_c$ | 대응 critical normal spacing | m | $\lambda_c$가 정해진 뒤 DEFINITION |

## 2. 수직 spacing distribution

유한 $N$에서

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta(a-a_i(t)).
$$

열역학적 극한 상태는

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t)}.
$$

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $P_N(a,t)$ | finite empirical spacing density | 유한계의 normal-spacing distribution | m$^{-1}$ |
| $P(a,t)$ | thermodynamic-limit spacing density | 중심 구조상태 | m$^{-1}$ |
| $\bar a$ | $\int aP(a,t)\,da$ | 평균 수직 spacing | m |
| $\operatorname{Var}(a)$ | $\int(a-\bar a)^2P\,da$ | 수직 spacing variance | m$^2$ |
| $Q_c(t)$ | $\int_{a_c}^{\infty}P(a,t)\,da$ | 순간 unstable-tail occupancy | dimensionless |

정규화는

$$
\int P(a,t)\,da=1
$$

이다.

## 3. spacing-space dynamics

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $v(a,t)$ | $\langle\dot a_i\mid a_i=a\rangle$ | spacing space의 conditional transport velocity | m/s |
| $c$ | $\dot a$ | spacing velocity | m/s |
| $F(a,c,t)$ | spacing과 spacing velocity의 joint density | phase-space lift | s/m$^2$ |
| $A(a,c,t)$ | $\langle\ddot a_i\mid a_i=a,\dot a_i=c\rangle$ | conditional spacing acceleration | m/s$^2$ |

정확한 운동학적 항등식은

$$
\boxed{\partial_tP+\partial_a(Pv)=0.}
$$

이고 phase-space form은

$$
\partial_tF+\partial_a(cF)+\partial_c(AF)=0
$$

이다.

## 4. 수직 loading 변수

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $\sigma_n(t)$ | applied normal stress | 활성 loading axis 방향 반복응력 | Pa |
| $\sigma_m$ | mean normal stress | cycle mean | Pa |
| $\sigma_a$ | normal stress amplitude | 반복응력 amplitude | Pa |
| $\epsilon_n(t)$ | normal strain | work-conjugate normal strain | dimensionless |
| $f$ | cyclic frequency | 실험 cycle per second | Hz |
| $\omega$ | $2\pi f$ | angular frequency | rad/s |
| $T$ | $2\pi/\omega$ | loading period | s |
| $N_{\rm cyc}$ | cycle index/count | 완료된 cycle 수 | dimensionless |

사인 수직시험에서는

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t)
$$

이다.

## 5. generalized Lennard-Jones baseline

활성 pair potential은

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
$$

이다.

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $v(r)$ | pair interaction energy | 미시 normal interaction baseline | J |
| $\varepsilon_{\rm LJ}$ | energy scale | LJ energy parameter | J 또는 eV |
| $\sigma_{\rm LJ}$ | length scale | LJ distance parameter | m |
| $m,n$ | generalized LJ exponents | repulsive/attractive shape exponent | dimensionless |
| $\phi(\lambda)$ | 활성 chain의 normalized LJ potential | dimensionless energy function | dimensionless |
| $\phi'(\lambda)$ | $\phi$의 derivative | normalized tensile force/stress coordinate | dimensionless |
| $\phi''(\lambda)$ | second derivative | local tangent stiffness | dimensionless |

현재 normalized form은

$$
\boxed{
\phi(\lambda)=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

이고

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

을 만족한다.

현재 지수는

$$
m=12.19,
\qquad
n=6
$$

이다.

## 6. 사슬 에너지와 force

최근접 이웃 활성 normal chain은

$$
V=\sum_i\phi(\lambda_i)
$$

이다. force는 미분으로 나오며 별도의 fatigue-force law를 삽입하지 않는다.

국부 이상화 stability condition은

$$
\boxed{\phi''(\lambda_c)=0}
$$

이고 현재 지수에서

$$
\boxed{\lambda_c\approx1.1077715386}
$$

이다.

## 7. 에너지 수지

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $E_{\rm int}$ | kinetic + LJ configurational energy | direct chain simulation의 내부에너지 | J 또는 nondimensional |
| $W_{\rm ext}$ | 외부 power의 시간적분 | 외부 기계적 일 | J 또는 nondimensional |
| $A_H$ | $\oint\sigma_n\,d\epsilon_n$ | cycle당 normal stress-strain hysteresis work | J/m$^3$ |

수치모델은 time-discretization error 범위에서 적절한 work-energy balance를 만족해야 한다.

## 8. pair-distance hierarchy

최근접 이웃을 넘어가면

$$
R_i^{(k)}=a_i+a_{i+1}+\cdots+a_{i+k-1}
$$

이고 $P_k(r,t)$를 $R_i^{(k)}$의 density라 하면

$$
\mathcal U(t)=\sum_{k=1}^{\infty}\int v(r)P_k(r,t)\,dr
$$

이다.

$P_k\approx P^{*k}$는 **CONTROLLED APPROXIMATION**이지 항등식이 아니다.

## 9. 시간척도 변수

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $t_0$ | 현재 1D scaling에서 $\sqrt{m_{\rm Al}a_0/(EA_0)}$ | atomic mechanical time scale | s |
| $\omega^*$ | $\omega t_0$ | dimensionless angular frequency | dimensionless |
| $m_{\rm Al}$ | Al atomic mass | inertial scale | kg |
| $E$ | 현재 mapping의 reference Young's modulus | stress scale | Pa |
| $A_0$ | effective 1D-to-stress area | force/stress conversion scale | m$^2$ |

실험주파수를 주장하려면 $\omega$와 $\omega^*$ mapping을 명시해야 한다.

## 10. 균열개시 변수

| 기호 | 정의 | 의미 | 단위 |
|---|---|---|---|
| $\tau_c$ | 정해진 normal mechanical-instability condition의 first time | crack-initiation first-passage time | s |
| $S(t)$ | $P(\tau_c>t)$ | survival probability | dimensionless |
| $F_{\rm ci}(t)$ | $1-S(t)$ | cumulative initiation probability | dimensionless |
| $h(t)$ | $-\dot S/S$ | 정의 가능한 경우 initiation hazard | s$^{-1}$ |

first-passage 논증 없이 $Q_c(t)$와 $F_{\rm ci}(t)$를 동일시하지 않는다.
