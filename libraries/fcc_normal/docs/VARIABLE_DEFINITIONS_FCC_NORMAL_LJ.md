# Variable Definitions — FCC Normal-LJ Extension

This document defines every symbol newly introduced by the three-dimensional FCC normal-deformation generalized-Lennard-Jones extension.

The existing spacing/probability variables remain defined in `VARIABLE_DEFINITIONS_NORMAL_LJ.md`.

## Classification labels

- **EXACT / IDENTITY** — exact under the explicitly stated microscopic model.
- **DEFINITION** — chosen mathematical definition.
- **ASSUMPTION** — modeling assumption.
- **CONTROLLED APPROXIMATION** — simplification that must be checked.
- **EMPIRICAL INPUT** — externally supplied material quantity.

## 1. FCC geometry

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $a_{\rm lat}$ | conventional FCC cubic lattice constant | reference unit-cell edge length | m | EMPIRICAL INPUT |
| $\Omega_0$ | $a_{\rm lat}^3/4$ | reference atomic volume in FCC | m$^3$/atom | EXACT / IDENTITY under FCC geometry |
| $\mathbf R$ | reference FCC Bravais vector | vector from a reference atom to another lattice site | m | DEFINITION |
| $(i,j,k)$ | integer lattice-vector labels with $i+j+k$ even | FCC parity representation | dimensionless | DEFINITION |
| $\mathbf F$ | homogeneous deformation gradient | maps reference lattice vectors to deformed vectors | dimensionless | DEFINITION |
| $\lambda_n$ | axial / normal stretch | [001] stretch imposed along the active loading direction | dimensionless | DEFINITION |
| $\lambda_t$ | common transverse stretch | lateral contraction/expansion under [001] loading | dimensionless | DEFINITION |
| $R_{\rm cut}$ | spherical lattice-sum cutoff | finite radius used to approximate the infinite FCC sum | multiples of $a_{\rm lat}$ | CONTROLLED APPROXIMATION |

For the active [001] homogeneous normal calculation,

$$
\boxed{\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n).}
$$

The FCC reference vectors are represented by

$$
\boxed{
\mathbf R=\frac{a_{\rm lat}}2(i,j,k),
\qquad i+j+k\ \text{even}.
}
$$

## 2. FCC pair-potential energy

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $U(\mathbf F)$ | $\frac12\sum_{\mathbf R\ne0}v(|\mathbf F\mathbf R|)$ | configurational pair energy per atom under homogeneous deformation | J/atom | EXACT under infinite stated pair model |
| $q$ | $\sigma_{\rm LJ}/a_{\rm lat}$ | dimensionless LJ length ratio fixed by zero-pressure equilibrium | dimensionless | DEFINITION / calibration output |
| $A_m$ | $\frac12\sum_{\mathbf R\ne0}|\mathbf R/a_{\rm lat}|^{-m}$ | repulsive FCC lattice sum | dimensionless | EXACT in infinite-sum limit |
| $A_n$ | $\frac12\sum_{\mathbf R\ne0}|\mathbf R/a_{\rm lat}|^{-n}$ | attractive FCC lattice sum | dimensionless | EXACT in infinite-sum limit |
| $E_{\rm coh}$ | positive energy required to separate the solid into isolated atoms | cohesive / binding-energy magnitude | J/atom or eV/atom | EMPIRICAL INPUT when calibrated |

For the generalized LJ law,

$$
q^{m-n}=\frac{nA_n}{mA_m}
$$

follows from the zero-pressure isotropic equilibrium condition.

## 3. Normal stress and transverse relaxation

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $P_{33}$ | $\Omega_0^{-1}\partial U/\partial\lambda_n$ | axial nominal/engineering stress under the chosen stretch measure | Pa | DEFINITION under current convention |
| $E_{[001]}$ | small-strain slope of relaxed $P_{33}$ versus axial strain | [001] directional Young modulus | Pa | DEFINITION |
| $\nu_{[001]}$ | $-d\lambda_t/d\lambda_n$ at the reference state | [001] transverse contraction ratio | dimensionless | DEFINITION |
| $C_{11},C_{12},C_{44}$ | cubic elastic stiffness constants | second derivatives of energy density with respect to the relevant strain modes | Pa | DEFINITION |

The active free-transverse homogeneous condition is

$$
\boxed{\frac{\partial U}{\partial\lambda_t}=0.}
$$

For a cubic crystal,

$$
\boxed{
E_{[001]}
=\frac{(C_{11}-C_{12})(C_{11}+2C_{12})}{C_{11}+C_{12}}
}
$$

and

$$
\boxed{\nu_{[001]}=\frac{C_{12}}{C_{11}+C_{12}}.}
$$

## 4. Central-pair structural constraint

| Symbol / relation | Meaning | Classification |
|---|---|---|
| $C_{12}=C_{44}$ | Cauchy relation for a cubic central-pair crystal at zero pressure | EXACT structural property of the model class under the stated conditions |
| $\Delta_C=C_{12}-C_{44}$ | numerical diagnostic for recovery of the Cauchy relation | DEFINITION |

A nonzero experimental difference between $C_{12}$ and $C_{44}$ is evidence that a central pair potential is not a complete 3D interatomic energy model.

## 5. Calibration variables

| Symbol / code name | Meaning | Unit | Classification |
|---|---|---|---|
| `epsilon_cohesive` | $\varepsilon_{\rm LJ}$ obtained by fitting $E_{\rm coh}$ | J or eV | calibration output |
| `epsilon_normal` | $\varepsilon_{\rm LJ}$ obtained by fitting $E_{[001]}$ | J or eV | calibration output |
| `sigma_lj_m` | physical $\sigma_{\rm LJ}=q a_{\rm lat}$ | m | calibration output |
| `cutoff_lattice_constants` | finite spherical cutoff radius divided by $a_{\rm lat}$ | dimensionless | numerical parameter |
| `difference_step` | centered finite-difference step for stress derivative | dimensionless stretch | numerical parameter |
| `strain_step` | small strain used for elastic finite differences | dimensionless | numerical parameter |

## 6. Current directional empirical inputs

The current validation script uses the following external reference values:

$$
a_{\rm lat}=4.0495\ \text{Å},
$$

$$
E_{\rm coh}=3.43\ \text{eV/atom},
$$

$$
C_{11}=107\ \text{GPa},\qquad
C_{12}=61\ \text{GPa},\qquad
C_{44}=29\ \text{GPa}.
$$

These imply

$$
E_{[001]}\approx62.7024\ \text{GPa},
$$

$$
\nu_{[001]}\approx0.363095.
$$

They are **EMPIRICAL INPUTS**, not predictions of the current LJ model.

## 7. Current FCC-LJ outputs

For $(m,n)=(12.19,6)$ and the converged reference calculation:

| Quantity | Current value | Meaning |
|---|---:|---|
| $\sigma_{\rm LJ}$ | $2.62721$ Å | zero-pressure LJ length parameter from FCC lattice geometry |
| $\varepsilon_{\rm LJ}^{(E)}$ | $0.445621$ eV | energy scale fitted to $E_{[001]}$ |
| $\varepsilon_{\rm LJ}^{(\rm coh)}$ | $1.56683$ eV | energy scale fitted to cohesive energy |
| $E_{\rm coh}^{(E\text{-fit})}$ | $0.9755$ eV/atom | cohesive-energy prediction after normal elastic fitting |
| $E_{[001]}^{(\rm coh\text{-fit})}$ | $220.466$ GPa | normal modulus prediction after cohesion fitting |
| $\sigma_{\rm ideal}^{(E\text{-fit})}$ | $9.045$ GPa | maximum relaxed [001] engineering stress in the reference calculation |

These numerical values are model outputs, not universal constants of aluminum.

---

# 한국어 번역 — FCC Normal-LJ 확장 변수정의

이 문서는 3차원 FCC 수직변형 generalized-Lennard-Jones 확장에서 새로 도입된 모든 기호를 정의한다.

기존 spacing/probability 변수는 `VARIABLE_DEFINITIONS_NORMAL_LJ.md`에 계속 정의한다.

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 미시모델 아래에서 정확함.
- **DEFINITION** — 선택한 수학적 정의.
- **ASSUMPTION** — 모델링 가정.
- **CONTROLLED APPROXIMATION** — 검증해야 하는 단순화.
- **EMPIRICAL INPUT** — 외부에서 주어진 재료값.

## 1. FCC geometry

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $a_{\rm lat}$ | conventional FCC cubic lattice constant | 기준 unit-cell edge length | m | EMPIRICAL INPUT |
| $\Omega_0$ | $a_{\rm lat}^3/4$ | FCC 원자 하나당 기준부피 | m$^3$/atom | EXACT / IDENTITY under FCC geometry |
| $\mathbf R$ | reference FCC Bravais vector | 기준원자에서 다른 lattice site까지의 vector | m | DEFINITION |
| $(i,j,k)$ | $i+j+k$가 짝수인 integer lattice label | FCC parity representation | dimensionless | DEFINITION |
| $\mathbf F$ | homogeneous deformation gradient | reference lattice vector를 deformed vector로 mapping | dimensionless | DEFINITION |
| $\lambda_n$ | axial / normal stretch | active [001] loading direction의 수직 stretch | dimensionless | DEFINITION |
| $\lambda_t$ | common transverse stretch | [001] loading에 따른 횡수축/팽창 | dimensionless | DEFINITION |
| $R_{\rm cut}$ | spherical lattice-sum cutoff | infinite FCC sum을 유한하게 근사하는 반경 | multiples of $a_{\rm lat}$ | CONTROLLED APPROXIMATION |

활성 [001] 균질 normal calculation에서는

$$
\boxed{\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n)}
$$

을 사용한다.

FCC 기준 vector는

$$
\boxed{
\mathbf R=\frac{a_{\rm lat}}2(i,j,k),
\qquad i+j+k\ \text{even}
}
$$

으로 표현한다.

## 2. FCC pair-potential energy

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $U(\mathbf F)$ | $\frac12\sum_{\mathbf R\ne0}v(|\mathbf F\mathbf R|)$ | homogeneous deformation에서 원자 하나당 configurational pair energy | J/atom | EXACT under infinite stated pair model |
| $q$ | $\sigma_{\rm LJ}/a_{\rm lat}$ | zero-pressure equilibrium으로 결정되는 dimensionless LJ length ratio | dimensionless | DEFINITION / calibration output |
| $A_m$ | $\frac12\sum|\mathbf R/a_{\rm lat}|^{-m}$ | repulsive FCC lattice sum | dimensionless | EXACT in infinite-sum limit |
| $A_n$ | $\frac12\sum|\mathbf R/a_{\rm lat}|^{-n}$ | attractive FCC lattice sum | dimensionless | EXACT in infinite-sum limit |
| $E_{\rm coh}$ | solid를 isolated atom으로 분리하는 데 필요한 positive energy | cohesive/binding-energy magnitude | J/atom or eV/atom | EMPIRICAL INPUT when calibrated |

zero-pressure isotropic equilibrium condition에서

$$
q^{m-n}=\frac{nA_n}{mA_m}
$$

가 나온다.

## 3. Normal stress와 transverse relaxation

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $P_{33}$ | $\Omega_0^{-1}\partial U/\partial\lambda_n$ | 선택한 stretch measure의 axial nominal/engineering stress | Pa | DEFINITION |
| $E_{[001]}$ | relaxed $P_{33}$와 axial strain의 small-strain slope | [001] directional Young modulus | Pa | DEFINITION |
| $\nu_{[001]}$ | reference state에서 $-d\lambda_t/d\lambda_n$ | [001] 횡수축비 | dimensionless | DEFINITION |
| $C_{11},C_{12},C_{44}$ | cubic elastic stiffness constants | 관련 strain mode에 대한 energy-density 2차미분 | Pa | DEFINITION |

활성 free-transverse homogeneous condition은

$$
\boxed{\frac{\partial U}{\partial\lambda_t}=0}
$$

이다.

cubic crystal에서

$$
\boxed{
E_{[001]}
=\frac{(C_{11}-C_{12})(C_{11}+2C_{12})}{C_{11}+C_{12}}
}
$$

이고

$$
\boxed{\nu_{[001]}=\frac{C_{12}}{C_{11}+C_{12}}}
$$

이다.

## 4. Central-pair structural constraint

| Symbol / relation | Meaning | Classification |
|---|---|---|
| $C_{12}=C_{44}$ | zero pressure cubic central-pair crystal의 Cauchy relation | stated condition에서 model class의 EXACT structural property |
| $\Delta_C=C_{12}-C_{44}$ | Cauchy relation이 수치적으로 복원되는지 확인하는 diagnostic | DEFINITION |

experimental $C_{12}$와 $C_{44}$의 차이가 크다는 것은 central pair potential이 complete 3D interatomic energy가 아님을 뜻한다.

## 5. Calibration 변수

| Symbol / code name | Meaning | Unit | Classification |
|---|---|---|---|
| `epsilon_cohesive` | $E_{\rm coh}$를 맞춘 $\varepsilon_{\rm LJ}$ | J or eV | calibration output |
| `epsilon_normal` | $E_{[001]}$를 맞춘 $\varepsilon_{\rm LJ}$ | J or eV | calibration output |
| `sigma_lj_m` | physical $\sigma_{\rm LJ}=q a_{\rm lat}$ | m | calibration output |
| `cutoff_lattice_constants` | spherical cutoff radius / $a_{\rm lat}$ | dimensionless | numerical parameter |
| `difference_step` | stress derivative용 centered finite-difference step | dimensionless stretch | numerical parameter |
| `strain_step` | elastic finite difference에 사용하는 small strain | dimensionless | numerical parameter |

## 6. 현재 방향성 empirical input

현재 validation script는

$$
a_{\rm lat}=4.0495\ \text{Å},
$$

$$
E_{\rm coh}=3.43\ \text{eV/atom},
$$

$$
C_{11}=107\ \text{GPa},\qquad
C_{12}=61\ \text{GPa},\qquad
C_{44}=29\ \text{GPa}
$$

를 외부 reference로 사용한다.

이 값으로부터

$$
E_{[001]}\approx62.7024\ \text{GPa},
$$

$$
\nu_{[001]}\approx0.363095
$$

가 나온다.

이들은 현재 LJ model의 prediction이 아니라 **EMPIRICAL INPUT**이다.

## 7. 현재 FCC-LJ output

$(m,n)=(12.19,6)$와 converged reference calculation에서

| Quantity | Current value | Meaning |
|---|---:|---|
| $\sigma_{\rm LJ}$ | $2.62721$ Å | FCC lattice geometry에서 나온 zero-pressure LJ length parameter |
| $\varepsilon_{\rm LJ}^{(E)}$ | $0.445621$ eV | $E_{[001]}$에 맞춘 energy scale |
| $\varepsilon_{\rm LJ}^{(\rm coh)}$ | $1.56683$ eV | cohesive energy에 맞춘 energy scale |
| $E_{\rm coh}^{(E\text{-fit})}$ | $0.9755$ eV/atom | normal elastic fitting 후 cohesive-energy prediction |
| $E_{[001]}^{(\rm coh\text{-fit})}$ | $220.466$ GPa | cohesion fitting 후 normal modulus prediction |
| $\sigma_{\rm ideal}^{(E\text{-fit})}$ | $9.045$ GPa | reference calculation의 relaxed [001] maximum engineering stress |

이 수치들은 Al의 universal constant가 아니라 현재 model output이다.
