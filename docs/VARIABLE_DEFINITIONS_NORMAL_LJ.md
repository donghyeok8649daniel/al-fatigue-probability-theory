# Variable definitions — normal generalized-LJ chain

This file defines the symbols introduced by the mainline normal-opening generalized-LJ simulation. General project symbols remain in `docs/VARIABLE_DEFINITIONS.md`.

## Normal-chain microscopic variables

| Symbol | Definition | Physical meaning | Unit / scale | Classification |
|---|---|---|---|---|
| $x_i(t)$ | coordinate of atom/site $i$ | 1D normal atomic position | normalized length | DEFINITION |
| $\lambda_i(t)$ | $x_{i+1}-x_i$ | local normal bond spacing divided by equilibrium spacing | dimensionless | DEFINITION |
| $\lambda$ | continuous stretch coordinate | generic argument of the normalized LJ energy | dimensionless | DEFINITION |
| $N_b$ | number of bonds | finite number of local spacing samples | dimensionless | DEFINITION |
| $P_N(\lambda,t)$ | $N_b^{-1}\sum_i\delta(\lambda-\lambda_i)$ | finite empirical normal-spacing density | dimensionless-density convention | DEFINITION |
| $\langle\lambda\rangle$ | mean of $\lambda_i$ | mean normal spacing | dimensionless | DEFINITION |
| $\operatorname{Var}(\lambda)$ | variance of $\lambda_i$ | normal-spacing broadening measure | dimensionless$^2$ | DEFINITION |
| $\lambda_{\max}$ | $\max_i\lambda_i$ | largest instantaneous local normal opening | dimensionless | DEFINITION |

## Normalized generalized-LJ energy

The normalized pair energy is

$$
\boxed{
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}.
}
$$

It is normalized so that

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

| Symbol | Definition | Physical meaning | Classification |
|---|---|---|---|
| $\phi(\lambda)$ | normalized generalized-LJ energy | local normal pair-energy shape used by the reduced chain | DEFINITION / NORMALIZATION |
| $\phi'(\lambda)$ | $d\phi/d\lambda$ | dimensionless tensile generalized force | EXACT under the stated model |
| $\phi''(\lambda)$ | $d^2\phi/d\lambda^2$ | local tangent stiffness | EXACT under the stated model |
| $m$ | repulsive exponent | short-range repulsive shape parameter | EMPIRICAL / previous calibration input |
| $n$ | attractive exponent | attractive-tail shape parameter | EMPIRICAL / previous calibration input |
| $\lambda_c$ | solution of $\phi''(\lambda_c)=0$ | local ideal normal stability-loss stretch | DERIVED |
| $f_c$ | $\phi'(\lambda_c)$ | dimensionless static ideal normal-instability force | DERIVED |

For the current values $m=12.19$, $n=6$,

$$
\lambda_c=1.1077715386,
$$

$$
f_c=0.03703426967.
$$

## External normal loading

| Symbol | Definition | Physical meaning | Unit / scale | Classification |
|---|---|---|---|---|
| $f(t)$ | normalized external end force | cyclic normal generalized force applied to the right boundary | dimensionless | DEFINITION |
| $f_a$ | force amplitude | amplitude of $f(t)$ | dimensionless | DEFINITION |
| $f_m$ | mean force | mean value of $f(t)$ | dimensionless | DEFINITION |
| $\omega^*$ | dimensionless angular frequency | forcing frequency measured on the atomic normalized time scale | dimensionless | DEFINITION |
| $T^*$ | $2\pi/\omega^*$ | dimensionless loading period | dimensionless | DEFINITION |
| $N_{\rm cyc}$ | cycle index | loading cycle count | dimensionless | DEFINITION |

Under the earlier 1D stress mapping,

$$
\boxed{
f=\frac{\sigma_n}{E}.
}
$$

Therefore a 100 MPa amplitude with $E=69$ GPa gives

$$
f_a=1.44927536\times10^{-3}.
$$

## Time-scale variables

The normalized time scale used for dimensional interpretation is

$$
\boxed{
t_0=\sqrt{\frac{M a_0}{EA_0}}.
}
$$

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $t_0$ | $\sqrt{M a_0/(EA_0)}$ | atomic mechanical time scale of the normalized coordinate | s | DERIVED from stated inputs |
| $M$ | atomic mass used in the time mapping | Al atomic mass in the current dimensional estimate | kg | EMPIRICAL INPUT |
| $a_0$ | equilibrium/reference spacing | dimensional length used to map normalized spacing to physical spacing | m | EMPIRICAL / previous calibration input |
| $A_0$ | reference effective area | area used in the earlier 1D stress mapping | m$^2$ | previous calibration output |
| $E$ | Young's modulus | small-strain normal stiffness used for stress normalization | Pa | EMPIRICAL INPUT |
| $f_{\rm phys}$ | physical cyclic frequency | laboratory loading frequency | Hz | DEFINITION |

The mappings are

$$
\boxed{
\omega^*=2\pi f_{\rm phys}t_0,
}
$$

and

$$
\boxed{
f_{\rm phys}=\frac{\omega^*}{2\pi t_0}.
}
$$

For the current dimensional values,

$$
t_0\approx5.55\times10^{-14}\ \mathrm{s}.
$$

Thus $\omega^*=0.02$ is an atomic-scale dynamical test, whereas 20 Hz corresponds to $\omega^*\approx6.97\times10^{-12}$.

## Instability event

The first idealized local normal-instability event is defined by

$$
\boxed{
\tau_{\rm inst}
=
\inf\{t:\max_i\lambda_i(t)\ge\lambda_c\}.
}
$$

This event is a mechanical stability diagnostic. It is **not automatically equal to macroscopic crack initiation, fatigue life, or cumulative failure probability**.

The numerical parameter `runaway_spacing` in the code is only a post-instability integration stop. It is not a physical fracture threshold.

---

# 한국어 번역 — 수직 generalized-LJ chain 변수정의

이 파일은 메인 수직-opening generalized-LJ simulation에서 새로 도입된 기호를 정의한다. 프로젝트 전체 공통기호는 `docs/VARIABLE_DEFINITIONS.md`에서 관리한다.

## 수직 chain 미시변수

| 기호 | 정의 | 물리적 의미 | 단위 / 척도 | 분류 |
|---|---|---|---|---|
| $x_i(t)$ | atom/site $i$의 좌표 | 1D 수직 원자위치 | normalized length | DEFINITION |
| $\lambda_i(t)$ | $x_{i+1}-x_i$ | 평형간격으로 나눈 국부 수직 bond spacing | dimensionless | DEFINITION |
| $\lambda$ | 연속 stretch 좌표 | normalized LJ energy의 일반 argument | dimensionless | DEFINITION |
| $N_b$ | bond 수 | finite local spacing sample 수 | dimensionless | DEFINITION |
| $P_N(\lambda,t)$ | $N_b^{-1}\sum_i\delta(\lambda-\lambda_i)$ | finite empirical normal-spacing density | dimensionless-density convention | DEFINITION |
| $\langle\lambda\rangle$ | $\lambda_i$의 평균 | 평균 수직 spacing | dimensionless | DEFINITION |
| $\operatorname{Var}(\lambda)$ | $\lambda_i$의 분산 | normal-spacing broadening 지표 | dimensionless$^2$ | DEFINITION |
| $\lambda_{\max}$ | $\max_i\lambda_i$ | 순간적으로 가장 큰 국부 수직 opening | dimensionless | DEFINITION |

## normalized generalized-LJ energy

normalized pair energy는

$$
\boxed{
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

이다.

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

이 되도록 normalization했다.

| 기호 | 정의 | 물리적 의미 | 분류 |
|---|---|---|---|
| $\phi(\lambda)$ | normalized generalized-LJ energy | reduced chain에서 사용하는 국부 수직 pair-energy shape | DEFINITION / NORMALIZATION |
| $\phi'(\lambda)$ | $d\phi/d\lambda$ | 무차원 tensile generalized force | stated model 아래 EXACT |
| $\phi''(\lambda)$ | $d^2\phi/d\lambda^2$ | 국부 tangent stiffness | stated model 아래 EXACT |
| $m$ | repulsive exponent | short-range repulsive shape parameter | EMPIRICAL / previous calibration input |
| $n$ | attractive exponent | attractive-tail shape parameter | EMPIRICAL / previous calibration input |
| $\lambda_c$ | $\phi''(\lambda_c)=0$의 해 | 이상적인 국부 normal stability-loss stretch | DERIVED |
| $f_c$ | $\phi'(\lambda_c)$ | 무차원 static ideal normal-instability force | DERIVED |

현재 $m=12.19$, $n=6$이면

$$
\lambda_c=1.1077715386,
$$

$$
f_c=0.03703426967
$$

이다.

## 외부 수직하중

| 기호 | 정의 | 물리적 의미 | 단위 / 척도 | 분류 |
|---|---|---|---|---|
| $f(t)$ | normalized external end force | 오른쪽 경계에 가하는 반복 수직 generalized force | dimensionless | DEFINITION |
| $f_a$ | force amplitude | $f(t)$의 진폭 | dimensionless | DEFINITION |
| $f_m$ | mean force | $f(t)$의 평균 | dimensionless | DEFINITION |
| $\omega^*$ | dimensionless angular frequency | atomic normalized time scale 기준 forcing frequency | dimensionless | DEFINITION |
| $T^*$ | $2\pi/\omega^*$ | 무차원 loading period | dimensionless | DEFINITION |
| $N_{\rm cyc}$ | cycle index | loading cycle count | dimensionless | DEFINITION |

기존 1D stress mapping에서는

$$
\boxed{
f=\frac{\sigma_n}{E}
}
$$

이다.

따라서 $E=69$ GPa에서 100 MPa amplitude는

$$
f_a=1.44927536\times10^{-3}
$$

이다.

## 시간척도 변수

차원복원을 위한 normalized time scale은

$$
\boxed{
t_0=\sqrt{\frac{M a_0}{EA_0}}
}
$$

이다.

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $t_0$ | $\sqrt{M a_0/(EA_0)}$ | normalized coordinate의 atomic mechanical time scale | s | stated input에서 DERIVED |
| $M$ | time mapping에 사용한 atomic mass | 현재 dimensional estimate에서 Al atomic mass | kg | EMPIRICAL INPUT |
| $a_0$ | equilibrium/reference spacing | normalized spacing을 physical spacing으로 바꾸는 길이척도 | m | EMPIRICAL / previous calibration input |
| $A_0$ | reference effective area | 기존 1D stress mapping의 면적 | m$^2$ | previous calibration output |
| $E$ | Young's modulus | stress normalization에 사용하는 small-strain normal stiffness | Pa | EMPIRICAL INPUT |
| $f_{\rm phys}$ | physical cyclic frequency | 실험 loading frequency | Hz | DEFINITION |

mapping은

$$
\boxed{
\omega^*=2\pi f_{\rm phys}t_0
}
$$

및

$$
\boxed{
f_{\rm phys}=\frac{\omega^*}{2\pi t_0}
}
$$

이다.

현재 dimensional value를 사용하면

$$
t_0\approx5.55\times10^{-14}\ \mathrm{s}
$$

이다.

따라서 $\omega^*=0.02$는 atomic-scale dynamics test이고, 20 Hz는 $\omega^*\approx6.97\times10^{-12}$에 해당한다.

## instability event

첫 이상화된 국부 normal-instability event는

$$
\boxed{
\tau_{\rm inst}
=
\inf\{t:\max_i\lambda_i(t)\ge\lambda_c\}
}
$$

로 정의한다.

이 event는 mechanical stability diagnostic이며, **macroscopic crack initiation, fatigue life, cumulative failure probability와 자동으로 동일한 것이 아니다.**

코드의 `runaway_spacing`은 instability 이후 무한히 계산하지 않기 위한 numerical stop일 뿐 physical fracture threshold가 아니다.
