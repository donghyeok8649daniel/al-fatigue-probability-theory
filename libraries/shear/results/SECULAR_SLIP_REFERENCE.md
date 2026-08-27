# Secular Slip Reference Run

This file records the current nondimensional proof-of-principle reference values for the conservative nonlinear slip-bath model.

## Parameters

$$
M=m=k=k_c=b=1,
$$

$$
\Delta_\gamma=0.1,\qquad \omega=0.2,
$$

with a smooth two-cycle force ramp and no viscous damping.

The reference direct integration used 800 bath sites and velocity-Verlet integration with

$$
\Delta t=0.01.
$$

## Cycle-end states

### $F_a=0.34$

Last six cycle-end values:

$$
-0.0240,\,-0.0239,\,-0.0240,\,-0.0239,\,-0.0240,\,-0.0240.
$$

Interpretation: bounded intrawell response; no secular cycle drift.

### $F_a=0.40$

Last six cycle-end values:

$$
-1.9650,\,-1.9648,\,-1.9650,\,-1.9649,\,-1.9650,\,-1.9649.
$$

Interpretation: finite transient interwell relocation followed by a periodic state.

### $F_a=0.50$

Last six cycle-end values:

$$
-5.8529,\,-6.8542,\,-7.8523,\,-8.8538,\,-9.8519,\,-10.8534.
$$

The late-cycle increment is approximately

$$
\boxed{\Delta s_{\rm cycle}\approx-1.00.}
$$

Representative late-cycle work is approximately

$$
\boxed{\oint F\,ds\approx2.994}
$$

in the nondimensional units of this model.

## Energy balance

For the $F_a=0.50$ reference run, the relative error in

$$
E_{\rm int}(t)-E_{\rm int}(0)=\int_0^tF\dot s\,dt
$$

was approximately

$$
\boxed{1.8\times10^{-7}}.
$$

## Spacing-like distribution diagnostic

For the finite bath, define local relative-displacement samples

$$
q_0=s-u_1,
\qquad
q_j=u_{j+1}-u_j.
$$

In the $F_a=0.50$ run, the variance over these samples increased from approximately

$$
1.88\times10^{-3}
$$

at cycle 2 to

$$
3.53\times10^{-2}
$$

at cycle 11.

This demonstrates redistribution of deformation into unresolved lattice modes. **It must not yet be interpreted as a thermodynamic-limit fatigue prediction for $P(a,t)$**, because the statistic depends on the finite observation domain and contains propagating phonon strain as well as structural change.

## Interpretation boundary

The strongest valid claim from this reference run is:

$$
\boxed{
\text{conservative microscopic dynamics}
\;\Rightarrow\;
\text{hysteresis + inter-basin cycle-state change is possible}
}
$$

for a nonlinear periodic non-affine coordinate coupled to a lattice bath.

The run does **not** establish an Al S–N curve, a crack-initiation life, or a low-stress fatigue threshold.

---

# 한국어 번역 — 장기 slip 기준 계산

이 파일에는 보존적인 비선형 slip-bath 모델의 현재 무차원 원리증명 기준값을 기록한다.

## 파라미터

$$
M=m=k=k_c=b=1,
$$

$$
\Delta_\gamma=0.1,\qquad \omega=0.2
$$

를 사용했고, 외력은 두 cycle 동안 매끄럽게 ramp시켰다. 전체계에는 점성 damping을 넣지 않았다.

기준 직접적분에서는 bath site 800개와 velocity-Verlet 적분을 사용했으며

$$
\Delta t=0.01
$$

이다.

## Cycle 끝의 상태

### $F_a=0.34$

마지막 여섯 cycle의 끝값은

$$
-0.0240,\,-0.0239,\,-0.0240,\,-0.0239,\,-0.0240,\,-0.0240
$$

이다.

해석: 하나의 well 내부에 갇힌 bounded response이며 cycle이 증가해도 secular drift가 없다.

### $F_a=0.40$

마지막 여섯 cycle의 끝값은

$$
-1.9650,\,-1.9648,\,-1.9650,\,-1.9649,\,-1.9650,\,-1.9649
$$

이다.

해석: 초기에는 유한한 interwell relocation이 발생하지만 이후 새로운 periodic state에 도달한다.

### $F_a=0.50$

마지막 여섯 cycle의 끝값은

$$
-5.8529,\,-6.8542,\,-7.8523,\,-8.8538,\,-9.8519,\,-10.8534
$$

이다.

후반 cycle의 증가량은 대략

$$
\boxed{\Delta s_{\rm cycle}\approx-1.00}
$$

이다.

후반 cycle에서 대표적인 한 cycle의 일은 이 모델의 무차원 단위로 대략

$$
\boxed{\oint F\,ds\approx2.994}
$$

이다.

## 에너지 수지

$F_a=0.50$ 기준 계산에서

$$
E_{\rm int}(t)-E_{\rm int}(0)=\int_0^tF\dot s\,dt
$$

관계의 상대오차는 약

$$
\boxed{1.8\times10^{-7}}
$$

이었다.

## spacing과 유사한 분포 진단량

유한 bath에서 다음과 같은 국부 상대변위 sample을 정의한다.

$$
q_0=s-u_1,
\qquad
q_j=u_{j+1}-u_j.
$$

$F_a=0.50$ 계산에서 이 sample들의 분산은 cycle 2의 약

$$
1.88\times10^{-3}
$$

에서 cycle 11의

$$
3.53\times10^{-2}
$$

까지 증가했다.

이는 deformation이 unresolved lattice mode로 재분배되고 있음을 보여준다. 그러나 **이를 아직 열역학적 극한의 $P(a,t)$에 대한 피로 예측으로 해석해서는 안 된다.** 이 통계량은 유한한 관찰영역에 의존하며 구조변화뿐 아니라 propagating phonon strain도 포함하기 때문이다.

## 해석의 경계

이 기준 계산으로부터 할 수 있는 가장 강한 주장은

$$
\boxed{
\text{conservative microscopic dynamics}
\;\Rightarrow\;
\text{hysteresis + inter-basin cycle-state change is possible}
}
$$

이라는 것이다. 즉 비선형 주기형 non-affine coordinate가 lattice bath에 연결된 경우 보존적인 미시역학에서도 히스테리시스와 basin 간 cycle-state 변화가 가능하다.

하지만 이 계산은 Al의 S–N curve, crack-initiation life, 또는 low-stress fatigue threshold를 확립한 것이 아니다.
