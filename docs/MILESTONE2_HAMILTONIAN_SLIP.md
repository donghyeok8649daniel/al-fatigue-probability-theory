# Auxiliary Research Note — Hamiltonian non-affine slip proof of cycle-state evolution

## Status

**Auxiliary proof-of-principle only. This is not the main physical mechanism of the project.**

The main project direction is cyclic normal stress, normal interatomic spacing $a_i(t)$, the state density $P(a,t)$, and normal-opening crack initiation. See `MILESTONE2_NORMAL_DEFORMATION.md`.

This slip model is retained because it answers one narrower mathematical question:

> Can a fully conservative microscopic system with a nonlinear structural coordinate produce cycle-to-cycle state evolution without an empirical damage law?

For this auxiliary model, the answer is yes.

## 1. Auxiliary microscopic model

Introduce one resolved non-affine slip coordinate $s$ coupled to a long harmonic lattice bath with coordinates $u_j$:

$$
H(t)=\frac{P_s^2}{2M}+V_\gamma(s)
+\frac{k_c}{2}(s-u_1)^2
+\sum_{j=1}^{\infty}\frac{p_j^2}{2m}
+\frac{k}{2}\sum_{j=1}^{\infty}(u_{j+1}-u_j)^2
-F(t)s.
$$

The full system contains no viscous damping term.

The current periodic slip landscape is

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right].
$$

This one-harmonic landscape is a **CONTROLLED APPROXIMATION** and is not used as the main energy model for normal fatigue.

## 2. Exact energy balance under the stated model

Define the internal energy excluding external loading potential. The equations of motion give

$$
\boxed{
\frac{dE_{\rm int}}{dt}=F(t)\dot s(t).
}
$$

Therefore

$$
\boxed{
A_H=\oint F\,ds
}
$$

is transferred into unresolved lattice modes and/or retained in a changed structural state rather than being introduced through phenomenological damping.

## 3. Auxiliary cycle-map result

At cycle endpoints $t_N=NT$, define

$$
s_N=s(NT).
$$

The nondimensional reference simulation gives three regimes:

| $F_a$ | Long-time behavior | Interpretation |
|---:|---|---|
| 0.34 | $s_N\approx-0.024$ | bounded intra-basin periodic response |
| 0.40 | finite relocation, then periodic | transient structural relocation |
| 0.50 | approximately one period of drift per cycle | running inter-basin state |

For $F_a=0.50$,

$$
\boxed{
s_{N+1}-s_N\approx-1.
}
$$

The global energy-balance relative error is approximately

$$
\boxed{
1.8\times10^{-7}.
}
$$

Thus the running state is not explained by numerical energy loss in the reference calculation.

## 4. What this result means

The strongest valid conclusion is

$$
\boxed{
\text{conservative microscopic dynamics}
\rightarrow
\text{nonlinear cycle-state evolution is possible}.
}
$$

This supports the search for an analogous mechanism in the **normal-spacing sector**.

It does **not** establish that shear slip is the dominant fatigue mechanism in this project.

It does **not** replace the main target

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

under cyclic normal stress.

## 5. Why this model is kept

The model remains useful as a falsification/reference problem because it demonstrates that one does not need to insert an empirical damage evolution law merely to obtain a nontrivial cycle map.

Future mainline work should translate this lesson into a normal-deformation model based on fixed interatomic mechanics, preferably using the generalized Lennard-Jones pair-potential baseline and exact spacing/correlation dynamics.

---

# 한국어 번역 — Hamiltonian 비아핀 slip 보조 원리증명

## 상태

**이 문서는 보조적인 원리증명이다. 프로젝트의 주 물리메커니즘이 아니다.**

프로젝트의 메인 방향은 반복 수직응력, 수직 원자간격 $a_i(t)$, 상태밀도 $P(a,t)$, 그리고 수직 opening에 의한 균열개시다. 메인 이론은 `MILESTONE2_NORMAL_DEFORMATION.md`를 따른다.

이 slip 모델을 유지하는 이유는 다음의 좁은 수학적 질문에 답하기 위해서다.

> 경험적 damage law 없이 보존적인 미시계와 비선형 구조좌표만으로 cycle-to-cycle 상태진화가 가능한가?

이 보조모델에서는 가능했다.

## 1. 보조 미시모델

하나의 non-affine slip coordinate $s$를 긴 harmonic lattice bath $u_j$에 결합한다.

$$
H(t)=\frac{P_s^2}{2M}+V_\gamma(s)
+\frac{k_c}{2}(s-u_1)^2
+\sum_{j=1}^{\infty}\frac{p_j^2}{2m}
+\frac{k}{2}\sum_{j=1}^{\infty}(u_{j+1}-u_j)^2
-F(t)s.
$$

전체계에는 viscous damping term이 없다.

현재 주기 slip landscape는

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right]
$$

이다.

이 one-harmonic landscape는 **CONTROLLED APPROXIMATION**이며 normal fatigue의 메인 에너지모델로 사용하지 않는다.

## 2. 명시된 모델에서의 정확한 에너지수지

외부 loading potential을 제외한 내부에너지에 대해 운동방정식에서

$$
\boxed{
\frac{dE_{\rm int}}{dt}=F(t)\dot s(t)
}
$$

가 나온다.

따라서

$$
\boxed{
A_H=\oint F\,ds
}
$$

는 phenomenological damping으로 넣은 에너지가 아니라 unresolved lattice mode로 전달되거나 바뀐 구조상태에 남는 에너지다.

## 3. 보조 cycle-map 결과

cycle 끝 $t_N=NT$에서

$$
s_N=s(NT)
$$

를 정의한다.

무차원 기준 simulation에서는 세 영역이 나왔다.

| $F_a$ | 장기 거동 | 의미 |
|---:|---|---|
| 0.34 | $s_N\approx-0.024$ | 하나의 basin 안에서 주기응답 |
| 0.40 | 유한 relocation 후 주기상태 | transient 구조이동 |
| 0.50 | cycle마다 약 한 period씩 drift | running inter-basin state |

$F_a=0.50$에서는

$$
\boxed{
s_{N+1}-s_N\approx-1
}
$$

이고 전체 energy-balance 상대오차는 약

$$
\boxed{
1.8\times10^{-7}
}
$$

이다.

따라서 기준 계산에서 running state를 numerical energy loss로 설명하기 어렵다.

## 4. 이 결과의 의미

가장 강하게 말할 수 있는 결론은

$$
\boxed{
\text{보존적인 미시역학}
\rightarrow
\text{비선형 cycle-state evolution 가능}
}
$$

이라는 것이다.

이 결과는 **normal-spacing sector**에서 유사한 메커니즘을 찾을 수 있다는 가능성을 보여주는 보조증거다.

하지만 shear slip이 이 프로젝트의 주 피로메커니즘이라는 뜻은 아니다.

그리고 반복 수직응력에서의 메인 목표

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

를 해결한 것도 아니다.

## 5. 이 모델을 남기는 이유

경험적 damage evolution law를 넣지 않아도 nontrivial cycle map이 가능하다는 것을 보여주는 falsification/reference problem으로 가치가 있기 때문에 유지한다.

향후 메인 연구는 이 교훈을 고정된 interatomic mechanics, generalized Lennard-Jones pair-potential baseline, 그리고 정확한 spacing/correlation dynamics를 이용해 수직변형 모델로 옮겨야 한다.
