# 1D Normal-LJ Continuous-Time Energy-Feasibility Result

## Main result

For the active one-dimensional normal-spacing theory, normalization, mean spacing, and LJ configurational energy do **not** by themselves force a tensile crack tail.

The reason is exact: generalized-LJ repulsion diverges as $\lambda\to0^+$, so arbitrarily large energy can be stored in reverse compression while all probability remains at or below the tensile stability limit $\lambda_c$.

A third physical condition is therefore necessary.

The minimal direct condition used in the current theorem is

$$
\boxed{\lambda\ge\lambda_L(t)>0.}
$$

Under this condition and the crack-free condition

$$
\lambda\le\lambda_c,
$$

the LJ potential is convex and the complete admissible energy interval at fixed mean $\mu(t)$ is exactly

$$
\boxed{
\psi(\mu(t))
\le
\mathcal E(t)
\le
\mathcal E_{\rm safe}^{\max}(t),
}
$$

with

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_c).
}
$$

The upper bound is attained by a two-point endpoint distribution, so it is an exact extremum rather than a heuristic estimate.

Define

$$
M_E(t)=\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t).
$$

If a mechanically valid hard compression bound has been established and

$$
M_E(t)<0,
$$

then no crack-free distribution can satisfy normalization, mean, energy, and support simultaneously. Therefore some probability mass must lie beyond $\lambda_c$.

The continuous-time first-passage definition is

$$
\boxed{
\tau_E
=
\inf\{t\ge0:M_E(t)<0\}.
}
$$

This result uses no fitted probability family, cycle-evolution law, damage variable, damping coefficient, or 3D model.

## Numerical parameter sweep

The repository runner evaluates the exact formula for illustrative lower compression bounds. These values are **not Al material inputs**.

At $\mu=1$:

| illustrative $\lambda_L$ | dimensionless $\mathcal E_{\rm safe}^{\max}$ |
| ---: | ---: |
| 0.90 | 0.00704373 |
| 0.95 | 0.00215815 |
| 0.98 | 0.000650058 |
| 0.99 | 0.000296088 |

The strong dependence on the compression bound confirms that the remaining physical problem is not a minor numerical correction.

## Next step

The next active target is to derive or independently constrain

$$
\boxed{\lambda_L(t)}
$$

from the one-dimensional normal-LJ mechanics itself.

---

# 한국어 번역 — 1D Normal-LJ 연속시간 에너지 실현가능성 결과

## 핵심 결과

활성 1차원 수직-spacing 이론에서 정규화, 평균 spacing, LJ configurational energy만으로는 tensile crack tail을 강제로 만들 수 없다.

이유는 정확하다. generalized-LJ repulsion은 $\lambda\to0^+$에서 발산하므로 모든 확률질량을 tensile stability limit $\lambda_c$ 이하에 유지하면서도 reverse compression 쪽에 임의로 큰 에너지를 저장할 수 있다.

따라서 세 번째 물리조건이 필요하다.

현재 theorem에서 사용하는 가장 직접적인 최소조건은

$$
\boxed{\lambda\ge\lambda_L(t)>0}
$$

이다.

이 조건과 crack-free 조건

$$
\lambda\le\lambda_c
$$

아래에서는 LJ potential이 convex이고, 평균 $\mu(t)$가 주어졌을 때 가능한 전체 energy interval은 정확히

$$
\boxed{
\psi(\mu(t))
\le
\mathcal E(t)
\le
\mathcal E_{\rm safe}^{\max}(t)
}
$$

이다.

여기서

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_c)
}
$$

이다.

upper bound는 실제 endpoint two-point distribution에서 달성되므로 heuristic estimate가 아니라 정확한 extremum이다.

$$
M_E(t)=\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t)
$$

를 정의한다.

역학적으로 유효한 hard compression bound가 확보되어 있고

$$
M_E(t)<0
$$

이면 정규화, 평균, 에너지 및 support를 동시에 만족하는 crack-free distribution이 존재할 수 없다. 따라서 일부 확률질량은 반드시 $\lambda_c$를 넘어야 한다.

연속시간 first-passage 정의는

$$
\boxed{
\tau_E
=
\inf\{t\ge0:M_E(t)<0\}
}
$$

이다.

이 결과에는 fitted probability family, cycle-evolution law, damage variable, damping coefficient 또는 3D model이 들어가지 않는다.

## 수치 parameter sweep

repository runner는 illustrative lower compression bound에 대해 정확한 식을 계산한다. 이 값들은 **Al material input이 아니다.**

$\mu=1$에서:

| illustrative $\lambda_L$ | dimensionless $\mathcal E_{\rm safe}^{\max}$ |
| ---: | ---: |
| 0.90 | 0.00704373 |
| 0.95 | 0.00215815 |
| 0.98 | 0.000650058 |
| 0.99 | 0.000296088 |

compression bound에 대한 강한 의존성은 남아 있는 물리문제가 사소한 numerical correction이 아니라는 것을 보여준다.

## 다음 단계

다음 활성 목표는 1차원 normal-LJ mechanics 자체로부터

$$
\boxed{\lambda_L(t)}
$$

를 유도하거나 독립적으로 제약하는 것이다.
