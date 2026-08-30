# Milestone 16 — Probability, Energy Hysteresis, and 1D FEM Coupling

## Status

**CANDIDATE KINETIC EXTENSION / NUMERICALLY VERIFIED INTERFACE — NOT YET A CALIBRATED ALUMINUM FATIGUE-LIFE LAW**

This milestone connects four quantities that were previously discussed separately:

1. normalization of the layer-spacing density;
2. mean layer spacing;
3. mean layer-interaction energy;
4. mechanical energy hysteresis.

It also attaches this state to every element of the existing one-dimensional tensile FEM without changing the validated C continuum solver.

## 1. Common reduced variables

Let

$$
\lambda=\frac{a}{a_0},
\qquad
U(a)=E_0\phi(\lambda),
\qquad
E_0=EA_0a_0,
$$

and use the calibrated generalized-LJ layer potential

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}.
$$

The stress carried by the representative layer patch produces

$$
F=A_0\sigma,
\qquad
f=\frac{Fa_0}{E_0}=\frac{\sigma}{E}.
$$

Therefore the reduced force-biased energy is

$$
\boxed{w(\lambda,t)=\phi(\lambda)-f(t)\lambda.}
$$

## 2. Four coupled equations

### 2.1 Normalization

For the current conditional intact density,

$$
\boxed{
\int_{\lambda_{\min}}^{\lambda_{\max}}
p(\lambda,t)\,d\lambda=1.
}
$$

The finite computational interval is part of the intact-basin approximation. It must not be hidden as a global tensile Gibbs equilibrium.

### 2.2 Mean spacing

$$
\boxed{
\bar\lambda(t)=\int \lambda p(\lambda,t)\,d\lambda,
\qquad
\bar a(t)=a_0\bar\lambda(t).
}
$$

The spacing strain used in the energy loop is

$$
\bar\epsilon_a(t)=\bar\lambda(t)-1.
$$

### 2.3 Mean layer energy

With the equilibrium energy shifted to zero,

$$
\psi(\lambda)=\phi(\lambda)-\phi(1),
$$

the mean interaction energy per represented patch is

$$
\boxed{
\bar U(t)=E_0\int\psi(\lambda)p(\lambda,t)\,d\lambda.
}
$$

Because the patch reference volume is $A_0a_0$ and $E_0=EA_0a_0$, the corresponding energy density is

$$
\boxed{
u_{\rm LJ}(t)
=E\int\psi(\lambda)p(\lambda,t)\,d\lambda.
}
$$

For diagnostics, the code also evaluates the reduced nonequilibrium free energy

$$
\mathcal G[p;f]
=
\int[\psi(\lambda)-f(\lambda-1)]p\,d\lambda
+\chi^{-1}\int p\ln p\,d\lambda.
$$

### 2.4 Energy hysteresis

The mechanical work density accumulated along the mean-spacing path is

$$
\boxed{
w_h(t)
=
\int_0^t \sigma(s)\,d\bar\lambda(s)
=
\int_0^t\sigma(s)\dot{\bar\lambda}(s)\,ds.
}
$$

For load cycle $k$,

$$
\boxed{
H_k
=
\oint_k\sigma\,d\bar\lambda.
}
$$

$H_k$ has units J/m³. It is the area of the stress–mean-spacing-strain loop. It is not automatically equal to permanent damage energy; under the present reflecting intact-domain model it is rate-dependent dissipated path work.

## 3. Kinetic equation that closes the four quantities

Assume that omitted lattice degrees of freedom act as an isothermal Markov bath and that velocity relaxation is faster than the resolved spacing evolution. With

$$
\chi=\frac{E_0}{k_BT},
\qquad
\tau=\frac{t}{t_r},
$$

the candidate overdamped evolution equation is

$$
\boxed{
\frac{\partial p}{\partial\tau}
=
\frac{\partial}{\partial\lambda}
\left[
(\phi'(\lambda)-f(\tau))p
+\chi^{-1}\frac{\partial p}{\partial\lambda}
\right].
}
$$

Equivalently,

$$
\partial_\tau p=-\partial_\lambda J,
\qquad
J=(f-\phi')p-\chi^{-1}\partial_\lambda p.
$$

The current implementation imposes

$$
J(\lambda_{\min},t)=J(\lambda_{\max},t)=0.
$$

Consequences:

- total intact probability remains exactly normalized;
- finite $t_r$ makes loading and unloading distributions differ;
- $H_k$ can be nonzero;
- probability is not yet lost irreversibly into a broken state.

The next crack-initiation extension must replace the upper no-flux boundary by a physically derived first-passage/commitment rule and add a broken-state probability $q(t)$.

## 4. FEM interface

For element $e$, the C solver supplies

$$
\sigma_e(t),\qquad
\epsilon_e(t),\qquad
u_e(t).
$$

The probability post-processor uses only

$$
f_e(t)=\frac{\sigma_e(t)}{E}
$$

to solve

$$
p_e(\lambda,t).
$$

It then exports

$$
\bar a_e(t),\quad
\operatorname{Var}_e[\lambda],\quad
u_{{\rm LJ},e}(t),\quad
w_{h,e}(t),\quad
\int_{\lambda_c}^{\lambda_{\max}}p_e\,d\lambda.
$$

The last quantity is an instability-tail diagnostic, not yet a crack probability.

The FEM element length remains a numerical discretization length. It is not identified with the atomic spacing, $\ell_{\rm stat}$, or an independent statistical-cell size. A future statistical cell may contain multiple FEM integration points, or one FEM element may contain many correlated statistical cells.

## 5. Numerical method and verification

The probability equation uses a cell-centered Chang–Cooper finite-volume flux and backward Euler time stepping. The face drift is evaluated as a discrete potential difference so that the declared finite-volume Gibbs density is an exact zero-current state under constant load.

Regression tests verify:

- normalization to numerical precision;
- preservation of constant-load local equilibrium;
- positive closed-cycle loop area at finite relaxation time;
- conversion from mean stretch to physical mean spacing;
- rejection of invalid time histories.

The reproducible demo uses a uniform-area bar, so all element stresses and probability states are spatially uniform. This is intentional: it validates the coupling without inventing a stress concentration. A later variable-area/notched 1D geometry will produce a nonuniform probability map using the same interface.

---

# 마일스톤 16 — 확률분포·에너지 히스테리시스·1D FEM 결합

## 핵심 결과

현재 네 식은 다음 하나의 계산 흐름으로 연결된다.

$$
\sigma_e(t)
\longrightarrow
p_e(\lambda,t)
\longrightarrow
\left\{
\bar a_e,
\operatorname{Var}_e[\lambda],
u_{{\rm LJ},e},
w_{h,e},
p_{{\rm tail},e}
\right\}.
$$

정규화식은

$$
\int p_e(\lambda,t)d\lambda=1,
$$

평균거리식은

$$
\bar a_e(t)=a_0\int\lambda p_e(\lambda,t)d\lambda,
$$

원자층 평균에너지식은

$$
u_{{\rm LJ},e}(t)
=E\int[\phi(\lambda)-\phi(1)]p_e(\lambda,t)d\lambda,
$$

에너지 히스테리시스식은

$$
H_{e,k}=\oint_k\sigma_e\,d\bar\lambda_e
$$

이다.

이를 닫는 후보 동역학은

$$
\frac{\partial p_e}{\partial\tau}
=
\partial_\lambda
\left[(\phi'-\sigma_e/E)p_e+\chi^{-1}\partial_\lambda p_e\right]
$$

이다. 유한한 완화시간 $t_r$ 때문에 같은 응력에서도 loading과 unloading의 $p_e$가 달라지고 히스테리시스 면적이 생긴다.

다만 현재 상단경계는 no-flux이므로 비가역 균열확률은 아직 아니다. 지금 계산되는 $\lambda_c$ 이상 tail은 불안정 접근도 진단값이다. 다음 단계에서 LJ barrier에 근거한 first-passage flux와 broken-state probability $q_e(t)$를 붙여야 실제 initiation probability로 발전한다.

메시 시각화에서는 실제 C FEM의 node와 element 경계를 표시한다. 현재 예제는 균일 단면이므로 모든 요소의 응력이 동일하다. 따라서 균일한 색이 정상이며, 취약부 지도를 만들려면 다음 단계에서 1D variable-area/notched bar를 추가해야 한다.
