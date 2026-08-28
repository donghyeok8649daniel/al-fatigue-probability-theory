# Variable Definitions — Exact 1D Distribution Transport

## Classification labels

- **EXACT / IDENTITY** — exact under the stated finite-$M$ definitions and original nonlinear 1D layer-LJ mechanics.
- **DEFINITION** — a chosen mathematical definition.
- **CONTINUUM REPRESENTATION** — a smooth representation of an exact empirical-measure identity.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $M$ | number of represented layer spacings | finite represented system size | count | DEFINITION |
| $\lambda_i$ | $a_i/a_0$ | normalized spacing of represented layer pair $i$ | dimensionless | DEFINITION |
| $v_i$ | $\dot\lambda_i$ | normalized spacing velocity | inverse reduced time | DEFINITION |
| $a_i^{(\lambda)}$ | $\ddot\lambda_i$ | normalized spacing acceleration | inverse reduced time$^2$ | DEFINITION |
| $P_M(\lambda,t)$ | $M^{-1}\sum_i\delta(\lambda-\lambda_i)$ | finite-$M$ empirical spacing density | inverse stretch | DEFINITION |
| $J_M(\lambda,t)$ | $M^{-1}\sum_i v_i\delta(\lambda-\lambda_i)$ | spacing-space probability flux | inverse reduced time | DEFINITION |
| $F_M(\lambda,v,t)$ | $M^{-1}\sum_i\delta(\lambda-\lambda_i)\delta(v-v_i)$ | finite-$M$ spacing-velocity phase-space empirical measure | inverse stretch-velocity | DEFINITION |
| $G_M(\lambda,v,t)$ | $M^{-1}\sum_i a_i^{(\lambda)}\delta(\lambda-\lambda_i)\delta(v-v_i)$ | acceleration flux in phase space | corresponding flux unit | DEFINITION |
| $\bar v(\lambda,t)$ | $J/P$ | conditional mean spacing velocity | inverse reduced time | CONTINUUM REPRESENTATION |
| $\bar a(\lambda,v,t)$ | $G/F$ | conditional mean spacing acceleration | inverse reduced time$^2$ | CONTINUUM REPRESENTATION |
| $K(\lambda,t)$ | $\int v^2F\,dv$ | second spacing-velocity moment density | velocity$^2$/stretch | DEFINITION |
| $A_1(\lambda,t)$ | $\int G\,dv$ | one-point conditional acceleration source | inverse reduced time$^2$/stretch | DEFINITION |
| $P_2^+$ | central/right neighboring-spacing joint density | right-neighbor ordering information | inverse stretch$^2$ | DEFINITION |
| $P_2^-$ | central/left neighboring-spacing joint density | left-neighbor ordering information | inverse stretch$^2$ | DEFINITION |
| $M_r$ | $M^{-1}\sum_i\lambda_i^r$ | finite empirical raw moment of order $r$ | dimensionless | DEFINITION |

## Exact transport identities

$$
\boxed{\partial_tP_M+\partial_\lambda J_M=0}
$$

and

$$
\boxed{
\partial_tF_M
+\partial_\lambda(vF_M)
+\partial_vG_M=0.
}
$$

For a smooth continuum representation,

$$
\boxed{\partial_tP+\partial_\lambda(P\bar v)=0.}
$$

The exact finite-$M$ moment identities are

$$
\boxed{
\dot M_r=r\langle\lambda^{r-1}v\rangle
}
$$

and

$$
\boxed{
\ddot M_r
=r(r-1)\langle\lambda^{r-2}v^2\rangle
+r\langle\lambda^{r-1}a^{(\lambda)}\rangle.
}
$$

---

# 한국어 번역 — 정확한 1D Distribution Transport 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 finite-$M$ 정의와 원래 nonlinear 1D layer-LJ mechanics 아래 정확.
- **DEFINITION** — 선택한 수학적 정의.
- **CONTINUUM REPRESENTATION** — exact empirical-measure identity의 매끄러운 연속 표현.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $M$ | represented layer spacing 수 | finite represented system size | count | DEFINITION |
| $\lambda_i$ | $a_i/a_0$ | represented layer pair $i$의 normalized spacing | 무차원 | DEFINITION |
| $v_i$ | $\dot\lambda_i$ | normalized spacing velocity | reduced time의 역수 | DEFINITION |
| $a_i^{(\lambda)}$ | $\ddot\lambda_i$ | normalized spacing acceleration | reduced time$^{-2}$ | DEFINITION |
| $P_M(\lambda,t)$ | $M^{-1}\sum_i\delta(\lambda-\lambda_i)$ | finite-$M$ empirical spacing density | inverse stretch | DEFINITION |
| $J_M(\lambda,t)$ | $M^{-1}\sum_i v_i\delta(\lambda-\lambda_i)$ | spacing-space probability flux | inverse reduced time | DEFINITION |
| $F_M(\lambda,v,t)$ | $M^{-1}\sum_i\delta(\lambda-\lambda_i)\delta(v-v_i)$ | finite-$M$ spacing-velocity phase-space empirical measure | inverse stretch-velocity | DEFINITION |
| $G_M(\lambda,v,t)$ | $M^{-1}\sum_i a_i^{(\lambda)}\delta(\lambda-\lambda_i)\delta(v-v_i)$ | phase-space acceleration flux | 대응 flux 단위 | DEFINITION |
| $\bar v(\lambda,t)$ | $J/P$ | conditional mean spacing velocity | inverse reduced time | CONTINUUM REPRESENTATION |
| $\bar a(\lambda,v,t)$ | $G/F$ | conditional mean spacing acceleration | reduced time$^{-2}$ | CONTINUUM REPRESENTATION |
| $K(\lambda,t)$ | $\int v^2F\,dv$ | second spacing-velocity moment density | velocity$^2$/stretch | DEFINITION |
| $A_1(\lambda,t)$ | $\int G\,dv$ | one-point conditional acceleration source | reduced time$^{-2}$/stretch | DEFINITION |
| $P_2^+$ | central/right neighboring-spacing joint density | 오른쪽 이웃 ordering 정보 | inverse stretch$^2$ | DEFINITION |
| $P_2^-$ | central/left neighboring-spacing joint density | 왼쪽 이웃 ordering 정보 | inverse stretch$^2$ | DEFINITION |
| $M_r$ | $M^{-1}\sum_i\lambda_i^r$ | $r$차 finite empirical raw moment | 무차원 | DEFINITION |

## 정확한 transport identity

$$
\boxed{\partial_tP_M+\partial_\lambda J_M=0}
$$

및

$$
\boxed{
\partial_tF_M
+\partial_\lambda(vF_M)
+\partial_vG_M=0
}
$$

가 성립한다.

매끄러운 continuum representation에서는

$$
\boxed{\partial_tP+\partial_\lambda(P\bar v)=0}
$$

이다.

정확한 finite-$M$ moment identity는

$$
\boxed{
\dot M_r=r\langle\lambda^{r-1}v\rangle
}
$$

및

$$
\boxed{
\ddot M_r
=r(r-1)\langle\lambda^{r-2}v^2\rangle
+r\langle\lambda^{r-1}a^{(\lambda)}\rangle
}
$$

이다.
