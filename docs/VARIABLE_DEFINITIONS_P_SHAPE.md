# Variable Definitions — Exact Instantaneous $P$ Shape

## Classification labels

- **EXACT / IDENTITY** — exact under the stated smooth moment representation and original nonlinear 1D layer-LJ mechanics.
- **DEFINITION** — a mathematical definition.
- **SPECIAL-CASE ASSUMPTION** — used only to recover a restricted diagnostic form; not part of the active global model.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $P(\lambda,t)$ | $\int F(\lambda,v,t)\,dv$ | smooth one-point spacing density | inverse stretch | DEFINITION |
| $u(\lambda,t)$ | $\mathbb E[v\mid\lambda]$ | conditional mean spacing velocity | reduced time$^{-1}$ | DEFINITION |
| $\Theta(\lambda,t)$ | $\operatorname{Var}(v\mid\lambda)$ | conditional spacing-velocity variance | reduced time$^{-2}$ | DEFINITION |
| $K(\lambda,t)$ | $P(u^2+\Theta)$ | second velocity-moment density | reduced time$^{-2}$/stretch | EXACT / IDENTITY |
| $\bar a(\lambda,t)$ | $\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]$ | conditional spacing acceleration | reduced time$^{-2}$ | DEFINITION |
| $D_tu$ | $\partial_tu+u\partial_\lambda u$ | material acceleration of conditional mean flow in spacing space | reduced time$^{-2}$ | DEFINITION |
| $m_+(\lambda,t)$ | $\mathbb E[\phi'(\lambda_{i+1})\mid\lambda_i=\lambda]$ | right-neighbor conditional LJ force-gradient mean | reduced force unit | DEFINITION |
| $m_-(\lambda,t)$ | $\mathbb E[\phi'(\lambda_{i-1})\mid\lambda_i=\lambda]$ | left-neighbor conditional LJ force-gradient mean | reduced force unit | DEFINITION |
| $C(t)$ | normalization factor in the integrated shape law | fixes $\int P\,d\lambda=1$ | appropriate density unit | DEFINITION |
| $f_N(t)$ | spacing-independent conditional neighbor-force mean used only in the restricted special case | diagnostic force mean | reduced force unit | SPECIAL-CASE ASSUMPTION |

## Central exact shape identities

Where $P>0$ and $\Theta>0$,

$$
\boxed{
\partial_\lambda\ln P
=
\frac{\bar a-D_tu}{\Theta}
-
\partial_\lambda\ln\Theta
}
$$

and therefore

$$
\boxed{
P(\lambda,t)
=
\frac{C(t)}{\Theta(\lambda,t)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\bar a(s,t)-D_tu(s,t)}{\Theta(s,t)}\,ds
\right].
}
$$

For the original nonlinear interior layer-LJ equation,

$$
\boxed{
\bar a=m_++m_- -2\phi'(\lambda).
}
$$

No global Taylor expansion or finite harmonic representation is used.

---

# 한국어 번역 — 정확한 순간 $P$ 함수형 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시한 smooth moment representation과 원래 nonlinear 1D layer-LJ mechanics 아래 정확.
- **DEFINITION** — 수학적 정의.
- **SPECIAL-CASE ASSUMPTION** — 제한적 diagnostic form을 복원할 때만 사용하며 active global model에는 포함하지 않음.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $P(\lambda,t)$ | $\int F(\lambda,v,t)\,dv$ | smooth one-point spacing density | inverse stretch | DEFINITION |
| $u(\lambda,t)$ | $\mathbb E[v\mid\lambda]$ | conditional mean spacing velocity | reduced time$^{-1}$ | DEFINITION |
| $\Theta(\lambda,t)$ | $\operatorname{Var}(v\mid\lambda)$ | conditional spacing-velocity variance | reduced time$^{-2}$ | DEFINITION |
| $K(\lambda,t)$ | $P(u^2+\Theta)$ | second velocity-moment density | reduced time$^{-2}$/stretch | EXACT / IDENTITY |
| $\bar a(\lambda,t)$ | $\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]$ | conditional spacing acceleration | reduced time$^{-2}$ | DEFINITION |
| $D_tu$ | $\partial_tu+u\partial_\lambda u$ | spacing space에서 conditional mean flow의 material acceleration | reduced time$^{-2}$ | DEFINITION |
| $m_+(\lambda,t)$ | $\mathbb E[\phi'(\lambda_{i+1})\mid\lambda_i=\lambda]$ | 오른쪽 이웃 conditional LJ force-gradient mean | reduced force unit | DEFINITION |
| $m_-(\lambda,t)$ | $\mathbb E[\phi'(\lambda_{i-1})\mid\lambda_i=\lambda]$ | 왼쪽 이웃 conditional LJ force-gradient mean | reduced force unit | DEFINITION |
| $C(t)$ | integrated shape law의 normalization factor | $\int P\,d\lambda=1$을 맞춤 | appropriate density unit | DEFINITION |
| $f_N(t)$ | 제한적 special case에서만 사용하는 spacing-independent conditional neighbor-force mean | diagnostic force mean | reduced force unit | SPECIAL-CASE ASSUMPTION |

## 핵심 exact shape identity

$P>0$, $\Theta>0$인 곳에서

$$
\boxed{
\partial_\lambda\ln P
=
\frac{\bar a-D_tu}{\Theta}
-
\partial_\lambda\ln\Theta
}
$$

이고 따라서

$$
\boxed{
P(\lambda,t)
=
\frac{C(t)}{\Theta(\lambda,t)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\bar a(s,t)-D_tu(s,t)}{\Theta(s,t)}\,ds
\right]
}
$$

이다.

원래 nonlinear interior layer-LJ equation에서는

$$
\boxed{
\bar a=m_++m_- -2\phi'(\lambda)
}
$$

이다.

global Taylor expansion이나 finite harmonic representation은 사용하지 않는다.
