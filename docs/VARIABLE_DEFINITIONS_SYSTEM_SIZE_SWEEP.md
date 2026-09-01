# Variable Definitions — 1D Layer-LJ System-Size Closure Test

## Classification labels

- **EXACT / IDENTITY** — exact under the stated 1D reduced model.
- **DEFINITION** — chosen mathematical definition.
- **CONTROLLED NUMERICAL PROTOCOL** — numerical scaling introduced to compare systems consistently.
- **NUMERICAL DIAGNOSTIC** — computed comparison quantity.
- **EXPLORATORY NUMERICAL EXTRAPOLATION** — finite-data extrapolation with no theorem status.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $M$ | number of represented normal layer spacings | finite represented system size | dimensionless count | DEFINITION |
| $N_a$ | $M+1$ | number of represented layer nodes in the finite chain | dimensionless count | DEFINITION |
| $\omega(M)$ | $0.62/M$ in the current sweep | dimensionless loading angular frequency used to preserve dynamic similarity | dimensionless | CONTROLLED NUMERICAL PROTOCOL |
| $\chi$ | $\omega M$ | dynamic-similarity parameter | dimensionless | DEFINITION |
| $t_s$ | $2T=4\pi/\omega$ | phase-locked snapshot time | dimensionless time | DEFINITION |
| $t_s/M$ | $4\pi/(\omega M)$ | snapshot time normalized by chain-length scale | dimensionless | EXACT / IDENTITY |
| $D_{\rm KS}(M)$ | Kolmogorov distance between empirical and closure CDFs | one-point shape mismatch | dimensionless | NUMERICAL DIAGNOSTIC |
| $\gamma_{1,\rm sim}(M)$ | empirical third central moment divided by empirical variance$^{3/2}$ | deterministic spacing skewness | dimensionless | NUMERICAL DIAGNOSTIC |
| $\gamma_{1,\rm closure}(M)$ | closure third central moment divided by closure variance$^{3/2}$ | closure skewness at the same mean and energy | dimensionless | NUMERICAL DIAGNOSTIC |
| $C_1(t)$ | nearest-neighbor spacing covariance | first spatial-correlation diagnostic | dimensionless$^2$ in normalized spacing | DEFINITION |
| $C_k(t)$ | lag-$k$ spacing covariance | spatial correlation at separation $k$ | dimensionless$^2$ | DEFINITION |

The current dynamic-similarity condition is

$$
\boxed{
\chi=\omega M=0.62.
}
$$

Therefore

$$
\boxed{
\frac{t_s}{M}
=
\frac{4\pi}{\chi}
\approx20.26834.
}
$$

## Narrow-closure numerical quantities

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $\lambda_*$ | solution of $\alpha+\beta\psi'(\lambda_*)=0$ | mode used to center numerical quadrature in a sharply concentrated closure | dimensionless | DEFINITION from the closure exponent |
| $s_*$ | $[\beta\psi''(\lambda_*)]^{-1/2}$ | local width scale for numerical integration | dimensionless | DEFINITION |
| $\alpha$ | closure length multiplier | determined by mean-stretch constraint | dimensionless | closure variable |
| $\beta$ | closure energy multiplier | determined by mean-energy constraint | inverse normalized energy | closure variable |

The switch to mode-centered quadrature is a numerical resolution strategy only. It is not a material law and does not change the closure family.

---

# 한국어 번역 — 1D Layer-LJ System-Size Closure 시험 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 1D reduced model 아래 정확.
- **DEFINITION** — 선택한 수학적 정의.
- **CONTROLLED NUMERICAL PROTOCOL** — system을 일관되게 비교하기 위해 도입한 numerical scaling.
- **NUMERICAL DIAGNOSTIC** — 계산된 비교량.
- **EXPLORATORY NUMERICAL EXTRAPOLATION** — theorem status가 없는 finite-data extrapolation.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $M$ | represented normal layer spacing 수 | finite represented system size | 무차원 count | DEFINITION |
| $N_a$ | $M+1$ | finite chain의 represented layer node 수 | 무차원 count | DEFINITION |
| $\omega(M)$ | 현재 sweep에서 $0.62/M$ | dynamic similarity를 유지하기 위한 dimensionless loading angular frequency | 무차원 | CONTROLLED NUMERICAL PROTOCOL |
| $\chi$ | $\omega M$ | dynamic-similarity parameter | 무차원 | DEFINITION |
| $t_s$ | $2T=4\pi/\omega$ | phase-locked snapshot time | dimensionless time | DEFINITION |
| $t_s/M$ | $4\pi/(\omega M)$ | chain-length scale로 normalize한 snapshot time | 무차원 | EXACT / IDENTITY |
| $D_{\rm KS}(M)$ | empirical CDF와 closure CDF 사이 Kolmogorov distance | one-point shape mismatch | 무차원 | NUMERICAL DIAGNOSTIC |
| $\gamma_{1,\rm sim}(M)$ | empirical third central moment / empirical variance$^{3/2}$ | deterministic spacing skewness | 무차원 | NUMERICAL DIAGNOSTIC |
| $\gamma_{1,\rm closure}(M)$ | closure third central moment / closure variance$^{3/2}$ | 같은 mean과 energy를 가진 closure skewness | 무차원 | NUMERICAL DIAGNOSTIC |
| $C_1(t)$ | nearest-neighbor spacing covariance | 첫 spatial-correlation diagnostic | normalized spacing에서는 무차원$^2$ | DEFINITION |
| $C_k(t)$ | lag-$k$ spacing covariance | separation $k$의 spatial correlation | 무차원$^2$ | DEFINITION |

현재 dynamic-similarity 조건은

$$
\boxed{
\chi=\omega M=0.62
}
$$

이다.

따라서

$$
\boxed{
\frac{t_s}{M}
=
\frac{4\pi}{\chi}
\approx20.26834
}
$$

이다.

## Narrow-closure numerical quantity

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $\lambda_*$ | $\alpha+\beta\psi'(\lambda_*)=0$의 해 | sharply concentrated closure의 numerical quadrature 중심 mode | 무차원 | closure exponent에서의 DEFINITION |
| $s_*$ | $[\beta\psi''(\lambda_*)]^{-1/2}$ | numerical integration의 local width scale | 무차원 | DEFINITION |
| $\alpha$ | closure length multiplier | mean-stretch constraint로 결정 | 무차원 | closure variable |
| $\beta$ | closure energy multiplier | mean-energy constraint로 결정 | inverse normalized energy | closure variable |

mode-centered quadrature 전환은 numerical resolution strategy일 뿐이다. material law가 아니며 closure family 자체도 바꾸지 않는다.
