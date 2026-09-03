# Crack Initiation Definition / 균열개시 정의

이 문서는 활성 1D normal-only 이론에서 crack initiation을 first-passage 문제로 정의한다.

## 1. Local mechanical threshold / 국소 기계적 임계점

active generalized-LJ interaction은

$$
\phi(\lambda)=\frac{\lambda^{-m}}{m(m-n)}-\frac{\lambda^{-n}}{n(m-n)}
$$

이다.

local tangent stiffness loss를 operational initiation threshold로 사용한다.

$$
\phi''(\lambda_c)=0
$$

따라서

$$
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
$$

이다.

이 기준은 adopted interaction model 아래의 mechanical criterion이다. 실제 재료에서의 experimental crack-initiation definition과의 calibration은 별도 문제다.

## 2. Local first-passage time / 국소 최초통과 시간

각 represented spacing $i$에 대해

$$
\tau_i^c=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}
$$

로 정의한다.

이 정의는 instantaneous tail mass와 다르다.

$$
Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)\,d\lambda
$$

$Q_c$는 현재 시각의 nonabsorbing tail이며 cumulative first passage가 아니다.

## 3. Finite empirical survival / 유한 경험적 생존

local survival indicator는

$$
\chi_i(\tau)=I[\tau<\tau_i^c]
$$

이다.

empirical survival fraction은

$$
S_M(\tau)=\frac{1}{M}\sum_{i=1}^{M}\chi_i(\tau)
$$

이고 cumulative local first-passage fraction은

$$
F_{\mathrm{ci},M}^{\mathrm{local}}(\tau)=1-S_M(\tau)
$$

이다.

분포론적으로 first-passage event density는

$$
-\frac{dS_M}{d\tau}=\frac{1}{M}\sum_i\delta(\tau-\tau_i^c)
$$

로 쓸 수 있다.

## 4. Survivor phase-space subdensity / 생존 위상공간 부분밀도

survivor phase-space subdensity를

$$
F_b(\lambda,c,\tau)
$$

로 두고 intact domain을

$$
0<\lambda<\lambda_c
$$

로 둔다.

right boundary $\lambda=\lambda_c$에서 failed side로부터의 재유입을 막기 위해 incoming velocity $c<0$에 대해

$$
F_b(\lambda_c,c,\tau)=0
$$

을 둔다.

이 조건은 kinetic absorbing boundary다.

## 5. Escape flux / 탈출 플럭스

outgoing first-passage flux는

$$
j_{\mathrm{esc}}(\tau)=\int_0^\infty cF_b(\lambda_c^-,c,\tau)\,dc
$$

이다.

$$
j_{\mathrm{esc}}\ge0
$$

이다.

lower boundary loss가 없다고 하면 survivor mass는

$$
S(\tau)=\int_0^{\lambda_c}\int_{-\infty}^{\infty}F_b(\lambda,c,\tau)\,dc\,d\lambda
$$

이다.

따라서

$$
\frac{dS}{d\tau}=-j_{\mathrm{esc}}
$$

이다.

cumulative local initiation fraction은

$$
F_{\mathrm{ci}}^{\mathrm{local}}=1-S
$$

이다.

## 6. Hazard / 위험률

$S>0$일 때 nondimensional hazard는

$$
h_\tau=\frac{j_{\mathrm{esc}}}{S}
$$

이고

$$
h_\tau=-\frac{d}{d\tau}\ln S
$$

이다.

physical-time hazard는

$$
h_t=\frac{h_\tau}{t_0}
$$

이다.

## 7. Survivor marginal and conditional density / 생존 주변밀도와 조건부 밀도

survivor spacing subdensity는

$$
P_b(\lambda,\tau)=\int_{-\infty}^{\infty}F_b(\lambda,c,\tau)\,dc
$$

이다.

이 함수는 normalized density가 아니라 subdensity다.

$$
\int_0^{\lambda_c}P_b(\lambda,\tau)\,d\lambda=S(\tau)
$$

이다.

survivor-conditioned normalized density는

$$
\widehat P_b(\lambda,\tau)=\frac{P_b(\lambda,\tau)}{S(\tau)}
$$

이다.

survivor-conditioned mean spacing은

$$
\bar\lambda_{\mathrm{surv}}(\tau)=\frac{1}{S(\tau)}\int_0^{\lambda_c}\lambda P_b(\lambda,\tau)\,d\lambda
$$

이다.

survivor-conditioned configurational energy는

$$
\bar U_{\mathrm{surv}}(\tau)=\frac{U_{\mathrm{ref}}}{S(\tau)}\int_0^{\lambda_c}[\phi(\lambda)-\phi(1)]P_b(\lambda,\tau)\,d\lambda
$$

이다.

## 8. Exact full-flow local survival / 전체 흐름 국소 생존

full initial-state measure $\mu_0$를 사용하면 local survival을 trajectory path functional로 직접 쓸 수 있다.

$$
S_{\mathrm{local}}(\tau)=\frac{1}{M}\sum_i\int I\left[\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

따라서

$$
F_{\mathrm{ci}}^{\mathrm{local}}(\tau)=1-S_{\mathrm{local}}(\tau)
$$

이다.

## 9. Specimen first-initiation time / 시편 최초 균열개시 시간

한 realization에서 specimen first-initiation time은

$$
\tau_{\mathrm{spec}}^c=\min_i\tau_i^c
$$

이다.

이 값은 represented spacings 중 하나라도 처음 threshold를 통과하는 시각이다.

## 10. Specimen survival probability / 시편 생존확률

ensemble measure $\mu_0$가 선언되면

$$
S_{\mathrm{spec}}(\tau)=\int I\left[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

이다.

specimen crack-initiation cumulative probability는

$$
F_{\mathrm{ci}}^{\mathrm{spec}}(\tau)=1-S_{\mathrm{spec}}(\tau)
$$

이다.

## 11. Local fraction is not automatically specimen probability / 국소 비율과 시편 확률의 구분

$$
1-S_M
$$

은 한 deterministic realization 내부의 local spatial first-passage fraction일 수 있다. 이것을 독립 cell product로 specimen-to-specimen probability로 바꾸지 않는다.

specimen probability는 $\mu_0$와 spatial correlation 구조를 포함하는 ensemble definition이 필요하다.

## 12. Post-initiation scope / 개시 이후 범위

$$
\tau_{\mathrm{spec}}^c
$$

에 도달한 뒤 현재 intact pre-crack chain은 실제 crack propagation model이 아니다.

따라서 specimen first-initiation 이후의 추가 local crossings는 post-crack model이 추가되지 않는 한 mathematical diagnostic으로만 해석한다.

## 13. Status summary / 상태 요약

- $\phi''(\lambda_c)=0$: MODEL-based local criterion
- $\tau_i^c$: DEFINITION
- kinetic absorbing boundary: EXACT formulation for the declared survivor transport
- $S_{\mathrm{local}}$: EXACT once the full flow and $\mu_0$ are declared
- $S_{\mathrm{spec}}$: EXACT path-integral formula once $\mu_0$ is declared
- physical construction of $\mu_0$ and specimen correlation scale: OPEN
- post-initiation crack propagation: outside current model
