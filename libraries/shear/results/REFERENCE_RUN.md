# Reference Run — Rubin-chain mechanics-derived hysteresis

The following numbers were produced from the reference nondimensional case implemented in `theory/rubin_chain.py`.

## Parameters

- $M=m=1$
- $K_0=k=1$
- $F_a=0.1$
- $\omega=0.5$
- chain band edge $\omega_D=2$
- finite-chain simulation: 1200 masses, $\Delta t=0.02$, 60 periods, 5-period smooth ramp
- loop statistics: cycles 10 through 49, before a far-boundary reflection returns

## Analytic semi-infinite result

$$
Z(\omega)=0.875+0.4841229182759271i
$$

$$
Q_a=0.1
$$

$$
\phi=0.5053605102841573\ \mathrm{rad}=28.95502437185985^\circ
$$

$$
\boxed{A_H^{\mathrm{analytic}}=0.015209170034901047}
$$

## Full conservative finite-chain integration

$$
\boxed{\langle A_H^{\mathrm{numeric}}\rangle=0.015208839984912282}
$$

Cycle-to-cycle standard deviation:

$$
1.921149725428978\times10^{-7}
$$

Relative analytic-vs-numeric loop-area error:

$$
\boxed{2.1700723182610268\times10^{-5}}
$$

Final internal energy:

$$
E_{\mathrm{int}}=0.8644685639287875
$$

Integrated external work:

$$
W_{\mathrm{ext}}=0.8644577442919658
$$

Relative energy-balance error:

$$
\boxed{1.2516096816928344\times10^{-5}}
$$

## Interpretation

This run demonstrates that a nonzero loop in a **resolved coordinate** can arise from a fully conservative microscopic chain without adding a viscous damping coefficient. The loop area is energy transferred into propagating unresolved modes.

This is a Milestone-1 proof-of-principle, not yet a quantitative aluminum fatigue prediction and not yet Milestone 2 secular fatigue accumulation.

---

# 한국어 번역 — Rubin-chain 역학 유도형 히스테리시스 기준 계산

다음 수치들은 `theory/rubin_chain.py`에 구현된 기준 무차원 문제에서 얻은 값이다.

## 파라미터

- $M=m=1$
- $K_0=k=1$
- $F_a=0.1$
- $\omega=0.5$
- 사슬의 band edge: $\omega_D=2$
- 유한 사슬 simulation: 질량 1200개, $\Delta t=0.02$, 60주기, 5주기 smooth ramp
- loop 통계: far-boundary reflection이 되돌아오기 전인 10번째 cycle부터 49번째 cycle까지 사용

## 준무한 사슬 해석 결과

$$
Z(\omega)=0.875+0.4841229182759271i
$$

$$
Q_a=0.1
$$

$$
\phi=0.5053605102841573\ \mathrm{rad}=28.95502437185985^\circ
$$

그리고 히스테리시스 면적은

$$
\boxed{A_H^{\mathrm{analytic}}=0.015209170034901047}
$$

이다.

## 전체 보존 유한사슬 직접 적분

수치적으로 얻은 평균 loop area는

$$
\boxed{\langle A_H^{\mathrm{numeric}}\rangle=0.015208839984912282}
$$

이다.

cycle 간 표준편차는

$$
1.921149725428978\times10^{-7}
$$

이며, 해석값과 수치값의 loop-area 상대오차는

$$
\boxed{2.1700723182610268\times10^{-5}}
$$

이다.

최종 내부에너지는

$$
E_{\mathrm{int}}=0.8644685639287875
$$

이고, 적분된 외부 일은

$$
W_{\mathrm{ext}}=0.8644577442919658
$$

이다.

따라서 에너지 수지의 상대오차는

$$
\boxed{1.2516096816928344\times10^{-5}}
$$

이다.

## 해석

이 계산은 점성 damping coefficient를 추가하지 않은 완전 보존 미시사슬에서도 **관심 좌표만 관찰하면** 0이 아닌 히스테리시스 loop가 생길 수 있음을 보여준다. loop area는 propagating unresolved mode로 전달된 에너지에 해당한다.

이 결과는 Milestone 1의 원리 증명이다. 아직 Al의 정량적 피로 예측도 아니고, Milestone 2의 장기적인 fatigue accumulation을 달성한 것도 아니다.
