# Auxiliary Symbol Index / 보조 수학기호 Index — 1D active derivations

This file covers auxiliary symbols that appear in detailed derivations but are not primary state variables.  
이 파일은 상세 유도에서 실제로 사용되지만 주 상태변수는 아닌 보조 기호를 정의한다.

The same mandatory rule applies: equation definition + English/Korean term + mathematical meaning + physical meaning + unit/scaling + status + dependencies.

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\nu$ | **not an active mean-rate symbol** | reserved Greek nu | 예약된 그리스문자 nu | not used as a state variable in the active theory | no active physical meaning | — | RESERVED | any legacy $\nu=\mathbb E[c\mid\lambda]$ is a typography error; use $u$ |
| $q^*$ | $q(\tau_L)=q(\tau_U)=q^*$ | matched force level | 동일비교 하중수준 | scalar forcing value used to pair two times | same instantaneous applied load on loading/unloading branches | 1 | DEFINITION | $q,\tau_L,\tau_U$ |
| $\tau_L$ | $q(\tau_L)=q^*,\;\dot q(\tau_L)>0$ | loading-branch comparison time | 하중증가 구간 비교시각 | time selected on increasing branch | instant at which load reaches $q^*$ while increasing | 1 | DEFINITION | $q,q^*$ |
| $\tau_U$ | $q(\tau_U)=q^*,\;\dot q(\tau_U)<0$ | unloading-branch comparison time | 하중감소 구간 비교시각 | time selected on decreasing branch | instant at which load reaches $q^*$ while decreasing | 1 | DEFINITION | $q,q^*$ |
| $\mathcal R_2$ | $\mathcal R_2(\tau)=\{P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)\}$ | reduced history descriptor | 축약 이력상태 기술자 | ordered set of retained one-point fields | reduced microscopic state used to compare load paths | mixed functional object | DEFINITION | $P,u,\Theta$ |
| $\mathcal S$ | hypothetical $\mathcal R_2=\mathcal S[q]$ | memoryless state map | 무기억 상태사상 | hypothetical single-valued map from instantaneous load to reduced state | constitutive map ruled out when same-load states differ | map | CONDITIONAL/HYPOTHETICAL | $q,\mathcal R_2$ |
| $Q_c$ | $Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)d\lambda$ | instantaneous critical tail mass | 순간 임계 초과확률질량 | upper-tail integral of nonabsorbing $P$ | instantaneous fraction of local spacings above threshold, allowed to return | 1 | DEFINITION | $P,\lambda_c$ |
| $S_M$ | $S_M=\frac1M\sum_i\chi_i$ | finite empirical survivor fraction | 유한계 경험적 생존비율 | average of local survival indicators | fraction of represented spacings not yet first-passed | 1 | EXACT finite definition | $M,\chi_i$ |
| $F_{{\rm ci},M}^{\rm local}$ | $1-S_M$ | finite local initiation fraction | 유한계 국소 균열개시 비율 | complement of $S_M$ | fraction of represented spacings already first-passed | 1 | DEFINITION | $S_M$ |
| $h_\tau$ | $h_\tau=j_{\rm esc}/S=-d\ln S/d\tau$ | nondimensional hazard | 무차원 위험률 | event rate per nondimensional time conditional on survival | local initiation rate in $\tau$ units | 1/$\tau$ | DEFINITION | $j_{\rm esc},S$ |
| $h_t$ | $h_t=h_\tau/t_0$ | physical-time hazard | 물리시간 위험률 | dimensional hazard | local initiation rate per second | s$^{-1}$ | DEFINITION | $h_\tau,t_0$ |
| $\bar\lambda_{\rm surv}$ | $\bar\lambda_{\rm surv}=S^{-1}\int_0^{\lambda_c}\lambda P_b d\lambda$ | survivor-conditioned mean spacing | 생존조건부 평균간격 | conditional first moment under $\widehat P_b$ | mean spacing among intact states only | 1 | DEFINITION | $S,P_b,\lambda_c$ |
| $\bar U_{\rm surv}$ | $\bar U_{\rm surv}=\frac{U_{\rm ref}}S\int_0^{\lambda_c}\Delta\phi P_b d\lambda$ | survivor-conditioned configurational energy | 생존조건부 배치에너지 | conditional expectation of $\Delta\phi$ over survivors | mean recoverable energy among intact states | J | DEFINITION | $U_{\rm ref},S,\Delta\phi,P_b$ |
| $V_{\rm phys}$ | $V_{\rm phys}=U_{\rm ref}V^*$ | physical chain configurational energy | 물리 사슬 배치에너지 | dimensionalized total potential | total recoverable potential energy of represented chain | J | DEFINITION | $U_{\rm ref},V^*$ |
| $W_{\rm ext}^{\rm cyc}$ | $W_{\rm ext}^{\rm cyc}=\int_{\rm cycle}q\dot x_M d\tau$ in normalized units | cycle external work | 사이클 외부일 | line/time integral of external power over one cycle | work delivered by prescribed end load in one cycle | normalized energy | DEFINITION | $q,x_M$ |
| $\Delta E_{\rm mech}^{\rm cyc}$ | $E_{\rm mech}^*(\tau_{n+1})-E_{\rm mech}^*(\tau_n)$ | cycle mechanical-energy change | 사이클 기계에너지 변화 | difference of mechanical energy across cycle endpoints | recoverable energy stored/released over the cycle | normalized energy | DEFINITION | $E_{\rm mech}^*$ |
| $D_{\rm irr}^{\rm cyc}$ | $\int_{\rm cycle}\dot D_{\rm irr}^*d\tau$ | cycle irreversible dissipation | 사이클 비가역 소산 | cycle integral of nonnegative dissipation rate | irreversible energy lost in one cycle | normalized energy | OPEN physically | $\dot D_{\rm irr}^*$ |
| $D_{\rm irr}^*(\tau;\Gamma_0)$ | $\int_{\tau_0}^{\tau}\dot D_{\rm irr}^*(\Phi_{s,\tau_0}^q(\Gamma_0),s)ds$ | trajectory irreversible functional | 궤적 비가역 소산함수 | path functional on one full-state trajectory | accumulated irreversible energy along one microscopic realization | normalized energy | OPEN physically | $\dot D_{\rm irr}^*,\Phi^q,\Gamma_0$ |
| $\lambda_{\rm ph}(s)$ | $d\lambda_{\rm ph}/ds=c_{\rm ph}$ | projected phase-space characteristic spacing | 투영 위상공간 특성 간격 | coordinate of a characteristic of the projected kinetic PDE | spacing coordinate carried along projected phase-space flow | 1 | DEFINITION | $c_{\rm ph},A$ |
| $c_{\rm ph}(s)$ | $dc_{\rm ph}/ds=A(\lambda_{\rm ph},c_{\rm ph},s)$ | projected phase-space characteristic rate | 투영 위상공간 특성 간격속도 | rate coordinate of projected phase-space characteristic | spacing-rate coordinate carried along projected flow | 1/$\tau$ | DEFINITION | $A,\lambda_{\rm ph}$ |
| $F_0$ | $F_0(\lambda,c)=F(\lambda,c,\tau_0)$ | initial phase-space density | 초기 위상공간 밀도 | initial condition for projected kinetic transport | initial joint spacing/rate population | phase-space density | DEFINITION | $F,\tau_0$ |
| $P_0$ | $P_0(\lambda)=P(\lambda,\tau_0)$ | initial spacing density | 초기 간격밀도 | initial marginal density | initial spacing population | inverse $\lambda$ | DEFINITION | $P,\tau_0$ |
| $u_0$ | $u_0(\lambda)=u(\lambda,\tau_0)$ | initial conditional mean rate | 초기 조건부 평균 간격속도 | initial conditional first moment | initial mean opening/closing rate at each spacing | 1/$\tau$ | DEFINITION | $u,\tau_0$ |
| $\Theta_0$ | $\Theta_0(\lambda)=\Theta(\lambda,\tau_0)$ | initial conditional variance | 초기 조건부 분산 | initial conditional second central moment | initial unresolved rate spread | 1/$\tau^2$ | DEFINITION | $\Theta,\tau_0$ |
| $\Gamma_0^*$ | deterministic initial state used when $\mu_0=\delta_{\Gamma_0^*}$ | single deterministic initial state | 단일 결정론적 초기상태 | atom of the initial measure | one prescribed microscopic realization | mixed state | DEFINITION | $\Gamma_0,\mu_0$ |

## Auxiliary exact relations / 보조 정확식

Same-force history test:

$$
\boxed{
q(\tau_L)=q(\tau_U)=q^*,
\qquad
\dot q(\tau_L)>0,
\qquad
\dot q(\tau_U)<0
}
$$

A memoryless load-only description would require

$$
\boxed{
\mathcal R_2(\tau)=\mathcal S[q(\tau)]
}
$$

and is contradicted for a trajectory when

$$
\boxed{
\mathcal R_2(\tau_L)\neq\mathcal R_2(\tau_U)
}
$$

Cycle energy balance with a future irreversible mechanism:

$$
\boxed{
W_{\rm ext}^{\rm cyc}
=\Delta E_{\rm mech}^{\rm cyc}+D_{\rm irr}^{\rm cyc}
}
$$

Projected phase-space characteristic density:

$$
\boxed{
F(\lambda_{\rm ph}(\tau),c_{\rm ph}(\tau),\tau)
=F_0(\lambda_{\rm ph}(\tau_0),c_{\rm ph}(\tau_0))
\exp\left[-\int_{\tau_0}^{\tau}\partial_cA(\lambda_{\rm ph}(s),c_{\rm ph}(s),s)ds\right]
}
$$
