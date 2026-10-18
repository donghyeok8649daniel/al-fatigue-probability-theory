# 지배방정식과 적용 예시

## 1. 기준과 단위

이 문서는 [N=1 UI 기준](../app/solver_adapter.py), [원래 TwoRowLJ](../solver_v1/model.py), [보존형 PDE 코드](../solver_v1/probability_pde_2d.py)의 표기를 연결합니다. 일반 $2N$ 이론은 [canonical manuscript](../paper/fatigue_probability_main_v4.tex)의 correlated joint law를 사용하며 product closure를 가정하지 않습니다. 과거 N=3 stochastic screening은 현재 N=1 PDE/UI 기준과 별개입니다.

같은 공통 길이 척도를 사용할 때

$$q_{\rm phys}=L_0q^*,\quad G_{\rm phys}=E_0G^*,\quad
t_{\rm phys}=t_0t^*,\quad \theta^*=\frac{k_BT}{E_0},$$

$$P_{N,\rm phys}=\frac{P_N^*}{L_0^{2N}},\qquad
M_{i,\rm phys}=\frac{L_0^2}{t_0E_0}M_i^*.$$

$M_{\rm phys}$의 단위는 $\mathrm{m^2/(J\,s)}$입니다. L0/E0만으로 t0를 결정할 수 없습니다. 선언된 kinetic reference L0/E0는 TwoRowLJ 실제 보정값으로 자동 전용하지 않습니다. 현재 기본 kinetic 파일은 미보정입니다.

## 2. 에너지와 하중

N=1에서 원래 tilt는

$$G^*(a^*,s^*,t^*)=
W^*(a^*,s^*)-f^*(t^*)[(a^*-a_0^*)+\chi s^*],$$

$$\varepsilon_{\rm ax}=\frac{(a^*-a_0^*)+\chi s^*}{a_0^*},\quad
f^*=\kappa_{\rm ax}^*\frac{\sigma_{\rm Pa}}{E_{x,\rm Pa}},\quad
\kappa_{\rm ax}^*=\frac{a_0^*}{c^TH_0^{*-1}c},\quad c=(1,\chi)^T.$$

$W^*$는 무한 Poisson/Bessel two-row 에너지입니다. N>1에서는 upper-cell 상호작용을 보존합니다. 별도 normal-chain 에너지를 중복 가산하지 않습니다.

물리적 대표계 이론의 $\mathcal G_N=\mathcal U_N-V_c\sigma\varepsilon$를 해석하려면 같은 에너지 counting과 실제 기하/work mapping이 필요합니다. stress ratio mapping 자체는 실제 Vc를 측정하지 않습니다. $V_{\rm eff}$와 $V_{\rm geom}$을 동일 값으로 채우지 않습니다.

## 3. Smoluchowski 지배방정식

상수 positive diagonal mobility를 쓰는 원래 방정식은

$$\partial_{t^*}P_N^*=-\nabla_*\cdot J^*,\qquad
J^*=-M^*(P_N^*\nabla_*G_N^*+\theta^*\nabla_*P_N^*).$$

이는 correlated joint density의 직접 진화입니다. 원래 모형의 확률을 Monte Carlo 빈도나 독립 one-cell product로 교체하지 않습니다. 아래 FVM 설명은 실제 구현된 N=1의 $q^*=(a^*,s^*)$에 한정합니다.

## 4. FVM에서 density와 raw 질량

$$\nu_{B,k}=\int_{B_k}dq^*,\qquad
p_k=\int_{B_k}\rho^*dq^*,\qquad m=\sum_kp_k.$$

원래 uniform grid에서는 $\nu_B=\Delta a^*\Delta s^*$이고 코드의 cell_volume은 이 구성공간 적분측도입니다. 재료 부피가 아닙니다.

한 내부 face의 left-to-right density flux는

$$J_{LR}=\frac{D_i^*}{\Delta q_i^*}
[B_{\rm SG}(\psi_{LR})\rho_L-B_{\rm SG}(-\psi_{LR})\rho_R],$$

$$D_i^*=M_i^*\theta^*,\quad
\psi_{LR}=\frac{G_R^*-G_L^*}{\theta^*},\quad
B_{\rm SG}(z)=\frac{z}{e^z-1},\quad B_{\rm SG}(0)=1.$$

이미 적분된 raw mass의 face rate는

$$T_{LR}=k_{LR}p_L-k_{RL}p_R,\quad
k_{LR}=\frac{D_i^*}{(\Delta q_i^*)^2}B_{\rm SG}(\psi_{LR}),\quad
k_{RL}=\frac{D_i^*}{(\Delta q_i^*)^2}B_{\rm SG}(-\psi_{LR}).$$

a-face에서는 $T_{LR}=\Delta s^*J_{LR}$, s-face에서는 $T_{LR}=\Delta a^*J_{LR}$입니다. raw rate에 transverse factor나 nu_B를 다시 곱하지 않습니다. 전체 내부 face를 조립하면 $\dot p=L_pp$이고 내부 conservative block의 열 합은 0입니다.

density와 raw mass의 일반 연산자 연결은 같은 partition/order/measure에서 $p=D_B\rho$와 $L_p=D_BL_\rho D_B^{-1}$를 함께 결속해야 합니다. 위 uniform-grid 코드를 임의 nonuniform grid 검증으로 확장하지 않습니다.

## 5. 원래 opening 사건과 확률

원래 mechanical opening basin의 binary keep mask를 $K_k\in\{0,1\}$로 쓰면 한 stage에서

$$p_k^+=K_kp_k^-,\qquad
\Delta A_{\rm open}=\sum_k(1-K_k)p_k^-.$$

rho 또는 p를 재정규화하여 제거 질량을 survivor에 되돌리지 않습니다. source-supported normalized initial cohort와 원래 opening-only absorber에서 $P_{\rm init,local}=A_{\rm open}$, $S_{\rm local}=1-A_{\rm open}$입니다. 실제 density 적분 m과의 차이 $R_{\rm mass}=m+A_{\rm open}-1$는 수치 잔차로 별도 기록합니다. 무조건 $1-m$을 물리 균열확률로 채택하지 않습니다.

같은 survivor의 관측량은 $m>0$에서 $\langle X\rangle_{\rm surv}=\sum p_kX_k/m$입니다. m=0에서 조건부 평균은 정의되지 않습니다. specimen aggregation에는 실제 상관/집단 법칙이 추가로 필요합니다.

## 6. 적용 예시 하나: 원래 N=1 모형의 내부 a-face

이 예시는 전체 격자의 한 face 기여를 뽑아 계산하는 기호 예시입니다. 새로운 시편·하중·grid·P0·시간 값을 지정하지 않습니다.

### 확인된 native 기준

| 입력 | 값 | 출처 |
|---|---:|---|
| N | 1 | app/solver_adapter.py canonical_model_params |
| b*, epsilon_LJ*, sigma_LJ* | 1, 1, 0.82 | ModelParams 상속 |
| theta*, Ma*, Ms*, chi | 0.02, 1, 0.05, 0.20 | canonical_model_params |
| a0*, kappa_ax* | 0.7713438268704838, 86.29296488740997 | AGENTS.md의 기록된 N=1 기준 |

일반 ModelParams의 N=3, theta*=.009, Ms*=.15, chi=.40은 이 예시에 섞지 않습니다. 위 값은 native mechanism reference이며 실제 Al 물성/온도/mobility가 아닙니다.

### 주어진 source 입력을 식에 연결

동일 stage에서 주어질 $r_\sigma=\sigma_{\rm Pa}/E_{x,\rm Pa}$, $a_L^*,a_R^*,s_0^*$, $\Delta a^*>0$, raw masses $p_L,p_R$를 그대로 기호로 둡니다. 같은 s-face level에서

$$f^*=86.29296488740997\,r_\sigma,$$

$$\Delta G^*=W^*(a_R^*,s_0^*)-W^*(a_L^*,s_0^*)-f^*\Delta a^*,$$

$$\psi=\frac{\Delta G^*}{0.02},\quad D_a^*=1\times0.02=0.02,$$

$$T_{LR}=\frac{0.02}{(\Delta a^*)^2}
[B_{\rm SG}(\psi)p_L-B_{\rm SG}(-\psi)p_R].$$

이 face의 기여는 $\dot p_L=-T_{LR}$, $\dot p_R=T_{LR}$입니다. 두 raw 질량의 합에는 이 내부 face가 기여하지 않습니다. 적법한 h*가 주어질 때 explicit 기여는 $p_L-h^*T_{LR}$와 $p_R+h^*T_{LR}$이며, 실제 step은 전체 grid outgoing rates/CFL과 다른 face를 함께 사용해야 합니다. 이 두 항만으로 전체 solver 결과를 산출한 것은 아닙니다.

그 다음 실제 원래 source가 공급한 binary opening mask를 5절의 식에 적용해 제거 raw mass를 계산합니다. 임의 fraction이나 임계 spacing을 넣지 않습니다.

### 수치 결과를 얻기 위해 남은 원자료

실제 load/modulus, 전체 grid·raw initial law·same-stage energy, mechanical mask, 허용 step/source order가 필요합니다. 물리 s/Hz·m/J/Pa 해석에는 같은 모형의 energy/length/work/kinetic mapping이 추가로 필요합니다. 이 예시는 없는 값을 채우지 않으므로 실제 오차·피로수명·시편 확률 수치는 보고하지 않습니다.
