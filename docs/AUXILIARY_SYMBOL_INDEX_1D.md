# Auxiliary Symbol Index / 보조기호 Index

이 문서는 활성 1D 이론의 보조기호를 정의한다. exact microscopic layer와 reduced laboratory layer의 기호를 구분한다.

## 1. Loading-history symbols / 하중 이력 기호

### $q^*$

- English: matched load value
- 한국어: 동일 비교 하중값
- Mathematical definition: $q(\tau_L)=q(\tau_U)=q^*$
- Physical meaning: loading과 unloading 상태를 같은 외부 하중에서 비교하기 위한 값
- Unit/scaling: dimensionless
- Status: DEFINITION
- Dependencies: $q,\tau_L,\tau_U$

### $\tau_L$

- English: loading-branch time
- 한국어: loading branch 시각
- Mathematical definition: $q(\tau_L)=q^*$ and $\dot q(\tau_L)>0$
- Unit/scaling: nondimensional time
- Status: DEFINITION

### $\tau_U$

- English: unloading-branch time
- 한국어: unloading branch 시각
- Mathematical definition: $q(\tau_U)=q^*$ and $\dot q(\tau_U)<0$
- Unit/scaling: nondimensional time
- Status: DEFINITION

### $\mathcal R_2(\tau)$

$$
\mathcal R_2(\tau)=\{P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)\}
$$

- English: second-order reduced descriptor
- 한국어: 2차 축약 기술자
- Physical meaning: same-load non-retracing을 비교하기 위한 exact-chain reduced state
- Status: DEFINITION

## 2. Density-shape symbols / 밀도 형상 기호

### $\lambda_*$

- English: reference spacing point
- 한국어: 기준 간격점
- Mathematical definition: arbitrary reference point where $P>0$ and $\Theta>0$
- Physical meaning: integration origin only
- Unit/scaling: dimensionless
- Status: DEFINITION

### $\eta$

- English: spacing integration variable
- 한국어: 간격 적분변수
- Mathematical meaning: dummy spacing variable
- Unit/scaling: dimensionless
- Status: DEFINITION

### $\mathcal N_P(\tau)$

$$
\int_0^\infty P(\lambda,\tau)d\lambda=1
$$

- English: density normalization factor
- 한국어: 확률밀도 정규화 계수
- Mathematical meaning: positive integration factor fixed by total probability mass
- Status: DEFINITION

## 3. Initial-value symbols / 초기값 기호

### $P_0(\lambda)$

$$
P_0(\lambda)=P(\lambda,\tau_0)
$$

- English: initial spacing density
- 한국어: 초기 간격밀도
- Physical meaning in the laboratory closure: structural/prestress spacing density at the declared reference phase; not the instantaneous thermal-displacement PDF
- Unit/scaling: dimensionless density in $\lambda$
- Status: DEFINITION / ACTIVE REDUCED INPUT

### $u_0(\lambda)$

$$
u_0(\lambda)=u(\lambda,\tau_0)
$$

- English: initial conditional mean spacing rate
- 한국어: 초기 조건부 평균 간격속도
- Status: EXACT-LAYER DEFINITION

### $\Theta_0(\lambda)$

$$
\Theta_0(\lambda)=\Theta(\lambda,\tau_0)
$$

- English: initial conditional spacing-rate variance
- 한국어: 초기 조건부 간격속도 분산
- Status: EXACT-LAYER DEFINITION

### $\tau_0$

- English: initial nondimensional time
- 한국어: 초기 무차원 시간
- Mathematical meaning: lower temporal bound for exact microscopic initial data
- Status: DEFINITION

## 4. Exact-layer characteristic symbols / 정확층 특성곡선 기호

### $s$

- English: dummy time variable
- 한국어: 시간 적분변수
- Status: DEFINITION

### $r$

- English: inner dummy time variable
- 한국어: 내부 시간 적분변수
- Status: DEFINITION

### $X(s;\alpha)$

$$
\frac{dX}{ds}=u(X(s),s),
\qquad
X(\tau_0)=\alpha
$$

- English: exact reduced mean-flow characteristic
- 한국어: 정확 축약 평균흐름 특성곡선
- Physical meaning: mean probability-transport path in spacing space
- Status: DEFINITION

### $\alpha$

- English: characteristic label
- 한국어: 특성곡선 초기 라벨
- Mathematical definition: $X(\tau_0;\alpha)=\alpha$
- Status: DEFINITION

### $\mathcal I_u$

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)dr
$$

- English: accumulated mean-rate gradient
- 한국어: 누적 평균속도 구배
- Status: DEFINITION

### $S_\Theta$

$$
S_\Theta=2\Psi-\frac{1}{P}\partial_\lambda(PC_3)
$$

- English: variance source
- 한국어: 분산 소스항
- Status: DEFINITION

## 5. First-passage symbols / 최초통과 기호

### $Q_c(\tau)$

$$
Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)d\lambda
$$

- English: instantaneous nonabsorbing tail mass
- 한국어: 순간 비흡수 꼬리질량
- Physical meaning: current tail only; not cumulative first passage
- Status: DEFINITION

### $\chi_i(\tau)$

$$
\chi_i(\tau)=I[\tau<\tau_i^c]
$$

- English: local survival indicator
- 한국어: 국소 생존 지시함수
- Status: DEFINITION

### $S_M(\tau)$

$$
S_M(\tau)=\frac{1}{M}\sum_i\chi_i(\tau)
$$

- English: empirical local survival fraction
- 한국어: 경험적 국소 생존비율
- Status: DEFINITION

### $F_{\mathrm{ci},M}^{\mathrm{local}}$

$$
F_{\mathrm{ci},M}^{\mathrm{local}}=1-S_M
$$

- English: cumulative local first-passage fraction
- 한국어: 누적 국소 최초통과 비율
- Status: DEFINITION

### $h_\tau$

$$
h_\tau=\frac{j_{\mathrm{esc}}}{S}
=-\frac{d}{d\tau}\ln S
$$

- English: nondimensional initiation hazard
- 한국어: 무차원 균열개시 위험률
- Status: EXACT when $S>0$ and differentiable

### $h_t$

$$
h_t=\frac{h_\tau}{t_0}
$$

- English: physical-time initiation hazard
- 한국어: 물리시간 균열개시 위험률
- Unit: 1/s
- Status: DEFINITION

### $\widehat P_b$

$$
\widehat P_b(\lambda,t)=\frac{P_b(\lambda,t)}{S(t)}
$$

- English: survivor-conditioned density
- 한국어: 생존조건부 밀도
- Status: DEFINITION for $S>0$

### $\tau_{\mathrm{spec}}^c$

$$
\tau_{\mathrm{spec}}^c=\min_i\tau_i^c
$$

- English: specimen first-initiation time
- 한국어: 시편 최초 균열개시 시간
- Status: DEFINITION for one realization

### $F_{\mathrm{ci}}^{\mathrm{spec}}$

$$
F_{\mathrm{ci}}^{\mathrm{spec}}(\tau)=1-S_{\mathrm{spec}}(\tau)
$$

- English: specimen crack-initiation cumulative probability
- 한국어: 시편 균열개시 누적확률
- Status: DEFINITION

## 6. Energy-history symbols / 에너지 이력 기호

### $D_{\mathrm{irr}}$

$$
D_{\mathrm{irr}}(t)=\int_0^t\dot D_{\mathrm{irr}}(t')dt'
$$

- English: accumulated irreversible dissipation
- 한국어: 누적 비가역 소산
- Unit: J
- Status: DEFINITION; zero in the conservative exact baseline

### $r_j^{\mathrm{irr}}$

$$
\dot D_{\mathrm{irr}}^*=-\sum_jr_j^{\mathrm{irr}}\dot x_j
$$

- English: future irreversible node force
- 한국어: 미래 비가역 노드 힘
- Unit/scaling: normalized force in the displayed equation
- Status: OPEN

### $W_{\mathrm{ext}}^{\mathrm{cyc}}$

$$
W_{\mathrm{ext}}^{\mathrm{cyc}}
=\Delta E_{\mathrm{mech}}^{\mathrm{cyc}}+D_{\mathrm{irr}}^{\mathrm{cyc}}
$$

- English: external work over one cycle
- 한국어: 한 사이클 외력 일
- Status: CONDITIONAL on a physically valid irreversible mechanism

## 7. Reduced laboratory thermal-first-passage symbols / 실험실 축약 열 최초통과 기호

### $q_{\mathrm{ref}}$

$$
q_{\mathrm{ref}}=\frac{\sigma_{\mathrm{ref}}}{E}
$$

- English: reference reduced normal traction
- 한국어: 기준 무차원 normal traction
- Mathematical meaning: external reduced load at the phase where $P_0$ is defined
- Physical meaning: reference load used to interpret each $\lambda_0$ as a local structural equilibrium
- Unit/scaling: dimensionless
- Status: ACTIVE REDUCED INPUT
- Dependencies: $\sigma_{\mathrm{ref}},E$

### $q_r(\lambda_0)$

$$
q_r(\lambda_0)=\phi'(\lambda_0)-q_{\mathrm{ref}}
$$

- English: residual conjugate bias
- 한국어: 잔류 켤레 바이어스
- Mathematical meaning: local bias that makes $\lambda_0$ an equilibrium at the reference load
- Physical meaning: minimal local-prestress embedding of structural $P_0$; not a reconstruction of finite-chain neighbour ordering
- Unit/scaling: dimensionless
- Status: ACTIVE REDUCED CLOSURE
- Dependencies: $\lambda_0,\phi,q_{\mathrm{ref}}$

### $\Lambda(\lambda_0,t)$

$$
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
\qquad
\phi''[\Lambda(\lambda_0,t)]>0
$$

- English: quasistatic stable-branch spacing map
- 한국어: 준정적 안정가지 간격 사상
- Mathematical meaning: stable solution carrying the reference spacing label $\lambda_0$ under the applied load history
- Physical meaning: reversible structural spacing response in the laboratory-frequency reduced model
- Unit/scaling: dimensionless
- Status: ACTIVE REDUCED STATE MAP
- Dependencies: $\lambda_0,q(t),q_{\mathrm{ref}},\phi$

### $q_c$

$$
q_c=\phi'(\lambda_c)
$$

- English: maximum stable reduced traction
- 한국어: 최대 안정 무차원 traction
- Physical meaning: traction at the declared tangent-stiffness-loss point
- Unit/scaling: dimensionless
- Status: DERIVED MODEL THRESHOLD
- Dependencies: $\lambda_c,\phi$

### $A_c$

- English: characteristic cohesive area
- 한국어: 특성 응집 면적
- Mathematical meaning: area multiplying the local normal-instability energy scale
- Physical meaning: coherently participating area for one local normal first-passage event
- Unit: m^2
- Status: OPEN CALIBRATION PARAMETER
- Dependencies: independent experimental/higher-fidelity calibration
- Warning: not $A_0$, not a FEM element area, and not automatically an independent statistical-cell area

### $\Delta\psi_c(\lambda)$

$$
\Delta\psi_c(\lambda)
=
[\phi(\lambda_c)-\phi'(\lambda)\lambda_c]
-[\phi(\lambda)-\phi'(\lambda)\lambda]
$$

- English: dimensionless climb to the operational instability boundary
- 한국어: 균열개시 경계까지의 무차원 에너지 상승량
- Mathematical meaning: effective-energy difference from the stable state to $\lambda_c$ at the same local conjugate traction
- Unit/scaling: dimensionless
- Status: DERIVED REDUCED QUANTITY
- Dependencies: $\phi,\lambda,\lambda_c$

### $\Delta G_c(\lambda)$

$$
\Delta G_c(\lambda)=EA_ca_0\Delta\psi_c(\lambda)
$$

- English: characteristic-domain normal-instability barrier
- 한국어: 특성영역 normal 불안정 활성화 장벽
- Physical meaning: energy cost for the characteristic cohesive domain to reach the operational absorbing boundary
- Unit: J
- Status: DERIVED WITH SYMBOLIC $A_c$
- Dependencies: $E,A_c,a_0,\Delta\psi_c$

### $\nu_s(\lambda)$

$$
\nu_s(\lambda)=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
$$

- English: local transition-state attempt frequency
- 한국어: 국소 전이상태 시도 주파수
- Physical meaning: small-oscillation frequency inherited from the normal potential and inertial scale
- Unit: 1/s
- Status: DERIVED REDUCED PREFACTOR
- Dependencies: $\phi'',t_0$
- Warning: $\nu_s$ is not the conditional mean spacing-rate field $u$

### $k_c(\lambda,T;A_c)$

$$
k_c(\lambda,T;A_c)
=\nu_s(\lambda)
\exp\left[-\frac{\Delta G_c(\lambda)}{k_BT}\right]
$$

- English: local thermal normal-instability first-passage rate
- 한국어: 국소 열적 normal 불안정 최초통과율
- Mathematical meaning: rare-event sink rate used in the reduced survivor transport
- Physical meaning: positive-flux transition-state approximation for a locally re-equilibrated intact well
- Unit: 1/s
- Status: ACTIVE REDUCED LAW under the declared rare-event and local-equilibration assumptions
- Dependencies: $\nu_s,\Delta G_c,k_B,T,A_c$

### $W(\lambda_0,t)$

$$
W(\lambda_0,t)
=\exp\left[-\int_{t_0}^{t}k_c(\Lambda(\lambda_0,s),T;A_c)ds\right]
$$

- English: characteristic survivor weight
- 한국어: 특성곡선 생존 가중치
- Mathematical meaning: surviving mass fraction carried by one initial spacing label
- Unit/scaling: dimensionless
- Status: EXACT SOLUTION FACTOR of the declared reduced PDE
- Dependencies: $k_c,\Lambda$

### $\mathcal H_c(\lambda_0)$

$$
\mathcal H_c(\lambda_0)
=\int_0^{T_f}k_c[\Lambda(\lambda_0,t),T;A_c]dt
$$

- English: one-cycle integrated local hazard
- 한국어: 한 사이클 누적 국소 hazard
- Physical meaning: cumulative rare-event exposure of one structural reference label during one loading period
- Unit/scaling: dimensionless
- Status: DERIVED REDUCED CYCLE QUANTITY
- Dependencies: $k_c,\Lambda,T_f$

### $T_f$

$$
T_f=\frac1f
$$

- English: physical fatigue-loading period
- 한국어: 물리적 피로하중 주기
- Unit: s
- Status: DEFINITION
- Dependencies: loading frequency $f$

### $f$

- English: fatigue-loading frequency
- 한국어: 피로하중 주파수
- Unit: Hz
- Status: LOADING INPUT
- Dependencies: prescribed waveform

### $\mathcal H_f$

$$
\mathcal H_f
=\frac1f\int_0^1 k_c[\Lambda_*(\theta),T;A_c]d\theta
$$

- English: one-cycle hazard at loading frequency $f$
- 한국어: 주파수 $f$에서의 한 사이클 누적 hazard
- Mathematical meaning: cycle hazard for a fixed phase-shaped waveform in the strict quasistatic fast-equilibration limit
- Physical meaning: exposes the elapsed-time-controlled scaling of the fast thermal renewal hypothesis
- Unit/scaling: dimensionless
- Status: DERIVED NO-GO QUANTITY
- Dependencies: $f,k_c,\Lambda_*,T,A_c$
- Consequence: $\mathcal H_f\propto1/f$

### $N_{50}$

$$
S_{N_{50}}=\frac12
$$

- English: local median cycle count
- 한국어: 국소 중앙 생존 cycle 수
- Physical meaning: cycle count at 50 percent local survivor probability in the declared reduced regime
- Unit/scaling: cycles
- Status: DERIVED SURVIVAL OBSERVABLE
- Dependencies: $S_N,\mathcal H_c,P_0$
- Strict fast-equilibration consequence: $N_{50}\propto f$

### $S_N$

$$
S_N=\int_{\mathrm{stable}}P_0(\lambda_0)e^{-N\mathcal H_c(\lambda_0)}d\lambda_0
$$

- English: local survival after $N$ repeated cycles
- 한국어: $N$회 반복하중 후 국소 생존확률
- Unit/scaling: dimensionless
- Status: DERIVED REDUCED SURVIVAL
- Dependencies: $P_0,N,\mathcal H_c$
