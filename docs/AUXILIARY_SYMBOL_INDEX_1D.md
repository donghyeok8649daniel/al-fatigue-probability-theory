# Auxiliary Symbol Index / 보조기호 Index

이 문서는 활성 1D 이론의 보조기호를 정의한다. 수식은 렌더러 호환성을 위해 단순한 문법만 사용한다.

## 1. Loading-history symbols / 하중 이력 기호

### $q^*$

- English: matched load value
- 한국어: 동일 비교 하중값
- Mathematical definition: scalar satisfying $q(\tau_L)=q(\tau_U)=q^*$
- Physical definition: loading과 unloading 상태를 같은 외부 하중에서 비교하기 위한 값
- Status: DEFINITION

### $\tau_L$

- English: loading-branch time
- 한국어: loading branch 시각
- Mathematical definition: time satisfying $q(\tau_L)=q^*$ and $\dot q(\tau_L)>0$
- Status: DEFINITION

### $\tau_U$

- English: unloading-branch time
- 한국어: unloading branch 시각
- Mathematical definition: time satisfying $q(\tau_U)=q^*$ and $\dot q(\tau_U)<0$
- Status: DEFINITION

### $\mathcal R_2(\tau)$

$$
\mathcal R_2(\tau)=\{P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)\}
$$

- English: second-order reduced descriptor
- 한국어: 2차 축약 기술자
- Physical definition: same-load non-retracing을 비교하기 위한 reduced state
- Status: DEFINITION

## 2. Density-shape symbols / 밀도 형상 기호

### $\lambda_*$

- English: reference spacing point
- 한국어: 기준 간격점
- Mathematical definition: arbitrary reference point in a smooth interval where $P>0$ and $\Theta>0$
- Physical definition: integration origin only; no independent physical state
- Status: DEFINITION

### $\eta$

- English: spacing integration variable
- 한국어: 간격 적분변수
- Mathematical definition: dummy variable in the density-shape integral
- Status: DEFINITION

### $\mathcal N_P(\tau)$

$$
\int_0^\infty P(\lambda,\tau)\,d\lambda=1
$$

- English: density normalization factor
- 한국어: 확률밀도 정규화 계수
- Mathematical definition: positive integration factor fixed by total probability mass
- Status: DEFINITION

## 3. Initial-value symbols / 초기값 기호

### $P_0(\lambda)$

$$
P_0(\lambda)=P(\lambda,\tau_0)
$$

- English: initial spacing density
- 한국어: 초기 간격밀도
- Status: DEFINITION

### $u_0(\lambda)$

$$
u_0(\lambda)=u(\lambda,\tau_0)
$$

- English: initial conditional mean spacing rate
- 한국어: 초기 조건부 평균 간격속도
- Status: DEFINITION

### $\Theta_0(\lambda)$

$$
\Theta_0(\lambda)=\Theta(\lambda,\tau_0)
$$

- English: initial conditional spacing-rate variance
- 한국어: 초기 조건부 간격속도 분산
- Status: DEFINITION

### $\tau_0$

- English: initial nondimensional time
- 한국어: 초기 무차원 시간
- Mathematical definition: lower temporal bound for initial data and integral representations
- Status: DEFINITION

## 4. Characteristic symbols / 특성곡선 기호

### $s$

- English: dummy time variable
- 한국어: 시간 적분변수
- Mathematical definition: integration variable in $[\tau_0,\tau]$
- Status: DEFINITION

### $r$

- English: inner dummy time variable
- 한국어: 내부 시간 적분변수
- Mathematical definition: integration variable used inside $\mathcal I_u$
- Status: DEFINITION

### $X(s;\alpha)$

$$
\frac{dX}{ds}=u(X(s),s)
$$

$$
X(\tau_0)=\alpha
$$

- English: characteristic curve
- 한국어: 특성곡선
- Mathematical definition: integral curve of the reduced mean-flow field
- Physical definition: mean probability-transport path in spacing space
- Status: DEFINITION

### $\alpha$

- English: characteristic label
- 한국어: 특성곡선 초기 라벨
- Mathematical definition: initial coordinate satisfying $X(\tau_0;\alpha)=\alpha$
- Status: DEFINITION

### $\mathcal I_u$

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
$$

- English: accumulated mean-rate gradient
- 한국어: 누적 평균속도 구배
- Physical definition: cumulative local compression or dilation along a characteristic
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
Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)\,d\lambda
$$

- English: instantaneous nonabsorbing tail mass
- 한국어: 순간 비흡수 꼬리질량
- Mathematical definition: current probability mass above $\lambda_c$
- Physical definition: instantaneous tail only; not cumulative first passage
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
$$

$$
h_\tau=-\frac{d}{d\tau}\ln S
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
\widehat P_b(\lambda,\tau)=\frac{P_b(\lambda,\tau)}{S(\tau)}
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
D_{\mathrm{irr}}(t)=\int_0^t\dot D_{\mathrm{irr}}(t')\,dt'
$$

- English: accumulated irreversible dissipation
- 한국어: 누적 비가역 소산
- Unit: J
- Status: DEFINITION

### $r_j^{\mathrm{irr}}$

미래에 실제 microscopic irreversible node force가 유도될 경우 사용할 기호다.

$$
\dot D_{\mathrm{irr}}^*=-\sum_jr_j^{\mathrm{irr}}\dot x_j
$$

비가역 소산으로 해석하려면

$$
\dot D_{\mathrm{irr}}^*\ge0
$$

이어야 한다.

- English: irreversible node force
- 한국어: 비가역 노드 힘
- Status: OPEN

### $W_{\mathrm{ext}}^{\mathrm{cyc}}$

$$
W_{\mathrm{ext}}^{\mathrm{cyc}}=\Delta E_{\mathrm{mech}}^{\mathrm{cyc}}+D_{\mathrm{irr}}^{\mathrm{cyc}}
$$

- English: external work over one cycle
- 한국어: 한 사이클 외력 일
- Status: CONDITIONAL on a physically valid irreversible mechanism
