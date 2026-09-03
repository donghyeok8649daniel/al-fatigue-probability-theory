# Variable and Mathematical Symbol Index / 수학 기호·변수 Index

이 문서는 활성 1D $P$-$u$-$\Theta$ 이론의 기준 영·한 기호 사전이다.

각 기호는 수식적 정의, 수학적 의미, 물리적 의미, 단위 또는 스케일, 상태, 선행 정의를 함께 기록한다.

## 1. Status vocabulary / 상태 분류

- MODEL: 채택한 축약 물리모델 또는 구성 선택
- DEFINITION: 표기와 수학적 정의
- EXACT: 명시한 모델과 조건 아래 정확식
- CONDITIONAL: 추가 조건 아래 정확식
- OPEN: 수학적 자리는 있으나 물리법칙 또는 calibration이 미완성

## 2. Time, length, force, and energy scales / 시간·길이·힘·에너지 스케일

### $t$

- English: physical time
- 한국어: 물리시간
- Mathematical definition: independent real time variable
- Physical definition: laboratory elapsed time
- Unit: s
- Status: DEFINITION

### $a_0$

- English: equilibrium reference spacing
- 한국어: 평형 기준간격
- Mathematical definition: positive reference length
- Physical definition: spacing represented by $\lambda=1$
- Unit: m
- Status: input or calibration

### $t_0$

$$
t_0=\sqrt{\frac{m_a a_0}{EA_0}}
$$

- English: microscopic mechanical time scale
- 한국어: 미시 기계 시간척도
- Mathematical definition: positive time normalization constant
- Physical definition: inertia-stiffness time of the reduced chain
- Unit: s
- Status: DEFINITION under current calibration
- Dependencies: $m_a,a_0,E,A_0$

### $\tau$

$$
\tau=\frac{t}{t_0}
$$

- English: nondimensional time
- 한국어: 무차원 시간
- Unit: 1
- Status: DEFINITION
- Dependencies: $t,t_0$

### $\lambda$

$$
\lambda=\frac{a}{a_0}
$$

- English: normalized spacing coordinate
- 한국어: 무차원 간격 좌표
- Mathematical definition: positive dimensionless coordinate
- Physical definition: local normal separation relative to $a_0$
- Unit: 1
- Status: DEFINITION
- Dependencies: $a,a_0$

### $F_{\mathrm{ref}}$

$$
F_{\mathrm{ref}}=EA_0
$$

- English: reference force
- 한국어: 기준힘
- Unit: N
- Status: DEFINITION

### $U_{\mathrm{ref}}$

$$
U_{\mathrm{ref}}=EA_0a_0
$$

- English: reference energy
- 한국어: 기준에너지
- Unit: J
- Status: DEFINITION

### $q(\tau)$

$$
q(\tau)=\frac{F_{\mathrm{ext}}(t)}{EA_0}=\frac{\sigma_n(t)}{E}
$$

- English: nondimensional end load
- 한국어: 무차원 끝단 하중
- Mathematical definition: scalar forcing in the normalized boundary equation
- Physical definition: applied normal stress written in chain force units
- Unit: 1
- Status: DEFINITION under current bridge

## 3. Microscopic chain / 미시 사슬

### $M$

- English: number of represented spacings
- 한국어: 대표 간격 수
- Mathematical definition: positive integer
- Physical definition: number of nearest-neighbour gaps in the finite chain
- Unit: count
- Status: DEFINITION

### $x_j(\tau)$

- English: normalized node position
- 한국어: 무차원 노드 위치
- Mathematical definition: scalar microscopic coordinate
- Physical definition: node position along the loading axis in $a_0$ units
- Unit: 1
- Status: DEFINITION

### $\lambda_i(\tau)$

$$
\lambda_i=x_i-x_{i-1}
$$

- English: normalized local spacing
- 한국어: 무차원 국소 간격
- Mathematical definition: nearest-neighbour difference coordinate
- Physical definition: local normal opening between adjacent represented nodes
- Unit: 1
- Status: DEFINITION

### $a_i(t)$

$$
a_i=a_0\lambda_i
$$

- English: physical local spacing
- 한국어: 물리적 국소 간격
- Unit: m
- Status: DEFINITION

### $c_i(\tau)$

$$
c_i=\frac{d\lambda_i}{d\tau}
$$

- English: spacing rate
- 한국어: 간격 변화율
- Mathematical definition: first derivative of $\lambda_i$
- Physical definition: local opening or closing rate
- Unit: 1 per $\tau$
- Status: DEFINITION

### $\phi(\lambda)$

$$
\phi(\lambda)=\frac{\lambda^{-m}}{m(m-n)}-\frac{\lambda^{-n}}{n(m-n)}
$$

- English: normalized generalized-LJ energy
- 한국어: 무차원 generalized-LJ 에너지
- Mathematical definition: scalar potential for $\lambda>0$
- Physical definition: adopted nearest-neighbour normal interaction energy
- Unit: 1
- Status: MODEL

### $V^*$

$$
V^*=\sum_{i=1}^{M}\phi(\lambda_i)
$$

- English: normalized configurational energy
- 한국어: 무차원 배치에너지
- Physical definition: recoverable nearest-neighbour potential energy of the chain
- Unit: 1
- Status: EXACT under MODEL

### $T^*$

$$
T^*=\frac{1}{2}\sum_{j=1}^{M}\dot x_j^2
$$

- English: normalized kinetic energy
- 한국어: 무차원 운동에너지
- Status: EXACT

### $E_{\mathrm{mech}}^*$

$$
E_{\mathrm{mech}}^*=T^*+V^*
$$

- English: normalized mechanical energy
- 한국어: 무차원 기계에너지
- Status: EXACT

### $G_\lambda$

$$
G_\lambda=L^TL
$$

$$
(G_\lambda)_{k\ell}=M-\max(k,\ell)+1
$$

- English: spacing-coordinate mass metric
- 한국어: 간격좌표 질량 메트릭
- Mathematical definition: symmetric positive-definite metric induced by cumulative node positions
- Physical definition: accounts for shared node inertia in spacing coordinates
- Status: EXACT kinematics

## 4. Probability and phase space / 확률과 위상공간

### $F_M(\lambda,c,\tau)$

$$
F_M(\lambda,c,\tau)=\frac{1}{M}\sum_i\delta(\lambda-\lambda_i)\delta(c-c_i)
$$

- English: empirical phase-space measure
- 한국어: 경험적 위상공간 측도
- Mathematical definition: finite normalized Dirac sum on $(\lambda,c)$
- Physical definition: mechanically generated population of spacing and rate states
- Status: DEFINITION

### $P_M(\lambda,\tau)$

$$
P_M(\lambda,\tau)=\int F_M(\lambda,c,\tau)\,dc
$$

- English: empirical spacing measure
- 한국어: 경험적 간격 측도
- Status: DEFINITION

### $F(\lambda,c,\tau)$

- English: smooth phase-space density
- 한국어: 연속 위상공간 밀도
- Mathematical definition: normalized joint density or smooth representation of $F_M$
- Physical definition: joint population of spacing and spacing-rate states
- Status: DEFINITION

### $P(\lambda,\tau)$

$$
P(\lambda,\tau)=\int F(\lambda,c,\tau)\,dc
$$

- English: spacing marginal density
- 한국어: 간격 주변확률밀도
- Mathematical definition: marginal of $F$ over $c$
- Physical definition: probability or spatial fraction density of local spacing states
- Unit: inverse $\lambda$
- Status: DEFINITION

### $A(\lambda,c,\tau)$

$$
A(\lambda,c,\tau)=\mathrm{E}[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c,\tau]
$$

- English: conditional phase-space acceleration
- 한국어: 조건부 위상공간 가속도
- Mathematical definition: conditional expectation at fixed spacing and rate
- Physical definition: mean microscopic spacing acceleration at a phase-space point
- Status: DEFINITION

## 5. Conditional moments / 조건부 모멘트

### $u(\lambda,\tau)$

$$
u(\lambda,\tau)=\mathrm{E}[c\mid\lambda,\tau]
$$

- English: conditional mean spacing rate
- 한국어: 조건부 평균 간격속도
- Mathematical definition: first conditional moment of $c$
- Physical definition: mean opening or closing rate at fixed spacing
- Status: DEFINITION

### $\Theta(\lambda,\tau)$

$$
\Theta(\lambda,\tau)=\mathrm{Var}(c\mid\lambda,\tau)
$$

$$
\Theta(\lambda,\tau)=\mathrm{E}[(c-u)^2\mid\lambda,\tau]
$$

- English: conditional spacing-rate variance
- 한국어: 조건부 간격속도 분산
- Mathematical definition: second conditional central moment
- Physical definition: unresolved rate spread at fixed spacing
- Status: DEFINITION

### $C_3(\lambda,\tau)$

$$
C_3(\lambda,\tau)=\mathrm{E}[(c-u)^3\mid\lambda,\tau]
$$

- English: third conditional central moment
- 한국어: 3차 조건부 중심모멘트
- Status: DEFINITION

### $\mathcal A(\lambda,\tau)$

$$
\mathcal A(\lambda,\tau)=\mathrm{E}[\ddot\lambda_i\mid\lambda_i=\lambda,\tau]
$$

- English: one-point conditional acceleration
- 한국어: 1점 조건부 가속도
- Status: DEFINITION

### $\Psi(\lambda,\tau)$

$$
\Psi(\lambda,\tau)=\mathrm{Cov}(c,\ddot\lambda\mid\lambda,\tau)
$$

$$
\Psi(\lambda,\tau)=\mathrm{E}[(c-u)\ddot\lambda\mid\lambda,\tau]
$$

- English: spacing-rate and acceleration covariance
- 한국어: 간격속도-가속도 공분산
- Mathematical definition: conditional covariance source in the variance balance
- Physical definition: coupling between unresolved rate fluctuations and acceleration
- Status: DEFINITION

### $D_\tau$

$$
D_\tau=\partial_\tau+u\partial_\lambda
$$

- English: material derivative in spacing space
- 한국어: 간격공간 물질미분
- Status: DEFINITION

### $R_r$

$$
R_r(\lambda,\tau)=\int c^rF\,dc
$$

- English: raw rate-moment density
- 한국어: 원시 속도모멘트 밀도
- Status: DEFINITION

### $B_r$

$$
B_r(\lambda,\tau)=\int c^{r-1}AF\,dc
$$

- English: acceleration-moment source
- 한국어: 가속도 모멘트 소스
- Status: DEFINITION

## 6. Neighbour correlation objects / 이웃 상관 객체

### $P_2^+$ and $P_2^-$

- English: ordered central-neighbour joint spacing densities
- 한국어: 순서가 있는 중심-이웃 결합 간격밀도
- Mathematical definition: two-spacing joint densities with central marginal $P$
- Physical definition: retain left/right neighbour correlation
- Status: DEFINITION

### $F_2^+$ and $F_2^-$

- English: central-rate-neighbour joint densities
- 한국어: 중심속도-이웃 결합밀도
- Physical definition: source objects for $\Psi$
- Status: DEFINITION

### $m_+$ and $m_-$

$$
m_+=\frac{1}{P}\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)\,d\lambda'
$$

$$
m_-=\frac{1}{P}\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)\,d\lambda'
$$

- English: conditional neighbour-force means
- 한국어: 조건부 이웃 힘 평균
- Status: EXACT definitions from the joint densities

## 7. Integral representation / 적분 표현

### $\Gamma$

$$
\Gamma=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
$$

- English: full microscopic chain state
- 한국어: 전체 미시 사슬 상태
- Status: DEFINITION

### $\Gamma_0$

- English: initial full state
- 한국어: 초기 전체상태
- Status: DEFINITION

### $\Phi^q_{\tau,\tau_0}$

$$
\Gamma(\tau)=\Phi^q_{\tau,\tau_0}(\Gamma_0)
$$

- English: deterministic flow map
- 한국어: 결정론적 흐름 사상
- Status: EXACT under the microscopic ODE

### $\mu_0$

$$
\int\mu_0(d\Gamma_0)=1
$$

- English: initial full-state measure
- 한국어: 초기 전체상태 측도
- Mathematical definition: normalized probability measure over admissible $\Gamma_0$
- Physical definition: ensemble of initial microscopic or specimen realizations
- Status: DEFINITION; physical choice OPEN

### $\Lambda_i$

$$
\Lambda_i(\tau;\Gamma_0)=x_i(\tau;\Gamma_0)-x_{i-1}(\tau;\Gamma_0)
$$

- English: trajectory spacing projection
- 한국어: 궤적 간격 투영
- Status: EXACT projection

### $C_i$

$$
C_i=\frac{d\Lambda_i}{d\tau}
$$

- English: trajectory spacing rate
- 한국어: 궤적 간격속도
- Status: EXACT

### $A_i$

$$
A_i=\frac{d^2\Lambda_i}{d\tau^2}
$$

- English: trajectory spacing acceleration
- 한국어: 궤적 간격가속도
- Status: EXACT

### $X(s;\alpha)$

$$
\frac{dX}{ds}=u(X(s),s)
$$

$$
X(\tau_0)=\alpha
$$

- English: spacing-space characteristic
- 한국어: 간격공간 특성곡선
- Status: DEFINITION

### $\mathcal I_u$

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
$$

- English: accumulated velocity-gradient integral
- 한국어: 누적 속도구배 적분
- Status: DEFINITION

### $S_\Theta$

$$
S_\Theta=2\Psi-\frac{1}{P}\partial_\lambda(PC_3)
$$

- English: variance source in characteristic form
- 한국어: 특성곡선형 분산 소스
- Status: DEFINITION

## 8. G1-G4 observables / G1-G4 관측량

### $\bar\lambda$

$$
\bar\lambda(\tau)=\int_0^\infty\lambda P(\lambda,\tau)\,d\lambda
$$

- English: mean normalized spacing
- 한국어: 평균 무차원 간격
- Status: G1 DEFINITION

### $\bar a$

$$
\bar a(t)=a_0\bar\lambda(t/t_0)
$$

- English: mean physical spacing
- 한국어: 평균 물리간격
- Unit: m
- Status: G1 DEFINITION

### $\bar U$

$$
\bar U=U_{\mathrm{ref}}\int[\phi(\lambda)-\phi(1)]P(\lambda,\tau)\,d\lambda
$$

- English: mean intrinsic configurational energy
- 한국어: 평균 고유 배치에너지
- Unit: J
- Status: G2 DEFINITION and EXACT under the active energy model

### $\dot D_{\mathrm{irr}}$

$$
\dot D_{\mathrm{irr}}\ge0
$$

- English: irreversible dissipation power
- 한국어: 비가역 소산율
- Status: G3 OPEN physically

### $E_{\mathrm{hyst}}$

$$
E_{\mathrm{hyst}}(t)=\int_0^t\dot D_{\mathrm{irr}}(t')\,dt'
$$

- English: accumulated irreversible history energy
- 한국어: 누적 비가역 이력에너지
- Status: G3 DEFINITION

### $\lambda_c$

$$
\phi''(\lambda_c)=0
$$

- English: local mechanical stability threshold
- 한국어: 국소 기계적 안정성 임계간격
- Status: G4 MODEL-based criterion

### $\tau_i^c$

$$
\tau_i^c=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}
$$

- English: local first-passage time
- 한국어: 국소 최초통과 시간
- Status: G4 DEFINITION

### $F_b$

- English: survivor phase-space subdensity
- 한국어: 생존 위상공간 부분밀도
- Mathematical definition: subdensity on $0<\lambda<\lambda_c$ with kinetic absorbing boundary
- Status: DEFINITION

### $S(\tau)$

$$
S(\tau)=\int_0^{\lambda_c}\int_{-\infty}^{\infty}F_b(\lambda,c,\tau)\,dc\,d\lambda
$$

- English: local survival mass
- 한국어: 국소 생존질량
- Status: G4 DEFINITION

### $j_{\mathrm{esc}}$

$$
j_{\mathrm{esc}}=\int_0^\infty cF_b(\lambda_c^-,c,\tau)\,dc
$$

- English: escape flux
- 한국어: 탈출 플럭스
- Status: EXACT under the absorbing formulation

### $S_{\mathrm{spec}}$

$$
S_{\mathrm{spec}}(\tau)=\int I\left[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

- English: specimen survival probability
- 한국어: 시편 생존확률
- Mathematical definition: path-functional expectation over the full-state measure
- Physical definition: probability that no represented local spacing has first-passed by $\tau$
- Status: EXACT once $\mu_0$ is declared; physical $\mu_0$ remains OPEN

## 9. Interpretation rules / 해석 규칙

1. $P$는 mechanics에서 생성되며 named PDF를 가정하지 않는다.
2. $\Theta$는 fitted damage parameter가 아니라 conditional spacing-rate variance다.
3. $\frac12(u^2+\Theta)$는 전체 chain kinetic energy가 아니다.
4. exact variance balance에는 $2\Psi$가 포함된다.
5. $(P,u,\Theta)$는 history-bearing reduced descriptor이지만 autonomous closed state는 아니다.
6. 같은 응력에서 non-retracing이 나타나는 것만으로 $\dot D_{\mathrm{irr}}>0$를 뜻하지 않는다.
7. reduced PDE hierarchy가 닫히지 않아도 full deterministic flow의 exact integral representation은 존재한다.
8. local first-passage fraction은 자동으로 specimen-to-specimen probability가 아니다.
