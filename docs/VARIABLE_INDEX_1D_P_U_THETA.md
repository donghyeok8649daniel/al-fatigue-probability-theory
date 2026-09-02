# Variable & Mathematical Symbol Index / 수학 기호·변수 Index — active 1D $P$–$u$–$\Theta$ theory

This file is the **authoritative bilingual symbol dictionary** for the active theory.  
이 파일은 현재 이론의 **기준 영·한 수학기호 사전**이다.

A symbol is not considered fully defined unless this index contains its defining equation, mathematical meaning, physical meaning, unit/scaling, status, and dependencies.  
기호는 이 Index에 정의식, 수학적 의미, 물리적 의미, 단위/스케일, 상태, 선행 의존성이 모두 기록되어야 완전히 정의된 것으로 본다.

## 0. Mandatory definition rule / 기호 정의 강제 규칙

Every future symbol MUST be added here with the following fields.  
앞으로 새 기호를 만들 때는 반드시 아래 항목을 동시에 추가한다.

1. **Symbol / 기호** — exact LaTeX glyph.
2. **Equation definition / 수식적 정의** — defining relation, operator identity, or integral.
3. **English term / 영문 명칭**.
4. **Korean term / 한국어 명칭**.
5. **Mathematical definition / 수학적 정의** — scalar/vector/tensor/measure/operator role, domain, conditioning, normalization, etc.
6. **Physical definition / 물리적 정의** — represented physical quantity/process.
7. **Unit or scaling / 단위 또는 스케일**.
8. **Status / 상태** — MODEL, DEFINITION, EXACT, CONDITIONAL, OPEN.
9. **Dependencies / 선행 정의** — symbols/equations required before this symbol is meaningful.

If a defining equation exists, prose alone is not accepted as a definition.  
수식 정의가 존재하면 문장 설명만으로 정의하지 않는다.

The same glyph must not silently carry two unrelated meanings.  
동일한 기호가 서로 다른 수학적 객체를 암묵적으로 뜻하지 않도록 한다.

## 1. Status vocabulary / 식·기호 상태 분류

| Status | English definition | 한국어 정의 |
|---|---|---|
| **MODEL** | adopted reduced physical model or constitutive choice | 채택한 축약 물리모델 또는 구성 선택 |
| **DEFINITION** | mathematical definition chosen by notation | 표기법으로 선택한 수학적 정의 |
| **EXACT** | exact identity under the stated model/measure/boundary conditions | 명시된 모델·측도·경계조건 아래 정확히 성립 |
| **CONDITIONAL** | exact only under an additional explicitly stated condition | 추가 조건이 만족될 때만 정확히 성립 |
| **OPEN** | mathematical slot exists but the required physical law/calibration is not yet derived | 수학적 자리는 있으나 필요한 물리법칙/보정이 아직 미정 |

## 2. Physical scales and loading / 물리 스케일과 하중

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $t$ | independent variable | physical time | 물리시간 | real time coordinate | laboratory elapsed time | s | DEFINITION | none |
| $t_0$ | $\displaystyle t_0=\sqrt{\frac{m_a a_0}{EA_0}}$ | atomic mechanical time scale | 원자 기계 시간척도 | positive time normalization constant | characteristic inertia/stiffness time of the reduced chain | s | DEFINITION under calibration | $m_a,a_0,E,A_0$ |
| $\tau$ | $\displaystyle \tau=t/t_0$ | nondimensional time | 무차원 시간 | dimensionless independent variable | physical time measured in atomic mechanical units | 1 | DEFINITION | $t,t_0$ |
| $m_a$ | input parameter | represented atomic/repeat mass | 대표 원자/반복질량 | positive scalar parameter | inertia assigned to one represented microscopic repeat | kg | input / calibration | none |
| $a_0$ | equilibrium reference | equilibrium spacing | 평형 기준간격 | positive length scale | equilibrium normal spacing represented by $\lambda=1$ | m | input / calibration | none |
| $E$ | input parameter | reference Young modulus | 기준 영률 | positive stress scale | macroscopic elastic calibration scale | Pa | empirical input | none |
| $A_0$ | calibration parameter | effective reference area | 유효 기준면적 | positive area mapping scalar | converts 1D normalized force to physical force/stress | m$^2$ | OPEN calibration quantity | $E$ for force mapping |
| $U_{\rm ref}$ | $U_{\rm ref}=EA_0a_0$ | reference energy | 기준에너지 | energy normalization constant | physical energy represented by unit normalized bond-energy scale | J | DEFINITION | $E,A_0,a_0$ |
| $F_{\rm ref}$ | $F_{\rm ref}=EA_0$ | reference force | 기준힘 | force normalization constant | physical force represented by unit normalized end force | N | DEFINITION | $E,A_0$ |
| $\sigma_n(t)$ | prescribed function | applied normal stress | 가해진 수직응력 | scalar loading history | external tensile/compressive stress | Pa | input history | $t$ |
| $\sigma_m$ | $\displaystyle \sigma_m=\frac{\sigma_{\max}+\sigma_{\min}}2$ when cycle extrema are used | mean stress | 평균응력 | cycle mean scalar | mean normal fatigue stress | Pa | DEFINITION | $\sigma_n$ |
| $\sigma_a$ | $\displaystyle \sigma_a=\frac{\sigma_{\max}-\sigma_{\min}}2$ | stress amplitude | 응력진폭 | half-range scalar | cyclic normal stress amplitude | Pa | DEFINITION | $\sigma_n$ |
| $f$ | cycles per second | loading frequency | 하중주파수 | positive scalar frequency | laboratory fatigue cycling frequency | Hz | input | $t$ |
| $\omega$ | $\omega=2\pi f$ | physical angular frequency | 물리 각주파수 | angular rate | angular frequency of applied loading | rad/s | DEFINITION | $f$ |
| $\omega^*$ | $\omega^*=\omega t_0=2\pi f t_0$ | nondimensional angular frequency | 무차원 각주파수 | dimensionless frequency | loading frequency relative to microscopic time scale | 1 | DEFINITION | $f,t_0$ |
| $F_{\rm ext}(t)$ | prescribed end force | external end force | 외부 끝단 힘 | scalar boundary load | tensile force applied to the loaded end | N | input history | $t$ |
| $q(\tau)$ | $\displaystyle q=\frac{F_{\rm ext}}{EA_0}=\frac{\sigma_n}{E}$ | nondimensional end force | 무차원 끝단 하중 | scalar boundary forcing in normalized ODE | applied normal stress written in chain force units | 1 | DEFINITION under current bridge | $F_{\rm ext},E,A_0,\sigma_n$ |

## 3. Microscopic chain geometry and dynamics / 미시 사슬 기하와 동역학

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $M$ | positive integer | number of represented spacings | 대표 간격 수 | size of the finite spacing population | number of nearest-neighbour gaps in the reduced chain | count | DEFINITION | none |
| $i,j,k,\ell$ | integer indices | microscopic indices | 미시 인덱스 | discrete node/spacing indices | labels nodes, spacings, or matrix entries | 1 | DEFINITION | $M$ |
| $x_j(\tau)$ | state coordinate | normalized node position | 무차원 노드 위치 | scalar component of microscopic configuration vector | position of node/atom $j$ along the loading axis, measured in $a_0$ units | 1 | DEFINITION | $\tau,a_0$ |
| $\lambda_i$ | $\lambda_i=x_i-x_{i-1}$ | normalized spacing | 무차원 층간거리 | scalar nearest-neighbour spacing coordinate | local normal opening/stretch between adjacent represented nodes | 1 | DEFINITION | $x_i,x_{i-1}$ |
| $a_i$ | $a_i=a_0\lambda_i$ | physical spacing | 물리적 층간거리 | dimensional version of $\lambda_i$ | physical local inter-layer/inter-node separation | m | DEFINITION | $a_0,\lambda_i$ |
| $\boldsymbol\lambda$ | $\boldsymbol\lambda=(\lambda_1,\ldots,\lambda_M)^T$ | spacing vector | 간격 벡터 | vector in $\mathbb R^M$ | complete normal-spacing configuration of the chain | 1 | DEFINITION | $\lambda_i,M$ |
| $c_i$ | $c_i=\dot\lambda_i=d\lambda_i/d\tau$ | spacing rate | 간격 변화율 | first derivative of spacing coordinate | local opening/closing rate | 1 per $\tau$ | DEFINITION | $\lambda_i,\tau$ |
| $\boldsymbol c$ | $\boldsymbol c=\dot{\boldsymbol\lambda}$ | spacing-rate vector | 간격 변화율 벡터 | vector in $\mathbb R^M$ | all local opening/closing rates | 1 per $\tau$ | DEFINITION | $\boldsymbol\lambda$ |
| $\ddot\lambda_i$ | $d^2\lambda_i/d\tau^2$ | spacing acceleration | 간격 가속도 | second derivative of spacing coordinate | local acceleration of normal opening/closing | 1 per $\tau^2$ | DEFINITION | $\lambda_i,\tau$ |
| $m$ | model exponent, $m>n>1$ | repulsive exponent | 반발 지수 | positive scalar exponent | controls short-range LJ repulsion | 1 | MODEL parameter | $n$ |
| $n$ | model exponent, $1<n<m$ | attractive exponent | 인력 지수 | positive scalar exponent | controls long-range LJ attraction | 1 | MODEL parameter | $m$ |
| $\phi(\lambda)$ | $\displaystyle \phi=\frac{\lambda^{-m}}{m(m-n)}-\frac{\lambda^{-n}}{n(m-n)}$ | normalized generalized-LJ energy | 무차원 generalized-LJ 에너지 | scalar potential function on $\lambda>0$ | adopted nearest-neighbour normal interaction energy | 1 | MODEL | $m,n,\lambda$ |
| $\phi'(\lambda)$ | $\displaystyle \phi'=\frac{\lambda^{-n-1}-\lambda^{-m-1}}{m-n}$ | normalized bond-force coordinate | 무차원 결합 힘 좌표 | first derivative $d\phi/d\lambda$ | force-like response conjugate to spacing | 1 | EXACT derivative | $\phi$ |
| $\phi''(\lambda)$ | $\displaystyle \phi''=\frac{(m+1)\lambda^{-m-2}-(n+1)\lambda^{-n-2}}{m-n}$ | tangent stiffness | 접선강성 | second derivative $d^2\phi/d\lambda^2$ | local normal tangent stiffness of one represented interaction | 1 | EXACT derivative | $\phi$ |
| $V^*$ | $V^*=\sum_{i=1}^M\phi(\lambda_i)$ | normalized configurational energy | 무차원 배치에너지 | scalar function on $\mathbb R_+^M$ | recoverable nearest-neighbour potential energy of the chain | 1 | EXACT under MODEL | $\phi,\lambda_i$ |
| $T^*$ | $\displaystyle T^*=\frac12\sum_{j=1}^{M}\dot x_j^2$ | normalized kinetic energy | 무차원 운동에너지 | positive quadratic form in node velocities | kinetic energy of moving represented nodes | 1 | EXACT | $x_j$ |
| $E_{\rm mech}^*$ | $E_{\rm mech}^*=T^*+V^*$ | normalized mechanical energy | 무차원 기계에너지 | sum of kinetic and configurational energy | recoverable total mechanical energy | 1 | EXACT | $T^*,V^*$ |
| $\mathbf L$ | $L_{jk}=1$ for $k\le j$, else $0$ | cumulative-sum matrix | 누적합 행렬 | lower-triangular linear map with $\boldsymbol x=\mathbf L\boldsymbol\lambda$ | converts local spacings to node positions | 1 | DEFINITION | $M$ |
| $\mathbf G_\lambda$ | $\mathbf G_\lambda=\mathbf L^T\mathbf L$ | spacing-coordinate mass metric | 간격좌표 질량 메트릭 | symmetric positive-definite metric matrix | accounts for shared node inertia when using spacing rates as generalized velocities | 1 | EXACT kinematics | $\mathbf L$ |
| $(G_\lambda)_{k\ell}$ | $M-\max(k,\ell)+1$ | metric component | 질량 메트릭 성분 | matrix entry of $\mathbf G_\lambda$ | cross-inertial weight between spacing rates $c_k,c_\ell$ | 1 | EXACT | $M,k,\ell$ |

Microscopic equations / 미시 운동식:

$$
\boxed{
\ddot x_j=\phi'(\lambda_{j+1})-\phi'(\lambda_j),
\qquad j=1,\ldots,M-1
}
$$

$$
\boxed{
\ddot x_M=-\phi'(\lambda_M)+q(\tau)
}
$$

$$
\boxed{
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
}
$$

for bulk spacings.

## 4. Empirical probability and phase-space state / 경험적 확률·위상공간 상태

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\mathsf I$ | $\Pr(\mathsf I=i)=1/M$ | uniformly sampled spacing index | 균등 표본 간격 인덱스 | discrete random index on $\{1,\ldots,M\}$ | spatial counting device that turns one deterministic chain into a one-point empirical distribution | 1 | DEFINITION | $M$ |
| $c$ | phase-space coordinate | spacing-rate coordinate | 간격속도 좌표 | continuous variable conjugate to the rate dimension of $F$ | possible local spacing rate | 1 per $\tau$ | DEFINITION | $c_i$ |
| $\lambda$ | phase-space coordinate | spacing coordinate | 간격 좌표 | continuous variable on $\lambda>0$ | possible local normalized spacing | 1 | DEFINITION | $\lambda_i$ |
| $\delta(\cdot)$ | Dirac distribution | Dirac delta | 디랙 델타 | distribution satisfying $\int f(x)\delta(x-x_0)dx=f(x_0)$ | represents discrete microscopic states as an empirical measure | inverse unit of argument | DEFINITION | none |
| $F_M$ | $\displaystyle F_M=\frac1M\sum_i\delta(\lambda-\lambda_i)\delta(c-c_i)$ | empirical phase-space measure | 경험적 위상공간 측도 | probability measure/distribution on $(\lambda,c)$ generated by finite states | exact one-point population of spacing/rate states in the represented chain | density in $(\lambda,c)$ | DEFINITION | $M,\lambda_i,c_i,\delta$ |
| $P_M$ | $\displaystyle P_M(\lambda,\tau)=\int F_Mdc=\frac1M\sum_i\delta(\lambda-\lambda_i)$ | empirical spacing measure | 경험적 간격 측도 | marginal measure of $F_M$ over $c$ | distribution of local spacings in the finite chain | inverse $\lambda$ | DEFINITION | $F_M$ |
| $\mathcal G_M$ | $\displaystyle \mathcal G_M=\frac1M\sum_i\ddot\lambda_i\delta(\lambda-\lambda_i)\delta(c-c_i)$ | empirical acceleration flux | 경험적 가속도 플럭스 | signed distribution carrying acceleration weight in phase space | microscopic acceleration content of each spacing/rate state | acceleration × density | DEFINITION | $\ddot\lambda_i,F_M$ |
| $F$ | smooth/coarse representation of $F_M$ | phase-space density | 위상공간 밀도 | normalized density/measure on $(\lambda,c)$ | joint population of spacing and spacing-rate states | density in $(\lambda,c)$ | DEFINITION | $F_M$ conceptually |
| $P$ | $\displaystyle P(\lambda,\tau)=\int F(\lambda,c,\tau)dc$ | spacing marginal density | 간격 주변확률밀도 | marginal density on $\lambda$ | probability/spatial fraction density of local spacing states | inverse $\lambda$ | DEFINITION | $F$ |
| $A(\lambda,c,\tau)$ | $\displaystyle A=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c]$ | conditional phase-space acceleration | 조건부 위상공간 가속도 | conditional expectation field on $(\lambda,c,\tau)$ | mean microscopic spacing acceleration at fixed spacing and rate | 1/$\tau^2$ | DEFINITION | $F$, microscopic dynamics |
| $J$ | $J=Pu$ | spacing-space probability current | 간격공간 확률류 | scalar flux in $\lambda$-space | transport rate of probability/spatial population across spacing states | density × rate | DEFINITION | $P,u$ |

Exact projected transport / 정확한 투영 수송식:

$$
\boxed{
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
}
$$

## 5. Conditional moments and hierarchy / 조건부 모멘트와 계층

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $u(\lambda,\tau)$ | $u=\mathbb E[c\mid\lambda]$ | conditional mean spacing rate | 조건부 평균 간격속도 | first conditional moment of $c$ | mean local opening/closing rate among states having spacing $\lambda$ | 1/$\tau$ | DEFINITION | $F,P$ |
| $\Theta(\lambda,\tau)$ | $\Theta=\mathbb E[(c-u)^2\mid\lambda]=\operatorname{Var}(c\mid\lambda)$ | conditional spacing-rate variance | 조건부 간격속도 분산 | second conditional central moment | unresolved spread of opening/closing rates at fixed spacing | 1/$\tau^2$ | DEFINITION | $F,P,u$ |
| $C_3(\lambda,\tau)$ | $C_3=\mathbb E[(c-u)^3\mid\lambda]$ | third conditional central moment | 3차 조건부 중심모멘트 | third conditional central moment of rate | skewness-carrying part of local spacing-rate population | 1/$\tau^3$ | DEFINITION | $F,P,u$ |
| $\mathcal A(\lambda,\tau)$ | $\mathcal A=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]$ | one-point conditional acceleration | 1점 조건부 가속도 | conditional mean of microscopic acceleration at fixed spacing | average tendency of states at spacing $\lambda$ to accelerate open/closed | 1/$\tau^2$ | DEFINITION | microscopic dynamics, $P$ |
| $\Psi(\lambda,\tau)$ | $\Psi=\operatorname{Cov}(c,\ddot\lambda\mid\lambda)=\mathbb E[(c-u)\ddot\lambda\mid\lambda]$ | spacing-rate/acceleration covariance | 간격속도–가속도 공분산 | conditional covariance source term | correlation between unresolved rate fluctuations and microscopic acceleration fluctuations | 1/$\tau^3$ | DEFINITION | $c,u,\ddot\lambda$ |
| $D_\tau$ | $D_\tau=\partial_\tau+u\partial_\lambda$ | material derivative in spacing space | 간격공간 물질미분 | derivative along mean spacing-space flow | rate of change following mean probability transport | 1/$\tau$ operator | DEFINITION | $u$ |
| $R_r$ | $\displaystyle R_r=\int c^rFdc$ | raw rate-moment density | 원시 속도모멘트 밀도 | $r$th raw moment density | retained/projection moment of spacing-rate population | density × rate$^r$ | DEFINITION | $F,r$ |
| $B_r$ | $\displaystyle B_r=\int c^{r-1}AFdc$ | acceleration-moment source | 가속도 모멘트 소스 | source in raw moment hierarchy | acceleration contribution to $r$th rate moment | density × rate$^{r-1}$ × acceleration | DEFINITION | $F,A,r$ |
| $r$ | $r=0,1,2,\ldots$ | moment order | 모멘트 차수 | nonnegative integer | labels hierarchy order | 1 | DEFINITION | none |

Exact hierarchy / 정확한 계층식:

$$
\boxed{
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r
}
$$

with

$$
R_0=P,
\qquad
R_1=Pu,
\qquad
R_2=P(u^2+\Theta),
\qquad
R_3=P(u^3+3u\Theta+C_3)
$$

Exact reduced equations / 정확한 축약식:

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0
}
$$

$$
\boxed{
D_\tau u=\mathcal A-\frac1P\partial_\lambda(P\Theta)
}
$$

$$
\boxed{
D_\tau\Theta+2\Theta\partial_\lambda u+\frac1P\partial_\lambda(PC_3)=2\Psi
}
$$

The zero-right-hand-side $\Theta$ equation is only CONDITIONAL on $\Psi=0$.  
$\Psi=0$을 두는 짧은 $\Theta$ 식은 추가 조건이 있을 때만 정확하다.

## 6. Neighbour correlation objects / 이웃 상관 객체

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\lambda'$ | integration variable | neighbour-spacing variable | 이웃 간격 적분변수 | dummy spacing coordinate | possible adjacent spacing | 1 | DEFINITION | $\lambda$ |
| $P_2^+(\lambda,\lambda',\tau)$ | joint density with central marginal $P$ | central/right neighbour joint density | 중심–오른쪽 이웃 결합밀도 | ordered two-spacing joint density | correlation between a spacing and its right neighbour | inverse spacing$^2$ | DEFINITION | $P$ |
| $P_2^-(\lambda,\lambda',\tau)$ | joint density with central marginal $P$ | central/left neighbour joint density | 중심–왼쪽 이웃 결합밀도 | ordered two-spacing joint density | correlation between a spacing and its left neighbour | inverse spacing$^2$ | DEFINITION | $P$ |
| $F_2^+(\lambda,c,\lambda',\tau)$ | ordered joint density | central-rate/right-neighbour density | 중심속도–오른쪽이웃 결합밀도 | joint density over central spacing, central rate, and right spacing | source for rate–acceleration covariance from right neighbour | density in $(\lambda,c,\lambda')$ | DEFINITION | $F$ conceptually |
| $F_2^-(\lambda,c,\lambda',\tau)$ | ordered joint density | central-rate/left-neighbour density | 중심속도–왼쪽이웃 결합밀도 | joint density over central spacing, central rate, and left spacing | source for rate–acceleration covariance from left neighbour | density in $(\lambda,c,\lambda')$ | DEFINITION | $F$ conceptually |
| $m_+$ | $\displaystyle m_+=\frac1P\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)d\lambda'$ | conditional right-neighbour force mean | 조건부 오른쪽 이웃 힘 평균 | conditional expectation of $\phi'(\lambda_{i+1})$ | mean right-neighbour force contribution at fixed central spacing | force coordinate | EXACT | $P,P_2^+,\phi'$ |
| $m_-$ | $\displaystyle m_-=\frac1P\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)d\lambda'$ | conditional left-neighbour force mean | 조건부 왼쪽 이웃 힘 평균 | conditional expectation of $\phi'(\lambda_{i-1})$ | mean left-neighbour force contribution at fixed central spacing | force coordinate | EXACT | $P,P_2^-,\phi'$ |

Bulk conditional acceleration:

$$
\boxed{
\mathcal A_{\rm bulk}=m_++m_- -2\phi'(\lambda)
}
$$

Bulk acceleration covariance:

$$
\boxed{
\Psi_{\rm bulk}
=\frac1P\iint(c-u)\phi'(\lambda')[F_2^++F_2^-]dc\,d\lambda'
}
$$

No independence assumption is used.

## 7. Density-shape and normalization symbols / 확률밀도 형상·정규화 기호

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\lambda_*$ | arbitrary point in a smooth positive-support interval | reference spacing point | 기준 간격점 | lower integration reference for log-density reconstruction | no independent physical content; fixes integration origin | 1 | DEFINITION | support of $P,\Theta$ |
| $\eta$ | dummy integration variable | spacing integration variable | 간격 적분변수 | local dummy variable inside $\lambda$-integral | no independent physical state | 1 | DEFINITION | $\lambda$ |
| $\mathcal N_P(\tau)$ | fixed by $\int P d\lambda=1$ | density normalization factor | 확률밀도 정규화 계수 | positive scalar integration constant | ensures total one-point probability mass equals one | inverse normalization as required | DEFINITION | $P,\Theta$ |

Exact shape relation:

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta
}
$$

$$
\boxed{
P(\lambda,\tau)
=\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}{\Theta(\eta,\tau)}d\eta
\right]
}
$$

The divided formula is valid only on smooth regions with $P>0$ and $\Theta>0$.

## 8. Full-flow integral representation / 전체 흐름 적분표현

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\Gamma$ | $\Gamma=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)$ | full microscopic chain state | 전체 미시 사슬 상태 | point in the finite-dimensional state space | complete node-position/velocity state needed to propagate the deterministic chain | mixed normalized state | DEFINITION | $x_j,\dot x_j,M$ |
| $\Gamma_0$ | $\Gamma(\tau_0)=\Gamma_0$ | initial full state | 초기 전체상태 | initial condition in full state space | microscopic initial condition of one realization | mixed normalized state | DEFINITION | $\Gamma,\tau_0$ |
| $\tau_0$ | initial time | initial nondimensional time | 초기 무차원 시간 | lower temporal bound | time at which initial state/measure is prescribed | 1 | DEFINITION | $\tau$ |
| $\Phi_{\tau,\tau_0}^{q}$ | $\Gamma(\tau)=\Phi_{\tau,\tau_0}^{q}(\Gamma_0)$ | deterministic flow map | 결정론적 흐름 사상 | map from initial full state to later full state under forcing $q$ | microscopic LJ evolution operator | map | EXACT under MODEL | microscopic ODE, $q,\Gamma_0$ |
| $\mu_0(d\Gamma_0)$ | $\int\mu_0(d\Gamma_0)=1$ | initial full-state measure | 초기 전체상태 측도 | normalized probability measure on admissible $\Gamma_0$ | distribution of initial microscopic/specimen realizations; may reduce to a Dirac measure | probability measure | DEFINITION; physical choice OPEN | $\Gamma_0$ |
| $\Lambda_i(\tau;\Gamma_0)$ | $\Lambda_i=x_i-x_{i-1}$ evaluated along $\Phi^q$ | trajectory spacing projection | 궤적 간격 투영 | scalar observable of the full flow | spacing history of represented bond/layer $i$ for a realization | 1 | EXACT projection | $\Phi^q,\Gamma_0$ |
| $C_i(\tau;\Gamma_0)$ | $C_i=d\Lambda_i/d\tau$ | trajectory spacing rate | 궤적 간격속도 | derivative of $\Lambda_i$ along the full flow | opening/closing rate history of spacing $i$ | 1/$\tau$ | EXACT | $\Lambda_i$ |
| $A_i(\tau;\Gamma_0)$ | $A_i=d^2\Lambda_i/d\tau^2$ | trajectory spacing acceleration | 궤적 간격가속도 | second derivative of $\Lambda_i$ along the full flow | acceleration history of spacing $i$ | 1/$\tau^2$ | EXACT | $\Lambda_i$ |

Exact push-forward:

$$
\boxed{
F(\lambda,c,\tau)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\tau;\Gamma_0)]
\delta[c-C_i(\tau;\Gamma_0)]\mu_0(d\Gamma_0)
}
$$

and therefore

$$
\boxed{
P(\lambda,\tau)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\tau;\Gamma_0)]\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
Pu
=\frac1M\sum_i\int C_i\delta[\lambda-\Lambda_i]\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
P(u^2+\Theta)
=\frac1M\sum_i\int C_i^2\delta[\lambda-\Lambda_i]\mu_0(d\Gamma_0)
}
$$

These formulas show that lack of an autonomous three-field closure is not lack of an exact integral representation.

## 9. Volterra and characteristic symbols / 볼테라·특성곡선 기호

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $s$ | dummy time variable inside time integrals | integration-time variable | 시간 적분변수 | dummy temporal variable | no separate physical state | 1 in nondimensional time | DEFINITION | $\tau$ |
| $X(s;\alpha)$ | $\displaystyle \frac{dX}{ds}=u(X,s),\;X(\tau_0)=\alpha$ | spacing-space characteristic | 간격공간 특성곡선 | characteristic curve of the continuity field $u$ | path followed by mean one-point probability transport in spacing space | 1 | DEFINITION / EXACT characteristic construction | $u,\tau_0$ |
| $\alpha$ | $X(\tau_0;\alpha)=\alpha$ | characteristic label | 특성곡선 초기 라벨 | initial coordinate labeling a characteristic | initial spacing-space point of mean transport path | 1 | DEFINITION | $X$ |
| $\mathcal I_u(s;\alpha)$ | $\displaystyle \mathcal I_u=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)dr$ | accumulated velocity-gradient integral | 누적 속도구배 적분 | scalar path integral along $X$ | cumulative local compression/dilation of spacing-space probability flow | 1 | DEFINITION | $u,X$ |
| $S_\Theta$ | $\displaystyle S_\Theta=2\Psi-\frac1P\partial_\lambda(PC_3)$ | $\Theta$ source | $\Theta$ 소스항 | inhomogeneous source in characteristic $\Theta$ ODE | generation/removal of spacing-rate variance from acceleration covariance and third-moment transport | 1/$\tau^3$ | DEFINITION | $\Psi,P,C_3$ |

Exact characteristic forms:

$$
\boxed{
P(X(\tau),\tau)=P_0(\alpha)e^{-\mathcal I_u(\tau;\alpha)}
}
$$

$$
\boxed{
\Theta(X(\tau),\tau)
=e^{-2\mathcal I_u(\tau;\alpha)}
\left[
\Theta_0(\alpha)
+\int_{\tau_0}^{\tau}e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)ds
\right]
}
$$

## 10. G1–G2 observables / G1–G2 관측량

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\bar\lambda$ | $\displaystyle \bar\lambda=\int_0^\infty\lambda P(\lambda,\tau)d\lambda$ | mean normalized spacing | 평균 무차원 간격 | first moment of spacing density | mean local normal spacing relative to $a_0$ | 1 | G1 DEFINITION | $P$ |
| $\bar a$ | $\bar a=a_0\bar\lambda$ | mean physical spacing | 평균 물리간격 | dimensional rescaling of $\bar\lambda$ | mean local inter-layer separation | m | G1 DEFINITION | $a_0,\bar\lambda$ |
| $\Delta\phi$ | $\Delta\phi(\lambda)=\phi(\lambda)-\phi(1)$ | reference-subtracted bond energy | 기준차감 결합에너지 | scalar energy difference from equilibrium state | recoverable configurational energy above equilibrium per represented spacing | 1 | DEFINITION | $\phi$ |
| $\bar U$ | $\displaystyle \bar U=U_{\rm ref}\int\Delta\phi P d\lambda$ | mean intrinsic configurational energy | 평균 고유 배치에너지 | expectation of reference-subtracted potential energy | mean recoverable microscopic configurational energy per represented spacing | J | G2 DEFINITION; EXACT under active energy model | $U_{\rm ref},\Delta\phi,P$ |

Equivalent full-flow forms:

$$
\boxed{
\bar a=\frac{a_0}{M}\sum_i\int\Lambda_i\,\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
\bar U=\frac{U_{\rm ref}}{M}\sum_i\int[\phi(\Lambda_i)-\phi(1)]\mu_0(d\Gamma_0)
}
$$

## 11. G3 irreversibility symbols / G3 비가역성 기호

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\dot D_{\rm irr}$ | $\dot D_{\rm irr}\ge0$ | irreversible dissipation power | 비가역 소산율 | nonnegative scalar power functional | rate at which mechanical energy is irreversibly removed from recoverable storage | W physically | G3 OPEN physically | future irreversible mechanism |
| $D_{\rm irr}$ | $\displaystyle D_{\rm irr}(t)=\int_0^t\dot D_{\rm irr}(t')dt'$ | accumulated irreversible dissipation | 누적 비가역 소산 | time integral of nonnegative dissipation rate | permanently dissipated energy | J | DEFINITION | $\dot D_{\rm irr}$ |
| $E_{\rm hyst}$ | $\displaystyle E_{\rm hyst}(t)=\int_0^t\dot D_{\rm irr}(t')dt'$ | hysteresis accumulation observable | 히스테리시스 누적에너지 | currently identified with accumulated irreversible dissipation | irreversible fatigue-history energy observable | J | G3 DEFINITION | $\dot D_{\rm irr}$ |
| $r_j^{\rm irr}$ | $\displaystyle \dot D_{\rm irr}^*=-\sum_jr_j^{\rm irr}\dot x_j$ if introduced | irreversible node force | 비가역 노드 힘 | generalized force component constrained to nonpositive mechanical power | future physical force producing irreversible loss | normalized force | OPEN | future mechanism, $x_j$ |

Current conservative baseline:

$$
\boxed{
\dot D_{\rm irr}=0,
\qquad
E_{\rm hyst}=0
}
$$

## 12. G4 first-passage and survival symbols / G4 최초통과·생존 기호

| Symbol | Equation definition | English term | 한국어 명칭 | Mathematical definition | Physical definition | Unit / scaling | Status | Dependencies |
|---|---|---|---|---|---|---|---|---|
| $\lambda_c$ | $\phi''(\lambda_c)=0$ | local mechanical stability threshold | 국소 기계적 안정성 임계간격 | first positive tangent-stiffness-loss root in the active potential branch | operational normal-opening initiation threshold | 1 | MODEL-based mechanical criterion | $\phi''$ |
| $\tau_i^c$ | $\displaystyle \tau_i^c=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}$ | local first-passage time | 국소 최초통과 시간 | first hitting time of failure set | first time spacing $i$ reaches local tangent-stiffness loss | 1 in nondimensional time | DEFINITION | $\lambda_i,\lambda_c$ |
| $\chi_i(\tau)$ | $\chi_i=\mathbf1_{\{\tau<\tau_i^c\}}$ | local survival indicator | 국소 생존 지시함수 | Bernoulli indicator of not-yet-hit state | whether spacing $i$ remains uninitiated | 0 or 1 | DEFINITION | $\tau_i^c$ |
| $F_b(\lambda,c,\tau)$ | survivor subdensity on $0<\lambda<\lambda_c$ | survivor phase-space subdensity | 생존 위상공간 부분밀도 | sub-probability density after absorbing first passage | intact spacing/rate population | phase-space subdensity | DEFINITION | $F,\lambda_c$ |
| $P_b(\lambda,\tau)$ | $P_b=\int F_bdc$ | survivor spacing subdensity | 생존 간격 부분밀도 | marginal subdensity of $F_b$ | intact spacing population | inverse $\lambda$ | DEFINITION | $F_b$ |
| $S(\tau)$ | $\displaystyle S=\int_0^{\lambda_c}P_b(\lambda,\tau)d\lambda$ | local survival mass | 국소 생존질량 | total sub-probability mass remaining unabsorbed | fraction/probability of represented local states not yet initiated | 1 | G4 DEFINITION | $P_b$ |
| $j_{\rm esc}$ | $\displaystyle j_{\rm esc}=\int_0^\infty cF_b(\lambda_c^-,c,\tau)dc$ | escape flux | 탈출 플럭스 | outward probability current through kinetic absorbing boundary | rate at which intact local states first cross the instability threshold | 1/$\tau$ | EXACT under absorbing formulation | $F_b,\lambda_c$ |
| $F_{\rm ci}^{\rm local}$ | $F_{\rm ci}^{\rm local}=1-S$ | cumulative local initiation fraction | 누적 국소 균열개시 비율 | complement of local survival | cumulative fraction of local states that first-passed | 1 | G4 DEFINITION | $S$ |
| $h$ | $\displaystyle h=\frac{j_{\rm esc}}S=-\frac{d}{d\tau}\ln S$ | initiation hazard | 균열개시 위험률 | conditional event rate given survival | instantaneous local initiation rate among survivors | 1/$\tau$ | DEFINITION / EXACT when smooth | $S,j_{\rm esc}$ |
| $\widehat P_b$ | $\widehat P_b=P_b/S$ | survivor-conditioned density | 생존조건부 밀도 | normalized conditional density given survival | spacing distribution among intact states only | inverse $\lambda$ | DEFINITION | $P_b,S$ |
| $S_{\rm local}$ | $\displaystyle \frac1M\sum_i\int\mathbf1[\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c]\mu_0(d\Gamma_0)$ | full-flow local survival | 전체흐름 국소 생존 | path-functional expectation averaged over spacing index and $\mu_0$ | probability/fraction that a local represented spacing has never crossed threshold | 1 | EXACT given $\mu_0$ | $\Lambda_i,\mu_0,\lambda_c$ |
| $\tau_{\rm spec}^c$ | $\tau_{\rm spec}^c=\min_i\tau_i^c$ for one realization | specimen first-initiation time | 시편 최초 균열개시 시간 | minimum first-hitting time over represented spacings | first local initiation anywhere in one represented specimen | 1 in nondimensional time | DEFINITION | $\tau_i^c$ |
| $S_{\rm spec}$ | $\displaystyle \int\mathbf1[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c]\mu_0(d\Gamma_0)$ | specimen survival probability | 시편 생존확률 | path-functional probability over the full-state ensemble | probability that no represented local spacing has initiated by $\tau$ | 1 | EXACT formula once $\mu_0$ is declared; physical $\mu_0$ OPEN | $\Lambda_i,\mu_0,\lambda_c$ |
| $F_{\rm ci}^{\rm spec}$ | $F_{\rm ci}^{\rm spec}=1-S_{\rm spec}$ | specimen crack-initiation CDF | 시편 균열개시 누적확률 | cumulative distribution of specimen first-initiation event | probability a specimen has initiated by time $\tau$ | 1 | DEFINITION | $S_{\rm spec}$ |
| $\mathbf1[\cdot]$ | indicator equals 1 if condition true, else 0 | indicator function | 지시함수 | measurable event indicator | counts whether a trajectory survives/fails | 0 or 1 | DEFINITION | event condition |

For the kinetic absorbing boundary:

$$
\boxed{
F_b(\lambda_c,c,\tau)=0\quad\text{for incoming }c<0
}
$$

and

$$
\boxed{
\dot S=-j_{\rm esc}
}
$$

## 13. Core mathematical term glossary / 핵심 수학용어 사전

| English term | 한국어 | Equation / formal definition | Mathematical meaning | Physical meaning |
|---|---|---|---|---|
| generalized Lennard–Jones potential | generalized 레너드–존스 퍼텐셜 | $\phi(\lambda)$ above | nonlinear scalar potential | reduced normal interaction law |
| normalized spacing | 무차원 간격 | $\lambda=a/a_0$ | dimensionless coordinate | local opening/stretch relative to equilibrium |
| spacing rate | 간격 변화율 | $c=d\lambda/d\tau$ | phase-space velocity coordinate | local opening/closing rate |
| empirical measure | 경험적 측도 | finite Dirac sum | probability measure generated from finite states | direct mechanical state distribution without assuming a named PDF |
| spatial counting measure | 공간 계수 측도 | $\Pr(\mathsf I=i)=1/M$ | uniform discrete measure over spacing index | samples local states inside one deterministic chain |
| phase space | 위상공간 | $(\lambda,c)$ | state space of coordinate and rate | local normal spacing plus its rate |
| phase-space density | 위상공간 밀도 | $F(\lambda,c,\tau)$ | joint density on $(\lambda,c)$ | population of opening/rate states |
| marginal density | 주변확률밀도 | $P=\int Fdc$ | projection of joint density | spacing distribution after discarding rate coordinate |
| conditional expectation | 조건부 기댓값 | $\mathbb E[Y\mid X=x]$ | mean on a conditional fiber/event | average microscopic quantity among states sharing a retained variable |
| conditional variance | 조건부 분산 | $\operatorname{Var}(c\mid\lambda)=\Theta$ | second conditional central moment | rate spread hidden by the spacing marginal |
| covariance | 공분산 | $\operatorname{Cov}(X,Y\mid Z)$ | conditional joint fluctuation moment | dynamical correlation between fluctuations |
| probability current | 확률류 | $J=Pu$ | flux in reduced coordinate space | transport of local-state population across spacing values |
| material derivative | 물질미분 | $D_\tau=\partial_\tau+u\partial_\lambda$ | derivative along reduced mean flow | time change following mean spacing-space transport |
| raw moment | 원시 모멘트 | $R_r=\int c^rFdc$ | moment about zero | rate-statistics hierarchy |
| central moment | 중심 모멘트 | $\mathbb E[(c-u)^r\mid\lambda]$ | moment about conditional mean | local spread/skewness of spacing rates |
| moment hierarchy | 모멘트 계층 | $\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r$ | lower moments depend on higher moments/sources | information-loss chain created by projection |
| closure | 폐쇄/클로저 | express higher statistics using retained fields | reduction assumption | needed only for an autonomous low-order solver |
| push-forward measure | 푸시포워드 측도 | $\mu_\tau=(\Phi^q_{\tau,\tau_0})_\#\mu_0$ conceptually | image measure under deterministic flow | evolution of an initial microscopic population through mechanics |
| deterministic flow map | 결정론적 흐름 사상 | $\Gamma(\tau)=\Phi^q_{\tau,\tau_0}(\Gamma_0)$ | solution operator of the finite ODE | maps initial microscopic state to later state |
| Volterra integral equation | 볼테라 적분방정식 | present value contains $\int_{\tau_0}^{\tau}(\cdots)ds$ | causal time-integral equation | explicit history form of exact balances |
| characteristic curve | 특성곡선 | $dX/ds=u(X,s)$ | curve tangent to the transport field | mean probability-transport path in spacing space |
| integrating factor | 적분인자 | $e^{2\mathcal I_u}$ in the $\Theta$ ODE | converts a first-order linear ODE to exact integral form | separates compression/dilation from variance source |
| first passage | 최초통과 | $\tau^c=\inf\{\tau:X(\tau)\in\partial\Omega_f\}$ | first hitting time | first local mechanical instability event |
| absorbing boundary | 흡수경계 | no inflow from failed side | boundary removing first-passed mass | prevents failed states from re-entering intact population |
| survival probability | 생존확률 | $S=\Pr(\tau^c>\tau)$ | tail probability of first-hitting time | probability/fraction not yet initiated |
| hazard | 위험률 | $h=-d\ln S/d\tau$ | conditional event rate | instantaneous initiation rate among survivors |
| history dependence | 이력의존성 | same $q$ with different $(P,u,\Theta)$ | non-single-valued reduced response to instantaneous forcing | loading/unloading microscopic states can differ at the same stress |
| irreversible dissipation | 비가역 소산 | $\dot D_{\rm irr}\ge0$ | monotone loss functional | energy permanently removed from recoverable mechanical storage |
| configurational energy | 배치에너지 | $\bar U=U_{\rm ref}\mathbb E[\Delta\phi]$ | expectation of potential-energy difference | mean recoverable microscopic structural energy |
| local first-passage fraction | 국소 최초통과 비율 | $1-S_{\rm local}$ | average event indicator over local states | fraction/probability of represented local spacings already initiated |
| specimen survival | 시편 생존확률 | $S_{\rm spec}$ path integral | probability of no first passage anywhere | probability entire represented specimen remains uninitiated |

## 14. Critical interpretation rules / 핵심 해석 규칙

1. $P$ is mechanically generated; it is not assumed Gaussian, Weibull, Boltzmann, or any other named family.
2. $\Theta$ is an exact conditional spacing-rate variance, not a fitted damage parameter.
3. $\frac12(u^2+\Theta)$ is a local spacing-rate quadratic moment, **not** the complete chain kinetic energy; the latter uses $\mathbf G_\lambda$ and cross-spacing correlations.
4. The exact $\Theta$ balance contains $2\Psi$. Omitting $\Psi$ requires an explicit condition/closure.
5. $(P,u,\Theta)$ is history-bearing but not an autonomous closed Markov state.
6. Non-retracing at the same applied stress does not by itself imply $\dot D_{\rm irr}>0$.
7. The reduced PDE hierarchy is not the only solution representation: $P,u,\Theta,C_3,\Psi$ all have exact full-flow projection integrals.
8. $1-S_{\rm local}$ is not automatically specimen-to-specimen crack probability.
9. $S_{\rm spec}$ has an exact path-integral formula once $\mu_0$ is declared; the **physical choice of $\mu_0$ and correlation scale is open**.
10. A long-range/zeta energy must not be substituted into G2 while retaining nearest-neighbour equations of motion and then called mechanically exact.
11. Numerical normalized frequency $\omega^*$ must be mapped through $t_0$ before any laboratory-frequency interpretation.
12. New symbols must obey the mandatory definition rule in Section 0.
