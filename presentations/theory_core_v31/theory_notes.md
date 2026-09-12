# LJ/Bessel 무한 격자와 확률 동역학

Al 피로 연구의 수학적·물리적 구조

## 01. LJ/Bessel 무한 격자와 확률 동역학

Al 피로 연구의 수학적·물리적 구조
에너지 지형, 구성 이동, 개구 최초통과

팀 세미나용 이론 정리. 기존 two-row 생산 기준과 full FCC / 활성 계면 / 공간 결함 연구 확장을 구분한다. 이 자료는 실험 항복강도·피로수명·물리 Hz의 검증 완료를 주장하지 않는다.

출처:

- [AGENTS.md](../../AGENTS.md)
- [solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md](../../solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md)

## 02. 서로 다른 모델 계층

TwoRowLJ: 무한 하부 원자열과 상부 좌표의 축약 기준
Full FCC: 삼각 원자평면과 ABC 무한 적층의 벌크 기준
활성 계면: 한 계면의 개구와 두 성분 registry 이동
공간 결함: 비균일 변위장과 전위 코어를 포함하는 연구 확장

더 상세한 기하학이 이전 보정 파라미터의 전이를 보장하지 않는다. 각 단계는 에너지 정규화, 경계조건, 관측량을 새로 확인한다. 기존 생산 PDE와 연구용 계면/전위 계산의 구현 수준을 혼동하지 않는다.

출처:

- [AGENTS.md](../../AGENTS.md)
- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)
- [solver_v1/VECTOR_FCC_CORE_DERIVATION.md](../../solver_v1/VECTOR_FCC_CORE_DERIVATION.md)

## 03. 집단좌표와 확률밀도

    q = (a,s),     P = P(a,s,t)

    s = bn + ξ,     n = floor(s/b + 1/2)

    −b/2 ≤ ξ < b/2

a는 정상 분리, s는 누적 registry 변위다.
ξ는 우물 내부 좌표이며 n은 주기적 registry 우물의 정수 지표다.

s가 한 주기만큼 이동하면 국소 에너지는 같아도 누적 translation은 달라진다. n≠0만으로 실제 Al 소성을 입증하지 않는다. 벡터 registry에서는 scalar floor 분해를 두 방향에 무심코 복사하지 않고 결정격자와 경로를 명시해야 한다.

출처:

- [solver_v1/CONFIGURATIONAL_PLASTICITY.md](../../solver_v1/CONFIGURATIONAL_PLASTICITY.md)

## 04. LJ 쌍 퍼텐셜과 두 원자열

    φ(r) = 4ε_{LJ} [(σ_{LJ}/r)^{12} − (σ_{LJ}/r)^{6}]

    r_{n}^{2} = a^{2} + [(n+1/2)b − s]^{2}

    W(a,s) = 4ε_{LJ}[σ_{LJ}^{12}S_{6} − σ_{LJ}^{6}S_{3}]

Sₚ는 모든 정수 n에 대한 rₙ⁻²ᵖ의 합이다.
σLJ는 쌍 퍼텐셜의 길이 매개변수이며 적용 응력 σ와 구별한다.

현재 two-row W는 한 상부 사이트가 무한 하부 원자열과 갖는 교차 상호작용이다. 기존 normal-chain 에너지를 중복 가산하지 않는다. LJ pair functional form은 환경 항을 확장하더라도 유지한다.

출처:

- [solver_v1/BESSEL_LATTICE_DERIVATION.md](../../solver_v1/BESSEL_LATTICE_DERIVATION.md)
- [solver_v1/lattice_bessel.py](../../solver_v1/lattice_bessel.py)

## 05. Gamma 적분에서 Gaussian 변환으로

    (d^{2}+|x|^{2})^{−p} = Γ(p)^{−1}∫_{0}^{∞} t^{p−1}e^{−t(d²+|x|²)}dt

    ∫_{ℝᴰ} e^{−t|x|²}e^{−ik·x}dx = (π/t)^{D/2}e^{−|k|²/(4t)}

    f̂_{p,D}(k) = [π^{D/2}/Γ(p)]∫_{0}^{∞} t^{p−1−D/2}e^{−d²t−|k|²/(4t)}dt

Fourier 규약은 exp(−ik·x)이며 역변환에 (2π)⁻ᴰ가 붙는다.
d>0, p>D/2에서 영 모드의 적분도 유한하다.

멱함수를 Gamma 적분으로 바꾸면 격자 차원의 역할이 Gaussian 적분의 t^(−D/2)로 명시된다. 이 변환이 row와 plane의 Bessel 차수가 달라지는 이유다. 위 식은 절대수렴 조건에서 적분 순서를 바꾼다.

출처:

- [solver_v1/BESSEL_LATTICE_DERIVATION.md](../../solver_v1/BESSEL_LATTICE_DERIVATION.md)
- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)
- [https://dlmf.nist.gov/10.32#E10](https://dlmf.nist.gov/10.32#E10)

## 06. 차원 D에 따른 Bessel 차수

    f̂_{p,D}(k) = [2π^{D/2}/Γ(p)] (|k|/2d)^{p−D/2}K_{p−D/2}(d|k|)

    f̂_{p,D}(0) = [π^{D/2}Γ(p−D/2)/Γ(p)]d^{D−2p}

D=1인 원자열은 Kₚ₋₁/₂, D=2인 원자평면은 Kₚ₋₁을 준다.
k=0은 별도로 계산해 작은 인수에서의 수치 상쇄를 피한다.

적분 항등식 ∫t^(ν−1)exp(−βt−γ/t)dt=2(γ/β)^(ν/2)Kν(2√βγ)를 적용한다. ν=p−D/2, β=d², γ=|k|²/4. Bessel은 새로운 경험식이 아니라 기존 무한 격자합의 정확한 reciprocal 표현이다.

출처:

- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)
- [https://dlmf.nist.gov/10.32#E10](https://dlmf.nist.gov/10.32#E10)

## 07. 원자열의 Poisson 급수

    S_{p}(a,s) = A_{p}(a) + ∑_{m≥1} B_{p,m}(a)cos θ_{m}

    A_{p} = √π Γ(p−1/2)a^{1−2p}/[bΓ(p)]

    B_{p,m} = [4√π/(bΓ(p))](πm/ab)^{p−1/2}K_{p−1/2}(2πma/b)

θₘ=2πm(1/2−s/b): 반 격자 staggered 위상을 보존한다.
무한 실공간 이웃 수를 임의 cutoff로 정하는 방식과 구별한다.

Poisson 공식의 1/b와 ±m 쌍의 2를 함께 반영해야 계수가 맞는다. m=0은 평균 정상 결합, m≠0은 registry corrugation이다. 실제 코드는 허용오차로 reciprocal 급수를 절단하는 반해석적 평가다.

출처:

- [solver_v1/BESSEL_LATTICE_DERIVATION.md](../../solver_v1/BESSEL_LATTICE_DERIVATION.md)
- [solver_v1/lattice_bessel.py](../../solver_v1/lattice_bessel.py)

## 08. 정상 개구와 registry 장벽의 결합

    K_{ν}(z) ~ √(π/2z)e^{−z}

    S_{p,s} = ∑_{m≥1}(2πm/b)B_{p,m}sin θ_{m}

    K′_{ν}(z) = −[K_{ν−1}(z)+K_{ν+1}(z)]/2

원자열의 corrugation은 대략 exp(−2πma/b)로 줄어든다.
W의 해석적 미분이 drift와 Hessian을 함께 결정한다.

큰 a에서 s 장벽이 약해지는 경향은 같은 LJ 격자합에서 나온다. 점근 지수만으로 정확한 장벽이나 총 절단오차를 결정하지 않는다. a,s 혼합 미분도 동일 계수와 위상을 미분해 얻는다.

출처:

- [solver_v1/BESSEL_LATTICE_DERIVATION.md](../../solver_v1/BESSEL_LATTICE_DERIVATION.md)

## 09. 지수형 환경 밀도의 무한합

    f_{ρ}(r)=ρ_{0}e^{−βρ(r/rₑ−1)}=C_{ρ}e^{−κρr}

    ĝ_{1D}(k)=2aκ_{ρ}K_{1}(aQ)/Q,     Q=√(κ_{ρ}^{2}+k^{2})

    ρ_{parallel}=2C_{ρ}/[exp(κ_{ρ}b)−1]

같은 원자열에서는 자기 자신을 제외한 기하급수를 사용한다.
교차 원자열은 위 Fourier 계수에 같은 Poisson 위상을 적용한다.

환경 밀도 ρ는 확률밀도 P, 질량 밀도와 다른 양이다. Cρ=ρ0 exp(βρ), κρ=βρ/rₑ. 1D transform은 exp(−κr)/r의 2K0(aQ)를 κ로 미분해 검증할 수 있다. ρ0와ρref의 비율에는 gauge가 있다.

출처:

- [solver_v1/ANALYTIC_LJ_EAM_HYBRID.md](../../solver_v1/ANALYTIC_LJ_EAM_HYBRID.md)

## 10. LJ–EAM의 원자별 에너지 계산

    E = (1/2)∑_{i}∑_{j≠i}φ_{LJ}(r_{ij}) + ∑_{i}F(x_{i})

    ρ_{i}=∑_{j≠i}f_{ρ}(r_{ij}),     F(x)=−A√x,     x_{i}=ρ_{i}/ρ_{ref}

    (F∘x)_{αβ}=F″x_{α}x_{β}+F′x_{αβ}

각 원자의 전체 환경 밀도를 먼저 합산한 뒤 F를 적용한다.
두 원자열 셀의 2F와 full FCC의 원자당 F는 정규화가 다르다.

마지막 식 α,β는 집단좌표 미분 지표이고 F′,F″는 dimensionless density x에 대한 미분이다. 원자 지표 i,j와 구별한다. 이 hybrid는 기존 Al EAM 파라미터화를 자동으로 계승하지 않는다. 쌍/embedding 분해의 gauge 때문에 계수의 물리 해석에는 식별성 검사가 필요하다.

출처:

- [solver_v1/ANALYTIC_LJ_EAM_HYBRID.md](../../solver_v1/ANALYTIC_LJ_EAM_HYBRID.md)
- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)

## 11. FCC(111) 평면과 ABC 적층

    b=a_{lat}/√2,     h_{111}=a_{lat}/√3

    a_{1}=b(1,0),     a_{2}=b(1/2,√3/2)

    A_{atomic}=√3 b^{2}/2,     τ=(a_{1}+a_{2})/3

3τ는 평면 격자 벡터이며 A, B, C 배열이 반복된다.
결정축을 회전할 때는 +τ 적층에 맞는 실제 cubic-axis 대응을 쓴다.

기하학 frame e1=(1,−1,0)/√2, e2=(1,1,−2)/√6, e3=(1,1,1)/√3은 직교 정규다. 현재 +τ,+h 생성 격자에 cubic 텐서를 붙일 때 plane_basis_in_stacked_cubic_axes의 rows(−e1,−e2,e3)를 사용한다. 생성 원자위치와 cubic label을 혼동하면 전단 성분 부호가 틀릴 수 있다.

출처:

- [solver_v1/fcc111_geometry.py](../../solver_v1/fcc111_geometry.py)
- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)
- [solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md](../../solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md)

## 12. 삼각 평면의 2D Poisson 합

    S_{p}^{2D}(d,δ)=∑_{R∈Λ}[d^{2}+|R+δ|^{2}]^{−p}

    S_{p}^{2D} = A_{atomic}^{−1}∑_{G∈Λ*}f̂_{p,2}(G)e^{iG·δ}

    b_{1}=(2π/b)(1,−1/√3),     b_{2}=(2π/b)(0,2/√3)

aᵢ·bⱼ=2πδᵢⱼ를 만족하는 reciprocal 격자를 쓴다.
±G를 묶으면 실수 cos(G·δ) corrugation을 얻는다.

reciprocal shell은 h²−hk+k²가 같은 벡터들이다. 첫 shell의 degeneracy는 6이며 완전한 shell을 포함해야 대칭이 보존된다. D=2 Fourier 계수는 앞서 유도한 K_(p−1) 식을 사용한다. 지수 밀도는 Q=√(κ²+|G|²)에 대해 ĝ₂D=2πκ exp(−dQ)(dQ+1)/Q³이다. exp(−κr)/r의 2π exp(−dQ)/Q를 −∂κ로 미분해 얻는다. 모든 평면·층의 밀도를 원자별로 합산한 뒤 embedding을 적용한다.

출처:

- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)
- [solver_v1/fcc111_lattice_sum.py](../../solver_v1/fcc111_lattice_sum.py)

## 13. 무한 적층과 ζ(4), ζ(10)

    S_{p}^{full}=(1/2)Z_{triangle}(p)+∑_{ℓ≥1}S_{p}^{2D}(ℓa,Δ_{ℓ})

    ∑_{ℓ≥1} [π(ℓa)^{2−2p}/(A_{atomic}(p−1))]

    = πa^{2−2p}ζ(2p−2)/[A_{atomic}(p−1)]

LJ의 p=3,6은 각각 ζ(4), ζ(10)을 준다.
원자열 적층의 ζ(5), ζ(11)은 서로 다른 차원 축약에서 나온다.

양·음 층의 대칭이 쌍 에너지의 1/2과 상쇄되어 양의 층 합만 남는다. 동일 평면 항은 self를 뺀 Ztriangle/2다. 비영 reciprocal 모드는 exp(−ℓa|G|)로 감쇠한다. homogeneous Δℓ=ℓ(τ+s e_s)를 local fault와 등치하지 않는다.

출처:

- [solver_v1/FCC111_FULL_STACK_DERIVATION.md](../../solver_v1/FCC111_FULL_STACK_DERIVATION.md)

## 14. 한 계면만 움직이는 기하학

    z_{ℓ}=ℓh  (ℓ≤0),     z_{ℓ}=ℓh+δa  (ℓ≥1)

    u_{ℓ}=0  (ℓ≤0),     u_{ℓ}=u  (ℓ≥1)

    a=h+δa,     u=(u_{1},u_{2}),     u=s e_{s}는 경로 축약

상·하 반결정 내부 간격은 벌크 값 h로 유지한다.
개구 a와 slip u는 같은 활성 계면의 상대 변위다.

전체 층 간격을 동시에 늘리는 homogeneous dilation과 국소 cleavage는 다르다. full registry는 두 성분이며 direct110와 Shockley112 경로는 명시적으로 선택한다. 직선 경로의 최대가 full-vector index-one saddle인지 별도로 확인한다.

출처:

- [solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md](../../solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md)
- [solver_v1/VECTOR_REGISTRY_AND_STRENGTH_AUDIT.md](../../solver_v1/VECTOR_REGISTRY_AND_STRENGTH_AUDIT.md)

## 15. 계면 쌍 에너지의 유한 차분

    ΔE_{pair}=∑_{k≥1}k[w(kh+δa,kτ+u)−w(kh,kτ)]

    D_{r}(x)=∑_{k≥1}k[(k+x)^{−r}−k^{−r}]

    =ζ(r−1,1+x)−xζ(r,1+x)−ζ(r−1)

x=δa/h, r=2p−2이며 w는 plane–plane LJ kernel이다.
k층 떨어진 상·하 평면 쌍은 절단면을 가로질러 k개 존재한다.

동일 반결정 내부 상호작용은 정확히 상쇄된다. 큰 무한 벌크 상수를 수치적으로 빼지 않는다. 평균 모드는 Hurwitz zeta로 합산하고 registry-dependent 꼬리는 기존 reciprocal 급수로 평가한다. 이 계수 k를 누락하면 cleavage 정규화가 달라진다.

출처:

- [solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md](../../solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md)

## 16. 계면 근처의 원자별 embedding 변화

    Δρ_{r}=∑_{k≥r+1}[ρ_{plane}(kh+δa,kτ+u)−ρ_{plane}(kh,kτ)]

    ΔE_{emb}=2∑_{r≥0}[F(ρ_{bulk}+Δρ_{r})−F(ρ_{bulk})]

    W_{int}=ΔE_{pair}+ΔE_{emb},     W_{int}(h,0)=0

r은 계면에서의 원자층 깊이이며 양쪽 원자가 factor 2를 준다.
깊은 층은 벌크 밀도로 복귀하므로 embedding 차분도 수렴한다.

EAM은 원자마다 F를 평가한다. 한 global crystal density에 F를 한 번 적용하거나 층별 F(ρplane)를 합산하면 다른 이론이 된다. 이 scalar 식에 angular 연구 확장을 더할 때도 각 원자의 전체 tensor 환경 합 뒤에 비선형 함수를 적용한다.

출처:

- [solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md](../../solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md)

## 17. 환경의 방향 정보와 대칭

    U_{ijk}=∑_{neighbors}w(r)r_{i}r_{j}r_{k},     v_{i}=∑_{j}U_{ijj}

    Q_{ijk}=U_{ijk}−(v_{i}δ_{jk}+v_{j}δ_{ik}+v_{k}δ_{ij})/5

    I_{3}=||Q||^{2},     Q=0  (반전 대칭 환경)

정상 벌크와 결함 원자는 같은 scalar 밀도라도 방향 환경이 다를 수 있다.
r의 곱은 Fourier 공간에서 i∂/∂G에 대응한다.

STF rank3는 trace가 0이며 반전에서 부호가 바뀐다. 반전 대칭 affine 벌크의 모든 r,−r 쌍은 상쇄한다. 결함과 표면의 odd 환경은 남을 수 있다. 지수밀도 Fourier transform의 G 미분으로 무한 moment 합과 좌표 미분을 계산해 LJ/Bessel 틀을 유지한다.

출처:

- [solver_v1/ANALYTIC_ANGULAR_ENVIRONMENT.md](../../solver_v1/ANALYTIC_ANGULAR_ENVIRONMENT.md)

## 18. 현재의 해석적 환경 확장 계열

    F(x)=−A√x+B(x−1)+C(x−1)^{2}

    E_{env,i}=F(x_{i})+D_{1}x_{i}^{z}I_{1,i}+D_{3}I_{3,i}+K_{3}I_{3,i}^{2}

    +D_{2}g(I_{2,i})+D_{E}g(I_{E,i}),     g(I)=I/(1+αI/N_{2})

각 I는 원자별 무한 환경 moment를 합한 뒤 만든 불변량이다.
현재 탐색은 이 기존 계열 내부에서 수행하며 생산 기본식은 보존한다.

I_E는 두 기존 radial quadrupole의 선형 조합으로 T2g 응답을 소거하고 Eg 응답을 분리한 연구 채널이다. normalization N2는 고정 amplitude gauge이며 specimen 면적/길이가 아니다. 채널별 정규화는 코드와 동일하게 적용한다. K3는 pristine harmonic 자료만으로 식별되지 않는다. 여기의 양의 계수·포화형 가정은 연구 가설이며 검증된 Al 퍼텐셜이라는 뜻이 아니다.

출처:

- [solver_v1/vector_material_calibration.py](../../solver_v1/vector_material_calibration.py)
- [solver_v1/symmetry_resolved_material.py](../../solver_v1/symmetry_resolved_material.py)
- [solver_v1/quadrupole_saturation.py](../../solver_v1/quadrupole_saturation.py)
- [solver_v1/coordination_screening.py](../../solver_v1/coordination_screening.py)

## 19. 계면에 작용하는 외부 일

    t=σ·n,     t_{n}=n·σ·n,     τ_{α}=e_{α}·σ·n

    G_{int}=W_{int}−A_{atomic}[t_{n}(a−h)+τ_{1}u_{1}+τ_{2}u_{2}]

물리 좌표 [m], traction [Pa], 면적 [m²]를 쓰면 외부 일은 [J]다.
결정 방향과 실험 하중 방향의 변환을 명시해야 한다.

이 식은 에너지/원자 계면 셀로 정규화한 활성 계면의 물리적 traction work다. 기존 two-row PDE의 f*=kappa sigma/E 매핑을 σA로 바꾸라는 뜻이 아니다. 동일한 외부 일을 서로 다른 에너지/좌표 정규화에 맞추어 유도해야 한다. A_atomic과 통계적 A_c는 완전히 분리한다.

출처:

- [solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md](../../solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md)
- [solver_v1/KINETICS_LOADING_AND_STRESS_AUDIT.md](../../solver_v1/KINETICS_LOADING_AND_STRESS_AUDIT.md)

## 20. 완화 축방향 접선 보정

    H_{0}δq=cδf,     c=(1,χ)^{T}

    δε=(c^{T}δq)/a_{0}=(c^{T}H_{0}^{−1}c)δf/a_{0}

    κ_{axial}=a_{0}/(c^{T}H_{0}^{−1}c),     f*=κ_{axial}σ/E

무주파수 평형 접선을 σ/E에 맞추는 기존 reduced 규약이다.
새 에너지 모델은 자기 a₀와 Hessian으로 κ를 다시 계산한다.

원래 TwoRowLJ의 a0=.7713438268704838, chi=.2에서 kappa=86.29296488740997을 보존한다. finite-frequency 응답을 따라가게 하려고 kappa를 조절하지 않는다. chi의 결정학 후보는 (a0/h)(e·m)(e·n)이며 물리 좌표 대응이 확인되어야 한다.

출처:

- [solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md](../../solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md)
- [solver_v1/CONFIGURATIONAL_PLASTICITY.md](../../solver_v1/CONFIGURATIONAL_PLASTICITY.md)

## 21. Smoluchowski 확률 방정식

    ∂_{t}P=−∇·J,     J=−M(P∇G+θ∇P)

    M=diag(M_{a},M_{s}),     θ=k_{B}T,     D=θM

에너지 기울기가 drift를, 열항이 확률의 확산을 정한다.
실제 계산은 확률밀도 자체를 결정론적으로 진화시킨다.

θ는 사용한 에너지 단위에서의 kBT다. 생산 모델의 일정 diagonal mobility는 모델 가정이며 원자적 모든 주파수에서의 정확한 운동방정식이라고 주장하지 않는다. 국소 평형과 빠른 운동량 소거를 확인해야 한다. reference MD는 그 검증 데이터 생성 경로이며 생산 PDE를 Newtonian MD로 교체하지 않는다.

출처:

- [solver_v1/probability_pde_2d.py](../../solver_v1/probability_pde_2d.py)
- [solver_v1/PHYSICAL_TIME_AND_MOBILITY.md](../../solver_v1/PHYSICAL_TIME_AND_MOBILITY.md)

## 22. 자유에너지와 비음수 소산

    ℱ[P]=∫PGdq+θ∫P ln P dq,     μ=δℱ/δP

    J=−MP∇μ

    dℱ/dt=∫P∂_{t}Gdq−∫J^{T}M^{−1}J/P dq

마지막 식은 고정 반사 경계와 일정 온도에서 성립한다.
소산이 양수여도 순수 정상 지연만 있을 수 있어 소성과 구별한다.

∇μ=∇G+θ∇lnP를 대입하고 부분적분하면 위 식을 얻는다. 흡수/이동 경계에서는 경계 항이 더 필요하다. 정상 주기, 질량 보존, 일관된 conjugate work일 때만 한 주기의 입력 일이 소산과 같다. 누적 소산 자체를 임의 균열 임계값으로 쓰지 않는다.

출처:

- [solver_v1/probability_pde_2d.py](../../solver_v1/probability_pde_2d.py)
- [solver_v1/CONFIGURATIONAL_PLASTICITY.md](../../solver_v1/CONFIGURATIONAL_PLASTICITY.md)
- [solver_v1/INTERNAL_STEP_WORK_V30.md](../../solver_v1/INTERNAL_STEP_WORK_V30.md)

## 23. 보존형 Scharfetter–Gummel flux

    B(z)=z/(e^{z}−1),     Δ=(G_{R}−G_{L})/θ

    J_{L→R}=(D/Δq)[B(Δ)P_{L}−B(−Δ)P_{R}]

    P_{R}/P_{L}=e^{−Δ}  에서  J_{L→R}=0

같은 내부면 flux를 두 셀에 반대 부호로 반영한다.
명시·암시 적분은 같은 공간 generator를 사용한다.

B(−Δ)=exp(Δ)B(Δ)이므로 이산 Gibbs 정지분포를 보존한다. Δ≈0에서 expm1 또는 급수로 B를 안정 계산한다. positivity, 질량잔차, dt/격자 수렴을 독립 검증한다. 작은 음수 보정은 개구 흡수로 합산하지 않는다.

출처:

- [solver_v1/probability_pde_2d.py](../../solver_v1/probability_pde_2d.py)
- [solver_v1/RESULT_FIELDS.md](../../solver_v1/RESULT_FIELDS.md)

## 24. 생존 앙상블의 변형률 분해

    ⟨X⟩_{S}=∫XPdq/M_{raw},     M_{raw}=∫_{intact}Pdq

    ε_{a}=⟨a−a_{0}⟩_{S}/a_{0},     ε_{ξ}=χ⟨ξ⟩_{S}/a_{0}

    ε_{p}=χb⟨n⟩_{S}/a_{0},     ε_{total}=ε_{a}+ε_{ξ}+ε_{p}

흡수된 상태를 생존 앙상블 평균에 포함하지 않는다.
χ⟨s⟩/a₀=εξ+εp를 매 저장 시점에 검사한다.

정확한 연속 해에서 Mraw=Slocal이다. 수치상 차이는 mass residual이며 조건부 평균 분모는 실제 surviving density의 적분을 쓴다. 생존 질량이 소진되면 조건부 평균은 정의되지 않는다. n/xi 반열린 경계 규약과 음수 s indexing을 보존한다.

출처:

- [solver_v1/CONFIGURATIONAL_PLASTICITY.md](../../solver_v1/CONFIGURATIONAL_PLASTICITY.md)
- [solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md](../../solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md)

## 25. 우물 점유와 순·총 구성 이동

    P_{n}=∫_{Ωₙ}Pdq,     Σ_{n}P_{n}=M_{raw}

    dP_{n}/dt=𝒥_{n−1/2}−𝒥_{n+1/2}−Ȧ_{n}

    𝒥=j^{+}−j^{−},     𝒥_{gross}=j^{+}+j^{−}

양의 flux는 증가하는 s 방향이며 Ȧₙ은 해당 우물의 개구 손실이다.
총 양방향 활동과 방향성 있는 순이동은 다른 관측량이다.

SG의 한쪽률 분해는 정확한 이산 진단이다. gross는 |net|이 아니다. 단, 격자 경계의 확산 recrossing으로 SG gross는 대략 1/Δs에 의존하므로 mesh-independent committed hopping count로 해석하지 않는다. 물리 전이 횟수에는 well core/committor 정의가 추가로 필요하다.

출처:

- [solver_v1/CONFIGURATIONAL_PLASTICITY.md](../../solver_v1/CONFIGURATIONAL_PLASTICITY.md)
- [AGENTS.md](../../AGENTS.md)

## 26. 소성 유동과 선택적 개구의 분리

    d(∑nP_{n})/dt=∑𝒥_{n+1/2}−∑nȦ_{n}

    ε̇_{p,flow}=[χb/(a_{0}S)]∑𝒥_{n+1/2}

    ε̇_{p,open}=[χb/(a_{0}S)][−∑nȦ_{n}+⟨n⟩_{S}∑Ȧ_{n}]

생존 registry 평균은 우물 간 이동과 선택적 개구 손실로 바뀐다.
잔류 소성은 unload와 충분한 zero-stress hold 후에도 수렴해야 한다.

연속 질량 균형 Sdot=−ΣȦn을 사용해 조건부 평균의 미분을 전개한 식이다. 개구 항을 plastic-flow로 부르면 선택적으로 살아남는 집단의 변화를 소성으로 오인한다. 반사 outer s boundary 가정이며 open outer boundary가 있으면 추가항을 명시한다.

출처:

- [solver_v1/CONFIGURATIONAL_PLASTICITY.md](../../solver_v1/CONFIGURATIONAL_PLASTICITY.md)
- [solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md](../../solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md)

## 27. 균열 개시의 최초통과 정의

    A_{abs}(t)=∫_{0}^{t}∫_{Γopen}J·n dΓdt

    P_{init,local}=A_{abs},     S_{local}=1−A_{abs}

    R_{mass}=M_{raw}+A_{abs}−1

개구 dividing surface의 실제 흡수와 수치적 보정을 구별한다.
격자·시간간격·flux 일치 오차보다 작은 신호는 미분해로 표시한다.

이동 dividing surface에서는 상대 경계 flux 또는 swept mass를 일관되게 계산한다. raw 1−Mraw는 반올림 잔차를 포함할 수 있으므로 canonical physical probability로 쓰지 않는다. Ωinit의 순간 점유나 Dacc 임계값은 현재 opening first-passage 정의와 다르다.

출처:

- [solver_v1/RESULT_FIELDS.md](../../solver_v1/RESULT_FIELDS.md)
- [solver_v1/PROBABILITY_AND_PLASTICITY_RESOLUTION.md](../../solver_v1/PROBABILITY_AND_PLASTICITY_RESOLUTION.md)

## 28. 빠른 a의 유한온도 축약

    P_{0}(a,s,t)=p_{s}(s,t)ρ_{eq}(a|s,f)

    ρ_{eq}=e^{−G/θ}/Z_{a},     ℱ_{eff}=−θ ln Z_{a}

    ∂_{t}p_{s}=∂_{s}[M_{s}(p_{s}∂_{s}ℱ_{eff}+θ∂_{s}p_{s})]

고정 반사 a 경계에서 빠른 연산자의 nullspace를 먼저 구한다.
∂sℱeff=⟨Gs⟩와 solvability 조건이 느린 방정식을 준다.

느린 시간에서 Ma/Ms→∞의 선도해다. 임의 Ma 조정이나 a* delta 가정으로 시작하지 않는다. 저온에서는 ℱeff≈G(a*,s)+(θ/2)ln[Gaa(a*,s)/(2πθ)]이며 integration measure의 단위를 고정해야 한다. 축약된 정상 변형률에는 a* 대신 conditional ⟨a⟩를 사용한다.

출처:

- [solver_v1/FAST_A_REDUCTION.md](../../solver_v1/FAST_A_REDUCTION.md)

## 29. 이동 경계와 QSD

    Z_{b}=∫_{L(s)}^{U(s)}e^{−G/θ}da

    ∂_{s}ℱ_{b}=⟨G_{s}⟩_{b}−θ[e^{−G(U)/θ}U′−e^{−G(L)/θ}L′]/Z_{b}

    ℒ_{a}ρ_{QSD}=−k_{open}ρ_{QSD}

truncated Gibbs와 absorbing-process QSD를 구별한다.
국소 sink −kopen ps는 고유모드 분리와 시간척도 검증이 필요하다.

Ga(U)=0이라는 saddle 조건은 Leibniz 경계항을 0으로 만들지 않는다. 흡수 연산자에는 정상적인 영 고유값 Gibbs 상태가 없다. 단순히 Ma→∞로 보내면 escape도 빨라지므로 희귀 탈출/혼합시간의 분리가 필요하다. 현재 reduced crack probability는 생산 검증을 통과하지 않았다.

출처:

- [solver_v1/FAST_A_REDUCTION.md](../../solver_v1/FAST_A_REDUCTION.md)

## 30. 선형 완화와 주파수 응답

    δq̇=−MH_{0}δq+Mcδf,     R=M^{1/2}H_{0}M^{1/2}

    λ_{i}=eig(R),     τ_{i}=1/λ_{i}

    𝒢(ω)=(κ/a_{0})c^{T}(iωI+MH_{0})^{−1}Mc

𝒢(0)=1은 정적 접선 보정에서 나온다.
유한 주파수의 진폭·위상은 M, Hessian, ω가 함께 결정한다.

R은 대칭이므로 positive definite H와M에 대해 양의 relaxation rates를 준다. Was가 0이 아니면 eigenvector를 순수 a/s로 가정하지 않는다. 정상 좌표의 phase lag 자체를 plastic hysteresis라고 부르지 않는다. χ와 전달함수𝒢의 기호도 구별한다.

출처:

- [solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md](../../solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md)

## 31. 길이·에너지·시간의 차원화

    q_{phys}=L_{0}q*,     G_{phys}=E_{0}G*,     t_{phys}=t_{0}t*

    M*=(t_{0}E_{0}/L_{0}^{2})M_{phys},     [M_{phys}]=m²/(J·s)

    t_{0}=M*L_{0}^{2}/(M_{phys}E_{0}),     f_{Hz}=f_{model}/t_{0}

P*=L₀²Pphys이며 θ*=kBT/E₀로 확률과 열항도 변환한다.
에너지·길이 보정만으로 overdamped 시간은 정해지지 않는다.

한 공통 t0에는 Ma와Ms의 물리/모델 비율이 일치해야 한다. 독립 이동도가 기존1:.05와 다르면 상대 이동도 자체를 증거에 따라 검증해야 한다. 원자 질량을 숨은 Newtonian 시간으로 넣지 않는다. 생산 Ma,phys/Ms,phys/t0는 미보정 상태다.

출처:

- [solver_v1/PHYSICAL_TIME_AND_MOBILITY.md](../../solver_v1/PHYSICAL_TIME_AND_MOBILITY.md)

## 32. 같은 집단좌표의 이동도 측정

    C=⟨δqδq^{T}⟩,     K=∫_{0}^{∞}⟨δq(t)δq(0)^{T}⟩dt

    H=θC^{−1},     Γ=θC^{−1}KC^{−1},     M=CK^{−1}C/θ

    χ_{v}(ω)=v^{T}(iωI+MH)^{−1}Mv

위 관계는 일치하는 좌표의 조화·평형·overdamped 가정에서 유도한다.
유한 대역 스펙트럼 S(ω)/2를 K(0)와 곧바로 등치하지 않는다.

B=MH, C(t)=exp(−Bt)C이면 K=B^(−1)C다. Gibbs covariance C=θH^(−1)에서 Γ식과M식을 얻는다. 실제 MD는 운동량과 기억 효과를 갖는다. 저주파 한계, 강제 응답/FDT, 온도, 부호, 진폭, 기록 길이, 결합 모드를 검사해야 생산 이동도를 승인할 수 있다.

출처:

- [solver_v1/MODAL_KINETIC_CALIBRATION_V25.md](../../solver_v1/MODAL_KINETIC_CALIBRATION_V25.md)
- [solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md](../../solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md)

## 33. 좌표 축약에서 생기는 기억 효과

    M_{scalar}=(v^{T}Cv)^{2}/[θ(v^{T}Kv)] ≤ v^{T}Mv

    χ_{v}(0)=v^{T}Cv/θ,     −Imχ_{v}/ω → v^{T}Kv/θ

    H̄=H/N_{p},     M̄=N_{p}M,     θ̄=θ/N_{p}

scalar 저주파 응답과 full-system 투영 이동도는 일반적으로 다르다.
원자수 정규화에서는 drift와 diffusion을 동시에 보존해야 한다.

Cauchy–Schwarz를 K^(−1/2)Cv와K^(1/2)v에 적용하면 부등식을 얻는다. 평면별 gap은 Σq_l=0 제약을 먼저 제거한 독립 subspace에서 역산한다. M만Np배하고θ를 유지하면 확산이Np배 바뀐다. 이는 실제 온도를 바꾸라는 뜻이 아니라 같은 평면 확률법칙의 에너지 정규화다. local atomic-cell PMF로의 전이는 별도 검증 대상이다.

출처:

- [solver_v1/collective_generator_bridge.py](../../solver_v1/collective_generator_bridge.py)
- [solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md](../../solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md)

## 34. 이상강도와 실제 소성의 공간 구조

    γ(d)=W_{int}(d)/A_{atomic}

    E[d]≈E_{elastic}[d]+∫γ(d(x))dA−W_{external}

    d=(opening,slip_{1},slip_{2})

균일 반결정 이동은 결함 없는 coherent 계면의 이상강도를 탐색한다.
국소 slip에는 주변의 탄성 변형과 전위 코어·source 기하학이 필요하다.

실제 항복을 이상강도에 곱하는 보정계수로 만들지 않는다. 같은 에너지에서 공간적 변형 비용을 유도하고 결함의 안정성·구동·유한 source를 검증한다. 지름, 원자 셀 면적, A_c를 자동으로 activation 크기로 쓰지 않는다. 여기의 continuum+misfit은 long-wave 근사로서 원자 코어에서 별도 검사해야 한다.

출처:

- [solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md](../../solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md)
- [solver_v1/MATERIAL_TO_SPECIMEN_STRENGTH.md](../../solver_v1/MATERIAL_TO_SPECIMEN_STRENGTH.md)

## 35. 반무한 벌크 소거와 비국소 kernel

    K(k)=|k|(Z_{+}^{−1}+Z_{−}^{−1})^{−1}

    E_{elastic}=(1/2)∫d̂^{†}K(k)d̂ d²k/(2π)²

    K_{jump}=(BH^{+}B^{†})^{−1}  (이산 harmonic 제약)

Z±는 같은 벌크 탄성 텐서에서 구한 표면 impedance다.
원자적 finite-q kernel에는 이미 포함된 local stiffness를 중복 가산하지 않는다.

두 half-space의 공통 변위를 최소화하면 첫 Schur kernel을 얻는다. H+는 translation zero modes를 처리한 pseudoinverse이고 B는 명시된 displacement-jump 제약이다. continuum K(0)=0과 atomistic local term의 차이를 구별한다. 원자열 코어에서 y,z도 움직이면 Bessel 영 모드의 algebraic tail을 보존해야 한다.

출처:

- [solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md](../../solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md)
- [solver_v1/DISCRETE_FCC_SCREW_DERIVATION.md](../../solver_v1/DISCRETE_FCC_SCREW_DERIVATION.md)
- [solver_v1/VECTOR_FCC_CORE_DERIVATION.md](../../solver_v1/VECTOR_FCC_CORE_DERIVATION.md)

## 36. 유한 pinned source의 정적 장벽

    G[curve]=∫γ_{line}(θ)dℓ−τb A_{swept}

    (γ_{line}+γ″_{line})κ_{curve}=τb

    τ_{c}=2γ_{line}(π/2)/(bL)

L은 실제 source/고정점 기하학이며 통계적 상관 면적과 다르다.
안정 minor arc와 index-one major arc의 에너지 차가 장벽을 준다.

이 식은 양의 anisotropic line stiffness를 갖는 planar local-line 근사 안에서의 유도다. 선 에너지는 같은 벌크/계면 kernel에서 연결하지만 core와 장애물의 물리값이 추가로 필요하다. finite barrier를 얻어도 물리 속도나 실험 항복이 자동 검증되지 않는다. ΔG의 응력 미분으로 activation 형상량이 나오며 임의 characteristic volume을 넣지 않는다.

출처:

- [solver_v1/FINITE_SOURCE_BARRIER_DERIVATION.md](../../solver_v1/FINITE_SOURCE_BARRIER_DERIVATION.md)

## 37. 국소 확률과 시편 확률

    N_{eff}=A_{stressed}/A_{c}

    log S_{spec}=N_{eff}log1p(−P_{local})

    P_{spec}=−expm1[log S_{spec}]

Ac는 통계적으로 독립적인 등가 개시 영역을 가정하는 외부 보정량이다.
국소 미분해 신호를 확대해 인증된 시편 확률로 제시하지 않는다.

Aatomic은 원자 에너지의 면적 정규화, Ac는 확률 상관 면적, Astressed는 시편 하중 면적이다. 세 양을 같은 것으로 쓰지 않는다. strain과well activity에는Neff를 곱하지 않는다. 분산을 맞춘 correlation area가 void/survival 확률의 독립영역 면적과 자동으로 같지 않다.

출처:

- [solver_v1/STATISTICAL_CORRELATION_AREA.md](../../solver_v1/STATISTICAL_CORRELATION_AREA.md)
- [app/specimen_probability.py](../../app/specimen_probability.py)

## 38. 보정과 식별성의 수학적 조건

    r_{k}(θ)=[y_{k}(θ)−y_{k,target}]/s_{k}

    J_{ki}=∂r_{k}/∂θ_{i},     J=UΣV^{T}

    ρ_{0}→cρ_{0},  ρ_{ref}→cρ_{ref}  에서  x=ρ/ρ_{ref} 불변

단위와 온도·방향·완화 조건이 일치하는 관측량만 비교한다.
bulk 접선만 맞춘 후보는 계면 장벽과 finite-q 응답도 검증해야 한다.

밀도 scale gauge를 고정하고 singular values/rank/상관을 보고한다. fixed radial shape에서 에너지 계수에 선형인 행렬을 이용해 exact anchors와minimax residual을 분리할 수 있다. fitting 자료, inspected development 자료, 적합에서 제외한 자료를 명시한다. 작은 training loss도 실험 Al 검증이나 유일한 파라미터의 증거는 아니다.

출처:

- [solver_v1/ALUMINUM_ANALYTIC_CALIBRATION.md](../../solver_v1/ALUMINUM_ANALYTIC_CALIBRATION.md)
- [solver_v1/CORE_INTERFACE_COMPATIBILITY_V23.md](../../solver_v1/CORE_INTERFACE_COMPATIBILITY_V23.md)
- [solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md](../../solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md)

## 39. 이론의 연결과 검증 경계

해석적 기반: LJ 무한합, Poisson 변환, Bessel 계수와 미분
확률 기반: 보존형 Smoluchowski, 생존 조건부 평균, 개구 최초통과
연구 확장: full FCC, 활성 계면, 환경 moment, 공간 전위와 source
남은 물리 검증: Al 재료 적합, 국소 kinetics, 실제 항복·피로, Ac

수학적 항등식, 수치 구현 검증, 재료 적합, 실험 검증은 별도 단계다. 생산 PDE는 model-time 상태이며 연구 MD의 ps를 그 Hz로 재명명하지 않는다. 개구 최초통과, 누적 소산, 순간 위험영역 점유는 서로 다른 관측량이다. 각 정의의 경계조건과 적용 범위를 확인한다.

출처:

- [AGENTS.md](../../AGENTS.md)
- [solver_v1/SOLVER_VALIDATION_GATES.md](../../solver_v1/SOLVER_VALIDATION_GATES.md)
- [solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md](../../solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md)

## 40. 발표 후 다시 볼 기호

a, s, u: 정상 분리, scalar registry, 두 성분 registry 변위
W, G, ℱeff: 내부 에너지, 외부 일 포함 에너지, 조건부 자유에너지
P, Pn, Aabs: intact 확률밀도, 우물 점유 질량, 누적 개구 흡수
M, θ, t₀: 집단좌표 이동도, 열 에너지, 물리 시간 환산 척도
Aatomic, Ac: 원자 에너지 정규화 면적, 통계적 상관 면적

수식과 자세한 유도는 같은 폴더의 theory_notes.md 및 각 슬라이드 출처의 저장소 문서/코드에 있다. 표기상 applied stress σ와LJ길이σLJ, axial projectionχ와complex susceptibilityχv, stacking shiftτ와resolved shear tractionτ는 각각 문맥과첨자로 구별한다.

출처:

- [AGENTS.md](../../AGENTS.md)
- [solver_v1/RESULT_FIELDS.md](../../solver_v1/RESULT_FIELDS.md)
