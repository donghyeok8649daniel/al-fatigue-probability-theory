# 순수 Si SW 베셀 에너지와 물성 보정의 계약

현재는 별도 정적 연구 후보다. 생산 energy registry/PDE/mobility/clock에
연결하지 않았다. 물성·최초 균열 개시·유한 온도 자유에너지·초/Hz·수명은
미검증이다. 결과와 실제 재현은 `results/silicon_sw_material_fit/README.md` 참조.

## 원래 SW와 추가한 최소 각도 형상

원자별 에너지와 전체 에너지는

```
E = sum_i E_i
E_i = (1/2) sum_j phi(r_ij)
    + lambda*epsilon sum_(j<k) g(r_ij)g(r_ik) f(mu_jik-c0)
phi(r) = A*epsilon*[B*(sigma/r)^p - (sigma/r)^q]*exp(sigma/(r-rc))
g(r) = exp(gamma*sigma/(r-rc)), rc=sigma*a
f(theta) = theta^2*(1-beta*theta)^2
```

거리항과 kernel은 원래 SW의 유한 지지를 유지한다. r>=rc에서0.
beta=0은 원래 SW다. epsilon을 고정하고 repulsion/attraction/angle의
세 양의 진폭을 적합한다. epsilon과 A/lambda의 중복 scale을 세지 않는다.

`f'=2 theta (1-beta theta)(1-2beta theta)`다.
theta의 물리 구간은[-1-c0,1-c0]. max(beta*theta)<=1/2이면 선호각 양쪽의
단조성과 유일 각도 영점을 유지한다. 원래 c0=-1/3에서 검사한 대칭 beta
구간은[-.375,.375]. 이 선언된 구간을 모든 새 c0에 무심코 전이하지 않는다.

## 정확한 방향 moment 표현

P(mu)=f(mu-c0)=sum_(m=0..4) c_m mu^m, T_m=sum_j g_j n_j^(tensor m),
D=sum_j g_j²라 두면 각도합은

```
sum_(j<k) g_j*g_k*P(n_j dot n_k)
  = (1/2)[sum_m c_m ||T_m||_F^2 - D*P(1)]
```

모든 이웃 면·basis의 T_m을 합한 다음 제곱한다. 이웃별 자기항 D*P(1)을
뺀다. 대칭 monomial 35개에 m!/(ax!ay!az!)의 multiplicity를 곱하여 full
tensor norm을 회복한다. 이는 beta 한 개의 고정 표현이며35개 자유계수가 아니다.
Diamond의 T3는0이 아니다. bulk 배경을 site별로 보존·차감한다.

평면 방향의 유한 Fourier mode(-4..4)와 radial factor를 분리한다.
각 mode의 2D Hankel 적분은 J_|m|(G*rho)를 사용한다. Poisson 역격자 합은
완전한 G/-G disk와 exp(iG·delta)를 갖는다. value,normal/lateral gradient,
전체 Hessian을 같은 식에서 미분한다. 허수 잔차를 검사한다.
같은 평면의 nonself 이웃은 직접 합하며 d=0 특이 변환을 피한다.

## 벌크 곡률 보존의 정확한 범위

tetrahedral4-neighbor 상태에서 mu=c0=-1/3이면 f,f',f''가 beta0과 같다.
9차원 cell state(exx,eyy,ezz,gamma_yz,gamma_xz,gamma_xy,ux,uy,uz)의
full Hessian이 보존된다. 압축으로 second neighbor가 들어오거나 c0가
바뀌면 이 명제가 적용되지 않는다. EOS, 힘, 응력과 H를 새로 계산한다.

tetrahedral4NN와 원래 영 각도 배경, 평형 조건에서 relaxed cubic constants는

```
C44 = 3*(C11+2*C12)*(C11-C12)/(7*C11+2*C12)
```

두 독립 harmonic stiffness의 결과다. 현재 공개 normal anchor를 넣으면
C44=65.267599GPa로 목표72.17693에 약9.6%부족하다. quartic 각도 형상은
이4NN harmonic 제약을 바꾸지 않는다. 독립 finite-q 검사에서도 차이가 남았다.
이 조건부 식을 모든 SW 선호각/추가 이웃/베셀 계열에 대한 불가능 증명으로
확대하지 않는다. 선호각 변경과 추가 이웃에서는 angular bulk 배경과 H를
다시 계산한다. 선호각 단독9개/추가이웃16개/초기공동19개/각도경계9개는
벌크E guard에 실패했다. 기존 학습 형상으로 거리 sigma=1.8A를 세밀히 검사한
3개 중2개는 네 anchor와 bulkE guard를 함께 통과했다.

현재 선택은 local_boundary_fit_results profile1이다. sigma=1.75A,
gamma=1.5, rc=3.77118A, c0=-.37357245686009144,
beta=.3640142880725576=.5/(1-c0). c0 변경의 비영 벌크 각도항을 유지한다.
공개 격자5.4610215037A와 C11/C12/C44=153.28991/56.25009/72.17693GPa를
동시에 맞춘다. bulk E RMSE .02533607eV/atom, Si111 F/E는
1.02103563eV/A/.17002476eV/atom. 양의각도경계와주변30개중12개가조건을
통과했고 같은학습/EOS손실로선택했다. bulk F 및 외부GB/aSi/SF와 진동수 오차가
남아 재료 채택이 아니다. 원래영각도 후보의 조건부C44 식을 새c0에 적용하지 않는다.

## 정적 하중과 단위

계면 state q=(opening,u1,u2), W는 per-interface-cell 초과 에너지(eV).
physical static work는 W-A_atomic(Tn*opening+tau1*u1+tau2*u2)다.
SI/conversion을 먼저 맞춘다.

```
1 eV/A^3 = 160.2176634 GPa
1 eV/A^2 = 16.02176634 J/m^2
```

Ac는 에너지·응력·이동도 정규화에 들어가지 않는다. 강체 평면의 full3x3H
국소 안정성 상실은 정해진 경로의 이상적 기계 진단이며 실제 최초 균열 개시
또는 피로수명과 다르다. 표면 이완/공간 균열핵/온도/시편 mechanics가 빠져 있다.

## 데이터와 선택

PW91 dia489/surface00110/surface11012만 학습511. 나머지1964는 손실 제외.
교환-상관 함수를 구별하고 같은PW91 frame548의 per-atom E를 고정 기준으로 쓴다.
자료군의 총 가중치, force/E/stress/EOS scale은 명시된 선택 가중치다.
불확실성/CI로 부르지 않는다. 외부28구조는 locked selection 후 검사했다.

원래SW의 벌크 상대E RMSE .0288054 eV/atom을 guard로 둔다.
fixed-shape3amp quadratic fit에서 active sets/scalar dual bracket을 검사한다.
beta가 변해도 원자 관측량은 beta의2차, joint loss와 guard는4차다.
정지점·실수 경계근·endpoint를 나열해 비교한다. floating-point root enumeration은
interval/global shape 최적성 증명과 다르다. 형상 dictionary 전체를 자유적합으로
채택하거나 heldout 실패를 loss에 돌려 넣지 않았다.

## 공통 확률 이론으로 넘어가는 gate

`partial_t P = div[M(P grad F + kBT grad P)]`의 재료 입력은 집단좌표에
대응하는 자유에너지여야 한다. 현재W는0K rigid static energy다. F(q,T),
좌표에 대응한 mobility·memory/low-frequency 응답과 물리clock은 미확보.
harmonicTHz는 Newtonian vibration이고 overdamped mobility 보정이 아니다.
도펀트·전하·농도·산화층을 이번순수Si적합이 설명한다고 보고하지 않는다.
