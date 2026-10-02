# Si(111) 원자 에너지의 푸리에–베셀 표현

2026-10-02. 순수 Si의 원래 Stillinger–Weber(SW) 에너지를 유지한 정적 연구 구현이다.
MACE를 SW로 치환해 재료를 승인하거나, 첫 균열 개시 계산을 완성한 결과가 아니다.
Al의 LJ/Bessel 에너지, 생산 확률식, UI, 이동도와 물리 Hz 설정은 바꾸지 않는다.

## 무엇을 베셀로 옮겼는가

Al에서 쓰던 Poisson 격자 합의 틀과 삼각 격자의 역격자 셸을 재사용했다.
원래 SW의 두 원자 함수와 결합각 함수를 그대로 두고, Si(111) 결정면의 합을
J0/J1/J2 Hankel 적분과 역격자 위상 합으로 계산한다. 반경 함수를 LJ나 지수 함수에
맞춰 바꾸지 않았으며 새 Si 매개변수를 적합하지 않았다.

Al의 역멱 함수는 변형 베셀 K의 닫힌식으로 변환되지만, SW의 유한 범위 지수
함수는 같은 닫힌식이 아니다. 이번 구현의 반경 적분은 수치 구적이다.
이를 K의 유한 합으로 정확히 바꿨다고 하거나 더 빠른 계산이라고 주장하지 않는다.

## 같은 결합각 에너지

중심 원자의 이웃마다 r, 단위방향 n, 원래 SW 가중치 w를 사용한다.

```text
rc = a sigma
w(r) = exp[gamma sigma/(r-rc)]  (r<rc), 0 otherwise
rho = sum_j w_j
v = sum_j w_j n_j
Q = sum_j w_j (n_j tensor n_j - I/3)
D = sum_j w_j^2
c = cos(theta0)

u_site = (1/2) sum_j phi2(r_j)
       + lambda epsilon [Q:Q/2 - c v.v + (1/6+c^2/2)rho^2
                         - (1-c)^2 D/2]
```

이웃의 모든 결정면을 더한 후에 제곱한다. 면마다 결합각 에너지를 계산한 후 더하면
서로 다른 면에 속한 두 이웃 사이 각도 항을 잃는다. D는 같은 이웃을 두 번 고르는
자기항을 제거한다. 이 항도 값·기울기·Hessian에 함께 들어간다.
`angular_cross_check.json`은 면별로 따로 제곱한 계산이 틀리는 실제 대조를 남긴다.

## 결정면 변환

Fourier convention은 f_hat(G)=integral f(x) exp(-i G.x) dx이며,

```text
sum_R f(R+delta,d) = (1/A_atomic_cell) sum_G f_hat(G,d) exp(i G.delta)
```

를 쓴다. rho_xy를 평면 반경, r=sqrt(rho_xy^2+d^2), t=rho_xy/r,
u=d/r, e=G/|G|로 두고 I_m[f]=2pi integral rho_xy f J_m(|G|rho_xy) d rho_xy라 쓰면

```text
pair_hat = I0[phi2]
rho_hat  = I0[w]
D_hat    = I0[w^2]
v_xy_hat = -i e I1[w t]
v_z_hat  = I0[w u]

Q_ab_hat = delta_ab {(I0[w t^2]+I2[w t^2])/2 - I0[w]/3}
           - e_a e_b I2[w t^2]             (a,b in plane)
Q_az_hat = Q_za_hat = -i e_a I1[w t u]
Q_zz_hat = I0[w (u^2-1/3)]
```

G=0의 횡방향 vector/cross 성분은 0이다. signed d를 유지한다.
적분 상한은 sqrt(rc^2-d^2)다. 원래 SW cutoff에서 함수와 필요한 미분이 0으로
사라지므로 d 미분에 움직이는 상한의 항이 추가되지 않는다. d 기울기·둘째미분은
적분 전 원래 함수의 forward jet으로 계산하고, 평면 미분은 iG와 -G tensor G로 계산한다.
값만 맞춘 뒤 힘을 별도로 근사한 모델이 아니다.

d=0의 두 원자 연속 변환에는 자기 원자의 비적분 특이점이 있어 이 경로를 차단한다.
같은 면의 비자기 이웃은 원래 SW cutoff 내에서 정확한 직접 합으로 포함한다.
변형으로 같은 면의 이웃이 cutoff 안에 들어오는 경우도 누락하지 않는다.
임의 real-space cutoff를 새로 넣은 것이 아니라 원래 SW 정의를 사용한다.

## 다이아몬드 적층과 계면

기존 FCC(111) 기하에서 h=a_lat/sqrt(3), tau=(a1+a2)/3을 사용한다.
Si의 두 basis는 z=(l+3b/4)h, lateral=l tau, b=0,1이다.
따라서 순차 면 간격은 3h/4와 h/4이며 shuffle/glide 절단면을 따로 둔다.
기존 +ABC 기하의 cubic 축 대응 규약을 보존한다.

좌표는 q=(추가 opening, lateral x, lateral y)이며 두 횡방향을 모두 남긴다.
상부 결정만 이동시킨 후 각 원자의 환경을 다시 합한다.
각 영향받은 면에는 primitive surface cell당 원자 하나가 있다.

```text
W_int(q) = sum_affected_sites [u_site(q) - u_site(bulk)]
```

양쪽 반결정의 모든 영향받은 중심을 포함한다. 이미 더한 두 쪽에 다시 factor 2를
곱하지 않는다. bulk moment가 모두 0이라고 가정하지 않고, 면별 에너지를 먼저
만든 후 site별 bulk 에너지를 뺀다. infinite total bulk를 크게 더했다 빼지 않는다.
W는 eV/primitive surface cell, gradient는 eV/Angstrom, Hessian은 eV/Angstrom^2다.
면적은 원자 에너지 환산용이며 통계적 상관 면적이나 메시 면적이 아니다.

## 실제 검증과 분해능

`test_silicon_bessel_reference.py`의 17개 시험이 통과했다. 검사에는 signed 거리의
15개 환경 성분, 격자 주기성·반전, d와 두 lateral 미분, 두 절단면의 에너지·힘·전체
3x3 Hessian, 기존 독립 bond/angle 공식의 평형 Hessian, 자기항과 cutoff 규약이 있다.
직접 대조는 원래 SW 명시적 이웃쌍/삼체 항을 사용한다.

처음 shell index32에서는 12개 실패/5개 통과였다. 허용오차를 완화하지 않고
역격자 분해능을 높였다. 비대칭 점의 곡률 오차는 shell32에서 약8.05e-5,
shell96에서 약4.21e-9 eV/Angstrom^2로 줄었다(nodes1024 고정).
별도 nodes256/512/1024 대조도 남겼다. 구현 기본값은 shell96/nodes512다.
셸 비교와 구적 대조는 시험한 상태의 수렴 근거이며 모든 상태의 엄밀한 tail bound가 아니다.
특히 매우 작은 면간 거리, 심한 압축이나 새로운 kernel에는 추가 수렴 검사가 필요하다.

12개 실제 계면 상태와 bulk의 독립 DiamondCell 대조, 환경 제곱의 교차항 대조는
`run_silicon_bessel_reference.py`로 별도 저장한다. 시험 수와 원자 상태 수를 합산하지 않는다.
최종 수치는 results/silicon_bessel_reference/summary.json과 완료 요약을 따른다.
새 DFT, MD, 자유에너지 표본 추출, 확률 PDE 계산은 수행하지 않았다.

## 재현

```text
python -m pytest solver_v1/test_silicon_bessel_reference.py -q -p no:cacheprovider
python -m solver_v1.run_silicon_bessel_reference --output results/silicon_bessel_reference_replay
```

출력 폴더는 새 경로여야 한다. 기존 실패·결과를 덮어쓰지 않는다.
Python/NumPy/SciPy/pytest가 필요하며 특정 개인 실행 파일의 절대 경로는 소스에 넣지 않는다.
원래 SW 매개변수는 기존 감사된 source_Si.sw와 SHA를 사용한다.

## 남은 연구

이 결과는 초기 내부 균열이 없는 360원자 시편의 첫 개시 계산이 아니다.
무한 주기 결정면의 강체 계면을 열고 미끄러뜨리는 정적 기준이다.
기존 SW의 실제 Si 힘·에너지 오차는 그대로 남는다. 변환의 일치가 재료 보정의 성공은 아니다.

MACE 전체를 이 SW 식으로 정확히 변환하지 않았다. Tersoff와도 별도다.
도펀트와 산화층은 현재 포함하지 않는다. 도핑에 의한 주기성 깨짐은 host의 베셀 합과
실제 국소 환경의 보정, 또는 명시된 supercell을 검토해야 하며 별도 교차원소 에너지가 필요하다.
원래 SW host만으로 Si/B/P/O 재료·전하·균열 장벽을 정하지 않는다.
이후에는 독립 재료 자료로 반경/각도 함수족을 검증하고, 실제 시편의 재배열·개시 경로,
유한온도 자유에너지와 기계적 이동도로 연결해야 한다.

## 1차 자료

- [LAMMPS의 원래 SW 식과 cutoff 정의](https://docs.lammps.org/pair_sw.html)
- [NIST DLMF의 정수 차수 Bessel 적분](https://dlmf.nist.gov/10.9#E2)
- 기존 저장 source_Si.sw 및 원자료 SHA, 직접 삼체 모멘트 감사, FCC 기하와 역격자 구현은
  source_manifest.json에 결합한다. 새로운 문헌 숫자를 만들어 적합하지 않았다.
