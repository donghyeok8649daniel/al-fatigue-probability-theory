# 공간 FVM 연결: 첫 연구 검증

이 모듈은 공간 FVM 역학 솔버나 새로운 생산 PDE가 아니다. 기존 연구용
resolved_tensor_tractions를 재사용하고, 공간 표면 quadrature의 확률 집계만
정의한다. LJ/Poisson/Bessel, 에너지·이동도·온도·생산 응력 mapping은 불변.

## 두 공간과 외부 일

공간 좌표 x에서 국소 내부 상태는 P(a,s,t|x)다. x의 표면 패치와 (a,s)의
확률 제어체적을 구별한다. 원자 상태의 SG 플럭스를 공간 확률 확산으로
옮겨 해석하지 않는다. 공간 전파항이나 강성 저하는 이 단계에 없다.

동일 좌표계의 대칭 Cauchy 응력 sigma와 단위 계면 법선 n, 접선 m에 대해

    t_n = n^T sigma n
    tau_1 = m^T sigma n
    tau_2 = (n cross m)^T sigma n.

입력 MPa는 출력 MPa다. 두 번째 전단 성분도 보존한다. 기존 scalar-s
PDE가 tau_2를 표현한다고 주장하지 않는다. 명시적 계면의 외부 일은
단위·원자셀 정규화를 맞춘 traction-displacement pairing이 필요하며,
기존 f*=kappa sigma/E를 이 연구 함수가 변경하지 않는다.

## 표면 집계

표면 개시 및 독립 등가 국소 영역 가정에 한해서

    log S_spec(t) = sum_j (Delta A_j/A_c) log1p(-P_j(t))
    P_spec(t) = -expm1(log S_spec(t)).

Delta A_j는 공간 계산 패치 면적, A_c는 외부 통계적 상관 면적이다.
원자 에너지 정규화 A_atomic_cell과 둘 다 다르다. 여기에는 임의의
물리적 A_c 기본값이 없다. 테스트 면적은 순수 합성 수치 검증이다.

패치를 동일한 P의 여러 하위 패치로 쪼개도 면적 합이 같으면 log S는
대수적으로 같다. 비균일한 P를 공간적으로 더 정확하게 분해할 때의 변화는
별도의 공간 quadrature 수렴 문제다. 내부 개시/체적 집계는 범위 밖이다.

## 인증의 한계

입력은 누적 실제 개구 흡수량이며 음수·1 초과·감소·NaN을 거부한다.
수치 floor는 검증 결과로 공급해야 하며 여기서 임의 생성하지 않는다.
양의 면적 패치 전부가 floor 초과 및 convergence_certified일 때만
numerically_resolved_extrapolation을 반환한다. 그 밖에는 NaN이며 별도
mathematical_extrapolation은 남긴다. 인증 bool은 외부 검증의 입력일 뿐
이 함수가 수렴 계산을 수행하거나 인증서를 검증하지 않는다.
확률 0인 패치도 보수적으로 미분해 처리하며 '정확한 부재' 증명은 별도다.
수치 인증은 재료·독립성·공간응력·물리시간 보정 인증이 아니다.

## 실제 검증 및 다음 gate

10 tests: 균일 패치 2/7/100 분할, 비균일 패치 분할, 영면적 제외,
noise amplification 차단, 누락 인증, 잘못된 확률 거부, 인장/순수 전단,
회전 공변성. 공간 응력 PDE나 원통 메쉬를 실제 풀었다는 뜻이 아니다.

다음 단계는 공간 역학의 균일 인장/전단 patch test와 응력장 수렴,
에너지 일관성 및 탄성 이중계산 방지, 국소 PDE 하중 이력 연결이다.
그 후 내부 격자/시간/공간의 독립 refinement를 수행한다. 독립 국소
진화에서 확률을 그대로 강성 감소율로 쓰거나 경험 손상식을 추가하지 않는다.
생산 물리 Ma/Ms/t0 및 Hz는 여전히 미보정이다.
