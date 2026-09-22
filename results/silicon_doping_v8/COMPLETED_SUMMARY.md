# Si 도핑 v8 검증 완료

As1.7의 측정 상한을80C로 수정했다. 이전85C 허용은 원문 Section III의 예외를
놓친 오류였다. 나머지시편85C/하한-40C, 기존탄성표 숫자와실온결과는 유지했다.
관련근거·유도는 solver_v1/SILICON_DOPING_VALIDATION_V8.md 참조.

## 실제 검증

| 검사 | 결과 |
|---|---|
| 집중 도핑·균열 검사 |72 PASS +6 subtests, 1.25s|
| 전체 solver |934 PASS +6 subtests, 3607.74s 기록|
| 전체 app |189 PASS +3 subtests, 3032.91s 기록|
| 데스크톱 smoke |exit0, 기존a0/kappa 유지|
| 독립 Fourier/Airy 탄성77조건 |H 상대차 최대2.89e-15|
| 합성 Fermi product27조건 |Hessian 절대차 최대2.56e-15|
| 고정평균N 곡률 step refinement |8.95e-6 →8.96e-10|
| 수정 continuum 재실행 |21조건/J42검사, 최대상대차1.78e-15|
| 별도폴더 재실행 |JSON/CSV7파일 byte일치|

시간은 각pytest 프로세스가 기록한 값이다. 실행 중 tool UTC clock과 기록시간의
증가폭이 달랐으며 원인은 미확인이다. 연속CPU시간이나 실제 경과라고 단정하지 않는다.

## 해석 범위

- 정규화된 sector 확률도 누락 전하상태가 있으면 평균N/자유에너지가 틀린다.
- 고정평균N의 Legendre 자유에너지와 정확히 고정된 canonical N을 구별했다.
- 동일한 E100/E110/E111을 가진 안정 tensor3개가 서로 다른 균열 H를 낸다.
  영률3개만으로 도핑별 full Cij나 K를 유일하게 맞추지 않는다.
- 이 수학검증은 실제도핑된Si 에너지/장벽/이동도 보정이 아니다.
  새DFT/MD/원자energy평가0회. Physical seconds/Hz는 계속사용불가다.

시작HEAD는4d1e3b984e3833f9800c4b09ba1af4b913ea3fc4, fresh origin과일치했다.
작업브랜치는silicon-wafer-research다. 생산Al/Si registry·PDE·UI 불변.
과거v7수치파일은 변경하지 않았다. 그 commit에서 과거manifest를재현할수있다.
현재소스는v8manifest로검사한다. 최종Git상태는git log와원격으로확인한다.
