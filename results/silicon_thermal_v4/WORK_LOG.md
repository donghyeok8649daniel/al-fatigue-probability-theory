# Si v4 — 진행 기록 (2026-09-21)

## 마감 상태

이하 내용은 초기 checkpoint다. 최종 완료/미완료 상태는 `COMPLETED_SUMMARY.md`와
각 `summary.json`/`interrupted_summary.json`을 따른다. 새 MD48회/1.74 ns 완료,
profile80/136 및 multibasin6/8에서 중단. 지정 기한 뒤 긴 시각 공백을 확인했고
남은 작업 프로세스를 종료했다. 전면8 경계 반사 record를 승인하지 않았다.
관련66개 검사, raw130 HMC/48 MD와 독립 triple39개 확인. 생산 gate는 불변이다.

사용자는 남은 연구를 4시간 동안 진행하라고 지시했다. 시작은 2026-09-20
19:11:55 UTC(한국 9/21 04:11:55), 작업 예산 끝은 23:11:55 UTC다.
브랜치는 silicon-wafer-research, 시작 local/origin HEAD db4bbbd663b2 일치,
fresh fetch 성공, clean. 생산/Al/UI/physical-clock 소스는 바꾸지 않는다.

## 진행 중인 계산과 판정

- Gaussian joint importance pilot 9조건×512개 완료. 300/600 K에서 ESS가
  낮아 단독 보정 채택을 거부했다. 같은 표준정규수를 여러 조건에 사용했으므로
  조건 간 표본을 독립이라고 세지 않는다.
- pCN pilot 36개 체인 완료. 600 K 혼합 불량. Harmonic-reference HMC를 추가했고
  독립 적분/상세평형/시간역전 검증을 통과했다. HMC도 물리 MD가 아니다.
- HMC pilot 36개 완료. 100/300 K에서 원래 국소 가지의 반복 일치는 양호하나,
  600 K에서는 다른 배치/box 경계 방문 때문에 global PMF 판정을 할 수 없다.
- 600 K seed103의 고정-gap quench에서 기존 가지보다 -2.213422554 eV 낮은
  배치를 발견했다. max atom displacement 1.35679 A, 원래 1 A bath box 밖으로
  이완하므로 범위/여러 가지 문제가 실제로 드러났다. 실제 Si의 특정 표면
  재구성을 인증한 것이 아니며 MCMC 전이를 물리 전이율로 세지 않는다.
- 해당 새 가지를 같은 q에서 추적 중이다. 초기 conditional Hessian 최소고유값
  약 0.10999 eV/A²로 양수. 자유에너지/장벽/기억의 후속 해석은 완료 후 기록한다.
- 100/300 K 원래 가지의 mean-force profile 17점×4체인×2온도를 계산 중이다.
  원래 가지의 국소 표본화와 전체 열평형을 구별한다.
- 전면 8반복(4031 bath DOF) 초기점 100/300 K 독립 HMC 체인 계산 중이다.

## 구현된 수학적 구분

q와 무관한 공통 orthonormal bath box에서 적분한다. 유한-range 원자 에너지의
무한 이탈 구성을 포함한 비정규화 적분을 Gaussian 근사로 덮지 않는다.
평균 힘의 분산 감소에는 face마다 0인 vector field를 쓴 integration-by-parts
항을 사용해 유한 영역의 경계항을 보존한다. 원자수로 kBT를 나누지 않는다.
Constrained NVE는 같은 box의 elastic reflection을 구현하고 반사 횟수를 기록한다.
Nonlinear frozen-q force correlation을 exact nonlinear memory라고 부르지 않는다.

Harmonic memory의 Lanczos 양의 spectral quadrature는 실제 Hessian으로만
계수를 정한다. 적합 매개변수나 마찰을 추가하지 않는다. 짧은 시간 오차가
작아도 5/20 ps 별도 창에서 오차가 커지는 것을 그대로 기록한다.

## 아직 최종 보고가 아님

새 테스트는 현재 thermal7, dynamics3, memory2 PASS. 최종 회귀, raw-state 검증,
최종 그림 확인, 전체 hash/인계/Git 정리는 진행 중이다. 이 문서는 checkpoint이며
재료/피로/도핑/생산 Hz 보정 완료 선언이 아니다.
