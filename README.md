# Al Fatigue Probability Theory

원자적 에너지에서 확률밀도와 균열 개시확률을 계산하는 연구다. 최신 확정 결과는 2026년 11월 18일 완료한 독립 cell별 2차원 FVM 검증이다. 8개 격자와 192개 하중 구간이 수치 기준을 통과했다.

## 바로 열기

| 자료 | 링크 |
|---|---|
| 현재 연구와 브랜치별 자료 | [연구 인덱스](docs/REPOSITORY_INDEX_KO.md) |
| 2차원 FVM 코드와 실행 방법 | [research/fvm-cell2d](research/fvm-cell2d/README.md) |
| 최종 PDF와 Word 및 코드 ZIP | [최종 검증 자료](deliverables/fvm-cell2d-20261118/README.md) |
| 원 지배방정식과 베셀 유도 | [paper/main.tex](paper/main.tex), [유도 PDF](research/source/slip_lattice_energy_mn_K_derivation_KR_v3_23pages.pdf) |
| 기존 프로그램 | [app](app/README.md), [solver_v1](solver_v1) |

## 현재 결과의 범위

2차원 FVM은 각 cell의 수직 분리 a와 상대 미끄럼 s의 공동분포를 전파한다. 대표계 내부 두 이동 원자의 동조 운동과 물리 cell 사이의 독립성을 가정한다. 원래 G4 에너지의 모든 항을 사용하며, 개시영역에 처음 진입한 확률을 흡수상태에 누적한다.

공간 상대차이의 최대값은 0.7509%로 1% 기준을 통과했다. 시간전파·외곽 세분·질량·음수 질량 검사도 통과했다. 실제 특성체적, 재료 이동도와 물리 시간의 보정은 남아 있다.

상관·기억 이론, NAVER 지원자료와 예산, 기존 원고와 시험기는 [연구 인덱스](docs/REPOSITORY_INDEX_KO.md)에서 해당 브랜치로 연결된다. 지원자료는 각 작성일의 계획을 담고 있다.

기존 이론·프로그램 설명은 [이전 README](README_20260930.md)에 보존했다.
