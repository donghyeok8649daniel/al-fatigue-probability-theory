# 연구 자료 찾아보기

2026년 11월 20일 원격 브랜치와 로컬 확정본을 대조해 정리했다. 현재 확정된 FVM 수치 결과의 기준일은 11월 18일이다.

| 자료 | 위치 |
|---|---|
| 2차원 FVM 실행 코드와 결과 | [research/fvm-cell2d](../research/fvm-cell2d/README.md) |
| 최종 PDF와 편집용 Word | [deliverables/fvm-cell2d-20261118](../deliverables/fvm-cell2d-20261118/README.md) |
| 기존 지배방정식과 베셀 유도 | [paper/main.tex](../paper/main.tex), [원자열 유도 PDF](../research/source/slip_lattice_energy_mn_K_derivation_KR_v3_23pages.pdf) |
| 기존 프로그램 | [app](../app/README.md), [solver_v1](../solver_v1) |
| 물리 시간과 이동도 조건 | [PHYSICAL_TIME_AND_MOBILITY.md](../solver_v1/PHYSICAL_TIME_AND_MOBILITY.md) |
| 재료·솔버 검증 조건 | [SOLVER_VALIDATION_GATES.md](../solver_v1/SOLVER_VALIDATION_GATES.md) |
| 상관·기억 이론과 비교검증 | [원격 자료](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/codex/crack-correlation-memory-20261104/outputs/sg-coarse-memory-20261104/README_KO.txt) |
| NAVER 지원자료와 예산 | [원격 자료](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/codex/crack-correlation-memory-20261104/outputs/d2sf-campus-20261105/README_KO.txt) |
| 기존 지원자료의 그래프 해설 | [원격 자료](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/codex/crack-correlation-memory-20261104/outputs/d2sf-campus-20261105/GRAPH_EXPLANATION_KO.txt) |

## 브랜치 위치

| 브랜치 | 자료 |
|---|---|
| [probability-pde-solver-v1](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/probability-pde-solver-v1) | 이번 정리와 2차원 FVM 확정 결과 |
| [codex/crack-correlation-memory-20261104](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/codex/crack-correlation-memory-20261104) | 11월 5일까지 게시한 상관·기억 이론과 지원자료 |
| [research/20261018-integration](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/research/20261018-integration) | 10월 18일 통합 원고 기록 |
| [research/20261018-bessel-mechanics](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/research/20261018-bessel-mechanics) | 베셀 역학 연구 |
| [research/20261018-surface-bulk](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/research/20261018-surface-bulk) | 표면·벌크 연구 |
| [research/20261018-fvm-coupling](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/research/20261018-fvm-coupling) | FVM 결합 연구 |
| [fatigue-tester](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/fatigue-tester) | 시험기 파일 |
| [silicon-wafer-research](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/silicon-wafer-research) | 실리콘 연구 |

기존 D2SF 문서는 해당 작성일의 기술·예산안이다. 현재 팀 구성이나 최종 제출 상태는 문서의 작성시점과 함께 확인한다. 2차원 FVM의 통과 판정은 대표계 수치 검증 범위이며, 실제 시편 특성체적·재료 이동도·물리 시계 보정에는 추가 근거가 필요하다.
