# v35 계산 완료 및 STL 면하중 표시 수정

## 결과

기존 에너지 가족의 **training 오차를 1.66037% 줄였다.** 새 검증 데이터와
eta<=1 기준은 계속 실패하므로 재료 보정 성공으로 채택하지 않는다.
생산 LJ/Bessel/PDE, kinetic JSON, 물리 시간/Hz 설정은 그대로다.
사용자가 추가로 요청한 STL 표시 문제는 별도 UI 수정으로 해결했다.

## 실제 수행한 calibration

7개 정확 anchor를 항상 만족하도록 종속계수를 선형 solve로 제거하고 shape
변화에 대한 implicit Jacobian을 적용했다. 기존 6shape/10coeff/sign/source/scale은 유지했다.
Rank 실패를 regularization하지 않으며 새 에너지 항을 추가하지 않았다.
각 실행의 정의, 모든 callback의 제약 위반, 계수 LP 인증을 함께 저장했다.

| 실행 | 완료 상태 | 최저 적격 training eta | 계산 시간 |
|---|---|---:|---:|
| v34 시작점 | 과거 12iter 한도 종료 | 11.1777577908 | 과거 실행 |
| local | 23iter 정상 종료 | 11.0777285870 | 3376.664s |
| full_saturation | 40iter 한도 종료, best는 endpoint 아님 | 10.9921662122 | 5072.299s |
| saturation_scan | 8개 alpha 독립 LP 완료 | 11.1777577908 | 139.822s |
| refined_continuation | 6iter 정상 종료 | 10.9921657010 | 1008.120s |

마지막 continuation의 추가 차이5.112e-7은 LP certificate tolerance2.198e-6보다
작다. 그 미세 차이를 분해된 재료 개선으로 해석하지 않는다. 사전에 정한
raw training eta 최소 규칙에 따라 후보를 선택한 후 제외 데이터 검증을 했다.
최종 k2=24/p=-1은 원래 경계다. 두 domain/두 active tolerance에서 계산한
feasible-cone의 최소 d_eta는-1.566e-7이며, global 최적성 인증은 아니다.

## 독립 검증

각 실행 전에 선언한 새 계면 상태10개와 새 q6개를 loss/선택에서 제외했다.
기존에 관찰한 excluded q4개는 별도 비교이며 전체 연구 이력의 맹검이 아니다.

| 최대 full-Hessian 상대 오차(Frobenius norm) | 시작점 | 최종 후보 |
|---|---:|---:|
| 새 계면10개 | 118.2114% | 118.0056% |
| 과거 excluded q4개 | 176.7758% | 166.6661% |
| 새 excluded q6개 | 170.1685% | 165.1135% |

최대값이 작아져도 일부 계면은 악화됐다. 예를 들어 state6은 약29%→55%,
state9는 약25%→45%다. Training의 대각 곡률4개는 source와 반대 부호이며
두 크기가 모두 원래 discrepancy scale을 넘는다. 재료 gate는 미통과다.

최종 후보의 analytic tolerance refinement H차이는1.785e-13,
독립 gradient/H 차분 오차는4.991e-10/1.330e-9이다. Bulk radius12→16의
최대 차이는 F-norm2.135e-4, operator norm1.466e-4이며 radius16 tail bound는
operator norm5.147e-4다. 계산한14개 q에서 최소 고유값1.4920이지만 전체
Brillouin-zone 안정성 증명은 아니다. Source Al99 정적 비교이며 새 MD가 아니다.

## STL 수정과 검증

`GeometryWorkflow.import_model`이 가져온 geometry를 저장하면서 mesh를
None으로 만들어 면하중 picker가 빈 화면을 그렸다. STL/삼각형 OBJ를
가져온 즉시 최초 표면 mesh로 사용하도록 연결했다. Generate는 추가 세분화에
사용할 수 있다. 다른 형상으로 바꾸면 기존 삼각형 하중 배정은 초기화된다.

ASCII STL, binary STL, OBJ에서 수정 전 실패3개를 실제 재현했고 수정 후
모두 통과했다. 이동한 직육면체와0.5 단위 스케일을 사용해 실제 렌더링,
940×480 창에서의 표시, 투영 클릭, 연결 면 선택, 면적4mm², 언어/시점/
입력/결과 보존, 명시적 세분화와 이전 하중 해제를 검증했다.
바탕화면 Al Fatigue UI의 launcher가 현재 작업 폴더를 가리키는 것도 확인했다.
기존에 열린 앱은 재실행해야 새 코드가 적용된다.

## 실행한 검사

| 검사 | 실제 결과 |
|---|---|
| 전체 solver | 810 passed, 2510.83s, exit0 |
| 관련 수학 검사 | 14 passed, 1.04s, exit0 |
| STL 수정 후 전체 app | 115 passed +3subtests, 84.31s, exit0 |
| STL 수정 후 desktop smoke | PASS, a0=.7713438268704838 / kappa=86.29296488740997 |
| 최종 Python compileall | PASS |
| 설치된 UI launcher --check | exit0 |
| 최종 validator 재생 | CSV2/JSON2 byte-identical; summary는 elapsed 제외 동일 |
| 비교 그림 | PNG 시각 검토 완료, SVG 함께 저장 |

전체 solver 검사 이후의 수정은 새 reporting/validation 보조 코드와 UI다.
관련 로그는 로컬 `.cache/anchored_calibration_v35_*.log`에 있다.
보고기의 v34 schema 차이(`final_shape` 없음)는 저장된 shape 좌표로 해결했다.
공용 source loader의 옛 straight-row 단위가 이번 계면/벌크 단위로 오인되지
않게 source provenance와 실제 observable units를 분리하고 재생 확인했다.

## 파일

- [수학·실행 절차](../../solver_v1/EXACT_ANCHOR_CALIBRATION_V35.md)
- [Training 선택과 모든 실행 비교](comparison/training_selection.json)
- [개별 training residual](comparison/training_residuals.csv)
- [독립 검증 요약](excluded_validation/summary.json)
- [계면 전체 값](excluded_validation/interface.json) / [finite-q 전체 값](excluded_validation/finite_q.json)
- [비교 그림](comparison/calibration_comparison.png)

시작 Git local/origin은470aa28e09ec44478e499c78d187f3e40e5d7397,
브랜치는 probability-pde-solver-v1이었다. 최종 commit과 원격 일치는
`git log` 및 fresh remote ref로 확인한다. 실패한 재료 gate와 미보정
kinetic/물리Hz/실제 항복/피로 수명 gate는 계속 분리한다.
