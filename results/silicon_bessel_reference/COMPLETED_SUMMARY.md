# Si 원자 에너지의 베셀 표현: 수행 결과

2026-10-02. 순수 Si의 원래 SW 에너지를 Si(111) 주기 결정면의 푸리에–베셀 합으로
계산하는 별도 정적 연구 구현을 만들었다. 두 원자 항과 결합각 항을 함께 포함한다.
Al의 Poisson 격자 합과 삼각 역격자 기하를 재사용하고, SW 반경 함수는 J0/J1/J2
Hankel 적분으로 그대로 변환했다. 반경 적분은 수치 구적이며 Al의 닫힌 K식과 다르다.

## 실제로 수행한 계산

Shuffle과 glide 절단면에서 추가 opening과 두 횡방향 좌표를 사용했다. 각 절단면마다
평형, 압축·미끄럼, 순수 미끄럼, opening·미끄럼 두 상태, 완전 분리 상태를 대조했다.
총 12개 정적 상태의 에너지·기울기·전체 3×3 Hessian을 저장했다. 기준은 같은 SW의
명시적 이웃쌍·삼체 직접 합이다. 새 DFT, MD 또는 확률 PDE 계산은 수행하지 않았다.

| 12개 상태의 최대 절댓값 차이 | 결과 |
|---|---:|
| 에너지 | 1.0693668173189508e-12 eV/primitive surface cell |
| 기울기 | 8.179323884860423e-11 eV/Å |
| Hessian | 1.0200196243204118e-8 eV/Å² |

Shell index96, 구적점512를 사용했다. 실행 시간은 152.5105초였다. 위 오차는 직접 합과의
표현 오차이며 실험·DFT와의 Si 물성 오차가 아니다. 기본값이 모든 변형에 대해 엄밀한
오차한계를 보장한다는 뜻도 아니다.

두 diamond basis의 bulk 에너지 합을 기존 독립 DiamondCell 직접 계산과 대조한
검사도 실행 내에서 통과했다. 명시적 bond/angle 공식으로 얻은 평형 shuffle Hessian은
별도 테스트에서 비교했다.

## 결합각 교차항과 수렴

한 중심 원자에 속한 모든 면의 환경 모멘트를 합한 뒤 결합각 에너지를 계산했다.
이웃의 자기항 제거를 값과 두 차수 미분에 모두 적용했다. 실제 교차항 대조에서
올바른 합과 명시적 삼체 계산은 -0.9993201803896681 eV로 같았다. 면별로 따로 제곱하면
-2.1529350138709615 eV가 되어 1.1536148334812935 eV의 교차항이 빠진다.
이 두 이웃 대조는 교차항의 필요성을 보여주는 입력이며 실제 균열 시편은 아니다.

처음 shell32/nodes256에서는 새 시험 17개 중 12개가 실패했다. 판정 허용오차는 바꾸지
않고 역격자 분해능을 높였다. 별도 비대칭 대조점에서 nodes1024를 고정하면 Hessian
오차는 shell32의 8.04818e-5에서 shell96의 4.20955e-9 eV/Å²로 감소했다.
16/32/48/64/96 셸과 256/512/1024 구적점의 7개 대조가 refinement.json에 있다.
이 결과는 시험한 점의 수렴 근거이며 일반적인 급수 꼬리의 증명은 아니다.

## 구현 검사

- 새 베셀 시험: 17 passed in 131.73s. 최종 96/512 설정에서 완료했다.
- 기존 Si 환경·FCC 기하 회귀: 17 passed in 4.14s.
- 저장된 원자료의 차이를 별도 스크립트로 다시 계산하고 소스 SHA256을 확인했다.
  saved_results_verification.json을 따른다. 이 검사는 새 원자 계산으로 세지 않는다.

새 시험에는 signed 거리의 15개 환경 성분, 격자 주기성과 반전, 거리와 두 횡방향의
미분, 두 절단면의 직접 에너지·힘·Hessian 대조, 독립 평형 Hessian, cutoff와 자기항
규약이 포함된다. 시험 수와 정적 원자 상태 수를 더해 실행 수로 보고하지 않는다.

## 아직 해결하지 않은 부분

이 구현은 무한 주기 결정면 사이의 강체 계면을 다루는 정적 기준이다. 초기 균열이
없는 유한 시편에서 첫 균열이 생기는 경로·안정된 개시 상태·자유에너지를 확보한
결과가 아니다. 실제 Si 이동도, 물리 시간, Hz, 개시 확률과 수명은 보정되지 않았다.

원래 SW의 재료 오차는 변환 후에도 같다. MACE 또는 Tersoff 전체를 변환한 결과가
아니며 도펀트·산화층·전하를 포함하지 않는다. 기존 연구를 SW로 대체하지 않았다.
희박한 SW 이웃 직접 계산보다 더 빠르다는 성능 결과도 없다.

Al·UI·생산 에너지·PDE·물리 Hz gate 파일은 변경하지 않았다. 새 구현은 생산 registry에
등록하지 않았다. 검증할 새 재료 함수족을 만들 때 이 표현을 기준으로 쓸 수 있다.

## 파일과 재현

- 유도식: ../../solver_v1/SILICON_BESSEL_REFERENCE.md
- 구현: ../../solver_v1/silicon_bessel_reference.py
- 실행기: ../../solver_v1/run_silicon_bessel_reference.py
- 원자료: states.json, refinement.json, angular_cross_check.json
- 결과: summary.json, source_manifest.json, saved_results_verification.json
- 검사 기록: test_results.json

저장소 루트에서 실행한다. Python/NumPy/SciPy/pytest가 필요하다.

```text
python -m pytest solver_v1/test_silicon_bessel_reference.py -q -p no:cacheprovider
python -m pytest solver_v1/test_silicon_environment_research.py solver_v1/test_fcc111_geometry.py -q -p no:cacheprovider
python -m solver_v1.run_silicon_bessel_reference --output results/silicon_bessel_reference_replay
python results/silicon_bessel_reference/verify_saved_results.py
```

실행기의 출력은 기존 폴더를 덮어쓰지 않는 새 경로여야 한다. source_manifest.json에는
원래 SW 매개변수와 계산에 쓰인 소스의 해시를 기록했다. 개인 실행 경로는 포함하지 않는다.
최초 Git 바이트 검사에서 기존 FCC 기하·격자 합 파일 두 개의 CRLF/LF 차이를 발견했다.
실제 내용을 정규화해 완전히 같은 소스임을 확인하고, 계산 시 raw SHA와 LF SHA 및
Git blob SHA를 함께 남겼다. 기존 참조 소스는 수정하지 않았다. 새 파일과 원자료의
바이트 해시는 별도 package_manifest.json으로 묶었다.
