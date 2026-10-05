# 재현 방법

Si 연구 브랜치 루트에서 NumPy, SciPy, pytest, Matplotlib 환경을 사용한다.
실제 원자 샘플링에는 별도로 ASE, PyTorch, MACE와 README의 SHA가 일치하는
기존 모델 파일이 필요하다. 실행 환경 버전은 EXECUTION.json을 따른다.
원자료와 분석은 새 출력 폴더로 실행하며 기존 결과를 덮어쓰지 않는다.

## 실제 계산을 하지 않는 검산

```powershell
python -m pytest solver_v1/test_silicon_crack_records.py solver_v1/test_silicon_thermal_identity.py solver_v1/test_silicon_initiation_probability.py solver_v1/test_silicon_thermal_validity.py solver_v1/test_silicon_device_ensemble.py solver_v1/test_silicon_thermal_research.py -q -p no:cacheprovider
python -m results.silicon_wafer_feasibility.audit_crack_thermal_records --root . --pilot results/silicon_crack_thermal_validation/nonlinear_pilot --output NEW_ANALYSIS
python -m results.silicon_wafer_feasibility.audit_crack_geometry --root . --output NEW_GEOMETRY
python -m results.silicon_wafer_feasibility.verify_crack_thermal_package --root .
```

첫 분석은 130개 기존 평가의 난수, 좌표, 수용·거절, 반복 상태, 총에너지의
열단위 변환과 경계 항등식을 재계산한다. 재검산 자체의 새 원자 호출은 0이다.
기하 분석은 기존 두 구조의 원자 좌표만 사용하며 새 원자 호출은 0이다.
완료되지 않은 pilot은 완성본으로 분석하지 않는다.

## 실제 새 에너지 평가를 다시 수행할 때

```powershell
python -m results.silicon_wafer_feasibility.sample_nonlinear_ensemble_v15 --model MODEL_FILE --output NEW_PILOT --max-seconds 1200
```

OMP_NUM_THREADS, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS는 1로 설정했다.
코드가 PyTorch 연산 스레드를 2, interop 스레드를 1로 지정한다.
기존 6 GB 시작 메모리 보호선을 유지한다. 130회 제한된 pilot의 완료는
실제 열평형·혼합·자유에너지 또는 균열 확률의 통과를 의미하지 않는다.
MC 인덱스와 계산 실행 시간은 생산 모델의 물리 초/Hz가 아니다.

## 결과와 소스 바이트

package_manifest.json은 자기 자신을 제외한 결과 파일과 실행·검산·시험 소스의
SHA256을 기록한다. 기존 Hessian과 모델은 복제하지 않고 기존 경로와 SHA로
연결한다. JSON 원자료의 줄 끝을 포함한 바이트는 로컬 .gitattributes로 보존한다.
Git commit 전 실제 index blob과 대조한다. 사용자·임시 절대 경로나 모델
가중치는 결과 패키지에 넣지 않는다.
