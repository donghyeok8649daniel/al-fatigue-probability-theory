# v10 재현

이 폴더는 B/P 농도·배치에 따른 **neutral MACE 구속 경로**의 실제 계산이다.
실측 강도/파손확률/Hz 예측이 아니다. 새 DFT/MD/parameter fit은 수행하지 않았다.

## 환경과 모델

정확한 라이브러리 버전은 `runtime.json`을 따른다. Python3.13.5, CPU/float64,
MACE0.3.16, torch2.12.1, e3nn0.4.4, ASE3.29.0을 사용했다.
MACE는 별도 환경에 설치했으며 생산 솔버의 필수 의존성으로 추가하지 않았다.
프로세스당 torch thread2, OPENBLAS thread1을 사용한다.

- 가중치: [공식 MACE-MP-0b3 medium](https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model)
- SHA256: `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`
- 가중치 파일 자체는 커밋하지 않았다. 새 실행마다 코드가 hash를 확인한다.
- 아래 `MODEL_PATH`는 위 가중치의 로컬 경로로 바꾼다. 모든 output은 **새 폴더**다.

## 1. 독립 수학 검증

저장소 루트에서:

```text
python -m pytest solver_v1/test_silicon_concentration_research.py solver_v1/test_silicon_doping_research.py solver_v1/test_silicon_doping_validation.py solver_v1/test_silicon_charge_dynamics.py -q -p no:cacheprovider
```

농도 단위, 미소 음수 주기 좌표, 일반 비직교 셀의 일 미분, 회전 불변성,
화학 배치의 확률 질량 보존 및 기존 전하/SG 축약을 검사한다.

## 2. 실제 원자 농도 비교

```text
python results/silicon_wafer_feasibility/run_concentration_cleavage_v10.py --model MODEL_PATH --output .cache/v10-replay/fixed_cell
```

216자리·동일 기준 부피·B/P1/2/4/8개·세 선택 배치와 순수Si, 총25개 실행이다.
원자1개에서 두 near-plane 배치의 초기 위치는 동일하므로 그 두 실행은 재실행
대조이며 서로 독립인 화학 배치로 세지 않는다. 두 원소를 합쳐23개 고유 초기 배치다.
각 bulk force 수렴 뒤17개 q에서 강체 분리 에너지/힘/응력을 저장한다.
셀을 자유롭게 무응력 부피로 바꾼 계산이 아니며, 모든 결과에 잔류응력을 기록한다.

## 3. 독립 힘·대칭·크기 검사

```text
python results/silicon_wafer_feasibility/validate_concentration_cleavage_v10.py --model MODEL_PATH --source .cache/v10-replay/fixed_cell --output .cache/v10-replay/independent --cases Si_00_pure B_08_interface_cluster P_08_interface_cluster
```

q±h의 직접 에너지 차분, 별도 force-only 경로, 회전·원자순서 변경,
동일 화학농도의 두께 복제 및 면적 extensivity를 검사한다. 면적 복제는 같은
패턴을 반복하는 항등 대조이며, 독립 도펀트 배치나 전자 장거리 상호작용의 수렴이 아니다.

## 4. 같은 국소 하중의 구속 장벽

```text
python results/silicon_wafer_feasibility/run_loaded_concentration_v10.py --model MODEL_PATH --source .cache/v10-replay/fixed_cell --output .cache/v10-replay/loaded --traction-MPa 100
```

기존 양의 q 표에 압축 q를 실제 추가 계산한다. 에너지 곡선의 보간만으로 정지점을
선언하지 않고 실제 E'(q)-AT를 bracketed root로 푼다. 정지점마다 독립 에너지 차분을
검사한다. 같은100MPa는 국소 traction이며, 균열 시편의 원방 응력100MPa가 아니다.
기존 실행은 순수Si/도핑24조건을 두 폴더에 나눴지만 위 명령은25개를 한 폴더에서 재실행한다.

## 5. 선택된 고농도 배치의 최저 곡률

```text
python results/silicon_wafer_feasibility/check_concentration_stationarity_v10.py --model MODEL_PATH --state .cache/v10-replay/fixed_cell/B_08_interface_cluster/opening_0.000.npz --output .cache/v10-replay/lowest_B8 --max-products 240
python results/silicon_wafer_feasibility/validate_concentration_modes_v10.py --model MODEL_PATH --modes .cache/v10-replay/lowest_B8/modes.npz --output .cache/v10-replay/lowest_B8_energy
python results/silicon_wafer_feasibility/check_concentration_stationarity_v10.py --model MODEL_PATH --state .cache/v10-replay/fixed_cell/P_08_interface_cluster/opening_0.000.npz --output .cache/v10-replay/lowest_P8 --max-products 500 --krylov-dimension 64
python results/silicon_wafer_feasibility/validate_concentration_modes_v10.py --model MODEL_PATH --modes .cache/v10-replay/lowest_P8/modes.npz --output .cache/v10-replay/lowest_P8_energy
```

병진3개를 명시적으로 분리한 Gamma Cartesian Hessian의 낮은 고유쌍을 찾고
3개 차분 간격에서 확인한다. 대형 전체 행렬을 만든 검사가 아니며 full-q/셀변형/
finite-T 안정성 인증도 아니다. `solver_converged`와 각 residual을 반드시 확인한다.
제한에 도달하면 exit2와 부분 결과를 보존하며, 완료로 읽지 않는다.
P의 최초 ncv28/240product 검사는 한도 내 미수렴으로 보존했다. P 재검사는
ncv64/500product를 사용하며 tolerance1e-5와 에너지 모델은 같다. 별도의 energy 검사는 계산된
고유벡터 방향으로 실제 에너지를 다시 평가해 곡률을 확인한다.

## 6. 원시 결과 재집계·그림

```text
python results/silicon_wafer_feasibility/summarize_concentration_v10.py --source .cache/v10-replay/fixed_cell --loaded .cache/v10-replay/loaded --output .cache/v10-replay/analysis
python results/silicon_wafer_feasibility/audit_concentration_coverage_v10.py --protocol .cache/v10-replay/fixed_cell/protocol.json --output .cache/v10-replay/coverage
```

저장 원자료의 재집계와 모든 MACE 계산의 재실행은 다르다. 후자는 위 실제 원자
명령을 다시 실행해야 한다. 재현성 검증에서 무엇을 실제 재실행했는지는 완료 요약을 따른다.

IID 화학 배치 bookkeeping의 농도는 local cell의 **평균** 농도다. 실제로 원자수를
고정한 특정 periodic cell과 다른 ensemble이다. 선택 배치의 가중치를 다시 합1로
만들지 않았고, 대칭으로 추가할 수 있는 배치 수도 임의로 부여하지 않았다.

## 7. 파일과 Git 바이트 검사

```text
python results/silicon_wafer_feasibility/package_concentration_v10.py
python results/silicon_wafer_feasibility/package_concentration_v10.py --git-index
```

두 번째 명령은 실제 Git index blob을 읽어 manifest와 비교한다. 원자료를 저장소에
정상 add한 뒤 사용한다. `--write`는 완료 실행을 검토하고 새 manifest를 만드는
패키징용 옵션이므로, 기존 결과를 검증할 때 사용하지 않는다.

## 실패와 시간 기록

`fixed_cell`은 최초의 부분 실행이다. ASE의 기본 wrap epsilon이 허용한 미소 음수
좌표를 입력 검사가 거부해 B1에서 중단했다. 원래 코드·원시 결과를 보존했다.
`fixed_cell_v2`는 명시적 modulo 처리 후 새로 실행한 농도 비교다.
`lowest_B8_energy_control`은 13회 평가와 CSV 저장 후 numpy.bool_의 JSON
직렬화 오류로 중단했다. 실패 코드를 보존하고 Python 기본형으로 변환한 뒤
`lowest_B8_energy_control_v2`에서 원자 평가부터 다시 실행했다.
`lowest_P8_cluster`는 240product(481force평가) 한도로 미수렴한 결과다.
원래 실행 코드를 보존했다. 더 큰 Krylov 공간의 별도 실행은
`lowest_P8_cluster_v2`이며 실제 수렴·결과는 그 summary.json을 확인한다.
로그/경과시간에는 호스트 작업 공백이 포함될 수 있으며 연속 CPU 시간으로 해석하지 않는다.
