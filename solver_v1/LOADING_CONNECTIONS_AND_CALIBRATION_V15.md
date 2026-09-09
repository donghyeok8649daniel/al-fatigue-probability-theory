# 전단·방향별 응력 연결 상태와 보정 범위 — v15

이 문서는 이름이나 화면 레이블이 아닌 실제 호출 경로를 확인한 기록이다.
새 연구 결과를 기존 production 결과로 소급해서 해석하지 않는다.

## 1. 지금 어디까지 연결되어 있는가

| 계층 | 실제 지원 | 실제 경로 / 한계 |
|---|---|---|
| 응력 기하학 | 대칭3×3 Cauchy응력에서 normal·전단2성분 | `interface_static_scenarios.resolved_tensor_tractions`; 투영 연산이며 PDE가 아님 |
| 벡터 정적 연구 | W_int(a,s1,s2), 독립 normal/전단2성분의 정지점·Hessian·장벽 | `vector_registry_audit.stationary_state`, `activation_stress_sensitivity.barrier_and_derivative`; 물리 MPa의 static 연구 |
| 주기 확률 연구 | P(a,s), 독립 normal/shear 진폭·위상·모델 주기 | `run_low_stress_cyclic_diagnostic.run_cycles` → `low_stress_cyclic_diagnostic.run_reflecting_cycles`; reflecting-box SG 가설검사 |
| production확률 PDE | P(a,s), scalar축응력 σ/E→κ | `app.solver_adapter.run_ui_analysis` → `run_probability_pde_2d`; 독립 전단 입력 없음 |
| desktop UI | scalar mean/amplitude, 에너지모델 선택, 모델시간 | `app.desktop_ui`; 방향 입력은 disabled, 독립 전단·공간 구조물 solver 없음 |

연구 SG의 `build_surface`는 `matched_v3/monotone_opening/fitted_candidates.json`
속 과거 후보를 읽는다. 이번 v15 벌크/계면 보정 후보가 그 PDE에 연결된 것이
아니다. 해당 reflecting 실험의 opening_probability는 null이며, 이를 개구
흡수 확률로 읽지 않는다. 원자적 per-cell energy를 국소 activation energy로
보는 가정도 별도 검증이 필요하다. **전단 투영, 벡터 정적 계산, 확률 상태
차원, 공간 구조물 차원은 서로 다른 기능이다.**

`production_capabilities()`의 실제 값:

    state_coordinates=[a,s], probability_state_dimension=2
    independent_shear_input=False
    orientation_input_active=False
    vector_registry_pde=False
    spatial_specimen_solver=False

## 2. Production의 에너지·응력·시간 경로

    app.desktop_ui
      -> app.solver_adapter.UIAnalysisConfig
      -> build_time_basis_model / build_energy_model
      -> selected model energy, gradients, Hessian
      -> scalar generalized force f*=kappa_axial sigma/E
      -> run_probability_pde_2d

기본값은 `two_row_lj_reference`다. hypothetical hybrid와 reduced Al-target
best-feasible hybrid는 별도 ID이며, full-FCC active-interface 후보는 이
선택기에 등록되지 않았다. `canonical_model_params()`는 M_a=1,M_s=.05,
kT=.02,chi=.2를 유지한다. MPa/GPa 변환은 E_GPa×1000으로 수행한다.
이 경로는 수직·전단 tensortraction work를 새로 주입한 경로가 아니다.

`physical_load_conversion()`과 결과 metadata가 active model/parameter
source/calibration status 및 위 capability값을 기록한다. 방향 입력을
화면에 보인다는 이유만으로 orientation이 적용됐다고 말하면 안 된다.

    t_phys=t0 t_model, f_Hz=f_model/t0,
    M*=t0 E0/L0² M_phys.

변환 도구는 있으나 검증된 production M_a,phys/M_s,phys/t0는 없다.
`UIAnalysisConfig.validate_time_basis()`가 calibrated object와 동일 energy
fingerprint를 검사한다. 현재 JSON에서는 physical모드를 사용할 수 없으며
cycles/model-time을 Hz라고 바꾸지 않는다.

## 3. 이번에 실제로 맞추고 검사한 것

- 기존114/62/32GPa의 독립 cubic탄성,3.36eV/atom과 FCC4.05Å 기준을 유지한
  정적 보정 및 계면 오차/식별성 재검사. LJ/무한 Poisson–Bessel를 보존한다.
- QP의 잘못된 near-active 제약면 처리와 Gram정밀도 손실을 재현·수정.
  같은 primal/KKT 기준을 유지하고 독립 파수·작은 strain에서도 검사한다.
- 실제 Al99 MD1024frame:576atoms/plane,12periodic class의 정확한 에너지/
  covariance 정규화. 일부 합산값이 맞더라도 개별mode와 block오차를 공개한다.
- 실제 공개 Al미세선692개 rawrelaxation trace,643개 selected apparent
  activation row, 문헌 line-drag/core-resistance 단위와 조건을 분리하여 감사.
- MPa normal/shear 정적 장벽 응력미분과 anisotropic 유한 pinned-line 장벽.
  후자는 pin span/core/outer geometry가 가설이며 실험항복/속도/피로 검증이 아니다.

상세 결과는 `INTERFACE_CALIBRATION_DEVELOPMENT_V15.md`,
`PUBLIC_ALUMINUM_CALIBRATION_EVIDENCE.md`,
`ZERO_FREQUENCY_MOBILITY_AUDIT.md`,
`FINITE_SOURCE_BARRIER_DERIVATION.md` 및 각 results 디렉터리에 저장한다.

## 4. 아직 채택하면 안 되는 연결

- 불안정/계면 오차가 큰 static후보 → productionPDE/UI.
- 공개 MD의 물리 ps, 문헌 dislocation B_line → reduced a/s의 t0.
- 균일 interface per-cell 장벽 → 유한 dislocation source의 전체 activation.
- apparentactivation area, A_atomic_cell, wire단면적 → specimen통계 A_c.
- 방향투영도구 → 3D확률 PDE 또는 mesh공간예측.

이 구분을 해소하는 것은 상수 이름이나 축 단위를 바꾸는 일이 아니다.
각 상태·좌표·에너지·kinetic normalization을 도출하고, 독립 검증으로
통과시킨 뒤 연결해야 한다. 이번 작업은 그 검사를 실제 자료로 진행하며,
기존 production물리·UI 기본값은 변경하지 않는다.

## 5. 설명에 그치지 않고 실행한 방향별 정적 검사

`run_tensor_traction_audit.py`는 보존된 v14 후보, strain-stable cross 후보,
검증용 Al99 source를 각각 실제로 평가했다. x/y/z 수직50MPa, xy/xz/yz
전단10MPa, 반대 부호 이축25MPa, 명시된 혼합 tensor, 계면 순수 수직50MPa와
독립 전단 두 방향4MPa의11개 시나리오를 사용했다. 각 시나리오는
0→.5→1→.5→0→-.5→-1→-.5→0의 정적 하중 경로다.

총297개 상태를70.86초에 계산했다. 각 상태의3좌표 힘 평형과 전체 Hessian의
양의 고유값을 확인했고, 최대 무하중 복귀 좌표 오차는3.72e-15 L0 이하다.
실제 응력 tensor, 투영된 세 traction, 변위, 잔차와 모델 해시를
`results/public_aluminum_validation_v15/cubic_tensor_static_scenarios/`에 저장했다.
기준 cubic frame의 n=(1,1,1)/sqrt(3), m=(-1,1,0)/sqrt(2),
m2=(-1,-1,2)/sqrt(6)를 명시하며, 임의로 보정된 시편 방향이라 하지 않는다.

이는 선택한 하중에서 **국소 정적 최소점이 연속적으로 복귀한다**는 결과다.
시간적분, 열적 전이, 주기 P(a,s), zero-stress 동적 hold를 실행한 결과는
아니므로 잔류소성/항복/피로 검증으로 해석하면 안 된다. 이 계산을 추가해도
production PDE/UI의 독립 전단 입력은 여전히 연결되지 않은 상태다.

## 6. 온도 스케일도 모델별로 분리

실제 `kinetic_calibration_workflow.build_time_basis_model`의 model-time
경로는 역사적 kT=.02를 사용한다. TwoRowLJ는 에너지 척도가 물리적으로
정해지지 않은 무차원 reference이므로 이를 그대로 Kelvin으로 부르면 안 된다.
반면 registry의 두 reduced hybrid 생성자는 명시적으로 `kT_ev=.02`를
받는다. **.02eV/kB=232.09036K**의 열에너지 등가값이지300K가 아니다.
이는 온도 단위 환산이며, 해당 collective-coordinate PMF의 열적 보정이
검증되었다는 뜻도 아니다. 현재 결과를 실온 Al 동역학으로 부르지 않는다.

유효한 physical calibration이 있을 때만 기존 physical 경로가
kT=kB*T_calibration/E0와 그 보정의 mobility ratio를 적용한다. 지금은 그러한
보정이 없으므로 기본 kT나 이동도를 조용히 바꾸지 않았다. 공개300K MD의
비교에는 실제 source 온도/원자수 정규화를 별도로 사용한다. 온도/시간/에너지
척도를 한 상수로 함께 맞추는 것은 허용되지 않는다.
