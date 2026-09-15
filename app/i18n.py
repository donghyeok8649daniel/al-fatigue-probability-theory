"""Small, dependency-free localization layer for the desktop application."""
from __future__ import annotations

from collections.abc import Callable


DEFAULT_LANGUAGE = "ko"
LANGUAGE_NAMES = {"ko": "한국어", "en": "English"}


# A row-oriented source makes duplicate keys detectable before the per-language
# dictionaries are built.  Keep mathematical symbols unchanged in translations.
_ROWS = (
    ("app.title", "알루미늄 피로 — 확률 PDE", "Al Fatigue — Probability PDE"),
    ("app.brand", "AL FATIGUE", "AL FATIGUE"),
    ("app.subtitle", "확률 PDE 직접 해석 · 축방향 하중 · 이론 항상 적용", "Direct probability PDE · axial loading · theory always active"),
    ("language.label", "언어 / Language:", "Language / 언어:"),
    ("language.ko", "한국어", "한국어"),
    ("language.en", "English", "English"),
    ("nav.study_tree", "해석 트리", "STUDY TREE"),
    ("tree.study", "축방향 피로 해석", "Axial fatigue study"),
    ("tree.pre", "전처리", "Pre"),
    ("tree.material_load", "재료 / 시간 설정", "Material / time settings"),
    ("tree.axial_direction", "축 방향", "Axial direction"),
    ("tree.solve", "해석", "Solve"),
    ("tree.probability_pde", "확률 PDE (항상 적용)", "Probability PDE (always active)"),
    ("tree.spatial", "공간 역학: 미연결", "Spatial mechanics: not connected"),
    ("tree.post", "후처리", "Post"),
    ("tab.model", "1  모델", "1  MODEL"),
    ("tab.mesh", "2  메시", "2  MESH"),
    ("tab.load", "3  면 하중", "3  FACE LOAD"),
    ("tab.pre", "4  조건 설정", "4  PRE"),
    ("tab.solve", "5  해석", "5  SOLVE"),
    ("tab.post", "6  후처리", "6  POST"),
    ("section.face_load", "메시 면 하중 / 응력 조건", "Mesh face load / stress conditions"),
    ("load.scope", "위 평균·진폭이 국소 축응력 PDE의 유일한 하중 입력입니다. 선택 면은 공간 해석 준비 데이터입니다. 현재 전단·임의 텐서와 체적 역학은 미연결입니다.", "The mean/amplitude above are the sole load inputs to the local axial PDE. Selected faces prepare spatial boundary data. Shear/custom tensors and volume mechanics are not connected yet."),
    ("load.pick_help", "왼쪽 클릭: 연결된 평면 선택 · Ctrl+클릭: 추가/제거 · 오른쪽 드래그: 회전 · 가운데 드래그: 이동 · 휠: 확대. STL/OBJ 곡면은 평면 조각으로 선택됩니다.", "Left click: connected planar face · Ctrl+click: toggle · Right drag: rotate · Middle drag: pan · Wheel: zoom. STL/OBJ curved surfaces select faceted patches."),
    ("load.region_custom", "직접 선택", "Picked faces"),
    ("load.next", "면 하중으로 →", "Face load →"),
    ("load.tensor_reset", "기본 텐서 식 복원", "Restore default tensor expression"),
    ("load.refine", "선택 면 메시 추가 세분화", "Refine selected faces"),
    ("load.refine_failed", "선택 면과 메시를 확인하세요. 최대 20,000면까지 지원합니다.", "Check the selected faces and mesh. Maximum 20,000 triangles."),
    ("load.show_map", "메시별 생존확률 외삽 지도", "Mesh survival extrapolation map"),
    ("map.title", "미인증 수학적 외삽 · 표면 메시", "Uncertified mathematical extrapolation · surface mesh"),
    ("map.assumption", "모든 표면에 같은 국소 P를 가정합니다. 각 삼각형은 S_i=(1-P)^(A_i/A_c). 색 차이는 삼각형 면적 차이이며 응력집중 예측이 아닙니다. 재메시하면 면별 값은 변하지만 전체 생존확률의 곱은 같습니다. 체적은 기하 추정값이며 하중의 유효체적은 미해석입니다.", "Uniform local P on every surface is assumed. Each triangle uses S_i=(1-P)^(A_i/A_c). Color variation reflects triangle areas, not stress concentration. Remeshing changes patch values but preserves their survival product. Volume is a geometric estimate; effective loaded volume has not been solved."),
    ("map.geometry", "표면적 {area} mm² · 폐곡면 체적 추정 {volume} mm³ (자기교차 검증 제외)", "Surface area {area} mm² · enclosed volume estimate {volume} mm³ (self-intersections not checked)"),
    ("map.volume_unknown", "계산 불가", "unavailable"),
    ("map.need_result", "먼저 국소 PDE를 해석하세요. 지도는 PDE를 재실행하지 않습니다.", "Run the local PDE first. The map does not rerun it."),
    ("map.need_area", "해석 탭에서 양의 특성 상관 면적 A_c [mm²]를 입력하세요.", "Enter a positive correlation area A_c [mm²] in Solve."),
    ("map.survival", "면별 생존 외삽 [-]", "Patch survival extrapolation [-]"),
    ("map.initiation", "면별 개시 외삽 [-]", "Patch initiation extrapolation [-]"),
    ("map.refresh", "현재 결과 갱신", "Refresh current result"),
    ("map.values", "미인증: 국소 P={p} · 면별 범위 {minimum}–{maximum} · 전체 표면 개시 외삽={total}", "Uncertified: local P={p} · patch range {minimum}–{maximum} · whole-surface initiation extrapolation={total}"),
    ("load.empty_selection", "메시 생성 후 면을 선택하세요.", "Generate a mesh and select faces first."),
    ("load.unsupported_solver", "현재 생산 솔버는 위의 축응력 평균·진폭만 지원합니다. 전단 또는 임의 텐서가 입력되어 해석을 중단했습니다. 전단을 0으로 하고 기본 텐서 식을 사용하세요.", "The production solver supports the axial mean/amplitude above. Analysis stopped because shear or a custom tensor was entered. Use zero shear and the default tensor expression."),
    ("load.face_region", "면 영역", "Face region"),
    ("load.region_all", "전체 표면", "All surface"),
    ("load.region_top", "상단 면", "Top face"),
    ("load.region_bottom", "하단 면", "Bottom face"),
    ("load.region_lateral", "측면", "Lateral faces"),
    ("load.normal_mean", "법선 평균 응력", "Normal mean stress"),
    ("load.normal_amplitude", "법선 응력 진폭", "Normal stress amplitude"),
    ("load.shear_mean", "전단 평균 응력", "Shear mean stress"),
    ("load.shear_amplitude", "전단 응력 진폭", "Shear stress amplitude"),
    ("load.apply", "면 하중 조건 저장", "Store face load condition"),
    ("load.selection_summary", "선택 면: {faces}개 · 면적: {area} mm²", "Selected faces: {faces} · area: {area} mm²"),
    ("load.not_prepared", "아직 면 하중 조건을 저장하지 않았습니다.", "No face load condition stored yet."),
    ("load.prepared", "{faces}개 면 선택 저장. 국소 PDE는 위 축응력을 사용합니다. 면의 공간 역학은 아직 미연결입니다.", "Stored {faces} selected faces. The local PDE uses the axial stress above; spatial face mechanics are not connected."),
    ("load.invalid", "응력 값은 유한한 숫자여야 합니다.", "Stress values must be finite numbers."),
    ("load.tensor_section", "3×3 시간 응력 텐서", "3×3 time-dependent stress tensor"),
    ("load.tensor_help", "행은 ;, 열은 , 로 구분합니다. 허용 변수: t, f, pi, normal_mean, normal_amp, shear_mean, shear_amp. 허용 함수: sin, cos. 기본은 정상응력 sin + 전단 cos입니다.", "Separate rows with ; and columns with ,. Allowed variables: t, f, pi, normal_mean, normal_amp, shear_mean, shear_amp. Allowed functions: sin, cos. Default uses normal sin + shear cos."),
    ("load.tensor_preview_button", "텐서 함수 검증 / t=0 미리보기", "Validate tensor function / preview at t=0"),
    ("load.tensor_result", "t=0 응력 텐서 [MPa]: {matrix}", "Stress tensor at t=0 [MPa]: {matrix}"),
    ("load.tensor_invalid", "3×3 텐서 식이 잘못되었습니다. 숫자·허용 변수·sin/cos만 사용하세요.", "Invalid 3×3 tensor expression. Use numbers, allowed variables, and sin/cos only."),
    ("geometry.radius", "원통 반경 [mm]", "Cylinder radius [mm]"),
    ("geometry.length", "원통 길이 [mm]", "Cylinder length [mm]"),
    ("geometry.scale", "가져오기 단위: 파일 1단위당 mm", "Import units: mm per file unit"),
    ("geometry.cylinder", "기본 원통 적용", "Apply cylinder"),
    ("geometry.default", "기본 원통", "Default cylinder"),
    ("geometry.import", "모델 가져오기 (STL / OBJ)", "Import model (STL / OBJ)"),
    ("geometry.formats", "파일이 없으면 기본 원통을 사용합니다. STL/삼각형 OBJ만 지원합니다. STEP/IGES는 CAD에서 STL로 내보내세요. 가져오기 단위를 반드시 확인하세요.", "Without a file, use the default cylinder. Supports STL / triangulated OBJ. Export STEP/IGES to STL in CAD. Verify import units."),
    ("geometry.summary", "형상: {source} · 정점 {vertices} · 삼각형 {faces}", "Geometry: {source} · vertices {vertices} · triangles {faces}"),
    ("geometry.target", "표면 메시 최대 변 길이 [mm]", "Surface mesh maximum edge [mm]"),
    ("geometry.solve_scope", "현재 Solve는 국소 P(a,s,t) 해석입니다. 가져온 형상의 체적 FVM/FEM 응력 해석과 공간 확률장은 아직 연결되지 않았습니다.", "Solve currently evolves local P(a,s,t). Volume FVM/FEM stress analysis and spatial probability fields on the imported geometry are not connected yet."),
    ("geometry.generate", "표면 메시 생성", "Generate surface mesh"),
    ("geometry.next_mesh", "다음: 메시 →", "Next: mesh →"),
    ("geometry.next_pre", "다음: 조건 설정 →", "Next: Pre →"),
    ("geometry.unmeshed", "형상 미리보기 — 아직 표면 메시를 생성하지 않았습니다.", "Geometry preview — surface mesh not generated yet."),
    ("geometry.mesh_summary", "삼각형 {faces} · 전체 표면적 {area} mm² · {closed}", "Triangles {faces} · total surface area {area} mm² · {closed}"),
    ("geometry.closed", "모든 변에 면 2개 (자가교차 검증 아님)", "Two faces per edge (not a self-intersection check)"),
    ("geometry.open", "열린 경계 또는 비다양체 변 있음", "Open or non-manifold edges"),
    ("geometry.scope", "표면 메시/형상 확인용입니다. 체적 FVM 역학은 미연결이며 Solve는 기존 국소 확률 PDE입니다. 표면적을 A_c 또는 유효 응력 작용 면적으로 자동 대입하지 않습니다. 공간 균열 확률 컬러맵은 아직 제공하지 않습니다.", "Surface geometry only. Volume FVM mechanics is not connected; Solve remains the local probability PDE. Surface area is not assigned to A_c or effective stressed area. No spatial crack-probability colormap yet."),
    ("geometry.invalid_mesh", "잘못된 형상 파일입니다. 유한 좌표와 비퇴화 삼각형을 확인하세요. 기존 형상은 유지됩니다.", "Invalid geometry. Check finite coordinates and nondegenerate triangles. Existing geometry is preserved."),
    ("geometry.positive_dimensions", "치수·단위 배율·메시 길이는 유한한 양수여야 합니다.", "Dimensions, unit scale and mesh size must be finite and positive."),
    ("geometry.mesh_limit", "경량 UI 한도(삼각형 20,000개 / 가져오기 10 MB)를 초과합니다. 더 큰 메시 길이나 단순화된 파일을 사용하세요.", "Lightweight UI limit exceeded (20,000 triangles / 10 MB import). Use a larger mesh size or simplified file."),
    ("geometry.triangles_only", "OBJ 면을 CAD에서 삼각형으로 변환한 후 가져오세요.", "Triangulate OBJ faces in CAD before import."),
    ("geometry.supported_formats", "STL 및 삼각형 OBJ만 지원합니다.", "Only STL and triangulated OBJ are supported."),
    ("section.axial_inputs", "재료 · 시간 · 출력 조건", "Material, time and output settings"),
    ("field.young_gpa", "영률", "Young's modulus"),
    ("field.stress_mean", "평균 응력", "Mean stress"),
    ("field.stress_amplitude", "응력 진폭", "Stress amplitude"),
    ("field.model_frequency", "모델 주파수", "Model frequency"),
    ("field.cycles", "하중 사이클", "Load cycles"),
    ("field.output_resolution", "출력 해상도", "Output resolution"),
    ("field.axial_direction", "방향 입력 (현재 미연결)", "Orientation (not connected)"),
    ("unit.gpa", "GPa", "GPa"),
    ("unit.mpa", "MPa", "MPa"),
    ("unit.cycles_model_time", "사이클 / 모델 시간", "cycles / model time"),
    ("unit.cycles", "사이클", "cycles"),
    ("unit.steps_cycle", "스텝/사이클", "steps/cycle"),
    ("unit.direction", "[x y z]", "[x y z]"),
    ("section.scope", "해석 범위", "Scope"),
    ("scope.text", "현재 UI: 단일 축응력 → κ·σ/E, 고정 χ=0.2, 상태 P(a,s,t).\n독립 전단·방향·벡터 registry/3D 시편 해석은 아직 연결되지 않았습니다.\n완전 격자의 이상 강도는 실제 시편의 항복·피로 강도가 아닙니다.", "Current UI: scalar axial stress → κ·σ/E, fixed χ=0.2, state P(a,s,t).\nIndependent shear, orientation, vector registry and 3D specimen solving are not connected.\nIdeal pristine-lattice strength is not specimen yield or fatigue strength."),
    ("section.probability_pde", "전체 2D 확률 PDE", "Full 2D probability PDE"),
    ("status.always_active", "선택 없이 항상 적용", "Always active"),
    ("section.spatial_backend", "공간 표시 해석 방식", "Spatial display backend"),
    ("section.grid_quality", "확률 격자 품질", "Probability grid quality"),
    ("quality.preview", "미리보기", "Preview"),
    ("quality.resolved", "정밀 해석", "Resolved"),
    ("quality.preview_option", "미리보기 (21 × 31, explicit)", "Preview (21 × 31, explicit)"),
    ("quality.resolved_option", "정밀 해석 (81 × 91, implicit)", "Resolved (81 × 91, implicit)"),
    ("solve.explanation", "P(a,s,t)를 직접 진화시킵니다.\n초기 상태: σ(t=0)의 Gibbs 분포.\nMonte Carlo 계수 없음.\n초/Hz는 동일 집단좌표의 동역학 보정이 있을 때만 사용합니다.", "P(a,s,t) is evolved directly.\nInitial state: Gibbs at σ(t=0).\nNo Monte Carlo counting.\nSeconds/Hz require matching collective-coordinate kinetic calibration."),
    ("button.run", "해석 실행", "RUN ANALYSIS"),
    ("button.stop", "해석 중지", "STOP ANALYSIS"),
    ("button.stopping", "중지 중…", "STOPPING…"),
    ("summary.ready", "준비 완료. 확률과 모든 변형률 성분은 solver_v1에서 계산됩니다.\n", "Ready. Probability and all strain components come from solver_v1.\n"),
    ("status.ready", "준비 · 확률 PDE 항상 적용", "Ready · probability PDE always active"),
    ("status.stopping", "현재 PDE 스텝이 끝난 뒤 중지합니다…", "Stopping after the current PDE step…"),
    ("status.solving", "백그라운드에서 확률 PDE 직접 해석 중…", "Solving direct probability PDE in background…"),
    ("status.live", "PDE 모델 시간 {time:.5g} · S={survival:.12g} · ε={strain:.6g}", "PDE model time {time:.5g} · S={survival:.12g} · ε={strain:.6g}"),
    ("status.live_physical", "PDE 물리 시간 {time:.5g} s · S={survival:.12g} · ε={strain:.6g}", "PDE physical time {time:.5g} s · S={survival:.12g} · ε={strain:.6g}"),
    ("status.complete", "완료 · 기록 {records}개 · S={survival:.12g}", "Complete · {records} records · S={survival:.12g}"),
    ("status.failed", "해석 실패", "Solve failed"),
    ("dialog.invalid_title", "입력값 오류", "Invalid input"),
    ("dialog.invalid_detail", "입력값을 확인하십시오:\n{detail}", "Check the input values:\n{detail}"),
    ("dialog.solver_title", "PDE 솔버 오류", "PDE solver error"),
    ("dialog.solver_detail", "PDE 해석 중 오류가 발생했습니다:\n{detail}", "The PDE solve failed:\n{detail}"),
    ("error.direction", "축 방향에는 0이 아닌 값 세 개가 필요합니다.", "Axial direction requires three values and must be nonzero."),
    ("summary.solving", "공간 표시: {backend}\n확률 해석 품질: {quality} ({na} × {ns}, {integrator})\n응력: {stress_min:.6g} ~ {stress_max:.6g} MPa\n무차원 응력: {reduced_min:.6g} ~ {reduced_max:.6g}\n모델 주기: {period:.8g}\n국소 pristine De_fast / De_slow: {de_fast:.6g} / {de_slow:.6g}\n완화 축방향 κ: {kappa:.12g}\n초기 Gibbs preload force: {preload:.12g}\n계산 중…\n", "Spatial display: {backend}\nProbability solve quality: {quality} ({na} × {ns}, {integrator})\nStress: {stress_min:.6g} to {stress_max:.6g} MPa\nReduced stress: {reduced_min:.6g} to {reduced_max:.6g}\nModel period: {period:.8g}\nLocal pristine De_fast / De_slow: {de_fast:.6g} / {de_slow:.6g}\nRelaxed axial κ: {kappa:.12g}\nInitial Gibbs preload force: {preload:.12g}\nComputing…\n"),
    ("summary.complete", "확률 원천: {source}\n해석 품질: {quality} {grid}, {integrator}\n초기 조건: {initial_condition}\n초기 응력: {initial_stress:.8g} MPa\n모델 주파수: {frequency:.8g} 사이클 / 모델 시간\n모델 주기: {period:.8g}\n국소 pristine tau_fast / tau_slow: {tau_fast:.8g} / {tau_slow:.8g}\n국소 pristine omega*tau_fast / omega*tau_slow: {de_fast:.8g} / {de_slow:.8g}\n국소 선형 |G| / 위상: {transfer:.8g} / {phase:.6g} deg\n완화 축방향 κ: {kappa:.12g}\n최종 생존 확률: {survival:.12g}\n최종 균열 개시 확률: {initiation:.12g}\n최대 |확률질량 보존 잔차|: {residual:.3e}\n시간 상태: {time_status}\n", "Probability source: {source}\nSolve quality: {quality} {grid}, {integrator}\nInitial condition: {initial_condition}\nInitial stress: {initial_stress:.8g} MPa\nModel frequency: {frequency:.8g} cycles / model time\nModel period: {period:.8g}\nLocal pristine tau_fast / tau_slow: {tau_fast:.8g} / {tau_slow:.8g}\nLocal pristine omega*tau_fast / omega*tau_slow: {de_fast:.8g} / {de_slow:.8g}\nLocal linear |G| / phase: {transfer:.8g} / {phase:.6g} deg\nRelaxed axial κ: {kappa:.12g}\nFinal survival probability: {survival:.12g}\nFinal initiation probability: {initiation:.12g}\nMax |mass-balance residual|: {residual:.3e}\nTime status: {time_status}\n"),
    ("model.probability_source", "N=1 Smoluchowski/Fokker–Planck 확률 PDE 직접 해석", "N=1 direct Smoluchowski/Fokker–Planck probability PDE"),
    ("model.initial_condition", "σ(t=0)의 조건부 Gibbs 분포", "conditional Gibbs at σ(t=0)"),
    ("model.time_warning", "모델 시간은 아직 실제 초 단위로 보정되지 않았습니다.", "Model time is not calibrated to physical seconds."),
    ("integrator.explicit", "explicit", "explicit"),
    ("integrator.implicit", "implicit", "implicit"),
    ("plot.stress", "적용 응력", "Applied stress"),
    ("plot.strain_components", "변형률 성분", "Strain components"),
    ("plot.total_strain", "총 축방향 변형률", "Total axial strain"),
    ("plot.normal_strain", "수직 변형률", "Normal strain"),
    ("plot.intrawell_strain", "우물 내부 변형률", "Intrawell strain"),
    ("plot.plastic_strain", "소성 변형률", "Plastic strain"),
    ("plot.survival", "생존 확률", "Survival probability"),
    ("plot.initiation_probability", "균열 개시 확률", "Initiation probability"),
    ("plot.first_passage_flux", "최초통과 플럭스", "First-passage flux"),
    ("plot.mass_diagnostics", "확률질량 수치진단", "Probability-mass diagnostics"),
    ("plot.empty", "PDE 결과를 보려면 해석을 실행하십시오.", "Run an analysis to view PDE results."),
    ("axis.time", "무차원 솔버 시간", "Dimensionless solver time"),
    ("axis.stress", "축방향 응력 [MPa]", "Axial stress [MPa]"),
    ("axis.strain", "축방향 변형률 [-]", "Axial strain [-]"),
    ("axis.normal_strain", "수직 변형률 [-]", "Normal strain [-]"),
    ("axis.intrawell_strain", "우물 내부 변형률 [-]", "Intrawell strain [-]"),
    ("axis.plastic_strain", "소성 변형률 [-]", "Plastic strain [-]"),
    ("axis.survival", "생존 확률 [-]", "Survival probability [-]"),
    ("axis.initiation_probability", "균열 개시 확률 [-]", "Initiation probability [-]"),
    ("axis.first_passage_flux", "최초통과 플럭스 [1/모델 시간]", "First-passage flux [1/model time]"),
    ("axis.probability_mass", "확률질량 [-]", "Probability mass [-]"),
    ("legend.normal", "수직", "Normal"),
    ("legend.intrawell", "우물 내부", "Intrawell"),
    ("legend.plastic", "소성", "Plastic"),
    ("legend.total", "총 변형률", "Total"),
    ("legend.raw_intact_mass", "원시 건전 확률질량 (수치)", "Raw intact mass (numerical)"),
    ("legend.raw_one_minus_survival", "1 - 원시 건전 확률질량 (수치)", "1 - raw intact mass (numerical)"),
    ("legend.absorbed_mass", "누적 흡수 확률질량", "Cumulative absorbed mass"),
    ("legend.mass_residual", "확률질량 보존 잔차", "Mass-balance residual"),
    ("legend.negative_correction", "최대 음의 질량 보정", "Maximum negative-mass correction"),
    ("legend.cumulative_repair", "누적 수치 질량 보정", "Accumulated numerical mass repair"),
    ("legend.flux_discrepancy", "플럭스-흡수 불일치", "Flux-absorption discrepancy"),
    ("nav.instructions", "휠: 확대/축소 · Shift: X축만 · Ctrl: Y축만 · 드래그: 이동 · 더블클릭/Home: 초기화", "Wheel: zoom · Shift: X only · Ctrl: Y only · drag: pan · double-click/Home: reset"),
    ("section.specimen_probability", "시편 확률 집계", "Specimen probability aggregation"),
    ("field.probability_scale", "확률 범위", "Probability scale"),
    ("option.local", "국소", "Local"),
    ("option.specimen", "시편", "Specimen"),
    ("field.correlation_area", "특성 상관 면적 A_c", "Characteristic correlation area A_c"),
    ("field.stressed_area", "유효 응력 작용 면적 A_stressed", "Effective stressed area A_stressed"),
    ("unit.mm2", "mm²", "mm²"),
    ("diagnostic.N_eff", "유효 독립 영역 수 N_eff", "Effective independent regions N_eff"),
    ("diagnostic.rare_event_floor", "국소 희귀사건 수치 바닥", "Local rare-event numerical floor"),
    ("diagnostic.plastic_floor", "소성 신호 수치 바닥", "Plastic-signal numerical floor"),
    ("diagnostic.well_population", "구성 우물 점유 확률", "Configurational-well population"),
    ("diagnostic.interwell_flux", "우물 간 확률 플럭스", "Interwell probability flux"),
    ("diagnostic.transition_onset", "모델 구성전이 개시 응력", "Model configurational-transition onset stress"),
    ("diagnostic.residual_plasticity", "잔류 소성 검증", "Residual plasticity validation"),
    ("diagnostic.configurational_barrier", "구성 장벽", "Configurational barrier"),
    ("diagnostic.opening_barrier", "개구 장벽", "Opening barrier"),
    ("section.mechanism_diagnostics", "국소 메커니즘 진단", "Local mechanism diagnostics"),
    ("diagnostic.barrier_scope", "최대 하중의 bound-basin 진단이며 실험적 항복 기준이 아닙니다.", "Bound-basin diagnostic at peak load; not an experimental yield criterion."),
    ("diagnostic.intrawell_motion", "구성 우물 내부 운동", "Intrawell configurational motion"),
    ("diagnostic.interwell_transfer", "우물 간 구성 이동", "Interwell configurational transfer"),
    ("diagnostic.model_plastic_flow", "모델 소성 유동", "Model plastic flow"),
    ("diagnostic.residual_plastic_deformation", "잔류 소성 변형", "Residual plastic deformation"),
    ("status.plastic_below_numerical_resolution", "소성 수치 해상도 이하", "Plasticity below numerical resolution"),
    ("status.below_numerical_resolution", "수치 해상도 이하", "Below numerical resolution"),
    ("status.resolved", "수치적으로 분해됨", "Numerically resolved"),
    ("status.requires_convergence", "수렴 검증 필요", "Convergence validation required"),
    ("status.area_required", "전처리의 시편 확률 항목에 유효한 면적을 입력하세요: A_c > 0, A_stressed ≥ 0 [mm²].", "Enter valid areas in Pre / specimen probability: A_c > 0, A_stressed ≥ 0 [mm²]."),
    ("status.local_result_required", "먼저 국소 PDE 해석을 실행하십시오.", "Run the local PDE analysis first."),
    ("status.specimen_live", "해석 중: 현재까지의 국소 흡수량으로 계산한 미인증 외삽입니다.", "Analysis running: uncertified extrapolation of opening mass accumulated so far."),
    ("status.independence_assumption", "독립·동등 국소 영역 가정이며 A_c는 외부 보정값입니다.", "Assumes independent equivalent local regions; A_c is externally calibrated."),
    ("plot.local_survival", "국소 생존 확률", "Local survival probability"),
    ("plot.local_initiation", "국소 균열 개시 확률", "Local initiation probability"),
    ("plot.specimen_survival", "시편 생존 확률", "Specimen survival probability"),
    ("plot.specimen_initiation", "시편 균열 개시 확률", "Specimen initiation probability"),
    ("plot.well_populations", "구성 우물 점유 확률", "Configurational-well populations"),
    ("plot.interwell_flux", "우물 간 확률 플럭스", "Interwell probability flux"),
    ("plot.plastic_flow", "모델 소성 유동", "Model plastic flow"),
    ("plot.registry_transfer", "누적 순 레지스트리 이동", "Accumulated net registry transfer"),
    ("plot.specimen_unavailable", "인증된 시편 확률은 아직 없습니다. ‘시편 확률 수학적 외삽값’을 선택하면 미인증 외삽을 볼 수 있습니다. 물리적 인증에는 국소 신호의 수렴·분해능 검증이 필요합니다.", "Certified specimen probability is unavailable. Select Mathematical specimen extrapolation to view the uncertified estimate. Local convergence and signal resolution are required for certification."),
    ("axis.interwell_flux", "확률 플럭스 [1/모델 시간]", "Probability flux [1/model time]"),
    ("axis.plastic_flow", "변형률률 [1/모델 시간]", "Strain rate [1/model time]"),
    ("axis.registry_transfer", "레지스트리 우물 지수 질량 [-]", "Registry well-index mass [-]"),
    ("legend.well_minus", "P_-1", "P_-1"),
    ("legend.well_zero", "P_0", "P_0"),
    ("legend.well_plus", "P_+1", "P_+1"),
    ("legend.interwell_net", "순 우물 간 플럭스", "Net interwell flux"),
    ("legend.interwell_gross", "총 양방향 플럭스", "Gross bidirectional flux"),
    ("legend.net_plastic_flow", "순 우물 간 소성 유동", "Net interwell plastic flow"),
    ("legend.gross_slip_activity", "총 구성 슬립 활동", "Gross configurational slip activity"),
    ("legend.selective_opening", "개구 선택손실 기여 (소성 아님)", "Selective opening-loss contribution (not plastic)"),
    ("legend.registry_transfer", "누적 순 우물 이동", "Accumulated net well transfer"),
    ("legend.registry_absorption", "흡수된 우물지수 모멘트", "Absorbed well-index moment"),
    ("section.load_presets", "하중 해석 프리셋", "Load interpretation presets"),
    ("field.load_preset", "편집 가능한 하중 프리셋", "Editable load preset"),
    ("preset.custom", "사용자 입력", "Custom input"),
    ("preset.small_signal", "소신호 검증", "Small-signal verification"),
    ("preset.moderate", "중간 비선형 메커니즘 탐색", "Moderate nonlinear mechanism probe"),
    ("preset.compression", "압축 대조군", "Compression control"),
    ("preset.extreme", "극한 메커니즘 스트레스 테스트", "Extreme mechanism stress test"),
    ("button.apply_preset", "프리셋 적용", "APPLY PRESET"),
    ("stress.context", "σ_min={sigma_min:.6g} MPa · σ_max={sigma_max:.6g} MPa · σ/E_min={reduced_min:.6g} · σ/E_max={reduced_max:.6g}\n{regime}", "σ_min={sigma_min:.6g} MPa · σ_max={sigma_max:.6g} MPa · σ/E_min={reduced_min:.6g} · σ/E_max={reduced_max:.6g}\n{regime}"),
    ("stress.regime.small", "소신호 검증 범위", "Small-signal verification regime"),
    ("stress.regime.moderate", "비선형 메커니즘 탐색 범위", "Nonlinear mechanism-probe regime"),
    ("stress.regime.compression", "인장 opening 대조용 전압축 하중", "All-compressive tensile-opening control"),
    ("stress.regime.extreme", "주의: 정량 보정된 순수 Al 피로하중이 아닌 극한 메커니즘 스트레스 테스트입니다.", "Caution: extreme mechanism stress test, not a quantitatively calibrated pure-Al fatigue load."),
    ("field.time_basis", "시간 기준", "Time basis"),
    ("option.model_time", "모델 시간", "Model time"),
    ("option.physical_time", "물리 시간", "Physical time"),
    ("button.load_kinetics", "동역학 보정 불러오기…", "Load kinetic calibration…"),
    ("dialog.kinetics_title", "집단좌표 동역학 보정", "Collective-coordinate kinetic calibration"),
    ("filetype.kinetics_json", "동역학 보정 JSON", "Kinetic calibration JSON"),
    ("error.kinetics_busy", "해석 중에는 동역학 보정을 변경할 수 없습니다.", "Kinetic calibration cannot be changed during a solve."),
    ("error.kinetics_invalid", "보정을 적용하지 않았습니다. 모델·좌표·단위·온도·출처를 확인하세요.\n{detail}", "Calibration was not applied. Check model, coordinates, units, temperature and source.\n{detail}"),
    ("error.time_frequency", "시간 기준 변경 전에 양의 유한 주파수를 입력하세요.", "Enter a finite positive frequency before changing time basis."),
    ("status.kinetics_loaded", "보정 파일을 확인했습니다. 기존 결과는 그대로이며 새 해석에 물리 시간을 선택할 수 있습니다.", "Calibration file checked. Existing results are unchanged; physical time is available for a new solve."),
    ("status.kinetics_details", "입력 보정: t0={t0:.6g} s, T={temperature:.6g} K\n출처: {source}\n파일 일치 검사이며 재료·피로 검증 인증은 아닙니다.", "Loaded calibration: t0={t0:.6g} s, T={temperature:.6g} K\nSource: {source}\nFile consistency is not material/fatigue validation."),
    ("field.frequency", "주파수", "Frequency"),
    ("unit.hz", "Hz", "Hz"),
    ("axis.time_seconds", "시간 [s]", "Time [s]"),
    ("axis.first_passage_flux_seconds", "최초통과 플럭스 [1/s]", "First-passage flux [1/s]"),
    ("axis.interwell_flux_seconds", "확률 플럭스 [1/s]", "Probability flux [1/s]"),
    ("axis.plastic_flow_seconds", "변형률률 [1/s]", "Strain rate [1/s]"),
    ("status.kinetic_uncalibrated", "집단좌표 이동도 미보정: 실제 초와 Hz는 아직 사용할 수 없습니다.", "Collective-coordinate mobility is uncalibrated: physical seconds and Hz are not yet available."),
    ("status.kinetic_calibrated", "동역학 이동도 보정됨: 물리 시간 변환이 활성화되었습니다.", "Kinetic mobility calibrated: physical-time conversion is enabled."),
    ("status.physical_time_unavailable", "유효한 동역학 이동도 보정이 없어 물리 시간을 선택할 수 없습니다.", "Physical time is unavailable without a valid kinetic-mobility calibration."),
    ("diagnostic.t0", "모델 시간척도 t0", "Model time scale t0"),
    ("diagnostic.physical_mobility_a", "물리 이동도 M_a", "Physical mobility M_a"),
    ("diagnostic.physical_mobility_s", "물리 이동도 M_s", "Physical mobility M_s"),
    ("diagnostic.calibration_source", "동역학 보정 출처", "Kinetic calibration source"),
    ("button.run_convergence", "수렴 검증 실행", "RUN CONVERGENCE CHECK"),
    ("button.convergence_running", "수렴 검증 중…", "CHECKING CONVERGENCE…"),
    ("status.convergence_complete", "수렴 검증 완료", "Convergence check complete"),
    ("status.preview_uncertified", "미리보기 결과는 수렴 검증되지 않았습니다.", "Preview results are not convergence-certified."),
    ("diagnostic.local_probability", "국소 물리 PDE 확률", "Local physical PDE probability"),
    ("diagnostic.specimen_extrapolation", "시편 수학적 외삽값 (미인증)", "Mathematical specimen extrapolation (uncertified)"),
    ("diagnostic.specimen_certified", "인증된 물리 시편 확률", "Certified physical specimen probability"),
    ("plot.specimen_extrapolation", "시편 확률 수학적 외삽값 (미인증)", "Mathematical specimen probability extrapolation (uncertified)"),
    ("plot.specimen_extrapolation_note", "이 값은 독립영역 수학적 외삽이며 수렴 인증된 물리 확률이 아닙니다.", "This is an independent-region mathematical extrapolation, not a convergence-certified physical probability."),
    ("section.cycle_diagnostics", "사이클별 최초통과 진단", "Per-cycle first-passage diagnostics"),
    ("cycle.table_header", "사이클 | 평균 ε | 진폭 | 평균 ε_a | 평균 ε_xi | Δε_p | ε_p | 순이동 | 총활동 | 흡수 | ΔG_s | ΔG_open | 최대플럭스 | 종료S", "Cycle | mean ε | amplitude | mean ε_a | mean ε_xi | Δε_p | ε_p | net transfer | gross activity | absorbed | ΔG_s | ΔG_open | peak flux | end S"),
    ("summary.model_time_basis", "시간 기준: 모델 시간\n집단좌표 이동도 미보정: 실제 초와 Hz는 아직 사용할 수 없습니다.\n", "Time basis: model time\nCollective-coordinate mobility is uncalibrated: physical seconds and Hz are not yet available.\n"),
    ("field.energy_model", "에너지 모델", "Energy model"),
    ("energy_model.lj_reference", "TwoRowLJ 기준 모델", "TwoRowLJ reference"),
    ("energy_model.hybrid_hypothetical", "해석적 LJ–EAM 가상 민감도 모델", "Analytic LJ–EAM hypothetical sensitivity"),
    ("energy_model.al_best_feasible", "Al 목표 최적가능 하이브리드", "Al-target best-feasible hybrid"),
    ("summary.energy_model", "활성 에너지 모델: {name}\nPython 클래스: {python_class}\n보정 상태: {status}\n매개변수 출처: {source}\na0 / b / chi / kappa: {a0:.10g} / {b:.8g} / {chi:.8g} / {kappa:.10g}\n", "Active energy model: {name}\nPython class: {python_class}\nCalibration status: {status}\nParameter source: {source}\na0 / b / chi / kappa: {a0:.10g} / {b:.8g} / {chi:.8g} / {kappa:.10g}\n"),
    ("plot.gross_registry_activity", "누적 총 레지스트리 활동", "Cumulative gross registry activity"),
    ("axis.gross_registry_activity", "누적 우물 횡단 확률질량 [-]", "Cumulative well-crossing probability mass [-]"),
    ("legend.registry_forward", "증가 s 방향", "Increasing-s activity"),
    ("legend.registry_backward", "감소 s 방향", "Decreasing-s activity"),
    ("legend.registry_gross", "총 양방향 활동", "Gross bidirectional activity"),
    ("diagnostic.plasticity_status", "소성 해상도 상태", "Plasticity resolution status"),
    ("diagnostic.max_plastic_strain", "최대 |소성 변형률|", "Max |plastic strain|"),
    ("diagnostic.final_plastic_strain", "최종 소성 변형률", "Final plastic strain"),
    ("diagnostic.cumulative_net_registry", "누적 순 레지스트리 이동", "Cumulative net registry transfer"),
    ("diagnostic.cumulative_gross_registry", "누적 총 레지스트리 활동", "Cumulative gross registry activity"),
    ("status.plasticity_unresolved", "수치적 또는 미해결", "Numerical or unresolved"),
    ("status.plasticity_resolved", "우물 간 이동 해상됨", "Resolved interwell transfer"),
    ("summary.physical_time_basis", "시간 기준: 물리 시간\n주파수: {frequency_hz:.8g} Hz\n주기: {period_seconds:.8g} s\n총 시간: {duration_seconds:.8g} s\nt0: {t0:.8g} s\ntau_fast / tau_slow: {tau_fast:.8g} / {tau_slow:.8g} s\nM_a: {mobility_a:.8g} m²/(J·s)\nM_s: {mobility_s:.8g} m²/(J·s)\n출처: {source}\n", "Time basis: physical time\nFrequency: {frequency_hz:.8g} Hz\nPeriod: {period_seconds:.8g} s\nTotal time: {duration_seconds:.8g} s\nt0: {t0:.8g} s\ntau_fast / tau_slow: {tau_fast:.8g} / {tau_slow:.8g} s\nM_a: {mobility_a:.8g} m²/(J·s)\nM_s: {mobility_s:.8g} m²/(J·s)\nSource: {source}\n"),
)


def _build_translations() -> dict[str, dict[str, str]]:
    keys = [row[0] for row in _ROWS]
    if len(keys) != len(set(keys)):
        duplicates = sorted({key for key in keys if keys.count(key) > 1})
        raise RuntimeError(f"duplicate translation keys: {duplicates}")
    return {
        "ko": {key: ko for key, ko, _en in _ROWS},
        "en": {key: en for key, _ko, en in _ROWS},
    }


TRANSLATIONS = _build_translations()


def validate_catalog() -> None:
    """Raise when a locale is missing or has a different translation key set."""

    if set(TRANSLATIONS) != set(LANGUAGE_NAMES):
        raise RuntimeError("translation locales and language selector differ")
    expected = set(TRANSLATIONS[DEFAULT_LANGUAGE])
    for language, messages in TRANSLATIONS.items():
        if set(messages) != expected:
            raise RuntimeError(f"translation key mismatch for {language}")
        empty = sorted(key for key, value in messages.items() if not value)
        if empty:
            raise RuntimeError(f"empty translations for {language}: {empty}")


def tr(key: str, language: str = DEFAULT_LANGUAGE, **values: object) -> str:
    """Translate one stable key and optionally apply named format arguments."""

    try:
        template = TRANSLATIONS[language][key]
    except KeyError as exc:
        raise KeyError(f"unknown translation {language}:{key}") from exc
    return template.format(**values) if values else template


class Localizer:
    """Runtime locale state; changing it invokes render-only listeners."""

    def __init__(self, language: str = DEFAULT_LANGUAGE) -> None:
        if language not in TRANSLATIONS:
            raise ValueError(f"unsupported language: {language}")
        self._language = language
        self._listeners: list[Callable[[], None]] = []

    @property
    def language(self) -> str:
        return self._language

    def text(self, key: str, **values: object) -> str:
        return tr(key, self._language, **values)

    def subscribe(self, listener: Callable[[], None]) -> None:
        self._listeners.append(listener)

    def set_language(self, language: str) -> bool:
        if language not in TRANSLATIONS:
            raise ValueError(f"unsupported language: {language}")
        if language == self._language:
            return False
        self._language = language
        for listener in tuple(self._listeners):
            listener()
        return True


FIELD_TEXT_KEYS = {
    "stress": "plot.stress",
    "strain_components": "plot.strain_components",
    "strain": "plot.total_strain",
    "normal_strain": "plot.normal_strain",
    "intrawell_strain": "plot.intrawell_strain",
    "plastic_strain": "plot.plastic_strain",
    "local_survival_probability": "plot.local_survival",
    "local_initiation_probability": "plot.local_initiation",
    "specimen_survival_probability": "plot.specimen_survival",
    "specimen_initiation_probability": "plot.specimen_initiation",
    "specimen_probability_extrapolation": "plot.specimen_extrapolation",
    "first_passage_flux": "plot.first_passage_flux",
    "well_populations": "plot.well_populations",
    "interwell_flux": "plot.interwell_flux",
    "plastic_flow": "plot.plastic_flow",
    "registry_transfer": "plot.registry_transfer",
    "gross_registry_activity": "plot.gross_registry_activity",
    "mass_diagnostics": "plot.mass_diagnostics",
}


def plot_strings(
    field: str,
    language: str = DEFAULT_LANGUAGE,
    *,
    time_basis: str = "model",
) -> dict[str, object]:
    """Return localized plot metadata without touching numerical arrays."""

    if field not in FIELD_TEXT_KEYS:
        raise KeyError(field)
    ylabel_key = {
        "stress": "axis.stress",
        "strain_components": "axis.strain",
        "strain": "axis.strain",
        "normal_strain": "axis.normal_strain",
        "intrawell_strain": "axis.intrawell_strain",
        "plastic_strain": "axis.plastic_strain",
        "local_survival_probability": "axis.survival",
        "local_initiation_probability": "axis.initiation_probability",
        "specimen_survival_probability": "axis.survival",
        "specimen_initiation_probability": "axis.initiation_probability",
        "specimen_probability_extrapolation": "axis.initiation_probability",
        "first_passage_flux": "axis.first_passage_flux",
        "well_populations": "axis.probability_mass",
        "interwell_flux": "axis.interwell_flux",
        "plastic_flow": "axis.plastic_flow",
        "registry_transfer": "axis.registry_transfer",
        "gross_registry_activity": "axis.gross_registry_activity",
        "mass_diagnostics": "axis.probability_mass",
    }[field]
    legend_keys = {
        "strain_components": (
            "legend.normal", "legend.intrawell", "legend.plastic", "legend.total"
        ),
        "mass_diagnostics": (
            "legend.raw_intact_mass",
            "legend.raw_one_minus_survival",
            "legend.absorbed_mass",
            "legend.mass_residual",
            "legend.negative_correction",
            "legend.cumulative_repair",
            "legend.flux_discrepancy",
        ),
        "plastic_flow": (
            "legend.net_plastic_flow",
            "legend.gross_slip_activity",
            "legend.selective_opening",
        ),
        "registry_transfer": (
            "legend.registry_transfer",
            "legend.registry_absorption",
        ),
        "gross_registry_activity": (
            "legend.registry_forward",
            "legend.registry_backward",
            "legend.registry_gross",
        ),
    }.get(field, ())
    if time_basis not in {"model", "physical"}:
        raise ValueError("time_basis must be 'model' or 'physical'")
    if time_basis == "physical":
        ylabel_key = {
            "first_passage_flux": "axis.first_passage_flux_seconds",
            "interwell_flux": "axis.interwell_flux_seconds",
            "plastic_flow": "axis.plastic_flow_seconds",
        }.get(field, ylabel_key)
    return {
        "title": tr(FIELD_TEXT_KEYS[field], language),
        "xlabel": tr(
            "axis.time_seconds" if time_basis == "physical" else "axis.time",
            language,
        ),
        "ylabel": tr(ylabel_key, language),
        "legend": tuple(tr(key, language) for key in legend_keys),
    }


validate_catalog()
