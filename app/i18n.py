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
    ("tree.material_load", "재료 / 하중", "Material / loading"),
    ("tree.axial_direction", "축 방향", "Axial direction"),
    ("tree.solve", "해석", "Solve"),
    ("tree.probability_pde", "확률 PDE (항상 적용)", "Probability PDE (always active)"),
    ("tree.spatial", "공간 표시: FVM / FEM", "Spatial display: FVM / FEM"),
    ("tree.post", "후처리", "Post"),
    ("tab.pre", "1  전처리", "1  PRE"),
    ("tab.solve", "2  해석", "2  SOLVE"),
    ("tab.post", "3  후처리", "3  POST"),
    ("section.axial_inputs", "물리 축방향 하중 입력", "Physical axial-load inputs"),
    ("field.young_gpa", "영률", "Young's modulus"),
    ("field.stress_mean", "평균 응력", "Mean stress"),
    ("field.stress_amplitude", "응력 진폭", "Stress amplitude"),
    ("field.model_frequency", "모델 주파수", "Model frequency"),
    ("field.cycles", "하중 사이클", "Load cycles"),
    ("field.output_resolution", "출력 해상도", "Output resolution"),
    ("field.axial_direction", "인장응력 방향", "Tensile-stress direction"),
    ("unit.gpa", "GPa", "GPa"),
    ("unit.mpa", "MPa", "MPa"),
    ("unit.cycles_model_time", "사이클 / 모델 시간", "cycles / model time"),
    ("unit.cycles", "사이클", "cycles"),
    ("unit.steps_cycle", "스텝/사이클", "steps/cycle"),
    ("unit.direction", "[x y z]", "[x y z]"),
    ("section.scope", "해석 범위", "Scope"),
    ("scope.text", "물리 응력은 E로 나눈 뒤 검증된 완화 축방향 κ로 매핑됩니다.\n현재 역학은 입력한 축에만 응력을 적용하며 횡방향 응력은 0입니다.", "Physical stress is divided by E and mapped with the verified relaxed axial κ.\nThe current mechanics applies stress only along the entered axis; transverse stress is zero."),
    ("section.probability_pde", "전체 2D 확률 PDE", "Full 2D probability PDE"),
    ("status.always_active", "선택 없이 항상 적용", "Always active"),
    ("section.spatial_backend", "공간 표시 해석 방식", "Spatial display backend"),
    ("section.grid_quality", "확률 격자 품질", "Probability grid quality"),
    ("quality.preview", "미리보기", "Preview"),
    ("quality.resolved", "정밀 해석", "Resolved"),
    ("quality.preview_option", "미리보기 (21 × 31, explicit)", "Preview (21 × 31, explicit)"),
    ("quality.resolved_option", "정밀 해석 (81 × 91, implicit)", "Resolved (81 × 91, implicit)"),
    ("solve.explanation", "P(a,s,t)를 직접 진화시킵니다.\n초기 상태: σ(t=0)의 Gibbs 분포.\nMonte Carlo 계수 없음.\n모델 시간은 아직 실제 초 단위로 보정되지 않았습니다.", "P(a,s,t) is evolved directly.\nInitial state: Gibbs at σ(t=0).\nNo Monte Carlo counting.\nModel time is not calibrated to physical seconds."),
    ("button.run", "해석 실행", "RUN ANALYSIS"),
    ("button.stop", "해석 중지", "STOP ANALYSIS"),
    ("button.stopping", "중지 중…", "STOPPING…"),
    ("summary.ready", "준비 완료. 확률과 모든 변형률 성분은 solver_v1에서 계산됩니다.\n", "Ready. Probability and all strain components come from solver_v1.\n"),
    ("status.ready", "준비 · 확률 PDE 항상 적용", "Ready · probability PDE always active"),
    ("status.stopping", "현재 PDE 스텝이 끝난 뒤 중지합니다…", "Stopping after the current PDE step…"),
    ("status.solving", "백그라운드에서 확률 PDE 직접 해석 중…", "Solving direct probability PDE in background…"),
    ("status.live", "PDE 모델 시간 {time:.5g} · S={survival:.12g} · ε={strain:.6g}", "PDE model time {time:.5g} · S={survival:.12g} · ε={strain:.6g}"),
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
    ("status.area_required", "시편 집계에는 두 면적 입력이 필요합니다.", "Both area inputs are required for specimen aggregation."),
    ("status.local_result_required", "먼저 국소 PDE 해석을 실행하십시오.", "Run the local PDE analysis first."),
    ("status.independence_assumption", "독립·동등 국소 영역 가정이며 A_c는 외부 보정값입니다.", "Assumes independent equivalent local regions; A_c is externally calibrated."),
    ("plot.local_survival", "국소 생존 확률", "Local survival probability"),
    ("plot.local_initiation", "국소 균열 개시 확률", "Local initiation probability"),
    ("plot.specimen_survival", "시편 생존 확률", "Specimen survival probability"),
    ("plot.specimen_initiation", "시편 균열 개시 확률", "Specimen initiation probability"),
    ("plot.well_populations", "구성 우물 점유 확률", "Configurational-well populations"),
    ("plot.interwell_flux", "우물 간 확률 플럭스", "Interwell probability flux"),
    ("plot.plastic_flow", "모델 소성 유동", "Model plastic flow"),
    ("plot.registry_transfer", "누적 순 레지스트리 이동", "Accumulated net registry transfer"),
    ("plot.specimen_unavailable", "국소 신호가 수치 해상도 이하이거나 면적 입력이 없어 시편 확률을 표시할 수 없습니다.", "Specimen probability is unavailable because the local signal is below resolution or area inputs are missing."),
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
    "first_passage_flux": "plot.first_passage_flux",
    "well_populations": "plot.well_populations",
    "interwell_flux": "plot.interwell_flux",
    "plastic_flow": "plot.plastic_flow",
    "registry_transfer": "plot.registry_transfer",
    "mass_diagnostics": "plot.mass_diagnostics",
}


def plot_strings(field: str, language: str = DEFAULT_LANGUAGE) -> dict[str, object]:
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
        "first_passage_flux": "axis.first_passage_flux",
        "well_populations": "axis.probability_mass",
        "interwell_flux": "axis.interwell_flux",
        "plastic_flow": "axis.plastic_flow",
        "registry_transfer": "axis.registry_transfer",
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
    }.get(field, ())
    return {
        "title": tr(FIELD_TEXT_KEYS[field], language),
        "xlabel": tr("axis.time", language),
        "ylabel": tr(ylabel_key, language),
        "legend": tuple(tr(key, language) for key in legend_keys),
    }


validate_catalog()
