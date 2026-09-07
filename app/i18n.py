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
    ("nav.instructions", "휠: 확대/축소 · Shift: X축만 · Ctrl: Y축만 · 드래그: 이동 · 더블클릭/Home: 초기화", "Wheel: zoom · Shift: X only · Ctrl: Y only · drag: pan · double-click/Home: reset"),
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
    "survival": "plot.survival",
    "initiation_probability": "plot.initiation_probability",
    "first_passage_flux": "plot.first_passage_flux",
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
        "survival": "axis.survival",
        "initiation_probability": "axis.initiation_probability",
        "first_passage_flux": "axis.first_passage_flux",
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
        ),
    }.get(field, ())
    return {
        "title": tr(FIELD_TEXT_KEYS[field], language),
        "xlabel": tr("axis.time", language),
        "ylabel": tr(ylabel_key, language),
        "legend": tuple(tr(key, language) for key in legend_keys),
    }


validate_catalog()
