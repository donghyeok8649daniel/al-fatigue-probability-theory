import re
from pathlib import Path

import numpy as np

from app.i18n import (
    DEFAULT_LANGUAGE,
    FIELD_TEXT_KEYS,
    LANGUAGE_NAMES,
    Localizer,
    TRANSLATIONS,
    plot_strings,
    tr,
    validate_catalog,
)


def test_translation_catalog_has_identical_complete_locale_keys() -> None:
    validate_catalog()
    assert DEFAULT_LANGUAGE == "ko"
    assert set(TRANSLATIONS) == {"ko", "en"}
    assert set(TRANSLATIONS["ko"]) == set(TRANSLATIONS["en"])
    assert len(TRANSLATIONS["ko"]) == len(set(TRANSLATIONS["ko"]))
    assert all(TRANSLATIONS[language][key] for language in LANGUAGE_NAMES for key in TRANSLATIONS[language])


def test_every_literal_translation_key_used_by_desktop_ui_exists() -> None:
    source = Path(__file__).with_name("desktop_ui.py").read_text(encoding="utf-8")
    literal_keys = set(re.findall(r'_tr\("([a-z0-9_.]+)"', source))
    literal_keys.update(FIELD_TEXT_KEYS.values())
    for key in literal_keys:
        assert key in TRANSLATIONS["ko"]
        assert key in TRANSLATIONS["en"]


def test_runtime_language_switch_is_render_only() -> None:
    localizer = Localizer()
    numerical_inputs = {
        "young_gpa": "69",
        "stress_mean_mpa": "-500",
        "analysis_quality": "resolved",
    }
    result = {
        "model_time": np.array([0.0, 0.01]),
        "initiation_probability": np.array([0.0, 2.0e-12]),
        "survival": np.array([1.0, 1.0 - 2.0e-12]),
    }
    selected_field = "local_initiation_probability"
    view_limits = {selected_field: ((0.0, 0.01), (0.0, 3.0e-12))}
    identities = {key: id(value) for key, value in result.items()}
    renders: list[dict[str, object]] = []
    analysis_runs = 0

    localizer.subscribe(
        lambda: renders.append(plot_strings(selected_field, localizer.language))
    )
    assert localizer.set_language("en")
    assert localizer.language == "en"
    assert len(renders) == 1
    assert analysis_runs == 0
    assert numerical_inputs == {
        "young_gpa": "69",
        "stress_mean_mpa": "-500",
        "analysis_quality": "resolved",
    }
    assert selected_field == "local_initiation_probability"
    assert view_limits[selected_field] == ((0.0, 0.01), (0.0, 3.0e-12))
    assert {key: id(value) for key, value in result.items()} == identities


def test_plot_titles_axes_and_legends_switch_language() -> None:
    english = plot_strings("local_initiation_probability", "en")
    korean = plot_strings("local_initiation_probability", "ko")
    assert english == {
        "title": "Local initiation probability",
        "xlabel": "Dimensionless solver time",
        "ylabel": "Initiation probability [-]",
        "legend": (),
    }
    assert korean == {
        "title": "국소 균열 개시 확률",
        "xlabel": "무차원 솔버 시간",
        "ylabel": "균열 개시 확률 [-]",
        "legend": (),
    }
    assert plot_strings("strain_components", "ko")["legend"] == (
        "수직", "우물 내부", "소성", "총 변형률"
    )
    assert tr("legend.mass_residual", "ko") == "확률질량 보존 잔차"
    assert "균열" not in tr("legend.mass_residual", "ko")


def test_probability_and_plasticity_terms_use_required_precise_korean() -> None:
    expected = {
        "plot.local_initiation": "국소 균열 개시 확률",
        "plot.specimen_initiation": "시편 균열 개시 확률",
        "field.correlation_area": "특성 상관 면적 A_c",
        "field.stressed_area": "유효 응력 작용 면적 A_stressed",
        "diagnostic.N_eff": "유효 독립 영역 수 N_eff",
        "status.below_numerical_resolution": "수치 해상도 이하",
        "diagnostic.well_population": "구성 우물 점유 확률",
        "diagnostic.interwell_flux": "우물 간 확률 플럭스",
        "diagnostic.transition_onset": "모델 구성전이 개시 응력",
        "diagnostic.residual_plasticity": "잔류 소성 검증",
    }
    for key, value in expected.items():
        assert tr(key, "ko") == value
        assert tr(key, "en")
