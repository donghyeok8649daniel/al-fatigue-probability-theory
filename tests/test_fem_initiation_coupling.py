# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: _uniform_element_history, test_fem_elements_receive_lambda_c_first_passage_fields
#   test_optional_initiation_channel_maps_by_element_identifier
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
from pathlib import Path

import numpy as np

from simulations.run_normal_lj_fem_probability import (
    _element_initiation_histories,
    write_initiation_element_history,
)
from simulations.fem_tension_app import initiation_snapshot


def _uniform_element_history():
    dtype = [("step", int), ("time_s", float), ("element", int),
             ("x_mid_m", float), ("stress_pa", float)]
    rows = []
    time = np.linspace(0, 0.5, 6)
    stress = 50e6 + 20e6*np.sin(2*np.pi*time/0.5)
    for element, x in ((0, 0.25), (1, 0.75)):
        rows.extend((step, t, element, x, value)
                    for step, (t, value) in enumerate(zip(time, stress)))
    return np.asarray(rows, dtype=dtype)


def test_fem_elements_receive_lambda_c_first_passage_fields(tmp_path: Path):
    elements = _uniform_element_history()
    histories = _element_initiation_histories(
        elements, 69e9, relaxation_time_s=1.0, inverse_temperature=80, grid_cells=60)
    assert histories.keys() == {0, 1}
    # Identical stress histories share one solve rather than treating element
    # count as a count of independent probability samples.
    assert histories[0] is histories[1]
    history = histories[0]
    assert np.all(np.diff(history.survival) <= 1e-12)
    assert np.allclose(history.initiation, 1-history.survival)
    assert np.all(history.tail_conditional == 0)

    output = tmp_path / "initiation_elements.csv"
    write_initiation_element_history(output, elements, histories, 1.0)
    header = output.read_text(encoding="utf-8").splitlines()[0]
    assert "survival" in header
    assert "initiation_probability" in header
    assert "hazard_per_s" in header


def test_optional_initiation_channel_maps_by_element_identifier():
    nodes = np.array([(2, 0.0, 0, 0.0, 0.0, 1.0), (2, 0.0, 1, 0.5, 0.0, 1.0),
                      (2, 0.0, 2, 1.0, 0.0, 1.0)],
        dtype=[("step", int), ("time_s", float), ("node", int), ("x_m", float),
               ("displacement_m", float), ("applied_stress_pa", float)])
    elements = np.array([(2, 0.0, 0, 0.25, 0.0, 1.0, 1.0),
                         (2, 0.0, 1, 0.75, 0.0, 1.0, 1.0)],
        dtype=[("step", int), ("time_s", float), ("element", int), ("x_mid_m", float),
               ("strain", float), ("stress_pa", float), ("applied_stress_pa", float)])
    initiation = np.array([(2, 1, 0.8, 0.2, 0.01), (2, 0, 0.9, 0.1, 0.02)],
        dtype=[("step", int), ("element", int), ("survival", float),
               ("initiation_probability", float), ("hazard_per_s", float)])
    snapshot = initiation_snapshot(nodes, elements, initiation, 2, "initiation")
    np.testing.assert_allclose(snapshot["scalar"], [0.1, 0.2])
    assert snapshot["field"] == "initiation"
