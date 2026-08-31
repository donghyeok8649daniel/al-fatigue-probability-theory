# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: test_miller_direction_normalization_and_zero_rejection
#   test_canonical_cubic_direction_formulas, test_isotropic_cubic_limit_is_orientation_independent
#   test_unstable_cubic_constants_are_rejected
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import math

import numpy as np
import pytest

from theory.cubic_normal_orientation import (
    CubicElasticConstants,
    directional_young_modulus,
    miller_unit_vector,
)


def test_miller_direction_normalization_and_zero_rejection():
    np.testing.assert_allclose(miller_unit_vector(1, 1, 0), [2**-0.5, 2**-0.5, 0])
    with pytest.raises(ValueError):
        miller_unit_vector(0, 0, 0)


def test_canonical_cubic_direction_formulas():
    constants = CubicElasticConstants(110e9, 60e9, 30e9)
    s11, s12, s44 = constants.compliances
    q = s11-s12-0.5*s44
    assert directional_young_modulus(constants, 1, 0, 0) == pytest.approx(1/s11)
    assert directional_young_modulus(constants, 1, 1, 0) == pytest.approx(1/(s11-q/2))
    assert directional_young_modulus(constants, 1, 1, 1) == pytest.approx(1/(s11-2*q/3))


def test_isotropic_cubic_limit_is_orientation_independent():
    # lambda=40 GPa, mu=30 GPa gives C11=lambda+2mu, C12=lambda, C44=mu.
    constants = CubicElasticConstants(100e9, 40e9, 30e9)
    values = [directional_young_modulus(constants, *direction)
              for direction in ((1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 3, 5))]
    assert max(values)-min(values) < 1e-4*values[0]


def test_unstable_cubic_constants_are_rejected():
    with pytest.raises(ValueError, match="stability"):
        directional_young_modulus(CubicElasticConstants(50e9, 60e9, 20e9), 1, 0, 0)
