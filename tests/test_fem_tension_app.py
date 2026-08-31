# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: TestTensionRunConfig, TestSolverCommand
# - 주요 함수/메서드: TestTensionRunConfig.test_default_unit_conversions_and_area
#   TestTensionRunConfig.test_invalid_geometry_is_rejected
#   TestTensionRunConfig.test_invalid_discrete_controls_are_rejected
#   TestTensionRunConfig.test_negative_mean_stress_is_allowed_but_negative_amplitude_is_not
#   TestSolverCommand._option_value, TestSolverCommand.test_command_contains_converted_inputs_numerically
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest
from pathlib import Path

from simulations.fem_tension_app import (
    TensionRunConfig,
    solver_command,
    validate_run_config,
)


class TestTensionRunConfig(unittest.TestCase):
    def test_default_unit_conversions_and_area(self) -> None:
        config = TensionRunConfig()
        self.assertAlmostEqual(config.length_m, 0.05)
        self.assertAlmostEqual(config.width_m, 0.01)
        self.assertAlmostEqual(config.thickness_m, 0.001)
        self.assertAlmostEqual(config.area_m2, 1.0e-5)
        self.assertAlmostEqual(config.young_pa, 69.0e9)

    def test_invalid_geometry_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(width_mm=0.0))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(thickness_mm=-1.0))

    def test_invalid_discrete_controls_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(elements=0))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(cycles=0))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(steps_per_cycle=1))

    def test_negative_mean_stress_is_allowed_but_negative_amplitude_is_not(self) -> None:
        validate_run_config(TensionRunConfig(stress_mean_mpa=-25.0, stress_amplitude_mpa=10.0))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(stress_amplitude_mpa=-1.0))

    def test_single_crystal_direction_is_required(self) -> None:
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(loading_h=0, loading_k=0, loading_l=0))

    def test_optional_cubic_constants_determine_axis_modulus(self) -> None:
        config_100 = TensionRunConfig(
            loading_h=1, loading_k=0, loading_l=0,
            cubic_c11_gpa=110, cubic_c12_gpa=60, cubic_c44_gpa=30)
        config_111 = TensionRunConfig(
            loading_h=1, loading_k=1, loading_l=1,
            cubic_c11_gpa=110, cubic_c12_gpa=60, cubic_c44_gpa=30)
        validate_run_config(config_100); validate_run_config(config_111)
        self.assertNotAlmostEqual(config_100.young_pa, config_111.young_pa)
        self.assertEqual(config_100.elastic_calibration_mode, "cubic_direction_projection")


class TestSolverCommand(unittest.TestCase):
    @staticmethod
    def _option_value(command: list[str], option: str) -> str:
        index = command.index(option)
        return command[index + 1]

    def test_command_contains_converted_inputs_numerically(self) -> None:
        config = TensionRunConfig(
            length_mm=60.0,
            width_mm=8.0,
            thickness_mm=2.0,
            young_gpa=70.0,
            elements=24,
            stress_mean_mpa=30.0,
            stress_amplitude_mpa=80.0,
            frequency_hz=5.0,
            cycles=3,
            steps_per_cycle=64,
        )
        command = solver_command(config, Path("solver"), Path("output"))

        self.assertEqual(int(self._option_value(command, "--elements")), 24)
        self.assertAlmostEqual(float(self._option_value(command, "--length-m")), 0.06)
        self.assertAlmostEqual(float(self._option_value(command, "--area-m2")), 1.6e-5)
        self.assertAlmostEqual(float(self._option_value(command, "--young-pa")), 70.0e9)
        self.assertAlmostEqual(float(self._option_value(command, "--stress-mean-mpa")), 30.0)
        self.assertAlmostEqual(float(self._option_value(command, "--stress-amplitude-mpa")), 80.0)
        self.assertAlmostEqual(float(self._option_value(command, "--frequency-hz")), 5.0)
        self.assertEqual(int(self._option_value(command, "--cycles")), 3)
        self.assertEqual(int(self._option_value(command, "--steps-per-cycle")), 64)
        self.assertEqual(command[-2:], ["--outdir", "output"])


if __name__ == "__main__":
    unittest.main()
