# === 한국어 파일 안내 시작 ===
# - 파일 역할: 통합 tensile FEM GUI의 입력 단위변환, 물리 입력 검증, C solver 명령행 변환을 검증한다.
# - 주요 클래스: TestTensionRunConfig, TestSolverCommand
# - 주요 함수/메서드: 기본 단위변환/면적 계산, 비정상 입력 거부, CLI 인자 생성 검증
# - 주의: 실제 C binary 실행과 2D/3D rendering은 GitHub Actions의 headless smoke 단계에서 별도로 검증한다.
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
