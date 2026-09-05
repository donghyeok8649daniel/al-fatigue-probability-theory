# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다.
# - 주요 클래스: TestTensionRunConfig, TestSolverCommand
# - 주요 함수/메서드: TestTensionRunConfig.test_default_unit_conversions_and_area
#   TestTensionRunConfig.test_invalid_geometry_is_rejected
#   TestTensionRunConfig.test_invalid_discrete_controls_are_rejected
#   TestTensionRunConfig.test_negative_mean_stress_is_allowed_but_negative_amplitude_is_not
#   TestTensionRunConfig.test_single_crystal_direction_is_required
#   TestTensionRunConfig.test_optional_cubic_constants_determine_axis_modulus
#   TestSolverCommand._option_value, TestSolverCommand.test_command_contains_converted_inputs_numerically
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from simulations.fem_tension_app import (
    TensionRunConfig,
    run_python_fem_solver,
    run_live_theory_solver,
    run_selected_solver,
    run_theory_spatial_solver,
    solver_command,
    theory_load_params,
    theory_solver_params,
    validate_run_config,
)
from simulations.fem_tension_ui import load_fem_history
from simulations.visualize_fem1d import load_numeric_csv


class TestTensionRunConfig(unittest.TestCase):
    def test_default_unit_conversions_and_area(self) -> None:
        config = TensionRunConfig()
        self.assertAlmostEqual(config.length_m, 0.05)
        self.assertAlmostEqual(config.width_m, 0.01)
        self.assertAlmostEqual(config.thickness_m, 0.001)
        self.assertAlmostEqual(config.area_m2, 1.0e-5)
        self.assertAlmostEqual(config.young_pa, 69.0e9)
        self.assertAlmostEqual(config.axial_stress_mpa(0.0), config.stress_mean_mpa)
        self.assertAlmostEqual(
            config.axial_stress_mpa(0.25 / config.frequency_hz),
            config.stress_max_mpa,
        )

    def test_circular_section_and_tensile_axis_are_independent_of_crystal_direction(self) -> None:
        config = TensionRunConfig(
            section_shape="circular",
            diameter_mm=6.0,
            loading_h=1,
            loading_k=1,
            loading_l=1,
            tensile_axis_x=0.0,
            tensile_axis_y=3.0,
            tensile_axis_z=4.0,
        )
        self.assertAlmostEqual(config.area_m2, np.pi * (0.006**2) / 4.0)
        np.testing.assert_allclose(config.tensile_unit_vector, [0.0, 0.6, 0.8])

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

    def test_invalid_poisson_ratio_and_tensile_axis_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(poisson_ratio=0.5))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(
                tensile_axis_x=0.0,
                tensile_axis_y=0.0,
                tensile_axis_z=0.0,
            ))

    def test_negative_mean_stress_is_allowed_but_negative_amplitude_is_not(self) -> None:
        validate_run_config(TensionRunConfig(stress_mean_mpa=-25.0, stress_amplitude_mpa=10.0))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(stress_amplitude_mpa=-1.0))
        with self.assertRaises(ValueError):
            validate_run_config(TensionRunConfig(theory_stress_scale_mpa=0.0))

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

    def test_axial_loading_exposes_stress_ratio(self) -> None:
        config = TensionRunConfig(
            loading_h=0,
            loading_k=0,
            loading_l=1,
            stress_mean_mpa=10.0,
            stress_amplitude_mpa=10.0,
        )
        self.assertAlmostEqual(config.stress_min_mpa, 0.0)
        self.assertAlmostEqual(config.stress_max_mpa, 20.0)
        self.assertAlmostEqual(config.stress_ratio, 0.0)


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

    def test_python_fem_fallback_preserves_uniform_bar_solution(self) -> None:
        config = TensionRunConfig(elements=4, cycles=1, steps_per_cycle=4)
        with TemporaryDirectory() as directory:
            completed = run_python_fem_solver(config, Path(directory))
            nodes, elements = load_fem_history(Path(directory))

        self.assertEqual(completed.returncode, 0)
        first_nodes = nodes[nodes["step"] == 0]
        first_elements = elements[elements["step"] == 0]
        np.testing.assert_allclose(first_elements["stress_pa"], 50.0e6, rtol=1.0e-12)
        np.testing.assert_allclose(first_elements["strain"], 50.0e6 / config.young_pa, rtol=1.0e-12)
        self.assertAlmostEqual(
            float(first_nodes[-1]["displacement_m"]),
            50.0e6 * config.length_m / config.young_pa,
        )

    def test_coupled_runner_always_writes_theory_and_selected_spatial_results(self) -> None:
        config = TensionRunConfig(elements=4, cycles=1, steps_per_cycle=4)
        with TemporaryDirectory() as directory:
            output = Path(directory)
            run_theory_spatial_solver(config, output, "FEM", auto_build=False)
            nodes, elements = load_fem_history(output)

            self.assertTrue((output / "theory" / "initiation_elements.csv").is_file())
            self.assertTrue((output / "spatial" / "elements.csv").is_file())
            self.assertTrue((output / "initiation_elements.csv").is_file())
            self.assertFalse((output / "sn_curve.csv").exists())
            self.assertFalse((output / "life_distribution.csv").exists())
            self.assertGreater(nodes.size, 0)
            self.assertGreater(elements.size, 0)

    def test_theory_load_uses_mean_amplitude_and_explicit_scale(self) -> None:
        config = TensionRunConfig(
            stress_mean_mpa=60.0,
            stress_amplitude_mpa=100.0,
            theory_stress_scale_mpa=40.0,
            steps_per_cycle=80,
        )
        load = theory_load_params(config)
        solver = theory_solver_params(config, load)

        self.assertAlmostEqual(load.force_min, 0.0)
        self.assertAlmostEqual(load.force_max, 4.0)
        self.assertAlmostEqual(load.value(2.5), 4.0)
        self.assertAlmostEqual(load.value(7.5), 0.0)
        self.assertLessEqual(solver.dt, 0.02)
        self.assertAlmostEqual(solver.dt * solver.record_stride, load.period / 80.0)

    def test_live_theory_stream_uses_native_solver_and_user_stop(self) -> None:
        records = []
        result = run_live_theory_solver(
            TensionRunConfig(
                stress_mean_mpa=50.0,
                stress_amplitude_mpa=85.0,
                cycles=1,
                steps_per_cycle=4,
            ),
            records.append,
            lambda: len(records) >= 4,
        )

        self.assertEqual(len(records), 4)
        self.assertTrue(all(records[i]["cycle"] < records[i + 1]["cycle"] for i in range(3)))
        self.assertAlmostEqual(records[1]["applied_stress_mpa"], 135.0)
        self.assertAlmostEqual(records[3]["applied_stress_mpa"], -35.0)
        self.assertEqual(records[3]["tensile_crack_drive_mpa"], 0.0)
        self.assertEqual(records[3]["force"], 0.0)
        self.assertEqual(result["time"].size, 0)

    def test_higher_stress_advances_first_passage(self) -> None:
        first_times = []
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for amplitude in (150.0, 200.0):
                output = root / str(int(amplitude))
                config = TensionRunConfig(
                    stress_mean_mpa=0.0,
                    stress_amplitude_mpa=amplitude,
                    theory_stress_scale_mpa=40.0,
                    frequency_hz=1.0,
                    cycles=1,
                    steps_per_cycle=50,
                )
                run_selected_solver(config, output, "Theory")
                initiation = load_numeric_csv(output / "initiation_elements.csv")
                crossed = initiation[initiation["initiation_probability"] > 0.0]
                self.assertGreater(crossed.size, 0)
                first_times.append(float(crossed[0]["time_s"]))

        self.assertLess(first_times[1], first_times[0])


if __name__ == "__main__":
    unittest.main()
