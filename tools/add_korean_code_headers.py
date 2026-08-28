#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""저장소의 코드 파일 상단에 한국어 파일 안내 헤더를 자동 생성한다.

Python 파일은 AST를 사용해 top-level 함수/클래스와 클래스 메서드를 추출하고,
C/H 파일은 선언/정의 패턴을 읽어 주요 심볼을 요약한다. 기존 헤더가 있으면
동일한 marker 구간만 교체하므로 반복 실행해도 중복되지 않는다.
"""
from __future__ import annotations

import argparse
import ast
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

PY_START = "# === 한국어 파일 안내 시작 ==="
PY_END = "# === 한국어 파일 안내 끝 ==="
C_START = "/* === 한국어 파일 안내 시작 ==="
C_END = "=== 한국어 파일 안내 끝 === */"

PURPOSES = {
    "theory/__init__.py": "활성 1D normal layer-LJ 이론 패키지의 진입점과 공개 모듈 범위를 설명한다.",
    "theory/normal_lj_chain.py": "calibrated generalized-LJ layer interaction을 사용하는 보존적 1D normal chain을 적분하고, 에너지·spacing·instability 진단값을 계산한다.",
    "theory/normal_lj_closure_validation.py": "deterministic layer-spacing snapshot과 mean/energy 기반 one-point distribution closure를 같은 조건에서 비교하는 검증 도구를 제공한다.",
    "theory/normal_lj_distribution.py": "fixed-length/fixed-energy ensemble에서 유도한 large-M one-point spacing closure와 수치 적분·moment 계산을 구현한다.",
    "theory/normal_lj_energy_feasibility.py": "1D layer-LJ potential의 convexity와 support constraint를 이용해 crack-free energy feasibility bound를 계산한다.",
    "theory/normal_lj_pushforward.py": "deterministic spacing field에서 one-point density를 얻는 push-forward 관련 계산 도구를 모은다. harmonic/Taylor 기반 항목은 active 전역 분포 가정으로 사용하지 않는다.",
    "theory/normal_lj_spatial_correlation.py": "layer-spacing 순서정보를 나타내는 C_k, rho_k 및 scaled correlation 진단량을 계산한다.",
    "theory/normal_lj_timescale.py": "1D reduced layer model의 무차원/물리 시간·주파수 변환과 scale separation 진단을 계산한다.",
    "simulations/generate_results.py": "활성 1D normal layer-LJ 결과 생성 스크립트를 순서대로 호출하는 통합 실행 진입점이다.",
    "simulations/run_normal_lj_chain.py": "보존적 1D layer-LJ chain을 실행하고 cycle history, energy balance, instability 관련 결과를 저장한다.",
    "simulations/run_normal_lj_closure_falsification.py": "동일한 mean/energy에서 deterministic spacing distribution과 two-moment closure를 직접 비교해 반증시험을 수행한다.",
    "simulations/run_normal_lj_closure_system_size.py": "omega*M을 고정한 동적 유사성 조건에서 represented system size에 따른 closure mismatch를 검사한다.",
    "simulations/run_normal_lj_distribution.py": "large-M spacing closure의 energy sweep을 실행하고 distribution, variance, tail 결과를 생성한다.",
    "simulations/run_normal_lj_energy_feasibility.py": "crack-free support constraint 아래 energy-feasibility bound 예제를 계산하고 결과를 저장한다.",
    "simulations/run_normal_lj_pushforward_clue.py": "push-forward 관련 이전 진단 결과를 재현하는 스크립트다. harmonic/Taylor 근사는 active 전역 분포 모델로 해석하지 않는다.",
    "simulations/run_normal_lj_spatial_correlation.py": "동적으로 matched된 여러 chain size에서 spacing spatial correlation을 계산하고 CSV/JSON/figure를 생성한다.",
    "simulations/run_normal_lj_timescale.py": "물리 주파수와 reduced-model time scale의 대응을 계산해 reference 결과를 생성한다.",
    "tests/test_normal_lj_chain.py": "1D layer-LJ chain의 equilibrium, instability, 보존성 및 기본 수치 동작을 회귀검증한다.",
    "tests/test_normal_lj_closure_falsification.py": "closure-vs-mechanics 비교 지표와 반증시험 계산이 재현되는지 검증한다.",
    "tests/test_normal_lj_closure_system_size.py": "system-size sweep의 dynamic-similarity 규칙과 수치 진단 함수를 검증한다.",
    "tests/test_normal_lj_distribution.py": "spacing closure의 normalization, moment recovery, energy relation 및 수치 안정성을 검증한다.",
    "tests/test_normal_lj_energy_feasibility.py": "safe-energy bound와 관련 convexity/feasibility identity를 검증한다.",
    "tests/test_normal_lj_pushforward.py": "push-forward 보조 계산의 수학적 identity를 검증한다. harmonic/Taylor 항목은 과거 진단용 테스트로만 유지한다.",
    "tests/test_normal_lj_spatial_correlation.py": "C_k, rho_k, permutation reference 등 spatial-correlation 계산을 검증한다.",
    "tests/test_normal_lj_timescale.py": "무차원·물리 시간/주파수 변환과 scale 계산을 검증한다.",
    "firmware/include/fatigue_controller.h": "피로시험기 제어기의 설정값, 상태, 파형/제어 인터페이스를 선언하는 공개 헤더다.",
    "firmware/include/fatigue_hal.h": "load cell, actuator, 시간, 안전정지 등 하드웨어 의존 기능을 제어 로직에서 분리하기 위한 HAL 인터페이스 헤더다.",
    "firmware/src/fatigue_controller.c": "fatigue_controller.h에 선언된 파형 생성, force reference 계산, 제어 상태 갱신 및 안전 로직을 구현한다.",
    "firmware/src/host_test.c": "실제 MCU 없이 host 환경에서 firmware 제어 로직의 기본 동작을 확인하는 테스트 실행 파일이다.",
    "firmware/src/main_loop_example.c": "MCU main loop에서 HAL과 fatigue controller를 연결하는 사용 예시를 제공한다.",
    "firmware/CMakeLists.txt": "host-side firmware scaffold와 테스트 실행 파일을 빌드하기 위한 CMake 설정이다.",
}


def default_purpose(path: Path) -> str:
    rel = path.relative_to(ROOT).as_posix()
    if rel.startswith("theory/"):
        return "활성 1D normal layer-LJ 이론 계산에 사용하는 Python 모듈이다."
    if rel.startswith("simulations/"):
        return "활성 이론을 실행해 재현 가능한 수치 결과를 생성하는 Python 스크립트다."
    if rel.startswith("tests/"):
        return "활성 1D normal layer-LJ 코드의 수학적·수치적 동작을 검증하는 회귀 테스트다."
    if rel.startswith("firmware/include/"):
        return "피로시험기 firmware에서 공유하는 인터페이스와 자료형을 선언하는 헤더다."
    if rel.startswith("firmware/src/"):
        return "피로시험기 firmware의 구현 또는 host-side 검증용 C 소스다."
    return "이 프로젝트에서 사용하는 코드/빌드 파일이다."


def python_symbols(text: str) -> tuple[list[str], list[str]]:
    tree = ast.parse(text)
    funcs: list[str] = []
    classes: list[str] = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs.append(node.name)
        elif isinstance(node, ast.ClassDef):
            classes.append(node.name)
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    funcs.append(f"{node.name}.{item.name}")
    return funcs, classes


def c_symbols(text: str) -> tuple[list[str], list[str]]:
    # 함수 prototype/definition을 보수적으로 추출한다. 매크로/제어문은 제외한다.
    func_re = re.compile(
        r"(?m)^\s*(?:static\s+)?(?:inline\s+)?(?:const\s+)?"
        r"[A-Za-z_][\w\s\*]*?\s+([A-Za-z_]\w*)\s*\([^;{}]*\)\s*(?:;|\{)"
    )
    banned = {"if", "for", "while", "switch", "return", "sizeof"}
    funcs = []
    for name in func_re.findall(text):
        if name not in banned and name not in funcs:
            funcs.append(name)

    type_re = re.compile(r"(?m)^\s*typedef\s+(?:struct|enum)\s*(?:\w+\s*)?\{[\s\S]*?\}\s*(\w+)\s*;", re.MULTILINE)
    types = []
    for name in type_re.findall(text):
        if name not in types:
            types.append(name)
    return funcs, types


def wrap_items(prefix: str, items: list[str], comment_prefix: str, width: int = 105) -> list[str]:
    if not items:
        return [f"{comment_prefix} - {prefix}: 없음 또는 외부 선언만 사용"]
    lines: list[str] = []
    current = f"{comment_prefix} - {prefix}: "
    for item in items:
        piece = item if current.endswith(": ") else f", {item}"
        if len(current) + len(piece) > width:
            lines.append(current)
            current = f"{comment_prefix}   {item}"
        else:
            current += piece
    lines.append(current)
    return lines


def strip_python_header(text: str) -> str:
    if text.startswith(PY_START):
        end = text.find(PY_END)
        if end >= 0:
            end += len(PY_END)
            while end < len(text) and text[end] in "\r\n":
                end += 1
            return text[end:]
    return text


def make_python_header(path: Path, text: str) -> str:
    funcs, classes = python_symbols(text)
    rel = path.relative_to(ROOT).as_posix()
    lines = [
        PY_START,
        f"# - 파일 역할: {PURPOSES.get(rel, default_purpose(path))}",
    ]
    lines += wrap_items("주요 클래스", classes, "#")
    lines += wrap_items("주요 함수/메서드", funcs, "#")
    lines += [
        "# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.",
        PY_END,
        "",
    ]
    return "\n".join(lines) + text


def strip_c_header(text: str) -> str:
    if text.startswith(C_START):
        end = text.find(C_END)
        if end >= 0:
            end += len(C_END)
            while end < len(text) and text[end] in "\r\n":
                end += 1
            return text[end:]
    return text


def make_c_header(path: Path, text: str) -> str:
    funcs, types = c_symbols(text)
    rel = path.relative_to(ROOT).as_posix()
    lines = [
        C_START,
        f"파일 역할: {PURPOSES.get(rel, default_purpose(path))}",
        f"주요 자료형: {', '.join(types) if types else '없음 또는 다른 헤더에서 정의'}",
        f"주요 함수: {', '.join(funcs) if funcs else '없음 또는 선언 없음'}",
        "주의: 이 안내는 코드 탐색용이며 안전한 실제 하드웨어 구동에는 별도의 HAL 구현과 limit/interlock 검증이 필요하다.",
        C_END,
        "",
    ]
    return "\n".join(lines) + text


def make_cmake_header(path: Path, text: str) -> str:
    rel = path.relative_to(ROOT).as_posix()
    marker_start = "# === 한국어 빌드 파일 안내 시작 ==="
    marker_end = "# === 한국어 빌드 파일 안내 끝 ==="
    if text.startswith(marker_start):
        end = text.find(marker_end)
        if end >= 0:
            end += len(marker_end)
            while end < len(text) and text[end] in "\r\n":
                end += 1
            text = text[end:]
    return (
        f"{marker_start}\n"
        f"# - 파일 역할: {PURPOSES.get(rel, default_purpose(path))}\n"
        "# - 주요 기능: firmware library/host executable의 source, include path, link 구성을 지정한다.\n"
        f"{marker_end}\n\n{text}"
    )


def candidate_files() -> list[Path]:
    files: list[Path] = []
    for directory in ("theory", "simulations", "tests"):
        files.extend(sorted((ROOT / directory).glob("*.py")))
    files.extend(sorted((ROOT / "firmware" / "include").glob("*.h")))
    files.extend(sorted((ROOT / "firmware" / "src").glob("*.c")))
    cmake = ROOT / "firmware" / "CMakeLists.txt"
    if cmake.exists():
        files.append(cmake)
    return files


def transformed(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".py":
        base = strip_python_header(text)
        return make_python_header(path, base)
    if path.suffix in {".c", ".h"}:
        base = strip_c_header(text)
        return make_c_header(path, base)
    if path.name == "CMakeLists.txt":
        return make_cmake_header(path, text)
    return text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="헤더가 최신인지 검사만 한다")
    args = parser.parse_args()

    changed: list[str] = []
    for path in candidate_files():
        new_text = transformed(path)
        old_text = path.read_text(encoding="utf-8")
        if new_text != old_text:
            changed.append(path.relative_to(ROOT).as_posix())
            if not args.check:
                path.write_text(new_text, encoding="utf-8")

    if args.check and changed:
        print("한국어 코드 헤더 갱신 필요:")
        for rel in changed:
            print(f" - {rel}")
        return 1

    if changed:
        print(f"한국어 코드 헤더 갱신: {len(changed)}개 파일")
        for rel in changed:
            print(f" - {rel}")
    else:
        print("모든 대상 파일의 한국어 코드 헤더가 최신 상태입니다.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
