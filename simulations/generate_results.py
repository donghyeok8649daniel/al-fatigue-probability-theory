# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 결과 생성 스크립트를 순서대로 호출하는 통합 실행 진입점이다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: main
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Generate all active 1D normal layer-LJ reference results.

Run from repository root:
    python -m simulations.generate_results

The active research path is strictly one-dimensional and continuous-time:
normal stress -> spacing field / probability state -> spatial correlation and
push-forward structure -> normal-opening feasibility. Archived FCC and shear
libraries are not imported by this workflow.
"""

from simulations.run_normal_lj_chain import main as run_normal_chain
from simulations.run_normal_lj_timescale import main as run_normal_timescale
from simulations.run_normal_lj_energy_feasibility import (
    main as run_normal_energy_feasibility,
)
from simulations.run_normal_lj_distribution import main as run_normal_distribution
from simulations.run_normal_lj_closure_falsification import (
    main as run_normal_closure_falsification,
)
from simulations.run_normal_lj_closure_system_size import (
    main as run_normal_closure_system_size,
)
from simulations.run_normal_lj_spatial_correlation import (
    main as run_normal_spatial_correlation,
)
from simulations.run_normal_lj_pushforward_clue import (
    main as run_normal_pushforward_clue,
)


def main() -> None:
    run_normal_chain()
    run_normal_timescale()
    run_normal_energy_feasibility()
    run_normal_distribution()
    run_normal_closure_falsification()
    run_normal_closure_system_size()
    run_normal_spatial_correlation()
    run_normal_pushforward_clue()


if __name__ == "__main__":
    main()
