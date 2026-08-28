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
