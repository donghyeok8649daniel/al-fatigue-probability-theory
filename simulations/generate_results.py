"""Generate all active 1D normal-LJ reference results.

Run from repository root:
    python -m simulations.generate_results

The active research path is strictly one-dimensional and continuous-time:
normal stress -> P(a,t) -> mean/energy constraints -> normal-opening feasibility.
Archived FCC and shear libraries are not imported by this workflow.
"""

from simulations.run_normal_lj_chain import main as run_normal_chain
from simulations.run_normal_lj_timescale import main as run_normal_timescale
from simulations.run_normal_lj_energy_feasibility import (
    main as run_normal_energy_feasibility,
)


def main() -> None:
    run_normal_chain()
    run_normal_timescale()
    run_normal_energy_feasibility()


if __name__ == "__main__":
    main()
