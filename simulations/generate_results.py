"""Generate all active normal-deformation reference results.

Run from repository root:
    python -m simulations.generate_results

The active research path is cyclic normal stress -> normal interatomic spacing
-> P(a,t) -> normal-opening instability. The active result workflow includes
both the reduced 1D normal chain and the 3D FCC normal pair-lattice validation.
Shear-oriented historical work remains isolated under libraries/shear/.
"""

from simulations.run_normal_lj_chain import main as run_normal_chain
from simulations.run_fcc_normal_lj import main as run_fcc_normal_lj


def main() -> None:
    run_normal_chain()
    run_fcc_normal_lj()


if __name__ == "__main__":
    main()
