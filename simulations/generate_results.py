"""Generate the active normal-deformation reference results only.

Run from repository root:
    python -m simulations.generate_results

The active research path is cyclic normal stress -> normal interatomic spacing
-> P(a,t) -> normal-opening instability. Non-normal proof-of-principle models
are not part of the active result workflow.
"""

from simulations.run_normal_lj_chain import main


if __name__ == "__main__":
    main()
