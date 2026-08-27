"""Generate active normal-deformation simulation results only.

Historical Rubin/shear/non-affine result generators are preserved under
libraries/shear/ and are intentionally excluded from the default workflow.
"""

from __future__ import annotations

from simulations.run_normal_lj_chain import main


if __name__ == "__main__":
    main()
