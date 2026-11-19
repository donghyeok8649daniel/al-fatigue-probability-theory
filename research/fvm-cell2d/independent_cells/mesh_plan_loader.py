"""Reuse the original single a/s axes without forming their intercell product."""
from pathlib import Path
import sys


def original_axes(a_level,s_level,outer):
    root=Path(__file__).resolve().parents[1]
    kernel=root/'outputs/fvm-integration/candidate/reference_kernel/outputs/raw-sg-convergence-20261108'
    sys.path.insert(0,str(kernel))
    from mesh_plan_v5 import mesh,validate
    validate()
    return mesh(a_level,'a',outer),mesh(s_level,'s')
