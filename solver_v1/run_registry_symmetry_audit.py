"""Trace the intrinsic-fault signal to its actual plane contributions."""
import math
import numpy as np

from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_geometry import SHOCKLEY_112
from .fcc111_lattice_sum import plane_power_sum_reciprocal, plane_exponential_sum_reciprocal
from .joint_fcc_interface_calibration import JointFit, build_joint_model, IDEAL_H
from .reference_eam_targets import MishinRigidFCCReference
from .run_full_fcc_calibration_audit import write_csv
from .run_matched_interface_study import ROOT, save_json
import json


def main():
    data=json.loads((ROOT/"fitted_candidates.json").read_text(encoding="utf-8"))["linear"]
    model=build_joint_model(JointFit("linear",data["decay"],np.array(data["coefficients"]),0.,None,None))
    # This decomposition uses the TARGET geometry, matching the fit observations.
    model.a0=IDEAL_H; interface=FCC111ActiveInterface(model,path_id=SHOCKLEY_112,tolerance=2e-12)
    rows=[]; slip=1/math.sqrt(3); density=model.density_params
    for k in range(1,9):
        base=interface._baseline_delta(k); shifted=interface._active_delta(k,slip)
        difference=0.
        for p,coefficient in ((6,4*model.p.epsilon*model.p.sigma_lj**12),
                              (3,-4*model.p.epsilon*model.p.sigma_lj**6)):
            left=plane_power_sum_reciprocal(k*IDEAL_H,base,p=p,geometry=model.geometry).value
            right=plane_power_sum_reciprocal(k*IDEAL_H,shifted,p=p,geometry=model.geometry).value
            difference+=k*coefficient*(right-left)
        left=plane_exponential_sum_reciprocal(k*IDEAL_H,base,kappa=density.kappa,
             amplitude=density.C_rho,geometry=model.geometry).value
        right=plane_exponential_sum_reciprocal(k*IDEAL_H,shifted,kappa=density.kappa,
             amplitude=density.C_rho,geometry=model.geometry).value
        rows.append(dict(model="analytic_scalar_joint",layer=k,pair_difference_ev_cell=difference,
                         normalized_plane_density_change=(right-left)/model.embedding.rho_ref))
    source=MishinRigidFCCReference()
    source_changes=[]; source_pair=0.
    for k in range(1,5):
        left=source.plane(k*source.h,source.geometry.abc_shift(k))
        right=source.plane(k*source.h,source.geometry.abc_shift(k)+source.geometry.tau)
        pair=k*(right[0]-left[0]); source_pair+=pair
        source_changes.append(right[1]-left[1])
        rows.append(dict(model="Mishin_source_rigid",layer=k,pair_difference_ev_cell=pair,
                         normalized_plane_density_change=(right[1]-left[1])/source._rho_bulk))
    depths=np.cumsum(source_changes[::-1])[::-1]
    embedding=float(2*np.sum(source.F(source._rho_bulk+depths)-source.F(source._rho_bulk)))
    value=interface.evaluate(IDEAL_H,slip)
    write_csv(ROOT/"intrinsic_fault_layer_decomposition.csv",rows)
    save_json(ROOT/"intrinsic_fault_mechanism.json",dict(
        first_plane_exact_radial_identity="w(h,tau)=w(h,2tau), because 3tau is a lattice translation",
        first_nontrivial_plane=2,shortest_reciprocal_attenuation=float(math.exp(-8*math.pi*math.sqrt(2)/3)),
        attenuation_is_not_a_complete_energy_bound=True,
        analytic_pair_ev_cell=value.pair_energy,analytic_embedding_ev_cell=value.embedding_energy,
        source_pair_ev_cell=source_pair,source_embedding_ev_cell=embedding,
        source_sum_ev_cell=source_pair+embedding,
        interpretation="scalar first-plane coordination is blind to this stacking change; radial farther-shell or angular information is required",
        no_change_to_LJ_pair=True))


if __name__=="__main__":
    main()
