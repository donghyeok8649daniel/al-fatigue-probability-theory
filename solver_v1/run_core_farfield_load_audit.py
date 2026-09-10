"""Executed nominal-vs-actual atomistic load/units check, static research only."""
import argparse
import hashlib
from pathlib import Path

from .core_farfield_load_audit import nominal_shear_audit
from .current_material_rows import CurrentMaterialScrewCore
from .isolated_screw_core import ScrewFarField
from .run_current_material_core import load_current_material
from .run_source_core_reference import load_source_material,build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('fresh traction-audit output required')
    rows=[];bindings={}
    for label,material,builder in [('unchanged',load_current_material(),CurrentMaterialScrewCore),
            ('joint_trial',load_current_material(args.trial),CurrentMaterialScrewCore),
            ('Mishin_target',load_source_material(),build_source_core)]:
        model,tensor,meta=material;bindings[label]=meta['parameter_sha256']
        far=ScrewFarField(tensor,model.geometry.b,(model.geometry.b*3**.5/12,model.h/2))
        for ring in (7,11,15):
            core=builder(model,far,free_radius=.85,ring=ring,burgers_sign=0,tolerance=2e-13)
            for stress in (-50.,-20.,-5.,0.,5.,20.,50.):
                result=nominal_shear_audit(core,stress,length_scale_m=meta['length_scale_m'])
                P=result['actual_traction_columns_MPa'];zero=result['zero_load_columns_MPa']
                rows.append(dict(model=label,ring=ring,nominal_shear_MPa=stress,
                    gamma_y=result['gamma'][0],gamma_z=result['gamma'][1],
                    actual_xy_MPa=P[0,0],actual_xz_MPa=P[0,1],xz_error_MPa=result['sigma_xz_error_MPa'],
                    actual_yz_MPa=P[1,1],actual_zz_MPa=P[2,1],zero_load_zz_MPa=zero[2,1],
                    normal_reaction_increment_MPa=P[2,1]-zero[2,1],
                    energy_change_eV_atom=result['actual_energy_change_eV_atom']))
    write_csv(args.out/'actual_farfield_tractions.csv',rows)
    save_json(args.out/'completion.json',dict(completed=True,parameter_sha256=bindings,
        actual_site_energy_differentiated=True,volume_preserving_affine_shear=True,
        exact_nonlinear_stress_control=False,load_input_altered=False,
        no_statistical_area=True,physical_Hz=False,
        csv_sha256=hashlib.sha256((args.out/'actual_farfield_tractions.csv').read_bytes()).hexdigest()))


if __name__=='__main__':main()
