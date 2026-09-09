"""Independent candidate validation after actual fitting, not re-labeling."""
import json
import time

import numpy as np

from .range_resolved_material import build_range_surface,RangeBulkValidation
from .run_low_stress_cyclic_diagnostic import ROOT,write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .vector_material_calibration import UNITS,IDEAL_H
from .vector_registry_audit import stationary_state


OUT=ROOT/'results/fcc111_active_interface/range_core_v12/material'


def main():
    began=time.perf_counter(); source,_,states=source_and_targets()
    study=json.loads((OUT/'calibration.json').read_bytes()); summary=[]; roots=[]; waves=[]; stresses=[]
    for name,fit in study['best'].items():
        model=build_range_surface(fit['scalar_decay'],fit['odd_decay'],fit['quadrupole_decay'],fit['coefficients'])
        info=dict(candidate=name,**dict(zip(['C11_GPa','C12_GPa','C44_GPa'],fit['cubic_GPa'])),
            fit_loss=fit['squared_loss'],heldout_normalized_rms=fit['heldout_normalized_rms'],adopted=False)
        for label,guess,index in [('perfect',[IDEAL_H,0,0],0),('fault',states['fault'],0),('saddle',states['saddle'],1)]:
            r=stationary_state(model,guess,expected_index=index); value=r['evaluation']
            row=dict(candidate=name,state=label,valid=r['valid'],a=r['q'][0],x=r['q'][1],y=r['q'][2],
                energy_J_m2=float(UNITS.energy_to_surface(value.energy)),
                force_residual_eV_L0=r['force_residual'],hessian_eigenvalues=str(r['eigenvalues'].tolist()))
            roots.append(row); info[label+'_valid']=r['valid']; info[label+'_energy_J_m2']=row['energy_J_m2']
        info['separation_40h_J_m2']=float(UNITS.energy_to_surface(model.evaluate((40*IDEAL_H,0,0)).energy))
        info['source_separation_40h_J_m2']=float(UNITS.energy_to_surface(source.evaluate((40*IDEAL_H,0,0)).energy))
        matrices={}
        for radius in (6.,8.,10.):
            validator=RangeBulkValidation(model,cutoff=radius); selected={}
            for path,end in [('Gamma_X',[1,0,0]),('Gamma_L',[.5,.5,.5]),('Gamma_K',[.75,.75,0])]:
                for fraction in (.02,.1,.25,.5,.75,1.):
                    key=path,fraction
                    q=validator.crystallographic_wavevector(np.array(end)*fraction)
                    result=validator.evaluate(q); selected[key]=result['matrix']
                    difference=float(np.max(abs(result['matrix']-matrices[key]))) if key in matrices else None
                    waves.append(dict(candidate=name,radius_over_L0=radius,path=path,fraction=fraction,
                        minimum_eigenvalue_eV_L0sq=result['eigenvalues'][0],
                        matrix_change_from_previous_radius=difference))
            matrices=selected
            print(f'{name}: finite-q radius {radius}, min={min(np.linalg.eigvalsh(m)[0] for m in matrices.values()):.6g}',flush=True)
        info['minimum_tested_finite_q_eigenvalue']=min(np.linalg.eigvalsh(m)[0] for m in matrices.values())
        q=np.array([IDEAL_H,0.,0.])
        for normal,shear in [(0,0),(10,5),(25,12.5),(50,25),(25,12.5),(0,0)]:
            force=UNITS.traction_mpa_to_force([normal,shear,0])
            r=stationary_state(model,q,force=force,expected_index=0)
            if r['valid']:
                q=r['q']
            stresses.append(dict(candidate=name,normal_MPa=normal,shear_MPa=shear,
                stable=r['valid'],a=r['q'][0],x=r['q'][1],y=r['q'][2],
                force_residual_eV_L0=r['force_residual'],minimum_hessian_eigenvalue=min(r['eigenvalues']),
                static_not_dynamic=True,residual_plasticity_validated=False))
        summary.append(info)
    write_csv(OUT/'validated_summary.csv',summary); write_csv(OUT/'stationary_states.csv',roots)
    write_csv(OUT/'finite_q_validation.csv',waves); write_csv(OUT/'static_load_unload.csv',stresses)
    save_json(OUT/'validation_status.json',dict(completed=True,elapsed_seconds=time.perf_counter()-began,
        material_accepted=False,finite_q_sampling_is_not_global_stability_proof=True,
        static_load_unload_is_not_residual_plasticity=True,
        decision='improved fit requires held-out/elastic compatibility; no automatic production promotion'))


if __name__=='__main__':
    main()
