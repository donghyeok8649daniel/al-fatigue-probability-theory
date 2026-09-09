"""Like-for-like stress/relaxation audit, never a material refit or PDE change."""
from pathlib import Path
import time
import numpy as np

from .fcc111_geometry import DIRECT_110, SHOCKLEY_112, registry_path
from .reference_eam_targets import MishinRigidFCCReference
from .run_low_stress_cyclic_diagnostic import build_surface, write_csv, save_json, ROOT
from .stress_scale_audit import normal_relaxed_registry
from .interface_static_scenarios import InterfaceUnits


OUT = ROOT / "results/kinetics_loading_audit"


def main():
    started = time.perf_counter(); summary=[]; comparisons=[]
    source = MishinRigidFCCReference()
    for path in (DIRECT_110, SHOCKLEY_112):
        surface, units, metadata = build_surface(path)
        b_angstrom = units.length_scale_m/1e-10
        def source_packed(a,s):
            v = source.interface_derivatives(a*b_angstrom,s*b_angstrom,path_id=path).copy()
            v[1:3] *= b_angstrom
            v[3:] *= b_angstrom**2
            return v
        for name, evaluate, h, period in (
            ("unchanged_analytic_research_candidate", surface.packed, surface.h, surface.period),
            ("Mishin_reference_same_interface", source_packed, source.h/b_angstrom,
             registry_path(path).period_over_b*source.geometry.b/b_angstrom),
        ):
            eval_units = units if name.startswith("unchanged") else InterfaceUnits(
                units.length_scale_m, source.geometry.atomic_cell_area*1e-20)
            for samples in (33,65):
                print(f"{path} {name} relaxed branch samples={samples}",flush=True)
                rows, peak = normal_relaxed_registry(evaluate,h=h,period=period,units=eval_units,samples=samples)
                for row in rows:
                    comparisons.append(dict(model=name,path=path,samples=samples,**row))
                summary.append(dict(model=name,path=path,samples=samples,
                    fixed_gap_sampled_peak_mpa=max(r["fixed_shear_mpa"] for r in rows),
                    normal_relaxed_first_peak_mpa=peak["shear_mpa"] if peak else None,
                    spinodal_s_reduced=peak["s_reduced"] if peak else None,
                    spinodal_a_reduced=peak["a_reduced"] if peak else None,
                    force_residual=max(r["normal_force_residual"] for r in rows),
                    schur_residual=peak["schur_curvature"] if peak else None,
                    status=peak["status"] if peak else "not found in sampled branch"))
                print(summary[-1],flush=True)
                write_csv(OUT/"normal_relaxation_curves.csv",comparisons)
                write_csv(OUT/"stress_scale_summary.csv",summary)
    check=[]
    for stress in (4.,15.,25.,50.,100.,2000.):
        force=float(units.traction_mpa_to_force(stress))
        independent=stress*1e6*units.atomic_cell_area_m2*units.length_scale_m/1.602176634e-19
        check.append(dict(traction_MPa=stress,generalized_force_eV=force,
            independent_SI_work_eV=independent,error_eV=force-independent,
            recovered_MPa=float(units.force_to_traction_mpa(force))))
    write_csv(OUT/"unit_audit.csv",check)
    save_json(OUT/"stress_audit_status.json",dict(elapsed_seconds=time.perf_counter()-started,
        parameter_metadata=metadata,source_sha256=source.sha256,refit=False,
        physical_time_available=False,macroscopic_yield_prediction=False,
        remaining="scalar path and rigid half-crystals; no vector/full atomic relaxation or nucleation"))


if __name__ == "__main__":
    main()
