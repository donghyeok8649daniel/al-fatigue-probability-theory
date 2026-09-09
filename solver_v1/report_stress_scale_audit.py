"""Finish first-branch comparisons; never compare a later peak to first yield.

Consumes the executed normal-relaxation curves, then verifies the first fixed-
gap maximum with the actual analytic Hessian on both sampling resolutions.
"""
import csv
from functools import lru_cache
import json

import numpy as np
from scipy.optimize import brentq

from .fcc111_geometry import DIRECT_110, SHOCKLEY_112
from .interface_static_scenarios import InterfaceUnits
from .reference_eam_targets import MishinRigidFCCReference
from .run_low_stress_cyclic_diagnostic import build_surface, write_csv, save_json
from .run_stress_scale_audit import OUT, ROOT


def main():
    with (OUT/"normal_relaxation_curves.csv").open(encoding="utf-8") as stream:
        curves=list(csv.DictReader(stream))
    with (OUT/"stress_scale_summary.csv").open(encoding="utf-8") as stream:
        summaries=list(csv.DictReader(stream))
    source=MishinRigidFCCReference(); rows=[]
    for path in (DIRECT_110,SHOCKLEY_112):
        surface,units,metadata=build_surface(path)
        b=units.length_scale_m/1e-10
        for name in ("unchanged_analytic_research_candidate","Mishin_reference_same_interface"):
            def evaluate(s):
                if name.startswith("unchanged"):
                    return surface.packed(surface.h,s)
                v=source.interface_derivatives(source.h,b*s,path_id=path).copy()
                v[1:3]*=b; v[3:]*=b*b
                return v
            evaluate=lru_cache(maxsize=512)(evaluate)
            u=units if name.startswith("unchanged") else InterfaceUnits(units.length_scale_m,
                source.geometry.atomic_cell_area*1e-20)
            for count in (33,65):
                selected=[r for r in curves if r["model"]==name and r["path"]==path
                          and int(r["samples"])==count]
                peak=None
                for i in range(1,len(selected)-1):
                    left,current,right=selected[i-1:i+2]
                    if (float(current["fixed_shear_mpa"])>float(left["fixed_shear_mpa"]) and
                            float(current["fixed_shear_mpa"])>float(right["fixed_shear_mpa"])):
                        sl,sr=float(left["s_reduced"]),float(right["s_reduced"])
                        if evaluate(sl)[5]<=0 or evaluate(sr)[5]>=0:
                            raise ValueError("sampled first fixed maximum lacks a stable curvature bracket")
                        # A tabulated-source interpolation can have several
                        # curvature zeros within ONE coarse peak bracket.
                        # brentq on that whole interval may return a later root.
                        # Subdivide, then take the first + -> - crossing.
                        refined=np.linspace(sl,sr,33)
                        curvature=[evaluate(float(x))[5] for x in refined]
                        brackets=[(x,y) for x,y,cx,cy in zip(refined[:-1],refined[1:],
                            curvature[:-1],curvature[1:]) if cx>0 and cy<0]
                        if not brackets:
                            raise ValueError("first curvature crossing unresolved inside peak bracket")
                        s=brentq(lambda x:evaluate(float(x))[5],*brackets[0],xtol=2e-12)
                        v=evaluate(s)
                        peak=dict(first_fixed_peak_s_reduced=s,
                            first_fixed_peak_shear_MPa=float(u.force_to_traction_mpa(v[2])),
                            fixed_curvature_residual=float(v[5]))
                        break
                if peak is None:
                    raise ValueError("first fixed-gap peak not resolved")
                summary=next(r for r in summaries if r["model"]==name and r["path"]==path
                             and int(r["samples"])==count)
                relaxed=float(summary["normal_relaxed_first_peak_mpa"])
                rows.append(dict(model=name,path=path,samples=count,**peak,
                    first_normal_relaxed_peak_MPa=relaxed,
                    relaxation_reduction_fraction=1-relaxed/peak["first_fixed_peak_shear_MPa"],
                    interpretation="ideal scalar-path interface instability; NOT yield/fatigue strength"))
                print(rows[-1],flush=True)
    write_csv(OUT/"first_branch_comparison.csv",rows)
    # Source elastic constants are an independently documented comparator, not
    # numbers substituted into the candidate to improve its response.
    elastic=json.loads((ROOT/"results/fcc111_active_interface/nonlocal_v5/material_metadata.json").read_text())
    elastic_rows=[]
    for key,target in (("C11_GPa",114.),("C12_GPa",62.),("C44_GPa",32.)):
        actual=elastic["elastic_constants_GPa"][key]
        elastic_rows.append(dict(quantity=key,candidate_GPa=actual,source_GPa=target,
            relative_error=(actual-target)/target,source="Mishin 1999; repository Al targets, 0 K"))
    write_csv(OUT/"elastic_fit_errors.csv",elastic_rows)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    for ax,path in zip(axes,(DIRECT_110,SHOCKLEY_112)):
        for name,label in (("unchanged_analytic_research_candidate","Analytic candidate"),
                           ("Mishin_reference_same_interface","Mishin reference")):
            r=next(r for r in rows if r["model"]==name and r["path"]==path and r["samples"]==65)
            # Stop at the FIRST instability; later parts of the registry curve
            # must not masquerade as the onset of the pristine branch.
            data=[v for v in curves if v["model"]==name and v["path"]==path and int(v["samples"])==65]
            summary=next(v for v in summaries if v["model"]==name and v["path"]==path and int(v["samples"])==65)
            for field,limit,style in (("fixed_shear_mpa",r["first_fixed_peak_s_reduced"],"--"),
                    ("relaxed_shear_mpa",float(summary["spinodal_s_reduced"]),"-")):
                d=[v for v in data if float(v["s_reduced"])<=limit]
                ax.plot([float(v["s_reduced"]) for v in d],[float(v[field]) for v in d],style,
                    label=label+(" fixed gap" if style=="--" else " normal relaxed"))
        ax.set(title=path,xlabel="Registry displacement / L0",ylabel="Shear traction [MPa]")
        ax.legend(fontsize=7)
    fig.suptitle("Ideal pristine-interface branches, not macroscopic yield")
    fig.savefig(OUT/"first_branch_traction.png",dpi=150)
    plt.close(fig)
    save_json(OUT/"interpretation.json",dict(
        first_branch_refinement_performed=True,curvature_bracket_subdivisions=32,
        coarse_source_peak_bracket_was_multiroot=True,unit_scale_error_detected=False,
        normal_relaxation_supported=True,material_fit_accepted=False,
        missing_macroscopic_mechanisms=["finite defect core / nucleation geometry", "vector and full atomic relaxation",
            "validated spatial ensemble coupling", "collective-coordinate kinetic data"],
        no_parameter_or_mobility_tuning=True,production_energy_unchanged=True,
        physical_seconds_and_hz_available_for_Al=False))


if __name__=="__main__":
    main()
