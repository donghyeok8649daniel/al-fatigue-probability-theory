"""Sample the actual two-coordinate research landscape, without PDE promotion."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .run_matched_interface_study import ROOT,surfaces_from_results,save_json
from .run_full_fcc_calibration_audit import write_csv
from .fcc111_geometry import DIRECT_110,SHOCKLEY_112
from .interface_static_scenarios import hessian_from_packed


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family",choices=("odd_moment","even_moment","monotone_opening"),default="odd_moment")
    folder=ROOT/parser.parse_args().family
    surfaces=surfaces_from_results(folder); rows=[]; summaries=[]
    for surface in surfaces:
        h=surface["h"]; period=surface["period"]; units=surface["units"]
        path=surface["path"]; model=surface["name"]
        fractions=np.linspace(0.,1. if path==DIRECT_110 else 1/3,21)
        openings=np.r_[np.linspace(.9,1.5,13),1.75,2.]
        energy=np.empty((len(openings),len(fractions)))
        for i,ratio in enumerate(openings):
            for j,fraction in enumerate(fractions):
                v=surface["evaluate"](float(ratio*h),float(fraction*period))
                energy[i,j]=units.energy_to_surface(v[0])
                rows.append(dict(candidate=model,path=path,a_over_h=ratio,s_over_period=fraction,
                    energy_ev_cell=v[0],energy_J_m2=energy[i,j],
                    normal_traction_GPa=float(units.force_to_traction(v[1])),
                    shear_traction_GPa=float(units.force_to_traction(v[2])),
                    min_hessian_eigenvalue=float(np.linalg.eigvalsh(hessian_from_packed(v))[0]),
                    coordinate_length_scale_m=units.length_scale_m))
        fig,axis=plt.subplots(figsize=(6,4))
        plot=axis.contourf(fractions,openings,energy,levels=24)
        fig.colorbar(plot,ax=axis,label="W_int [J/m2]")
        axis.set_xlabel("s / full registry repeat"); axis.set_ylabel("interface spacing / h")
        axis.set_title(f"{model}\n{path}; STATIC research only",fontsize=10)
        fig.tight_layout(); fig.savefig(folder/(model+"_"+path+"_surface.png"),dpi=135); plt.close(fig)
        summaries.append(dict(candidate=model,path=path,sampled_minimum_J_m2=float(np.min(energy)),
            sampled_maximum_J_m2=float(np.max(energy)),global_minimum_certified=False,
            interpretation="saddle regions may have negative Hessian; this is not instability of every intact state"))
        write_csv(folder/"coupled_energy_grid.csv",rows)
        print(model,path,"sampled minimum",float(np.min(energy)),flush=True)
    save_json(folder/"coupled_grid_summary.json",summaries)


if __name__=="__main__":
    main()
