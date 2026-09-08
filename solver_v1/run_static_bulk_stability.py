"""Run deterministic finite-q stability scenarios on research candidates."""
import numpy as np

from .run_matched_interface_study import ROOT,surfaces_from_results,save_json
from .run_full_fcc_calibration_audit import write_csv
from .static_bulk_stability import StaticBulkHessian


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extra-monotone-radius",type=float,
                        help="additional independent direct-neighbor refinement for final candidate")
    args=parser.parse_args()
    if args.extra_monotone_radius is not None and (not np.isfinite(args.extra_monotone_radius) or args.extra_monotone_radius<=12):
        raise ValueError("extra validation radius must be finite and exceed 12 L0")
    rows=[]; summaries=[]
    for folder in (ROOT,ROOT/"odd_moment",ROOT/"even_moment",ROOT/"monotone_opening"):
        if not (folder/"fitted_candidates.json").exists():
            continue
        for surface in surfaces_from_results(folder):
            if surface["bulk"] is None or surface["path"]!="direct_110":
                continue
            name=surface["name"]; D1=D2=D3=0.; decay=None
            if name!="linear":
                bound=surface["evaluate"].__self__
                D3=bound.amplitude_ev; D1=bound.vector_amplitude_ev
                D2=bound.quadrupole_amplitude_ev
                decay=bound.angular.kappa
            minima=[]
            radii=(5.,8.,12.)
            if name=="angular_monotone_opening" and args.extra_monotone_radius is not None:
                radii=(*radii,args.extra_monotone_radius)
            for cutoff in radii:
                evaluator=StaticBulkHessian(surface["bulk"],cutoff=cutoff,D3=D3,D1=D1,D2=D2,angular_decay=decay)
                assert np.max(abs(evaluator.evaluate(np.zeros(3))["matrix"]))==0
                values=[]
                for label,end in (("GX",[0,1,0]),("GL",[.5,.5,.5]),("GK",[.75,.75,0])):
                    for fraction in np.linspace(.05,1.,20):
                        q=evaluator.crystallographic_wavevector(np.array(end)*fraction)
                        value=evaluator.evaluate(q)
                        values.append(float(value["eigenvalues"][0]))
                        rows.append(dict(candidate=name,cutoff_L0=cutoff,neighbors=len(evaluator.R),
                            path=label,path_fraction=fraction,min_eigenvalue_eV_L0sq=values[-1],
                            middle_eigenvalue_eV_L0sq=float(value["eigenvalues"][1]),
                            max_eigenvalue_eV_L0sq=float(value["eigenvalues"][2]),
                            symmetry_error=float(np.max(abs(value["matrix"]-value["matrix"].T))),
                            status="static second variation; no mass, no frequency units"))
                minima.append(min(values))
            summaries.append(dict(candidate=name,radii_L0=radii,minima_by_radius=minima,
                sampled_nonnegative=bool(minima[-1]>=0),
                conclusion="sampled directions only, not proof for the entire Brillouin zone",
                physical_Hz_available=False))
            write_csv(ROOT/"finite_q_static_stability.csv",rows)
            print(name,minima,flush=True)
    save_json(ROOT/"finite_q_stability_summary.json",summaries)


if __name__=="__main__":
    main()
