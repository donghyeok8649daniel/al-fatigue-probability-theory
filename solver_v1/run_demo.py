from pathlib import Path
import csv, json
import numpy as np
import matplotlib.pyplot as plt
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.solver import LoadParams, SolverParams, run_ensemble

OUT = Path("solver_v1/output")
OUT.mkdir(parents=True, exist_ok=True)

model_p = ModelParams(kT=0.009, mobility_s=0.15, chi_axial_projection=0.40)
model = TwoRowLJ(model_p)
model._build_opening_table()
solver_p = SolverParams(n_trajectories=32, dt=0.02, seed=17, first_passage_stride=5)

summary = []
runs = {}
for fmax in (2.5, 3.2, 3.4, 3.6):
    out = run_ensemble(model_p, LoadParams(force_max=fmax, cycles=10), solver_p, model=model)
    runs[fmax] = out
    last = out["time"] >= out["time"][-1] - 10.0
    area = abs(float(np.trapz(out["force"][last], out["strain"][last])))
    fp = 1.0 - float(out["survival"][-1])
    plastic = float(np.nanmax(out["plastic_well_activity"]))
    barrier_min = float(np.nanmin(out["opening_barrier"]))
    summary.append({
        "force_max": fmax,
        "last_cycle_loop_area": area,
        "first_passage_fraction": fp,
        "max_mean_abs_well_index": plastic,
        "minimum_mean_opening_barrier": barrier_min,
    })

with (OUT/"summary.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=list(summary[0]))
    w.writeheader(); w.writerows(summary)
(OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")

fig = plt.figure(figsize=(7,5))
for fmax in (2.5, 3.4):
    out = runs[fmax]
    last = out["time"] >= out["time"][-1] - 10.0
    plt.plot(out["strain"][last], out["force"][last], label=f"Fmax={fmax}")
plt.xlabel("dimensionless axial strain")
plt.ylabel("dimensionless tensile force")
plt.title("Theory Core v1: cyclic hysteresis")
plt.legend(); plt.tight_layout()
plt.savefig(OUT/"hysteresis.png", dpi=170)
plt.close(fig)

fig = plt.figure(figsize=(7,5))
for fmax in (3.2, 3.4, 3.6):
    out = runs[fmax]
    plt.plot(out["time"]/10.0, 1.0-out["survival"], label=f"Fmax={fmax}")
plt.xlabel("cycles")
plt.ylabel("first-passage fraction")
plt.title("Theory Core v1: opening first passage")
plt.legend(); plt.tight_layout()
plt.savefig(OUT/"first_passage.png", dpi=170)
plt.close(fig)

print(json.dumps(summary, indent=2))
