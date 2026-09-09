"""Reproducible low-MPa probes. Nothing here registers a production model.

`grid` samples the unchanged analytic research surface. `cycles` executes a
noncanonical reflecting-box SG hypothesis test, NOT calibrated crack dynamics.
`static` compares MPa load/unload paths and energy/work scales. All fit parameters
are read from the previously saved static candidate without optimization.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .angular_environment_reference import AngularInterfaceResearchSurface
from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_geometry import DIRECT_110, SHOCKLEY_112
from .interface_static_scenarios import InterfaceUnits, continue_static_branch
from .joint_fcc_interface_calibration import JointFit, joint_equilibrium
from .low_stress_cyclic_diagnostic import (
    KB_EV_K, CyclicTractions, aligned_grid, periodic_energy_samples, run_reflecting_cycles,
)


ROOT = Path(__file__).resolve().parents[1]
FIT = ROOT/"results/fcc111_active_interface/matched_v3/monotone_opening/fitted_candidates.json"
OUT = ROOT/"results/fcc111_active_interface/low_stress_v4"


def write_csv(path, rows):
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, list(rows[0]) if rows else [])
        writer.writeheader(); writer.writerows(rows)


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def build_surface(path_id=DIRECT_110, tolerance=2e-10):
    raw = FIT.read_bytes()
    fit = json.loads(raw)["angular_monotone_opening"]
    c = dict(zip(fit["coefficient_order"], fit["coefficients"]))
    scalar = np.array([c[k] for k in ("u", "v", "A", "B", "C")])
    qfit = JointFit("saved_research_only", fit["scalar_decay"], scalar, fit["squared_loss"], None, None)
    bulk, stretch = joint_equilibrium(qfit)
    face = FCC111ActiveInterface(bulk, path_id=path_id, tolerance=tolerance/10)
    surface = AngularInterfaceResearchSurface(face, c["D3"], angular_decay=fit["angular_decay"],
        vector_amplitude_ev=c["D1"], quadrupole_amplitude_ev=c["D2"], tolerance=tolerance)
    length = 4.05/math.sqrt(2)*1e-10
    units = InterfaceUnits(length, bulk.geometry.atomic_cell_area*length**2)
    metadata = dict(model="angular_monotone_opening", path=path_id,
        parameter_source=FIT.relative_to(ROOT).as_posix(), parameter_sha256=hashlib.sha256(raw).hexdigest(),
        no_refit=True, research_only=True, length_scale_m=length,
        atomic_cell_area_m2=units.atomic_cell_area_m2, energy_unit="eV/atomic interface cell",
        h_reduced=surface.h, registry_period_reduced=surface.period, lattice_stretch=stretch,
        tolerance=tolerance, physical_hz=False, material_calibration_accepted=False)
    return surface, units, metadata


def grid_file(args):
    return OUT/f"grid_a{args.n_a}_s{args.cells_per_well}_upper{args.upper:g}.npz"


def prepare_grid(args):
    started = time.perf_counter()
    surface, units, metadata = build_surface()
    grid = aligned_grid(surface.h, surface.period, n_a=args.n_a,
        cells_per_well=args.cells_per_well, upper_over_h=args.upper)
    # Independently check translation/reflection before reusing exact symmetries.
    errors = [abs(surface.packed(a, s)[0]-surface.packed(a, s+surface.period)[0])
              for a, s in ((.93*surface.h, .173*surface.period), (1.2*surface.h, -.287*surface.period))]
    reflection = [abs(surface.packed(a, s)[0]-surface.packed(a, -s)[0])
        for a, s in ((.72*surface.h, .03), (.93*surface.h, .173), (1.2*surface.h, .287),
                     (args.upper*.99*surface.h, .49))]
    if max(errors+reflection) > 100*metadata["tolerance"]:
        raise ValueError("direct110 symmetry/translation verification failed")
    def progress(i, total):
        if i == 1 or i % 5 == 0 or i == total:
            print(f"analytic energy rows {i}/{total}; {time.perf_counter()-started:.1f}s", flush=True)
    energy = periodic_energy_samples(surface.packed, grid, surface.period, progress=progress,
                                     mirror_registry=True)
    target = grid_file(args)
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(target, a=grid.a, s=grid.s, da=grid.da, ds=grid.ds, energy=energy,
        h=surface.h, period=surface.period)
    metadata.update(grid_a=args.n_a, cells_per_well=args.cells_per_well, upper_over_h=args.upper,
        elapsed_seconds=time.perf_counter()-started, periodicity_abs_error_ev=max(errors),
        direct110_reflection_abs_error_ev=max(reflection),
        energy_min_ev=float(np.min(energy)), energy_max_ev=float(np.max(energy)),
        file_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    save_json(target.with_suffix(".json"), metadata)
    print(json.dumps(metadata, indent=2), flush=True)


def low_stress_static():
    started = time.perf_counter(); rows = []; scale_rows = []
    for path in (DIRECT_110, SHOCKLEY_112):
        surface, units, metadata = build_surface(path)
        # Repeat a complete quasistatic waveform to expose path retracing only.
        for amplitude in (10., 20., 30., 50.):
            print(f"static {path} axial amplitude {amplitude:g} MPa", flush=True)
            phase = np.linspace(0., 4*np.pi, 65)
            values = amplitude*.5/1000*np.sin(phase)
            loads = np.column_stack((values, values))
            history = continue_static_branch(surface.packed, h=surface.h, period=surface.period,
                units=units, loads_gpa=loads, max_traction_step=.002)
            rows.extend(dict(path=path, nominal_axial_amplitude_mpa=amplitude,
                phase_radians=float(p), **r) for p, r in zip(phase, history))
        value = surface.packed(surface.h, 0.)
        H = np.array([[value[3], value[4]], [value[4], value[5]]])
        rates = np.linalg.eigvalsh(np.diag([1., np.sqrt(.05)])@H@np.diag([1., np.sqrt(.05)]))
        for shear in (4., 5., 10., 15., 25., 50.):
            work = float(units.traction_to_force(shear/1000)*surface.period)
            scale_rows.append(dict(path=path, shear_amplitude_mpa=shear,
                work_over_full_registry_period_ev=work, kT_ev=KB_EV_K*293.15,
                work_over_kT=work/(KB_EV_K*293.15), tau_fast_model=float(1/rates[-1]),
                tau_slow_model=float(1/rates[0]), physical_seconds_available=False))
    write_csv(OUT/"quasistatic_mpa_cycles.csv", rows)
    write_csv(OUT/"work_and_relaxation_scales.csv", scale_rows)
    save_json(OUT/"static_execution.json", dict(elapsed_seconds=time.perf_counter()-started,
        protocol="two quasistatic cycles, R=-1, declared 45-degree example axis",
        no_time_evolution=True, no_residual_plasticity_claim=True))


def run_cycles(args):
    from .probability_pde_2d import Grid2D
    target = grid_file(args)
    with np.load(target) as data:
        grid = Grid2D(data["a"], data["s"], float(data["da"]), float(data["ds"]))
        energy = data["energy"].copy(); h = float(data["h"]); period = float(data["period"])
    metadata = json.loads(target.with_suffix(".json").read_text())
    if hashlib.sha256(target.read_bytes()).hexdigest() != metadata["file_sha256"]:
        raise ValueError("energy-grid checksum mismatch; do not evolve an altered table")
    if hashlib.sha256(FIT.read_bytes()).hexdigest() != metadata["parameter_sha256"]:
        raise ValueError("energy grid parameters do not match current saved fit")
    units = InterfaceUnits(metadata["length_scale_m"], metadata["atomic_cell_area_m2"])
    if args.wells != 3:
        # Add equivalent intact wells, not a larger coherent activation energy.
        grid = aligned_grid(h, period, n_a=args.n_a, cells_per_well=args.cells_per_well,
                            upper_over_h=args.upper, wells=args.wells)
        energy = np.tile(energy[:, :args.cells_per_well], (1, args.wells))
    cases = dict(zero=CyclicTractions(period_model=args.period),
        shear4=CyclicTractions(shear_amplitude_mpa=4., period_model=args.period),
        shear15=CyclicTractions(shear_amplitude_mpa=15., period_model=args.period),
        normal15=CyclicTractions(normal_amplitude_mpa=15., period_model=args.period),
        reversed30=CyclicTractions(normal_amplitude_mpa=15., shear_amplitude_mpa=15.,
                                   phase_s=np.pi, period_model=args.period),
        axial10=CyclicTractions.axial_45(10., period_model=args.period),
        axial20=CyclicTractions.axial_45(20., period_model=args.period),
        axial30=CyclicTractions.axial_45(30., period_model=args.period),
        axial50=CyclicTractions.axial_45(50., period_model=args.period),
        nonproportional30=CyclicTractions(normal_amplitude_mpa=15., shear_amplitude_mpa=15.,
                                         phase_s=np.pi/2, period_model=args.period))
    names = list(cases) if args.case == "all" else [args.case]
    for name in names:
        started = time.perf_counter()
        tag = f"{name}_p{args.period:g}_a{args.n_a}_s{args.cells_per_well}_w{args.wells}_dt{args.steps}_upper{args.upper:g}"
        if args.cycles != 16 or args.hold != 8:
            tag += f"_c{args.cycles}_h{args.hold}"
        if args.integrator != "implicit":
            tag += "_"+args.integrator
        if args.global_gibbs:
            tag += "_globalgibbs"
        print(f"RUN {tag}", flush=True)
        def progress(i, total):
            if i % 4 == 0 or i == total:
                print(f"  actual SG cycles+hold {i}/{total}; {time.perf_counter()-started:.1f}s", flush=True)
        result = run_reflecting_cycles(energy, grid, h=h, period=period, units=units,
            load=cases[name], cycles=args.cycles, steps_per_cycle=args.steps,
            hold_cycles=args.hold, central_only=not args.global_gibbs, progress=progress,
            integrator=args.integrator)
        destination = OUT/"cycles"/tag
        # Store full phases for the first/last driven and first/last hold cycle,
        # plus quarter-cycle samples elsewhere. ALL steps enter the audit and
        # per-cycle statistics; this only limits serialized plot data size.
        dense_cycles = {0, 1, args.cycles-1, args.cycles, args.cycles+1, args.cycles+args.hold}
        history = [r for i, r in enumerate(result["history"])
                   if r["cycle"] in dense_cycles or i % max(1, args.steps//4) == 0]
        write_csv(destination/"history.csv.gz", history)
        write_csv(destination/"per_cycle.csv", result["cycles"])
        # Preserve full precision at the states needed for the matched-control
        # and last-cycle comparisons. Every cycle's moments/balances are saved
        # separately. This is output sampling, never an evolution approximation.
        indices = sorted({0, 1, args.cycles-1, args.cycles,
                          min(args.cycles+1, args.cycles+args.hold), args.cycles+args.hold})
        np.savez_compressed(destination/"snapshots.npz", cycle_indices=indices,
                            same_phase=result["same_phase_mass"][indices])
        audit = dict(**result["audit"], parameter_sha256=metadata["parameter_sha256"],
            elapsed_seconds=time.perf_counter()-started, case=name,
            load=vars(cases[name]), grid_sha256=metadata["file_sha256"],
            history_sampling="dense first/last driven, first/last hold; quarter cycles elsewhere",
            snapshot_cycle_indices=indices,
            history_rows=len(history), integrated_steps=len(result["history"])-1)
        save_json(destination/"audit.json", audit)
        print(json.dumps(dict(case=tag, elapsed_seconds=audit["elapsed_seconds"],
            final_mean_n=result["cycles"][-1]["mean_well_index"],
            final_outside=result["cycles"][-1]["outside_central_mass"],
            mass_residual=audit["max_mass_residual"],
            well_balance=audit["max_well_balance_residual"], opening_probability=None)), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("static", "grid", "cycles"), required=True)
    parser.add_argument("--n-a", type=int, default=31)
    parser.add_argument("--cells-per-well", type=int, default=25)
    parser.add_argument("--upper", type=float, default=2.)
    parser.add_argument("--wells", type=int, default=3)
    parser.add_argument("--period", type=float, default=4.)
    parser.add_argument("--steps", type=int, default=128)
    parser.add_argument("--cycles", type=int, default=16)
    parser.add_argument("--hold", type=int, default=8)
    parser.add_argument("--case", default="all", choices=("all", "zero", "shear4", "shear15", "normal15", "reversed30", "axial10", "axial20", "axial30", "axial50", "nonproportional30"))
    parser.add_argument("--global-gibbs", action="store_true")
    parser.add_argument("--integrator", choices=("implicit", "explicit"), default="implicit")
    args = parser.parse_args()
    if args.phase == "grid": prepare_grid(args)
    elif args.phase == "static": low_stress_static()
    else: run_cycles(args)


if __name__ == "__main__":
    main()
