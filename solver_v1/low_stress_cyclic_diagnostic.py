"""Noncanonical, reflecting-domain thermalization falsification experiment.

This is NOT an active-interface crack solver or a material calibration. It
tests the *hypothesis* that a rigid infinite-interface energy per atomic cell
can be assigned a collective-coordinate Gibbs measure. That physical
normalization is unvalidated. No opening sink, crack probability, seconds,
or specimen area is defined here. Production energy/PDE registries are untouched.

The spatial generator is imported, unchanged, from the verified 2D SG solver.
Only an explicitly supplied energy array and a prescribed traction history
are used. The normal domain and outer registry boundaries are reflecting.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Callable

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from .interface_static_scenarios import InterfaceUnits
from .probability_pde_2d import Grid2D, _bernoulli, _sg_generator_2d


KB_EV_K = 8.617333262145e-5


@dataclass(frozen=True)
class CyclicTractions:
    """MPa are physical stress units; period is strictly MODEL time.

    The two components are local tractions, not interchangeable with nominal
    axial stress. phase_s radians allows a declared nonproportional probe.
    """
    normal_amplitude_mpa: float = 0.
    shear_amplitude_mpa: float = 0.
    normal_mean_mpa: float = 0.
    shear_mean_mpa: float = 0.
    phase_s: float = 0.
    period_model: float = 4.

    def __post_init__(self):
        values = np.array(list(vars(self).values()), dtype=float)
        if not np.all(np.isfinite(values)) or self.period_model <= 0:
            raise ValueError("finite tractions/phase and positive model period required")
        if min(self.normal_amplitude_mpa, self.shear_amplitude_mpa) < 0:
            raise ValueError("amplitudes must be nonnegative; use phase for sign")

    @classmethod
    def axial_45(cls, amplitude_mpa, **kwargs):
        """Explicit example e=(n+m)/sqrt(2); NOT a fitted specimen orientation."""
        return cls(normal_amplitude_mpa=amplitude_mpa/2,
                   shear_amplitude_mpa=amplitude_mpa/2, **kwargs)

    def phase_value_mpa(self, phase):
        return np.array([self.normal_mean_mpa+self.normal_amplitude_mpa*np.sin(phase),
                         self.shear_mean_mpa+self.shear_amplitude_mpa*np.sin(phase+self.phase_s)])


def aligned_grid(h, period, *, n_a=41, cells_per_well=33, wells=3,
                 lower_over_h=.70, upper_over_h=2.):
    """All registry-well faces are exactly represented; no boundary snapping."""
    if h <= 0 or period <= 0 or n_a < 5 or cells_per_well < 5:
        raise ValueError("positive geometry and at least five cells required")
    if wells < 3 or wells % 2 != 1 or cells_per_well % 2 != 1:
        raise ValueError("odd well count and odd cells per well required")
    if not 0 < lower_over_h < 1 < upper_over_h:
        raise ValueError("normal reflecting domain must contain pristine spacing")
    ae = np.linspace(lower_over_h*h, upper_over_h*h, n_a+1)
    se = np.linspace(-wells*period/2, wells*period/2, wells*cells_per_well+1)
    return Grid2D((ae[1:]+ae[:-1])/2, (se[1:]+se[:-1])/2,
                  float(ae[1]-ae[0]), float(se[1]-se[0]))


def periodic_energy_samples(evaluate: Callable, grid, period, *, progress=None, mirror_registry=False):
    """Sample the actual infinite-series energy, reuse proven translations.

    No interpolation, averaging, or finite-neighbor replacement. Optional exact
    direct110 reflection reuse requires independent verification by the caller.
    Repeated equal registry offsets are evaluated once per normal row. Returned
    table is the cell-centre FV representation, not a new material potential.
    """
    offsets = (grid.s+period/2) % period-period/2
    # Derive representatives from integer FV positions (no energy rounding).
    count = int(round(period/grid.ds))
    groups = np.arange(len(grid.s)) % count
    representatives = offsets[:count]
    output = np.empty((len(grid.a), len(grid.s)))
    for i, a in enumerate(grid.a):
        if mirror_registry:
            if count % 2 != 1:
                raise ValueError("odd periodic cell count required for reflection reuse")
            row = np.empty(count)
            for j in range(count//2, count):
                row[j] = float(evaluate(float(a), float(representatives[j]))[0])
                row[count-1-j] = row[j]
        else:
            row = np.array([float(evaluate(float(a), float(s))[0]) for s in representatives])
        output[i] = row[groups]
        if progress:
            progress(i+1, len(grid.a))
    return output


def registry_partition(grid, period):
    n = np.floor(grid.s/period+.5).astype(int)
    indices = np.flatnonzero(n[1:] != n[:-1])
    if np.any(np.diff(n)[indices] != 1):
        raise ValueError("all neighboring wells must be represented")
    faces = (grid.s[:-1]+grid.s[1:])[indices]/2
    err = np.max(np.abs(faces-(n[indices]+.5)*period))
    if err > 128*np.finfo(float).eps*max(1., np.max(np.abs(grid.s))):
        raise ValueError("registry-well boundaries do not align with grid faces")
    return n, np.unique(n), indices, float(err)


def gibbs_cell_mass(energy, kT, cell_wells, *, central_only=True):
    if kT <= 0 or np.any(~np.isfinite(energy)):
        raise ValueError("finite energies and positive thermal scale required")
    mass = np.exp(-(energy-np.min(energy))/kT)
    if central_only:
        mass[:, cell_wells != 0] = 0.
    return mass / np.sum(mass)


def interface_traffic(mass, energy, grid, kT, mobility_s, indices):
    """SG forward/backward *grid* traffic; sum is NOT continuum hop count.

    With cell mass p=P da ds, face rates D/ds^2 give integrated physical
    signed J. The gross sum diverges as 2D P_face/ds as ds -> 0: continuum
    Brownian recrossings cannot be identified with committed well transitions.
    """
    psi = (energy[:, 1:]-energy[:, :-1])/kT
    factor = mobility_s*kT/grid.ds**2
    plus = np.sum(factor*_bernoulli(psi)*mass[:, :-1], axis=0)[indices]
    minus = np.sum(factor*_bernoulli(-psi)*mass[:, 1:], axis=0)[indices]
    return plus, minus


def _moments(mass, grid, n, wells, h, period):
    total = float(np.sum(mass))
    populations = np.array([np.sum(mass[:, n == w]) for w in wells])
    mean_a = float(np.sum(mass*grid.a[:, None])/total)
    mean_s = float(np.sum(mass*grid.s[None, :])/total)
    mean_n = float(np.dot(wells, populations)/total)
    xi = grid.s-period*n
    mean_xi = float(np.sum(mass*xi[None, :])/total)
    return dict(raw_mass=total, mass_residual=total-1., mean_a=mean_a, mean_s=mean_s,
                normal_strain=(mean_a-h)/h, registry_shear=mean_s/h,
                intrawell_shear=mean_xi/h, registry_index_shear=period*mean_n/h,
                mean_well_index=mean_n, outside_central_mass=float(np.sum(populations[wells != 0])),
                outer_well_mass=float(populations[0]+populations[-1]),
                normal_edge_mass=float(np.sum(mass[[0,-1], :])),
                upper_normal_edge_mass=float(np.sum(mass[-1])),
                decomposition_residual=(mean_s-mean_xi-period*mean_n)/h,
                **{f"P_{w:+d}": float(p) for w, p in zip(wells, populations)})


def run_reflecting_cycles(base_energy, grid, *, h, period, units: InterfaceUnits,
                          load: CyclicTractions, cycles=12, steps_per_cycle=128,
                          hold_cycles=4, temperature_K=293.15, mobility_a=1.,
                          mobility_s=.05, central_only=True, progress=None,
                          integrator="implicit"):
    """Backward Euler + unchanged SG; no repair or mass renormalization.

    One periodic set of LU factors is reused exactly on later cycles. This
    changes cost only. Flux accounting uses p_(n+1), the same BE state as the
    solver. Hold is zero traction, not continued mean loading. No large negative
    values or numerical repairs are silently clipped.
    """
    kT = KB_EV_K*temperature_K
    if cycles < 1 or hold_cycles < 0 or steps_per_cycle < 8:
        raise ValueError("positive cycles, nonnegative hold, at least 8 steps/cycle")
    if min(temperature_K, mobility_a, mobility_s) <= 0:
        raise ValueError("positive temperature and mathematical mobilities required")
    if integrator not in ("implicit", "explicit"):
        raise ValueError("implicit or explicit SG diagnostic required")
    shape = (len(grid.a), len(grid.s))
    if np.shape(base_energy) != shape or not np.all(np.isfinite(base_energy)):
        raise ValueError("finite energy table with exact grid shape required")
    n, wells, indices, alignment = registry_partition(grid, period)
    # Only the existing energy-array SG interface is reused; no model registration.
    adapter = SimpleNamespace(p=SimpleNamespace(kT=kT, mobility_a=mobility_a, mobility_s=mobility_s))
    identity = sparse.eye(np.prod(shape), format="csc")
    dt = load.period_model/steps_per_cycle
    da_coord = grid.a[:, None]-h
    s_coord = grid.s[None, :]
    def operator(traction_mpa):
        forces = units.traction_to_force(np.asarray(traction_mpa)/1000.)
        energy = base_energy-forces[0]*da_coord-forces[1]*s_coord
        generator = _sg_generator_2d(energy, adapter, grid)
        if integrator == "implicit":
            solver = splu(identity-dt*generator)
        else:
            solver = generator
        return energy, solver
    mass = gibbs_cell_mass(base_energy, kT, n, central_only=central_only)
    initial = mass.copy()
    factors = []
    # Midpoint forcing gives a symmetric sample of the continuous sinusoid;
    # BE remains first order in time. No physical frequency is inferred.
    for j in range(steps_per_cycle):
        traction = load.phase_value_mpa(2*np.pi*(j+.5)/steps_per_cycle)
        factors.append((*operator(traction), traction))
    zero_energy, zero_lu = operator(np.zeros(2))
    signed = np.zeros(len(indices)); gross = signed.copy()
    maximum_balance = maximum_registry_balance = 0.
    minimum_mass = float(np.min(mass))
    history = []
    cycle_rows = []
    def row(t, cycle, segment, traction, plus, minus):
        moments = _moments(mass, grid, n, wells, h, period)
        return dict(time_model=t, cycle=cycle, segment=segment,
                    normal_traction_mpa=float(traction[0]), shear_traction_mpa=float(traction[1]),
                    signed_registry_rate=float(np.sum(plus-minus)),
                    gross_grid_traffic_rate=float(np.sum(plus+minus)),
                    cumulative_net_transfer=float(np.sum(signed)),
                    cumulative_gross_grid_traffic=float(np.sum(gross)),
                    cumulative_forward_grid_traffic=float(np.sum((gross+signed)/2)),
                    cumulative_backward_grid_traffic=float(np.sum((gross-signed)/2)),
                    **moments)
    p0, m0 = interface_traffic(mass, base_energy, grid, kT, mobility_s, indices)
    history.append(row(0., 0, "initial", np.zeros(2), p0, m0))
    same_phase = [mass.copy()]
    for cycle in range(cycles+hold_cycles):
        start = mass.copy(); before = history[-1]; first = len(history)
        segment = "cyclic" if cycle < cycles else "zero_stress_hold"
        for j in range(steps_per_cycle):
            if cycle < cycles:
                energy, lu, traction = factors[j]
            else:
                energy, lu, traction = zero_energy, zero_lu, np.zeros(2)
            old_pop = np.array([np.sum(mass[:, n == w]) for w in wells])
            if integrator == "implicit":
                updated = lu.solve(mass.ravel()).reshape(shape)
                plus, minus = interface_traffic(updated, energy, grid, kT, mobility_s, indices)
            else:
                # Positivity CFL of the SAME SG generator. Traction stays at
                # this phase's midpoint on the refinement substeps.
                substeps = max(1, int(np.ceil(dt*np.max(-lu.diagonal())/.45)))
                updated = mass.copy()
                plus = np.zeros(len(indices)); minus = plus.copy()
                for _ in range(substeps):
                    forward, backward = interface_traffic(updated, energy, grid, kT, mobility_s, indices)
                    plus += forward/substeps; minus += backward/substeps
                    updated += (dt/substeps)*(lu @ updated.ravel()).reshape(shape)
                    if np.min(updated) < 0:
                        raise FloatingPointError("explicit SG positivity failed")
            minimum_mass = min(minimum_mass, float(np.min(updated)))
            if np.min(updated) < 0:
                raise FloatingPointError("negative cell mass: diagnostic refuses numerical repair")
            mass = updated
            net = plus-minus
            new_pop = np.array([np.sum(mass[:, n == w]) for w in wells])
            balance = new_pop-old_pop-dt*(np.r_[0., net]-np.r_[net, 0.])
            maximum_balance = max(maximum_balance, float(np.max(np.abs(balance))))
            moment_balance = np.dot(wells, new_pop-old_pop)-dt*np.sum(net)
            maximum_registry_balance = max(maximum_registry_balance, abs(float(moment_balance)))
            signed += dt*net; gross += dt*(plus+minus)
            time = (cycle+(j+1)/steps_per_cycle)*load.period_model
            history.append(row(time, cycle+1, segment, traction, plus, minus))
        current = history[first:]
        end = history[-1]
        normal = np.array([r["normal_strain"] for r in current])
        registry = np.array([r["registry_shear"] for r in current])
        loads = np.array([[r["normal_traction_mpa"], r["shear_traction_mpa"]] for r in current])
        q = np.array([[r["mean_a"], r["mean_s"]] for r in [before]+current])
        # Discrete mean external work; reversible static cycle should give zero
        # under refinement. Finite-rate loop area is not irreversible plasticity.
        forces = units.traction_to_force(loads/1000.)
        work = float(np.sum(forces*np.diff(q, axis=0)))
        cycle_rows.append(dict(cycle=cycle+1, segment=segment,
            period_model=load.period_model, same_phase_L1=float(np.sum(np.abs(mass-start))),
            normal_mean=float(np.mean(normal)), normal_amplitude=float(np.ptp(normal)/2),
            registry_shear_amplitude=float(np.ptp(registry)/2),
            registry_index_increment=end["mean_well_index"]-before["mean_well_index"],
            cycle_net_transfer=end["cumulative_net_transfer"]-before["cumulative_net_transfer"],
            cycle_gross_grid_traffic=end["cumulative_gross_grid_traffic"]-before["cumulative_gross_grid_traffic"],
            external_work_ev_cell=work, **{k: end[k] for k in (
                "raw_mass", "mass_residual", "mean_well_index", "outside_central_mass",
                "normal_strain", "registry_shear", "intrawell_shear", "registry_index_shear",
                "cumulative_net_transfer", "cumulative_gross_grid_traffic", "outer_well_mass",
                "upper_normal_edge_mass")}, **{f"P_{w:+d}": end[f"P_{w:+d}"] for w in wells}))
        same_phase.append(mass.copy())
        if progress:
            progress(cycle+1, cycles+hold_cycles)
    audit = dict(boundaries="fixed reflecting on a and outer s", opening_probability=None,
        opening_status="NOT COMPUTED: no validated opening sink in this experiment",
        thermal_normalization="unvalidated one-atomic-cell collective-coordinate hypothesis",
        production_model_changed=False, material_calibration_accepted=False,
        physical_seconds_available=False, physical_hz_available=False,
        temperature_K=temperature_K, kT_ev=kT, mobility_a_model=mobility_a,
        mobility_s_model=mobility_s, initial_central_only=central_only,
        dt_model=dt, steps_per_cycle=steps_per_cycle, period_model=load.period_model,
        integrator=integrator,
        driven_cycles=cycles, zero_stress_hold_cycles=hold_cycles,
        grid_a=len(grid.a), grid_s=len(grid.s), domain_a_edges=[float(grid.a[0]-grid.da/2), float(grid.a[-1]+grid.da/2)],
        grid_da=grid.da, grid_ds=grid.ds, wells=len(wells),
        alignment_error=alignment, max_well_balance_residual=maximum_balance,
        max_registry_balance_residual=maximum_registry_balance,
        max_mass_residual=max(abs(r["mass_residual"]) for r in history),
        max_strain_decomposition_residual=max(abs(r["decomposition_residual"]) for r in history),
        minimum_cell_mass=minimum_mass, numerical_repair=0.,
        gross_activity_status="grid-dependent bidirectional SG traffic; not physical hop count")
    return dict(history=history, cycles=cycle_rows, audit=audit,
                initial_mass=initial, final_mass=mass, same_phase_mass=np.array(same_phase))
