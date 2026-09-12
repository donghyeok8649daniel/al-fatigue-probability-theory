# v30: internal-step conjugate-work resolution

Research reference MD only. Production LJ/Bessel energy, PDE, mobilities,
temperature convention and disabled physical seconds/Hz are unchanged.

## Why measure before extending the weak probe

v29 measured small phase lags, nonlinear strong-probe harmonics and a direct
power/energy discrepancy up to .061785 eV. Power was stored every25 fs while
the atomistic integrator stepped at1.25 or2.5 fs. Halving the integrator step
alone cannot remove a fixed output-sampling quadrature error. We first test
that explanation using the SAME forcing and restart, not retuned kinetics.

For q=(mean upper displacement-mean lower displacement).e, external energy
is -F(t)q. Equal-mass atoms receive +/-F e/N_plane, so

    P_ext(t)=sum_i f_ext,i.v_i=F(t)(v_upper-v_lower).e.

This is the conjugate external power, not energy difference relabeled as
dissipation. The internal energy includes pair/embedding plus kinetic energy,
but excludes the external potential. NVE implies dE_internal/dt=P_ext in
the continuous equations. Velocity-Verlet and quadrature have finite errors.

## All internal samples, bounded memory

Let the saved-frame stride be r, internal step dt, P_j the END-step power.
Native LAMMPS `fix ave/time 1 r r` supplies

    Pbar_k = (1/r) sum_(j=(k-1)r+1)^(kr) P_j.
    Qright_N = r dt sum_k Pbar_k.
    Wtrap_N = Qright_N - dt(P_N-P_0)/2.

This is exactly the trapezoid on ALL internal samples, not the trapezoid of
saved25 fs endpoints. Constant memory is used by the engine's one-block
average. Time0 is excluded from the right sum and reintroduced via the
endpoint correction. Native accumulation is opt-in; old sampling and forces
are preserved. It can be compared against explicitly saving every step.

The group velocity is the arithmetic plane mean because the source is
monatomic Al with equal masses. This implementation must not be generalized
to mixed masses without matching the actual collective-coordinate weights.

API semantics: [LAMMPS fix ave/time](https://docs.lammps.org/fix_ave_time.html)
defines the preceding Nrepeat samples ending on Nfreq. The installed engine
is also tested directly; API documentation alone is not numerical validation.

## Predeclared experiment

- Same source Al99 potential/restart hashes asv29, N6/864 atoms, NVE.
- Normal and direct110 at0.2 cycles/ps and positive4RMS reference force.
- dt2.5,1.25,.625 fs, each25 ps, saved25 fs with internal-step power.
- Each corresponding2 ps dense record saves every step independently.
- All12 records are new; aggregate162 ps. No selected sign/window.
- Dense record tests accumulation, not independent thermal replication.
- This25 ps audit does NOT measure low-frequency mobility; initial transients
  and onlyfive forcing cycles preclude that interpretation.

Compare sparse/dense trajectories and work over the identical first2 ps;
compare coarse-power, internal-power and integration-by-parts work with the
actual internal-energy change over25 ps. dt convergence is a separate axis
from frame-stride convergence. Signed work and all residuals are retained.

Runner: `python -m solver_v1.run_work_resolution_v30 run/analyze`.
All12 main records completed, plus6 additional legacy-default2 ps controls.
Aggregate main/legacy measurement174 ps (the1 ps API pilot is separate).

## Measured convergence

RMS of E(t)-E(0)-W(t) across ALL saved25 ps samples, eV:

| dt (fs) | Normal RMS | Slip RMS | Normal max | Slip max |
|---:|---:|---:|---:|---:|
| 2.5 | .00245401 | .00202913 | .00824767 | .00683047 |
| 1.25 | .000580752 | .000509438 | .00193533 | .00193965 |
| .625 | .000152057 | .000122884 | .000546090 | .000391963 |

RMS decreases aboutfourfold per timestep halving, consistent with second-order
integration/quadrature error in this tested regime. Final endpoint errors
alone are not monotone and were NOT used to declare convergence. Saved25 fs
power alone has max residuals .08935/.07354/.07984 eV normal and
.01544/.01312/.01346 eV slip: no analogous refinement trend.

The all-step/block accumulator matches dense independently recorded2 ps work
to6.59e-14 eV; trajectory differences<=1.94e-24 m. Six additional runs with
the original default (measurement OFF) give exactly the same sampled
coordinates and thermo arrays over2 ps as measurement ON. Thus the new
observer does not supply damping, alter force or repair energy numerically.

This resolves the specific coarse-output power quadrature limitation for
future reference measurements. It does not retroactively generate missing
internal-step powers for v29, and does not fix v29 nonlinear probe response,
thermal sampling, zero-frequency or local-coordinate mapping limitations.
Short25 ps work checks do not determine a mobility or a fatigue lifetime.
The long weak-probe/independent-restart experiment remains a subsequent task;
it was not run under this audit. Production physical seconds/Hz remain OFF.

Machine-readable output: results/work_resolution_v30, including source and
trajectory hashes, dense equivalence and legacy-default equivalence.

Actual regression: targeted16 passed2.69s; solver706 passed1824.40s;
app34 passed250.16s; desktop smoke passed. No tests removed or skipped.
