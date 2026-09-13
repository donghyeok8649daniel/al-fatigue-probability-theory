# v33 completed rank-four research controls — NOT adopted Al material

Five actual fixed-range studies,15 LPs including five baseline replays, and
an independent validation of the minimum-training-eta candidate are complete.
See `solver_v1/RANK_FOUR_ENVIRONMENT_V33.md` for the mathematical derivation,
bulk-background counting, numerical checks and exact results.

The best signed diagnostic is k4=8,D4=-.225446749 eV, eta12.450127>1. New
state Hessian error99.2024%, excluded-q error208.5574%. A smaller loss does
not establish Al compatibility or positive energy stability. No candidate
is promoted to a core, production PDE, kinetic calibration or UI model.

`definition.json` records the primary nonnegative-amplitude trial. The
separately named signed profile deliberately relaxes only D4's sign; its
stored `nonnegative` array is authoritative. These are not two combined
energy corrections. The radial ranges are declared research sensitivity
values, not physically identified Al parameters or fitted yield factors.

Source binding metadata originates from the common source loader. Actual
observable units: bulk energy eV/atom; interface energy eV/atomic interface
cell; displacement gradient eV/L0; Hessian eV/L0^2; cubic constants GPa.
L0=2.863782463805517e-10m. No dislocation-line length, statistical area,
physical mobility or seconds is inferred from these static values.

`joint_capacity_replay`: two additional actual LPs allow all five existing
rank4 ranges together. All15 earlier single-profile prediction vectors were
rebuilt with zero discrepancy. Elapsed284.663s. Nonnegative eta12.50729550;
signed eta11.55612833. The latter improves excluded-q error to31.7309% but
sets the repulsive LJ coefficient exactly0, so it is explicitly ineligible.
Both fail eta<=1. Conditional rank8 is not material validation. The initial
`joint_capacity` contains only the declaration before a list/tuple input
error; it is not a completed study. No earlier result is overwritten.

`nested_density_control`: separate same-quadratic-family closure control,
not another rank4 term. All13 LPs completed in81.584s. It mixes bulk-normalized
densities at the same decay, so lambda>0 is an equivalent quadratic envelope
and lambda0 is the exact old exponential endpoint. Baseline replay error0;
maximum mixture identity error2.00e-15. No sampled direction improves the
positive-LJ baseline eta12.50738498. See the quadratic-density derivation for
the equivalent width and the explicitly limited interpretation.
