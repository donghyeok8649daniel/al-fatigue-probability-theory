# Aluminum kinetic-calibration source status

No numerical kinetic target is currently accepted.

The repository contains static geometry, cohesion, elasticity, generalized
stacking-fault, surface, and vacancy data. It contains no trajectory or
published coefficient that maps without extra assumptions to the collective
interplanar coordinate $a$, the registry coordinate $s$, or their coupled
relaxation modes.

Bulk self-diffusion, phonon periods, fatigue lifetimes, and macroscopic yield
data are deliberately rejected as substitutes. Their variables and dynamical
reductions differ from this overdamped collective-coordinate model.

An acceptable source must record material/potential, orientation, temperature,
boundary conditions, collective-coordinate definition, sampling interval, and
uncertainty. Preferred raw data are $a(t),s(t)$ trajectories for a matrix
correlation fit. A collective-coordinate diffusivity is usable only when its
coordinate metric matches this model.

Until then `data/aluminum_kinetic_calibration.json` remains uncalibrated and
the desktop application cannot select physical seconds/Hz.

## 2026-09 source audit

Retrieved primary-source records for Olmsted, Hector, Curtin and Clifton,
*Atomistic simulations of dislocation mobility in Al, Ni and Al/Mg alloys*,
Modelling Simul. Mater. Sci. Eng.13 (2005), DOI10.1088/0965-0393/13/3/007,
author manuscript https://arxiv.org/abs/cond-mat/0412324 . This measures edge/
screw dislocation motion, not the current a/s reduced-cell coordinates.
It is relevant future kinetic evidence but supplies neither the required
normal mobility nor a proved mapping to the current registry density.
No coefficient from it was accepted or inserted in the default JSON.

The line-drag coefficient has Pa s units, whereas current cell friction is
J s/m^2. Conversion requires a derived coordinate/energy normalization; an
arbitrary line length or specimen A_c is prohibited. See
`KINETICS_LOADING_AND_STRESS_AUDIT.md` for the dimensional distinction.

The new CSV importer requires actual C_aa,C_as,C_sa,C_ss in m^2 and physical
lag seconds, temperature, same energy model, source and explicit measurement
tolerance. Synthetic exact correlations verify the software only. No suitable
Al covariance dataset was found in the repository or the reviewed source.
