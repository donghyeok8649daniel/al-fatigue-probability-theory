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
