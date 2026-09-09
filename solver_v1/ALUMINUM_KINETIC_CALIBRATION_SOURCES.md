# Aluminum kinetic-calibration source status

No numerical kinetic target is currently accepted for the production a/s PDE.
Actual published trajectory, line-drag and experimental relaxation data have
now been retrieved and analyzed; see the v15 update below. "No accepted cell
mobility" must not be paraphrased as "no physical kinetic data exist."

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
Al covariance dataset was found in the repository or the source reviewed at
that earlier stage. This historical statement is superseded by the v15 audit.

## v15: actual source data and why the production clock still stays disabled

`PUBLIC_ALUMINUM_CALIBRATION_EVIDENCE.md` contains the derivations, source
hashes and completed results. Public crystalline Al99 MD (Fransson/Erhart2023,
DOI10.5281/zenodo.10014454) now supplies physically timed plane-coordinate
covariances. Its resolved negative short-lag correlation eigenvalue is
incompatible with directly fitting the reversible harmonic overdamped a/s
model. The source's all-mobile periodic plane averages also differ from the
production coordinates. Equal-time harmonic/MD comparison separately checks
the576-atoms-per-plane energy normalization; it does not infer mobility.

Actual Olmsted2005 Table3 line-drag coefficients and Gorman1969 experimental
endpoints are now stored with Pa s units in
`data/aluminum_line_drag_references.json`. `line_drag_reference.py` derives
the conditional local-friction/profile metric needed to connect B_line to a
cell mobility, but marks it unvalidated and returns no opening mobility or t0.

Actual293K-analysis microwire relaxation/activation data (Verheyden et al.2018,
DOI10.1016/j.dib.2018.11.047) constrain specimen mechanisms, not direct a/s
diffusivity. Their apparent activation area is neither A_c nor atomic cell
area. Source exclusions, processed/raw stress offsets and a diameter mismatch
are retained. No experimental fatigue lifetime is fitted.

Accordingly, the production kinetic JSON is intentionally unchanged:
M_a,phys=null, M_s,phys=null, t0_seconds=null. Finding real kinetic data is
progress; pretending distinct coordinates share that clock is not.

## v15 low-frequency check beyond the short-lag rejection

`ZERO_FREQUENCY_MOBILITY_AUDIT.md` derives the coupled correlation-integral
route: Gamma0=kBT C0^-1 [integral C(t)dt] C0^-1. It can in principle describe
low-frequency friction even when short-time correlations oscillate. Actual
1024frame/512frame,4/8block and quadrature comparisons were run, not merely
proposed. Finite integrals change sign and remain below their empirical
block/quadrature/antisymmetry floor; no positive converged mobility is inferred.
The source NVT thermostat and all-mobile plane coordinate require independent
audits. A positive selected cutoff is not a production clock certificate.
