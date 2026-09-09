# Correlation-integral / low-frequency mobility audit

This is a **research diagnostic of actual MD**, not a replacement PDE.
Short-lag oscillations reject a direct harmonic overdamped Markov fit, but
do not by themselves disprove a lower-frequency Smoluchowski limit. We test
the latter separately instead of interpreting a phonon period as t0.

## 1. Exact linear-response identity and units

Let q be a specified physical collective coordinate [m],
C(t)=<q(t)q(0)^T> [m²], and C0=kBT H^-1 at equilibrium. For a harmonic
generalized Langevin projection with equilibrium initial conditions,

    m q_ddot + integral_0^t Gamma(t-u) q_dot(u) du + H q = noise,

the causal noise is orthogonal to q(0) and C_dot(0)=0 in the inertial case.
This equation is used to interpret MD, **not inserted into our production
Smoluchowski solver**. No atomic mass or damping value is chosen for Al.
For Laplace variable z,

    [m z² + z Gamma_tilde(z) + H] C_tilde(z)
        = [m z + Gamma_tilde(z)] C0.

Assuming decaying equilibrium correlations and finite zero-frequency memory,
define K=integral_0^infinity C(t)dt [m² s]. At z=0,

    H K = Gamma0 C0,
    Gamma0 = H K C0^-1 = kBT C0^-1 K C0^-1,
    M0 = Gamma0^-1 = C0 K^-1 C0/(kBT).

Thus [Gamma0]=J s/m² and [M0]=m²/(J s). Matrix ordering matters; H and the
friction need not commute. For a direct overdamped Ornstein-Uhlenbeck model,
C(t)=exp(-M H t)C0 gives K=(M H)^-1 C0 and exactly the same M0=M.
`test_zero_frequency_kinetics.py` verifies a noncommuting coupled example.
A separate deterministic, dimensionless damped-oscillator control has
negative C eigenvalues at short lag but positive K and the correct friction.
Its finite integration tail is checked, not ignored to loosen a tolerance.

Memory-based coarse graining is distinct from fitting one exponential to
each raw coordinate: see [Lei, Baker & Li (2016), PNAS,
DOI10.1073/pnas.1609587113](https://arxiv.org/abs/1606.02596). The matrix
identity above is explicitly derived here under its stated harmonic and
equilibrium assumptions; it is not a claim that all MD admits a constant M.

## 2. Finite-data estimator and safeguards

For measured lags, retain every finite K(T), not a hand-picked plateau.
Normalize with the measured positive C0:

    Q(T)=C0^-1/2 sym[K(T)] C0^-1/2  [s].

`correlation_integral_audit` reports all eigenvalues of Q(T), antisymmetric
sampling error, maximum time-block discrepancy, every-other-lag trapezoidal
change, and ||Q(T)-Q(T/2)||. Blocks are not independent confidence intervals.
Plane classes are correlated; their count is not a multiplier for sample
confidence. There is no clipping of negative integrals, eigenvalue repair,
or selected cutoff that automatically certifies M. The formal algebra helper
refuses a nonpositive integral and never returns a production calibration.

A convergent positive K would only give a low-frequency response for the
**same coordinate and PMF**. It would not establish the complete memory
kernel, opening escape kinetics, a Markov reduction over all loads, or the
production mobility ratio M_s/M_a. Bounded-coordinate MSD does not supply a
free-particle diffusivity by treating a saturation value as a diffusion rate.

## 3. Actual 300 K Al99 trajectory result

Source: [Fransson & Erhart, Zenodo10014454](https://zenodo.org/records/10014454),
periodic6912-atom crystal,576 atoms per(111) plane class,12 classes,25fs saved
spacing. The complete1024frame projected record spans25.575ps. Projection SHA:
`2037b0104c33573abc27532517610ff84f308e6ba55ded08a0b0b2bf834a7eb6`.
No cell/atom count is fitted. The full3.5GB upstream checksum is not claimed
verified from a compressed prefix. Source is NVT with1ps thermostat damping;
[the source workflow](https://dynasor.materialsmodeling.org/dev/get_started/workflow.html)
also cautions that thermostats influence measured dynamics. NVE/thermostat
sensitivity is a separate missing check, not a reason to invent a damping.

Actual outputs: `results/public_aluminum_validation_v15/zero_frequency_integrals/`.

| Study | Lag cutoff [ps] | Eigenvalues of Q [ps] | Empirical floor [ps] |
|---|---:|---|---:|
| first512frames,4blocks |1.55|-.0146321, .0000407, .0081562|.0193909|
| full1024frames,4blocks |3.15|-.0107767, -.0067566, -.0041907|.0382522|
| full1024frames,8blocks |1.55|-.0148557, -.0126256, -.0022011|.0587256|

None of the latter-half cutoffs has a smallest integral eigenvalue above its
empirical floor. The result is **unresolved**, not physically negative
friction. The short record cannot supply the converged zero-frequency
integral needed by the formula. Window changes and source thermostat effects
cannot be repaired by taking absolute values or choosing the first positive
integral. A later, longer record must be reported as a distinct study.

## 4. Exact remaining calibration requirement

Needed: longer stationary trajectories of matched opening/registry or
defect coordinates, declared boundary and thermal conditions, a compatible
finite-temperature PMF, coordinate/area normalization, and resolved memory/
low-frequency convergence. All-mobile periodic plane means are not rigid
half-crystal displacement coordinates. Even a successful source-plane M0
could not be copied into the reduced cell without deriving that projection.

Current production values remain M_a,phys=null, M_s,phys=null, t0=null.
The established conversion M*=t0 E0/L0² M_phys is unchanged. Physical PDE
seconds/Hz remain disabled. Actual MD picoseconds in these research plots
refer to the source record only.
