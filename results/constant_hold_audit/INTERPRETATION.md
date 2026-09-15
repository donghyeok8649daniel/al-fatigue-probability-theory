# Constant 1000 MPa hold: not a creep validation

These are new deterministic production-adapter runs with TwoRowLJ, E=69 GPa,
mean=1000 MPa, amplitude=0, kT*=0.02, Ma*=1, Ms*=0.05, chi=0.2,
three registry wells, a_upper=1.6 and total model time 0.4. They do not
reconstruct an unsaved user UI result whose model/settings are unknown.
No energy, mobility, temperature or probability rule was changed.

## Reproduction

`python -m app.run_constant_hold_audit` writes runs.json and summary.json.
The suite includes separate s, a and dt refinement, identical-duration controls
with different nominal frequency/cycles, and a zero-stress control. All are
implicit runs; explicit/implicit, domain and longer-time convergence are not
certified here. `elapsed_seconds` is program runtime, not model physical time.

| Run | a × s cells | max dt | final opening probability | final epsilon_p |
|---|---|---|---|---|
| reference | 41 × 60 | 0.0005 | 6.03448e-21 | 3.28347e-11 |
| s refinement | 41 × 90 | 0.0005 | 5.89022e-21 | 4.00628e-11 |
| a refinement | 61 × 60 | 0.0005 | 6.48299e-21 | 3.28993e-11 |
| dt refinement | 41 × 60 | 0.00025 | 6.05268e-21 | 3.28371e-11 |

## What is actually increasing?

The solver starts from a principal-well conditional Gibbs distribution at the
mean/preload force, not from an unloaded specimen subjected to a sudden step.
For the reference run total strain starts at 0.02287888840088428 and ends at
0.02287888840282664: change 1.94e-12. This does not demonstrate sustained
macroscopic creep strain. Grid changes shift the initial equilibrium quadrature
much more: final total strain is 0.02184710 for s refinement and 0.02299156
for a refinement. These are not experimentally calibrated creep curves.

Opening probability is accumulated removed mass, NOT 1 minus raw intact mass.
Finite diffusion/drift plus an absorbing opening boundary can yield an outward
flux under a constant load. A cumulative probability therefore need not become
constant when strain has become practically stationary.

Reference absorbed increments over successive 0.1-model-time intervals are
2.14515e-21, 1.37613e-21, 1.27023e-21, 1.24297e-21. The first interval is
larger, but later intervals still have computed outward flux. This is consistent
with depletion/adjustment of the conditional distribution's opening tail, not
proof of material creep damage or cyclic fatigue. No moving force boundary is
driven by the load here because amplitude is zero.

## Numerical resolution

Across the four runs the maximum absolute mass residual reaches 6.30607e-14,
far above the ~6e-21 absorbed signal. The maximum probability history change
from a refinement is 4.48511e-22. Under the required conservative floor
max(mass residual, repair, flux discrepancy, refinement differences), this
opening signal is BELOW NUMERICAL RESOLUTION. It must not become certified
specimen probability through area amplification.

This does NOT prove the stored absorption is merely subtraction roundoff:
it is computed boundary loss and could approximate a very small tail flux.
The present audit cannot certify its physical meaning at this scale. Stored
flux consistency is zero by the discrete accumulation bookkeeping, not an
independent continuum flux-quadrature proof.

The s-grid changes epsilon_p by 7.22801e-12 (~22% of reference). Net registry
transfer is 1.26634e-10 and gross SG traffic 1.40733e-10 in the reference.
The tiny partition moment increases while total strain barely changes, so
intrawell and registry contributions largely compensate. Neither nonzero
epsilon_p nor gross grid recrossing traffic establishes irreversible material
plasticity. No unload/hold or s-domain convergence was performed: residual
plasticity and physical creep remain unvalidated.

## Frequency control

f*=25, cycles=10 and f*=50, cycles=20 both give duration 0.4 with the same
record interval. Their reported final values agree exactly; complete array
comparison is stored/checked separately. At amplitude zero, changing nominal
frequency only changes the scheduled duration if cycle count is held fixed.
There is no oscillatory forcing and it must not be called a fatigue-frequency
effect. Physical seconds/Hz remain unavailable.

All 12 saved history arrays in the identical-duration pair were checked with
exact array equality. The zero-stress control has absorbed probability and
opening flux exactly zero, final epsilon_p=3.86270e-31, and gross SG traffic
1.18555e-16. Its nonzero normal strain baseline (0.00555868) is the finite-
temperature/discrete conditional mean relative to a0, not accumulated creep.
