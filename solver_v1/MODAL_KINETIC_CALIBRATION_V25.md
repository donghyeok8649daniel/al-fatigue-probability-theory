# Mode-resolved kinetic calibration: measured damping, not an invented clock

## Scope and completed checkpoint status

This work fits actual time correlations and independently generates reference
Al MD to test thermostat and timestep sensitivity. It is not a production
energy change or a replacement of the deterministic Smoluchowski equation.
The LJ infinite Poisson/Bessel pair and analytic environmental terms remain
unchanged. Al99 is used only as a kinetic **reference**. Its finite cutoff is
the published reference-potential convention, not a cutoff added to our LJ.

The first completed study uses the actual 8192-frame, 204.775 ps record in
[Fransson and Erhart (2023), Zenodo 10014454](https://zenodo.org/records/10014454).
The first 4096 frames are fitting data; the second 4096 are temporal holdout.
Refitting the second half is a sensitivity calculation, not replacement of
the first-half prediction. These halves and periodic planes are not claimed
independent replicas. Source provenance/checksums remain in v24.

Seven main reference MD records and their analyses are complete: five1ns
controls and two5ns controls. Separate thermal and sampling controls also
completed. All360 modal fit records succeeded and predicted the held-out
correlations better than their overdamped exponential controls; no fitted
parameter hit a bound. This is not production mobility certification.
The final comparison and limitations are in section14 below.

## 1. Exact normalization and mode definition

There are N=12 periodic (111) plane classes, Np=576 atoms per class. With
plane-mean displacement u_l and relative coordinate q_l=u_(l+1)-u_l,

    u_k = N^(-1/2) sum_l u_l exp(-2pi i kl/N),
    q_k = (exp(2pi i k/N)-1) u_k = B_k u_k.

The k=0 rigid translation vanishes exactly. Conjugate modes k,N-k represent
one real sine/cosine pair, not two independent MD samples. For k != N/2,
the real orthonormal coordinates are sqrt(2)Re(q_k), sqrt(2)Im(q_k).
E|q_k|^2 equals either variance only under phase-equivalent equilibrium;
finite-record phase anisotropy and cross-polarization must be audited.
The Nyquist mode k=N/2 is real and has no extra factor two.

For a harmonic per-atom Bloch Hessian H(k), the plane energy is

    E = Np/2 sum_k u_k^dagger H(k) u_k,
    C_q,k(0) = |B_k|^2 kBT [Np H(k)]^(-1).

This is the independently checked v24 static normalization. No statistical
correlation area is used. A collective plane coordinate is not one atomic
interface cell, nor a rigid semi-infinite half-crystal coordinate.

Under an invertible linear coordinate transformation on nonzero modes,
mobility transforms contravariantly:

    M_q(k) = B_k M_u(k) B_k^dagger,
    M_u(k) = M_q(k)/[4 sin^2(pi k/N)].

Even this exact transformation does not make M_u local or independent of k.
There is a DIFFERENT, exact energy-normalization identity: for the same
collective coordinate, write F_bar=F_plane/Np. Then

    qdot=-M_plane grad F_plane=-(Np M_plane) grad F_bar.

Thus M_bar=Np M_plane is exact for this declared rescaled collective PMF.
For the FULL stochastic generator the thermal term must also be normalized:

    L P = div M_plane[P grad F_plane+kBT grad P]
        = div M_bar[P grad F_bar+(kBT/Np) grad P].

This is energy bookkeeping, NOT a physical temperature of T/Np. Merely
changing F and M while leaving the numerical thermal-energy term kBT fixed
does NOT preserve the plane probability law: it amplifies its diffusion by Np.
The exact rescaling does not require statistically independent atoms. Interpreting M_bar as
a size-independent LOCAL atomic/interface-cell mobility, or replacing F_bar
by the production W, DOES require additional spatial/PMF equivalence tests.
Historical result keys containing `cell_extensivity_hypothesis` refer to this
unvalidated interpretation, not an uncertainty in the arithmetic factor Np.
No value is copied into the production a/s generator.

## 2. Damped position correlation, including the factor two

For a specified equilibrium linear mode the candidate correlation satisfies

    c'' + 2g c' + w0^2 c = 0, c(0)=1, c'(0)=0.

For wd^2=w0^2-g^2>0,

    c(t)=exp(-g t)[cos(wd t)+(g/wd)sin(wd t)].

This is the **position** correlation, not the velocity-current formula.
Its Laplace transform and integral are

    c_tilde(z)=(z+2g)/(z^2+2gz+w0^2),
    tau_integral=c_tilde(0)=2g/w0^2.

The envelope time 1/g is NOT tau_integral. For the convention in
[dynasor's DHO derivation](https://dynasor.materialsmodeling.org/dev/tutorials/dho_peak_fitting.html),
the damping coefficient Gamma=2g. Confusing Gamma and g gives a factor-two
error; confusing envelope decay and integral time can give a much larger one.
No atomic mass is chosen to set production time. Time units come from the
source MD timestamps. In a DHO interpretation m_eff=H/w0^2 is a diagnostic,
not permission to replace the production overdamped equation with inertia.

For the same coordinate and matched equilibrium PMF, H=kBT/C0 and

    Gamma0 = H tau_integral,
    M0 = C0/(kBT tau_integral).

This supplies a **conditional model-based modal estimate** only if the fitted
correlation describes the relevant zero-frequency memory. A good resonance
fit need not determine arbitrarily low-frequency tails. In matrix form,

    K = integral_0^infinity C(t)dt,
    Gamma0 = kBT C0^(-1) K C0^(-1),
    M0 = C0 K^(-1) C0/(kBT).

Do not silently replace coupled matrices by independently fitted raw axes.
See the existing ZERO_FREQUENCY_MOBILITY_AUDIT.md for the exact ordering and
[Lei, Baker and Li (2016)](https://arxiv.org/abs/1606.02596) for a general
data-driven memory approach. The present DHO is a tested low-complexity
candidate, not a proof of Markov dynamics or validated escape kinetics.

## 3. Deterministic fit protocol and actual first result

Six independent mode indices, three axes, and cutoffs 3.15/6.35/12.75 ps give
54 fits. Nine fixed starts per fit are saved. Log-positive bounds and least
squares prevent invalid parameters but **do not certify positive physical
friction**. All starts, bounds, Jacobian singular values, failures and
held-out errors are retained. Correlated lag errors are not treated as
independent likelihood samples; singular values are not confidence intervals.

All 54 first-half DHO fits have smaller heldout correlation RMS than a
single overdamped exponential control. Median RMS: .0752775 vs .3589432.
This resolves the need for an oscillatory source-coordinate description at
the measured short times, not the validity of the extrapolated low-frequency
mobility or the production PDE over barriers.

For k=1 and cutoff12.75 ps (first-half fit):

|Axis|g [ps^-1]|1/g [ps]|wd [rad/ps]|tau_integral [ps]|
|---|---:|---:|---:|---:|
|normal|.0965516|10.3572|15.2273|.000832770|
|direct110|.0743697|13.4463|7.04057|.00300028|

Second-half refits and weaker transverse modes show appreciable sensitivity.
No selected favorable window becomes the default. The largest normalized
off-diagonal modal covariance in the full record is .135725; local pooled
off-diagonals are much smaller (maximum .0219214). Small pooled coupling
alone does not certify every mode decoupled.

Machine-readable first-study outputs:
`results/modal_kinetic_calibration_v25/first_half_fit/`.

## 4. Independent reference MD controls (completed)

Engine: official serial LAMMPS22Jul2025 update4, shared-library version
20250722, EAM/alloy OpenMP. Installer SHA256 matches the official distribution:
47f1aeb0fcbeadcc211c9b1c258152c3545464e3381975a681bf558891a06986.
It was extracted into an ignored research cache, not installed system-wide.
The source Al99 potential MD5 is eb0f0b204ea40787274efcf8e44d6de2.

One independently generated 12^3 FCC crystal at a_lat=4.065 Angstrom is
equilibrated25 ps at nominal300 K with1 ps Nose-Hoover damping, initial
600 K velocities and fixed seed28459 following the source protocol. Atom
ordering/engine differ, so this is not claimed a bitwise published-record
reproduction. Actual atom coordinates and velocities in a saved restart are
shared by the matched12-plane comparisons. Six-plane size and lattice
controls have their own equilibrated reference states. The new NVT thermostat chain starts fresh; the
first25 ps of every sampled run are excluded by the analysis protocol.

Predeclared sampling:1000 ps,25 fs saved interval, same box/potential/state:

1. NVE,dt5 fs;
2. NVT300 K,damping1 ps,dt5 fs;
3. NVE,dt2.5 fs.

Actual mean temperature, energy drift/range, modal covariance, frequency,
damping, temporal halves and fitting-window changes must be compared.
Chaotic trajectories need not agree pointwise after timestep refinement;
their statistical observables and integration conservation are the checks.
Reference MD is not Monte Carlo fatigue-trajectory counting.

## 5. Production status and remaining projection

No production Ma,Ms,t0 has yet been certified. The exact dimensionalization
M*=t0 E0/L0^2 M_phys is unchanged. Source modal ps are physical MD timestamps;
production PDE still uses model time. Neither a fitted resonance lifetime nor
a scoped modal friction establishes rare opening escape, defect transport,
fatigue lifetime or the current assumed ratio Ms/Ma=.05.

Next decisions require the actual controls above, mode/phase uncertainty,
and a coordinate/PMF-consistent coarse-graining. Changing the static energy,
using an arbitrary area, or adjusting mobility to a desired fatigue result
is not part of this calibration.

## 6. Momentum-conserving local-drag hypothesis

A local drag on absolute plane displacement u would damp uniform translation.
That is not an intrinsic momentum-conserving bulk friction law. A different,
testable Rayleigh hypothesis acts on neighboring **relative velocities**:

    R = (Np eta/2) sum_l |qdot_l|^2,
    qdot=B udot,
    Gamma_u(k)=Np eta |B_k|^2.

Using the explicitly declared reference-MD mass m per atom gives

    2g_k=eta |B_k|^2/m,
    g_k=D |B_k|^2, D=eta/(2m).

This use of the reference MD's physical kinetic energy is not a hidden
atomic-mass clock for production: measured g is an independent required
input. The relative-coordinate mobility would be constant
M_q=1/(Np eta) on the zero-mean subspace only if this local Rayleigh model
is actually supported. Fit modes1–3 and predict the remaining modes; keep
per-mode coefficients and temporal/ensemble uncertainty. A failed hypothesis
is not repaired by choosing one favorable mode or inventing eta.

## 7. Independent experimental linewidth reference obtained

[Glensk et al. (2019), PRL123,235501](https://doi.org/10.1103/PhysRevLett.123.235501)
measured Al phonon linewidths with neutron scattering. Fig2(b), bottom-right,
contains transverse L-point q=(1/2,1/2,1/2) markers. The publisher PDF page3
and enlarged panel were visually inspected; marker centers were then taken
from PDF vector paths, not guessed from a search snippet.

|Approximate temperature from figure [K]|Gamma [THz, cycle-frequency]|g=pi Gamma [ps^-1]|Envelope 1/g [ps]|
|---:|---:|---:|---:|
|295|.271619|.853316|1.17190|
|700|.589945|1.85337|.539558|
|800|.753967|2.36866|.422180|
|900|.965610|3.03355|.329646|

Declared half-point plot-reading sensitivity is .00888 THz and7.11 K;
this is NOT the experiment's uncertainty or a confidence interval. The
room-temperature marker is near295 K, not silently declared exactly300 K.
PDF SHA256 and all marker/axis coordinates are saved with the extraction.

The paper's Eq2 uses cycle frequency nu, whereas dynasor uses angular
frequency. Consequently Gamma_PRL=g/pi, Gamma_dynasor=2g. In energy units,
Gamma_dynasor[meV]=2 hbar[meV ps]g[ps^-1]. Confusing THz, rad/ps and meV
would obscure a real reference-model discrepancy.

These experimental numbers are a held-out **finite-wavevector modal**
kinetic target. They are not a measurement of the zero-frequency a/s PMF
friction. Resonant damping measures memory near the phonon resonance; its
equality to zero-frequency friction needs an additional validated model.
No empirical t0 or production mobility is extracted from this one mode.

The separate Fransson et al.2021 dynasor study uses NVE sampling and fixed
4.05 Angstrom, unlike the public4.065 Angstrom NVT record. Its Al example
is an Al99 reference calculation, not an Al experiment. Comparisons must
retain this lattice/ensemble distinction.

## 8. Direct finite-band spectral mobility (independent of DHO extrapolation)

For the two-sided transform with cycle frequency nu,

    S_q(nu) = integral_(-infinity)^infinity C_q(t) exp(-2pi i nu t) dt,
    S_q(0)=2 K, K=integral_0^infinity C_q(t)dt

for reversible equilibrium correlations. Estimate S by two DPSS tapers with
NW=2 and3, temporal mean removed separately for every plane in each record.
For taper w and equivalent planes p, the exact normalization is

    S_hat_ij(nu) = dt/[N sum_t w_t^2]
                  sum_p FFT(w q_pi) FFT(w q_pj)^*.

The positive frequency bins are NOT doubled. A Parseval test explicitly
reconstructs the taper-weighted covariance with both frequency signs. The
finite-band proxy K_band=Re(mean_band S_hat)/2 is used in the full3x3 formula
Gamma_band=kBT C0^-1 K_band C0^-1, M_band=Gamma_band^-1.
These are **finite-band** estimates, not a demonstrated zero-frequency limit.
Positive semidefinite spectra arise by construction and do not certify them.
Imaginary cross spectra, all blocks and all chosen bands remain saved.

The initial protocol uses .02–.04,.04–.08,.08–.16 cycles/ps, full retained
record, halves and quarters, NW2/3. After seeing the first NVE full-record
normal estimate cluster around1.7–1.8e10 m²/(J s), an explicit FOLLOW-UP study
adds .005–.01 and .01–.02 cycles/ps. It is not falsely described as predeclared
before the first analysis. Bands below the taper half-bandwidth NW/T or with
fewer than three bins are rejected and recorded, not silently extrapolated.

Undamped-mode control: replace each measured mode by a sinusoid at each of
the two fitted temporal-half frequencies, keeping the measured variance.
The worst phase is calculated as the largest eigenvalue of the2x2 cosine/sine
spectral quadratic form; no random sampling. Orthogonal spatial modes add
with exact conjugate weights. For the first1ns NVE record, this worst-phase
window leakage is at most6.90e-5 of the measured low-band integral proxy.
Thus ordinary leakage from those undamped mode lines cannot explain the
observed low-band spectrum. This is NOT a complete statistical uncertainty or
a bound on unknown anharmonic spectral components.

## 9. Actual dimensional mobility calculations and their scope

The initial published8192frame DHO projection gives, across every prescribed
window and both temporal halves:

|Plane-gap coordinate|M [m²/(J s)]|Gamma [J s/m²]|
|---|---:|---:|
|normal|5.7283e10–8.1358e10|1.2291e-11–1.7457e-11|
|direct110 slip|1.1722e11–1.6776e11|5.9610e-12–8.5307e-12|

The slip/normal ratio ranges1.518–2.929. These are sensitivity ranges, not
confidence intervals. `mobility_projection/mobility.json` retains all values.
For unit interpretation,1 Å²/(eV ps)=6.241509074e10 m²/(J s).
The same plane-gap PMF can be normalized by576 atoms per represented plane.
Multiplication of M by576 then exactly preserves its dynamics under that
energy normalization ONLY with the matched thermal term kBT/576. Its identification with a LOCAL production cell is a
separately labeled **local extensivity/PMF hypothesis**, not a change of SI
units. A reduced two-atom cell must not be equated to this normalization
without its own counting and coordinate derivation.

An alternative harmonic normalization q_hat=sqrt(Np) delta_q gives
M_hat=Np M_plane while retaining kBT, but q_hat is no longer the measured
physical gap. The transformed potential is F_plane(q_hat/sqrt(Np)); outside
the harmonic neighborhood it is not the original local W(q_hat), and its
registry period is rescaled. This transformation cannot be used to manufacture
single-cell slip/opening kinetics from a plane-average trajectory.

The independently generated1ns NVE,12-plane,5fs record yields direct
full-record finite-band values1.7149e10–1.8297e10 normal and
4.4194e10–5.5313e10 slip in the initial bands. Its DHO values are instead
7.0588e10–7.2652e10 normal and1.3049e11–1.5226e11 slip. The separation is
larger than the DHO window variation: an accurate resonance correlation fit
does not by itself recover low-frequency mobility. The lower-band follow-up
widens the normal range to1.4655e10–1.8711e10 and slip to
4.4194e10–7.4670e10. No favorable value is selected as a clock.

## 10. Spatial locality is tested, not inferred from a plane average

For each nonzero Fourier mode, real sine/cosine coordinates are constructed
with sqrt(2) normalization (the Nyquist coordinate remains real). Their
spectra are symmetry-pooled and the full3x3 polarization response is retained.
This tests M_q(k), whereas a local plane marginal has already eliminated other
planes. Pooling does not certify finite-record phase-offdiagonal dynamics.

Even a hypothetical constant M_q on the zero-mean constrained subspace does
not make its scalar marginal mobility identical: if H_q is also constant,

    C_local=(N-1)/N C_mode,
    K_local=(N-1)/N K_mode,
    M_local=(N-1)/N M_mode.

For nonlocal H(k) these simple factors no longer suffice. Consequently matching
Np M_local across two sizes is useful but does not prove a constant local
Rayleigh coefficient or equivalence to a rigid active-interface PMF.
The measured first1ns NVE phase-pooled normal modal cell-hypothesis values
range roughly6.9e12–2.5e13 across modes/bands, with stronger slip variation.
Final timestep/ensemble comparisons must be included before a locality claim.

## 11. What these results can and cannot calibrate

Measured reference-MD modal damping and finite-band mobilities are now numeric,
dimensioned results, not missing-input placeholders. Their transfer to the
analytic production coordinate remains a separate model question:

1. Source/reference versus real-Al kinetics: compare the independent
   experimental L-point linewidth, not only the same-potential MD.
2. Resonance versus low-frequency memory: the DHO extrapolation mismatch is
   quantified by the direct spectrum, not repaired by changing static energy.
3. Local versus spatially coupled coordinates: retain mode dependence and
   the finite periodic constraint, rather than insert an arbitrary area.
4. Intrabasin versus defect/slip/opening kinetics: equilibrium phonons alone
   do not validate activated transport, fatigue life or absorbing flux.

The production M_a=1,M_s=.05, physical kinetic JSON and physical-Hz gate are
unchanged. A_c is absent from every mobility calculation. Atomic mass enters
only the explicitly identified reference-MD modal drag hypothesis, not an
invented overdamped time scale.

## 12. Why one measured phonon lifetime cannot identify the missing clock

The issue is an inverse problem, not merely a missing units conversion. For
a scalar passive exponential-memory example,

    Gamma_tilde(z)=gamma_infinity + A/(z+1/tau_b),
    Re Gamma_tilde(i omega_r)=gamma_infinity + A tau_b/(1+omega_r^2 tau_b^2),
    Gamma_tilde(0)=gamma_infinity + A tau_b.

For a measured resonance friction gamma_r, choose any positive A,tau_b with
A tau_b/(1+omega_r² tau_b²)<gamma_r, and set gamma_infinity to the positive
remainder. These kernels give the same resonance friction but different
zero-frequency mobility1/Gamma_tilde(0). They are a mathematical
identifiability example, NOT fitted Al parameters or a new production model.
The direct spectrum/DHO discrepancy supplies an actual reason to test memory
rather than assume away this freedom.

For a scalar equilibrium linear GLE, the correlation obeys

    m C''(t)+integral_0^t Gamma(t-u) C'(u)du+H C(t)=0,
    C'(0)=0,
    Gamma_tilde(z)=[(m z²+H) C_tilde(z)-m z C0]/[C0-z C_tilde(z)].

At z=0 this reduces exactly to Gamma0=H K/C0, with H=kBT/C0 in the
canonical harmonic PMF. The source MD's finite microcanonical/canonical
difference, anharmonic PMF, spatial projection and measured temperature must
still be assessed; the algebra alone is not an Al calibration certificate.

The follow-up thermal-channel test records plane velocities and both raw and
streaming-subtracted kinetic energy per atom. A complex linear spectral
predictor of q from the two adjacent plane kinetic energies is fit on one
temporal half and evaluated on the other, in both directions. A time-shifted
negative control is retained. Negative held-out explained fractions are not
clipped. Correlation with kinetic energy would suggest an omitted thermal
channel, but cannot by itself identify a heat-diffusion mode: potential energy,
energy current, causality and electron-mediated transport are not measured by
this test. No empirical heat correction is inserted into the PDE.

### Measured thermal contribution and sampling control

The completed N6, 1 ns NVE control gives a held-out normal spectral fraction
of 0.1946–0.4060 explained by adjacent-plane internal kinetic energies.
Time-shifted controls give negative fractions; direct110 slip does not show
a consistent improvement (approximately -0.0801 to +0.0392). Thus an omitted
thermal channel contributes to the normal low-frequency response, but this
is neither a complete memory model nor a causal heat-transport measurement.
The authoritative result is `thermal_channel_mass_checked/summary.json`.

An engine-level audit found that the Al99 EAM file overrides the preceding
`mass 1 26.98` command with 26.982 amu. Dynamics already used 26.982; only the
first kinetic-energy postprocessor used 26.98. Its relative energy mismatch
7.406052e-5 becomes 6.297410e-8 after the explicitly recorded correction.
The original arrays remain preserved. The current runner reads the mass from
the actual engine after pair initialization. Covariance/temperature/integral
mobility estimates do not use atomic mass and are unchanged by this repair.

A separate 250 ps trajectory saved every 5 fs was compared with its exact
25 fs decimation, using matched 225 ps records after the initial exclusion.
The decimated coordinates exactly reproduce the existing trajectory prefix.
Across the three initial bands and NW=2/3, the largest diagonal mobility
change is 2.410652e-4 (0.0241%). Therefore the large DHO/direct-spectrum
discrepancy is not explained by this sampling interval. This is a sampling
control, not a zero-frequency convergence certificate.

### External linewidth comparison

At the transverse L point, the five completed 1 ns reference controls give
linewidths 0.0713–0.1397 THz (the NVE subset 0.1080–0.1397 THz). The
Glensk et al. figure gives approximately 0.2716 THz for the ambient-temperature
experimental marker and 0.1861 THz for the separate 300 K DFT-MD marker.
These are source-MD/phonon frequency units, not the production PDE clock.
Figure-reading precision is not experimental uncertainty; temperature,
ensemble, fixed lattice and experimental broadening must remain explicit.
No linewidth ratio is multiplied into production mobility as a correction.

## 13. Specific remaining measurement, rather than an arbitrary conversion

To transfer these estimates into production, the kinetic measurement must use
the same collective coordinate and conjugate PMF as the selected analytic
energy model. A useful next experiment is a small imposed generalized force
on that coordinate, with momentum-balanced reactions on the surrounding
crystal, measuring the complex response over the low-frequency bands already
tested here. Compare this independently measured response with the covariance
prediction; do not fit a desired fatigue probability or strain amplitude.

For a constrained-coordinate force-correlation route, the projected random
force and orthogonal/constrained dynamics must be derived explicitly. The
unconstrained total force autocorrelation is not automatically the required
memory kernel. Check state dependence with opening and registry displacement,
finite-size dependence and a memory/Markov limit before choosing a constant M.
Neither experiment has been performed in this checkpoint. The existing
reference phonon estimates are useful kinetic data but not a substitute for
these coordinate/PMF and material checks.

## 14. Completed long-record comparison and decision

The N6 controls each contain200001 frames over5ns; after the common25ps
exclusion,4975ps remain. The dt2.5fs run completed in2933.73s wall time.
Its first1ns coordinates exactly reproduce the previous short run, and
adding thermal observations also leaves that trajectory exactly unchanged.
Nested prefixes are not independent replicas.

| Quantity | dt2.5fs | dt5fs |
|---|---:|---:|
| Mean T [K] |299.0682|298.1255|
| Energy range [eV/atom] |2.49533e-5|1.98929e-4|
| Energy slope [eV/atom/ps] |-1.15438e-9|2.64469e-8|
| Normal M, initial three bands [m²/(J s)] |6.70155e10–7.91291e10|6.95676e10–8.05442e10|
| direct110 M, initial three bands [m²/(J s)] |2.00431e11–2.27522e11|1.79906e11–2.38736e11|

These are full-record NW2/3 sensitivity ranges, not confidence intervals.
The timestep change of the taper-mean estimate in the initial bands is
0.86–4.32% normal and5.06–10.37% slip. Lower bands do not yield a uniformly
stable plateau: at.005–.01cycles/ps the slip discrepancy is51.75%; at
.001–.002cycles/ps it is30.45%. Both record-length and block variation remain.
The narrower taper-only range must not conceal this larger uncertainty.
The lowest band also cannot be independently estimated in the shorter blocks
under the stated resolution rule; invalid estimates remain in the JSON.

`long_record_comparison/paired_dt_summary.csv` contains every band and
half/quarter range, and `all_estimates.csv` retains both time steps and all
prefix lengths. The plot was visually inspected. `all_control_comparison`
adds thermodynamics, all seven modal summaries and external linewidth tests.
The5ns transverse L-point widths are0.11934–0.12130THz (dt2.5fs) and
0.11976–0.12020THz (dt5fs): resonance damping agrees much more closely than
the lowest-frequency mobility. This directly shows why resonance agreement
alone cannot certify the slow mobility.

**Decision:** dimensional reference-coordinate mobility has been calculated.
A unique constant zero-frequency production M_a/M_s has not been validated.
The primary remaining issues are memory/low-band uncertainty, spatial and
PMF normalization, and real-material transfer. No production physical clock
or physical-Hz enablement follows from this checkpoint.

Final executed tests: targeted12 passed1.42s; full solver673 passed906.20s;
app34 passed97.53s, no skips; desktop smoke passed. The reference MD runner
also passed a real mass-aware0.1ps runtime smoke. These verify implementation,
not experimental fatigue or local activated-transition kinetics.
