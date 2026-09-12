# Mobility calibration from conjugate impedance — v28

## Outcome and scope

We executed a deterministic **finite-low-band scalar drag calibration** on the
completed v25 reference-MD spectra and tested it against the completed v26/v27
conjugate forced responses. This is new parameter estimation and uncertainty
propagation, not a new MD trajectory, material fit, or production-clock update.

| Reference periodic plane-gap coordinate | Fitted drag [J s/m²] | Inverse drag [m²/(J s)] | Full control-fit range [m²/(J s)] |
|---|---:|---:|---:|
| normal | 1.380000849e-11 | 7.246372352e10 | 6.523118156e10–8.846382195e10 |
| direct110 slip | 4.585776422e-12 | 2.180655810e11 | 1.859261220e11–2.681768261e11 |

These are **Al99 reference-potential, periodic 144-atoms-per-plane coordinate
candidates**, not calibrated local production mobilities or measured aluminum
mobility. Their ratio is 3.009306870, not production's 0.05. Do not overwrite
production Ma=1, Ms=0.05 or enable a common physical clock with these numbers.
The ranges are observed sensitivity spans, not confidence intervals.

## 1. Same coordinate, same conjugate force

Keep the previously verified definition

    q_l = (mean(u_{l+1}) - mean(u_l)) dot e,
    G_ext = -F q_0,
    atomic forces = +/- F e/N_plane,
    chi = (c_cos - i c_sin)/F,
    q(t) = Re[chi F exp(i omega t)].

No area factor or density rescaling enters this force. There are six periodic
planes, N_plane=144, and sum_l q_l=0. Other plane gaps/components are not
independently clamped when one conjugate force is applied. Therefore

    1/chi_qq != (inverse(full susceptibility matrix))_qq

in general. The scalar inverse below includes the response of the eliminated
coordinates. It must not be called an element of the bare diagonal mobility.

## 2. What the measured phase actually determines

Define the scalar dynamic stiffness and its dissipative coefficient:

    Z(omega) = 1/chi(omega),
    Gamma_diss(omega) = Im Z(omega)/omega
                     = -Im chi/(omega |chi|²).

Only when Gamma_diss>0 define the point inverse-drag proxy

    M_diss(omega) = 1/Gamma_diss
                 = omega |chi|²/(-Im chi).

A positive noisy Im chi is preserved; abs/clip is not used to fabricate
positive friction. A finite-frequency inverse-drag proxy is not automatically
the zero-frequency mobility, nor the full complex inverse memory kernel.

For an inertial/memory *reference diagnostic* one can write

    Z(omega)=H - m omega² + i omega Gamma(omega).

Then Im Z/omega is Re Gamma, whereas the storage term may contain both bare
inertia and reactive memory. The saved diagnostic

    D_disp(omega) = [1/chi(0) - Re Z(omega)]/omega²

is **not a calibrated atomic mass**. Full-force v27 gives about 0.290 eV ps²/Å²
for both coordinates; this records storage dispersion, not a time scale to
insert into the overdamped production equation. Production remains inertialess.

Units with chi in Å²/eV and frequency in cycles/ps:

    omega=2 pi f,
    [Gamma_diss]=eV ps/Å²,
    1 Å²/(eV ps)=6.241509074460763e10 m²/(J s).

## 3. Exact propagation of an empirical complex-response disk

Given observed chi=z and an empirical disk |delta chi|<=r, if r<|z| then
the inverse map is exactly the disk

    center_Z = conjugate(z)/(|z|²-r²),
    radius_Z = r/(|z|²-r²).

Consequently

    Gamma_lower = (Im center_Z-radius_Z)/omega,
    Gamma_upper = (Im center_Z+radius_Z)/omega.

If Gamma_lower<=0<Gamma_upper, passive drag is feasible but **there is no
finite upper bound on M**. We save null for this unbounded/unavailable limit,
not a small positive clipped drag. If the input disk reaches the pole at
chi=0 we withhold finite inverse bounds; this is conservative also in the
tangent-disk special case. Feasibility of positive drag is tracked separately.

Two sensitivity sets are retained for every 180/360 ps response:

1. max |complex unforced lock-in displacement| / |F| over all matched windows
   and represented planes;
2. that radius plus observed timestep, force-amplitude, block and sign
   differences, added without an independence assumption.

The second sum can double-count correlated thermal errors and is deliberately
a conservative sensitivity construction, **not a rigorous probability bound**.
Neither disk is a statistical confidence interval. Shared restart, sites and
windows are not counted as independent replicas. v26 90 ps quarter fits remain
in their original results; no invented 90 ps null is assigned here. v26 has no
completed matching timestep study, explicitly marked incomplete rather than
claiming the zero numerical-error contribution used here is an error estimate.

Full-force finite-frequency point estimates and null-only inverse ranges:

| coordinate / reference f [cycles/ps] | dt [fs] | point M [m²/(J s)] | null-only M range |
|---|---:|---:|---:|
| normal / 2 | 2.5 | 1.866937317e11 | 9.460181161e10–6.975379462e12 |
| normal / 2 | 1.25 | 1.343060154e11 | 7.899714570e10–4.475436272e11 |
| slip / 1 | 2.5 | 2.271510529e11 | 1.301572350e11–8.893023484e11 |
| slip / 1 | 1.25 | 2.116723381e11 | 1.248123166e11–6.946742014e11 |

Including the observed controls removes a finite upper bound in **all four**
cases. At f=0.05 all full-record estimates are already unbounded under the
null-only construction. Thus detecting a signed phase is a weaker result than
accurately identifying its inverse dissipative coefficient.

## 4. Low-band calibration, objective and identifiability

For the scalar equilibrium coordinate, canonical FDT gives

    chi(0)=C0/(kBT),
    -Im chi(omega)=omega S(omega)/(2 kBT).

If a low-frequency limit exists and chi approaches chi(0),

    Gamma_0 = kBT K0/C0²,
    K0 = S(0)/2 = integral_0^infinity C(t) dt.

We substitute the measured **finite-band** K_band, not a certified K0.
No inference of zero-frequency convergence follows from positive PSD alone.
The scalar result C0²/(kBT K_band) differs from the corresponding diagonal
of C K^-1 C/(kBT) when cross-correlations are retained. Both values are saved.

Calibration bands are the existing broad bands [0.02,0.04), [0.04,0.08),
[0.08,0.16) cycles/ps. The central record is the full 5 ns, dt=2.5 fs,
NW=3, two-taper estimate, with its initial25 ps excluded. For each band,
the discrepancy scale s_k is the max-minus-min drag over the available
1/2/5 ns prefixes, dt=2.5/5 fs, NW=2/3, full/half/quarter records. This
mixes finite-record, numerical and thermal sensitivity; it is NOT a reported
experimental uncertainty. A small normalized loss can reflect broad scales.

Fit the physically positive one-parameter hypothesis:

    L(Gamma)=sum_k [(Gamma-Gamma_band,k)/s_k]²,
    Gamma_fit=sum_k Gamma_band,k/s_k² / sum_k 1/s_k².

This is the exact deterministic WLS minimum, not random search. A nonpositive
solution is rejected. We fit drag, not an average of noisy reciprocal values.
For p=Gamma/Gamma_reference the Jacobian has entries Gamma_reference/s_k,
rank1 with singular values 4.695598849 (normal), 3.640045340 (slip) for
the recorded reference=max |training drag|. This establishes identifiability
of the *one scalar hypothesis* only, not full matrix/PMF/material identifiability.
No covariance or parameter confidence interval is inferred from these scales.

The objectives are 0.1126388401 and 0.0340778481. We repeated the fit for
168 control groups (84 per coordinate), keeping the same discrepancy scales.
There are784 axis/band/control observations, not784 independent experiments.
Another196 matrix-band records fail their original bandwidth/bin-resolution
checks; all are retained with reasons in `unavailable_source_bands.csv` rather
than assigned a zero spectrum or silently treated as successful observations.

With full5 ns records only:

| dt [fs] / NW | normal M | slip M | units |
|---|---:|---:|---|
| 2.5 / 2 | 7.329516899e10 | 2.191082226e11 | m²/(J s) |
| 2.5 / 3 | 7.246372352e10 | 2.180655810e11 | m²/(J s) |
| 5 / 2 | 7.560017851e10 | 2.279826175e11 | m²/(J s) |
| 5 / 3 | 7.462107610e10 | 2.288950867e11 | m²/(J s) |

## 5. Lower-band and forced-response validation

The four lower bands are excluded from the loss. They are **post-fit validation
on previously inspected development data**, not blind held-out experiments.
For the central5 ns/NW3/dt2.5 estimate, fitted-drag relative errors include
-32.38% normal at[.001,.002), -34.60% slip at[.005,.01). Across every
represented control they reach90.78% and51.93%, respectively. Lower frequency
does not yet yield a tightly certified plateau. These failed tests are saved.

With the same measured static susceptibility and fitted drag, the *pure
overdamped* harmonic prediction is

    chi_OD(omega)=1/[1/chi(0)+i omega Gamma_fit].

At0.05 cycles/ps it predicts phase -0.03656° normal, -0.06142° slip. The
complex discrepancy is within the observed noise disk, but the loss remains
unresolved. Being compatible with a noisy result is not a positive precision
validation.

Reference temperature is approximately299 K; driven NVE records have the
previously reported finite heating. Canonical FDT versus the finite NVE
ensemble and source-potential versus real-material kinetics remain additional
uncertainties. This analysis does not fit a thermostat or absorb heating into
an invented production drag.

At the v27 high frequencies, full-force discrepancies are about10.06–10.08
observed null radii for normal and7.67–7.73 for slip, predominantly storage.
The measured storage is12.08%/15.86% above static susceptibility. A reversible
reflecting overdamped linear response of the **same conjugate coordinate** has

    chi(omega)=beta sum_k w_k lambda_k/(lambda_k+i omega), w_k>=0,
    Re chi(omega)<=chi(0).

Thus changing one positive M cannot fix that storage increase while retaining
the same static PMF. This does not prove absence of a lower-frequency Markov
limit or reject constant dissipative drag in a model retaining inertia/memory.
We do not add inertia to production or retune its Hessian to hide the mismatch.

## 6. Numerical verification and reproduction

New utilities verify units, harmonic impedance, exact inversion-disk extrema,
zero/noisy loss, scalar-versus-matrix distinction, deterministic constrained
fit, rejection/no-overwrite, and persistence of closed calibration gates.
The runner reproduces byte-identical outputs from the same completed inputs.

The separate raw-record control independently reconstructs the scalar
periodogram and covariance from both actual5 ns NPZ trajectories (dt2.5/5fs),
using source hashes and actual sample times. Maximum relative arithmetic
differences from the stored summaries are2.22e-16 and4.44e-16. This checks
data/units/implementation, not convergence or physical truth to that precision.
The plot is rendered and inspected; spans are explicitly not CIs.

    python -m solver_v1.run_impedance_calibration --out <fresh-directory>
    python -m solver_v1.run_impedance_controls verify --source <raw-MD-directory> \
      --saved-summary <matching-v25-summary.json> --out <fresh-verification.json>
    python -m solver_v1.run_impedance_controls plot --results <result-directory>

Results are in `results/impedance_calibration_v28`. Raw MD caches and previous
results are preserved. No new reference-MD execution is claimed in this stage.
See CURRENT_WORK_HANDOFF.md for the actual final test/Git results.

## 7. What has and has not been calibrated

Completed: finite-band scalar drag candidates with provenance, deterministic
loss/residuals,168 control fits, exact complex inverse-error propagation and
comparison to actual forced response. This provides numbers and falsifiable
tests rather than selecting one favorable phase or bandwidth.

Not completed: a precise zero-frequency plateau, low-frequency driven loss
resolution, isolated/local coordinate PMF correspondence, matrix mobility and
reference-potential-to-material kinetic validation. Plane normalization still
requires the matching kBT/N_plane when dividing its free energy by N_plane;
M alone cannot be transferred to a single-cell PDE at unchanged kBT.

Next useful measurement is a low-frequency conjugate response with *a priori*
duration/amplitude design and verified linearity/heating, followed by the same
coordinate/PMF spatial test. Do not use the high-frequency phase as a surrogate.

The existing equilibrium planning table at0.2 cycles/ps gives single-record
SNR3 durations13.235 ns normal and7.310 ns slip at one thermal-RMS force.
Under the stationary linear-noise assumptions, four times that force would
divide these durations by16 (about827/457 ps), but its linearity is **not yet
tested**. This is an external measurement design, not a mobility change or
an executed new trajectory. Increasing force does not remove the heat-budget
issue: with S the two-sided spectrum and L=-Im chi,

    T_required = 2 S r_SNR²/(F² L²),
    mean_power = omega F² L/2,
    expected_work_at_required_T = omega S r_SNR²/L
                               = 2 kBT r_SNR²  (canonical FDT).

Thus planned SNR3 has about18kBT=0.464 eV net work irrespective of amplitude
in this approximation. This is not a rigorous noise bound, nor a measured
heating prediction; finite NVE, transients and memory must still be checked.
It explains why simply amplifying the forcing is not by itself a kinetic
calibration solution. No forced higher-amplitude result is fabricated here.

LJ/Bessel, static calibration, full/reduced production PDE, Ma/Ms, physical
kinetic JSON, UI and A_c are unchanged. Physical production seconds and Hz
remain disabled. Reference-MD ps units are not production fatigue Hz.
