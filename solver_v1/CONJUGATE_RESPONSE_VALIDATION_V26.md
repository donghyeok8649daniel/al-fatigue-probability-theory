# Conjugate-force validation of the measured plane coordinate

## Scope

This repairs a missing validation path, not the production mobility by an
arbitrary multiplier. Al99 is an external reference potential only. The
analytic LJ/Poisson/Bessel energy, production PDE, static fit, specimen area
and physical-time gates are unchanged. Reference MD has its own physical
timestamps; these do not calibrate the production model clock.

The v25 finite-band plane mobility ratio is2.4692–3.1820 on matched full
records in the initial bands, compared with production's model ratio.05.
Naively transferring both coefficients would imply inconsistent clock ratios
49.38–63.64. That is not a license to tune the production ratio: the PMF,
coordinate normalization and Markov limit must first be matched.

## Exact generalized work

For fixed membership of adjacent periodic plane classes0 and1, let

    q = e · (mean_i_in_1 u_i - mean_i_in_0 u_i),
    G_ext = -F(t) q,
    F_i = +F(t)e/N1, i in plane1,
    F_i = -F(t)e/N0, i in plane0,
    F_i = 0 otherwise.

Thus sum_i F_i=0 and sum_i F_i·du_i=F dq exactly. The momentum balance does
not require equal atom counts. There is no extra atomic/statistical area in
this generalized force. In this calculation each periodic plane class has
144 atoms. This coordinate is not an isolated rigid interface: all other
atoms relax and sum_l q_l=0 on the periodic stack.

The actual engine applies the cosine force at EVERY integration step, using
[LAMMPS fix addforce](https://docs.lammps.org/fix_addforce.html) equal-style
time variables. The force is specified per atom in eV/Angstrom. Recorded
thermodynamic energy excludes the external time-dependent potential, so its
change is compared with integral F*qdot dt, NOT asserted to be conserved in
a driven run. NVE means no thermostat, not zero external work. Center-of-mass
velocity and work-energy residual remain independent checks.

## Response and covariance conventions

For F(t)=Re[Fhat exp(i omega t)] and qhat=chi(omega)Fhat, fit

    q(t)=offset+c_cos cos(omega t)+c_sin sin(omega t),
    chi=(c_cos-i c_sin)/Fhat.

Canonical equilibrium linear response for this same conjugate coordinate is

    R(t)=-C'(t)/(kBT),
    chi(omega)=[C(0)-i omega integral_0^infinity C(t)exp(-i omega t)dt]/kBT.

The integral of R truncated at T additionally contains
-C(T)exp(-i omega T)/(kBT). Both expressions and cutoff sensitivity are
saved; a nonzero endpoint cannot be quietly discarded. The conversion is
1 m²/J =16.02176634 Angstrom²/eV. A DHO synthetic identity verifies the sign
and prefactors independently of the MD fitting code.

FDT here is a canonical prediction tested against finite NVE records. Finite
microcanonical/canonical and anharmonic differences remain possible; fitting
chi does not prove all Markov, PMF or material assumptions. In particular,
small phase lag can be harder to resolve than the dominant in-phase elastic
response, and inaccurate Im(chi) precludes a reliable inverse-friction fit.

For an equilibrium scalar autocorrelation, the two-sided spectrum gives the
alternative exact identity Im chi(omega)=-omega S(omega)/(2kBT). A finite
band average is an approximation to S at the driven frequency, not an exact
point value. The comparison retains NW2/3, full/half records and bands
.04–.08,.04–.06,.045–.055cycles/ps around the.05cycles/ps drive. This avoids
silently trusting cancellation in one finite time-domain covariance integral.

For a scalar eliminated-coordinate response the dissipative impedance is
Re Gamma(omega)=Im[chi(omega)^(-1)]/omega. Its inverse is not the3x3 matrix
mobility measured in v25, nor automatically its diagonal element. At finite
frequency it is not a demonstrated zero-frequency coefficient. Cross-axis
responses are retained separately. A negative/noise-dominated dissipative
estimate must not be clipped or silently selected away.

## Declared protocol

The equilibrium input is the completed v25 six-plane,5ns,dt2.5fs NVE record.
At each axis use F_sigma=kBT/sqrt(C0), with matched force units. Amplitudes
are ±0.5F_sigma and ±F_sigma, chosen before seeing the response. They are
small-displacement mechanism probes, not calibrated fatigue/yield loading.
Each of eight runs lasts400ps at.05cycles/ps (20 cycles), with2.5fs integration
and25fs saved samples. The first40ps (two cycles) are excluded uniformly.
The same equilibrated restart is used; signed pairs are not independent
replicas. No random sampling of parameters or Monte Carlo counting is used.

Full, half and quarter response fits are saved with signed/amplitude controls.
Some quarter windows span a noninteger number of cycles; the fit includes an
offset and uses the actual time rather than a simplistic Fourier-bin pick.
The sign-paired chi mean removes deterministic even response in the linear
limit, while sign disagreement and block variation are retained, not treated
as independent-sample confidence intervals. Equilibrium FDT uses all planes
under their translation symmetry, with2,5,10,20,40ps integration cutoffs.

## Acceptance limits

Check actual completion, bounded small displacements, momentum/work balance,
signed and amplitude consistency, temporal sensitivity, and complex FDT
agreement. This is a first forced-response validation, not a completed
spatial/timestep/temperature convergence suite. Do not infer constant M from
the in-phase response alone. A failure remains a recorded result; production
Ma/Ms/t0 are not changed to make the reference agree.

## Planning the loss-component precision

For a stationary long record with two-sided noise spectrum S at the drive,
the asymptotic variance of a fitted quadrature is approximately2S/T_record.
Under unchanged weak-drive fluctuations, resolving the displacement loss
Fhat |Im chi| with planning signal/noise r therefore requires approximately

    T_record = 2 S r² / (Fhat² |Im chi|²).

This is a design estimate, not a measured confidence interval, nor a new
production time scale. Shared initial conditions and finite-window thermal
correlations prevent assuming independent replicas. With r=3, full-RMS
force, NW3 and.045–.055cycles/ps source band, the estimates are156ns normal
and97ns direct110. Thus400ps is a useful work/linear-response validation but
was not sufficient to certify the very small dissipative quadrature. An
increase of force cannot be used without rechecking linearity, and moving
the test to a higher frequency does not establish the zero-frequency limit.

## Independent engine phase benchmark

`run_force_phase_benchmark.py` uses a deliberately synthetic dimensionless
oscillator, m=K=1, drag=.002, omega=2pi*.05. It is NOT an Al model or a
mass-based production time scale. The exact susceptibility is
1/(K-m omega²+i omega drag). The engine uses the same every-step cosine
variable and sampled harmonic fit, with analytically initialized periodic
motion. At dt=.005,.0025,.00125 the phase errors are respectively
3.70e-10,1.12e-10,3.78e-11rad; complex response errors decrease from
6.37e-7 to1.55e-7. This rules out a gross sign/clock problem in this tested
forcing/analysis path, not all finite-temperature MD discretization errors.
The synthetic viscous force and spring are test-only; no damping was added
to the actual NVE Al reference runs. See the official
[viscous-force](https://docs.lammps.org/fix_viscous.html) convention.

## Completed results and disposition

All eight400ps runs completed. Mean temperatures are299.049–299.177K.
The table uses the signed-pair full-record response, with the longest40ps
covariance cutoff as an explicit reference; all other cutoffs remain saved.

| Axis | Force fraction | Re chi [A²/eV] | Relative Re error | Im chi [A²/eV] |
|---|---:|---:|---:|---:|
|normal|.5|.00236386477|+.2518%|+5.73899e-6|
|normal|1|.00235886141|+.0396%|+5.40434e-6|
|direct110|.5|.0119373361|+.1262%|+1.95909e-5|
|direct110|1|.0120639147|+1.1879%|-7.95602e-5|

Doubling force changes normalized Re chi by0.21% normal and1.06% slip.
This supports the conjugate normalization/near-linear in-phase response in
this limited study, not a global nonlinear or discretization validation.
Every full/half/quarter range of paired Im chi crosses zero. Signed-pair
disagreement is also retained. Thus the small dissipative quadrature is
**unresolved**, not negative friction and not a reliable mobility fit.

Maximum center-of-mass velocity change is1.25e-14Angstrom/ps. Maximum
work-energy residual is.012083eV using sampled F*qdot, or.005743eV using
integration by parts in F*q. These refer to the entire864-atom system, not
one cell. They remain visible. For comparison, the equilibrium-spectrum
loss predicts order.001–.002eV total dissipative work over20 cycles at the
full force scale. Therefore this work balance cannot certify such a small
physical dissipation. Negative finite-record work is not negative friction:
thermal exchange/noise and integration residual are still comparable.

The unforced current runner reproduces the old first five coordinate AND
thermodynamic frames exactly. The new force path is optional and does not
change the previous reference trajectory. No static or production parameter
was refitted; physical kinetic JSON and the production seconds/Hz gate remain
unchanged. This checkpoint repairs the missing conjugate-response test and
quantifies the measurement needed; it does not claim the overall mobility
calibration problem is solved.

Executed tests:18 targeted passed.74s;676 full solver passed968.59s;
34 app passed101.49s, no skips; desktop smoke passed. Optional actual engine
forced/unforced smoke and the three-step-size synthetic phase benchmark also
ran successfully. Scientific phase-resolution acceptance remains unpassed.
