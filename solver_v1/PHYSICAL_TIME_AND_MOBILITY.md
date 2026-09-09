# Physical time and reduced-coordinate mobility

## 1. Scope and current conclusion

The length/energy dimensionalization is complete, but physical seconds and
hertz are not currently available. No repository datum measures the mobility,
friction, diffusion, or relaxation of the same collective coordinates
$(a,s)$. The committed kinetic-calibration file is therefore explicitly
`uncalibrated`; it contains no fabricated mobility or time scale.

## 2. Implemented probability equation

The production $N=1$ finite-volume code evolves

$$
\frac{\partial P^*}{\partial t^*}=-\nabla_*\cdot\mathbf J^*,
$$

$$
J_i^*=-M_i^*\left(P^*\frac{\partial G^*}{\partial q_i^*}
+k_BT^*\frac{\partial P^*}{\partial q_i^*}\right),
\qquad(q_1^*,q_2^*)=(a^*,s^*).
$$

Scharfetter--Gummel uses $D_i^*=M_i^*k_BT^*$ and the energy jump divided by
$k_BT^*$. Explicit and backward-Euler modes use the same spatial generator.

## 3. Exact variable transformation

For the current common coordinate scale,

$$
a_{\mathrm{phys}}=L_0a^*,\quad
s_{\mathrm{phys}}=L_0s^*,\quad
G_{\mathrm{phys}}=E_0G^*,\quad
t_{\mathrm{phys}}=t_0t^*.
$$

The calibrated hybrid uses $L_0=b$, with
$b=2.8637824638\times10^{-10}$ m for the selected 4.05 angstrom FCC
reference. Its numerical energy parameters are in eV, so one convenient
conversion unit is $E_0=1\ \mathrm{eV}=1.602176634\times10^{-19}$ J. This
declaration does not provide a time scale. The older UI `TwoRowLJ` mechanism
surface remains a dimensionless reference unless its energy unit is declared
independently.

Temperature and force transform as

$$
k_BT^*=\frac{k_BT_{\mathrm{phys}}}{E_0},\qquad
F_{q,\mathrm{phys}}=\frac{E_0}{L_0}F_q^*.
$$

The stress input mapping remains

$$
\sigma\longrightarrow\frac{\sigma}{E}
\longrightarrow f^*=\kappa_{\mathrm{axial}}\frac{\sigma}{E}.
$$

No area is inserted in this constitutive map.

## 4. Probability-density scaling

Normalization is invariant:

$$
1=\int P_{\mathrm{phys}}\,da_{\mathrm{phys}}\,ds_{\mathrm{phys}}
=\int P^*\,da^*\,ds^*.
$$

Therefore, for two length coordinates,

$$
P_{\mathrm{phys}}=\frac{P^*}{L_0^2}.
$$

The one-coordinate reduced slow-$s$ density instead scales as
$P_{s,\mathrm{phys}}=P_s^*/L_0$.

## 5. Mobility dimensions and time scale

The physical overdamped law is

$$
\dot q_{i,\mathrm{phys}}=-M_{i,\mathrm{phys}}
\frac{\partial G_{\mathrm{phys}}}{\partial q_{i,\mathrm{phys}}},
$$

so

$$
[M_{i,\mathrm{phys}}]=\frac{\mathrm{m^2}}{\mathrm{J\,s}}.
$$

Substitution into either the overdamped law or Fokker--Planck equation gives

$$
M_i^*=\frac{t_0E_0}{L_0^2}M_{i,\mathrm{phys}},
\qquad
t_0=\frac{M_i^*L_0^2}{M_{i,\mathrm{phys}}E_0}.
$$

There is no omitted prefactor. Energy and length determine force and
stiffness, but overdamped time additionally requires mobility or friction.
Atomic mass cannot fill that gap without changing to inertial dynamics.

## 6. One common time and two mobilities

The current code uses

$$
M_a^*=1,\qquad M_s^*=0.05
$$

on one $t^*$ axis. Retaining the current generator requires

$$
\frac{M_{s,\mathrm{phys}}}{M_{a,\mathrm{phys}}}
=\frac{M_s^*}{M_a^*}=0.05.
$$

Independently measured mobilities with a different ratio would require an
evidence-backed update of the dimensionless ratio, not two incompatible
definitions of $t_0$.

## 7. Time and frequency conversion

Once a valid calibration exists,

$$
t_{\mathrm{phys}}=t_0t^*,\qquad t^*=t_{\mathrm{phys}}/t_0,
$$

$$
f_{\mathrm{Hz}}=f_{\mathrm{model}}/t_0,
\qquad f_{\mathrm{model}}=t_0f_{\mathrm{Hz}}.
$$

For $N$ cycles, $T_{\mathrm{phys}}=N/f_{\mathrm{Hz}}$ and the solver interval
is exactly $T^*=T_{\mathrm{phys}}/t_0$. These conversions do not define $t_0$.

## 8. Coupled relaxation calibration route

The reduced pristine dynamics are

$$
\delta\dot{\mathbf q}^*=-\mathbf M^*\mathbf H^*\delta\mathbf q^*.
$$

The dimensional Hessian and symmetric rate matrix are

$$
\mathbf H_{\mathrm{phys}}=\frac{E_0}{L_0^2}\mathbf H^*,
\qquad
\mathbf R_{\mathrm{phys}}=\mathbf M_{\mathrm{phys}}^{1/2}
\mathbf H_{\mathrm{phys}}\mathbf M_{\mathrm{phys}}^{1/2}.
$$

Its eigenvalues are physical inverse relaxation times. In a justified
decoupled limit, $M_{q,\mathrm{phys}}=1/(\tau_qH_{qq,\mathrm{phys}})$. For
coupled motion, rates and preferably mode shapes or the fitted relaxation
matrix are required. Two times alone can leave a coordinate-exchange branch
ambiguity. `physical_time.py` reports that limitation.

## 9. MD correlation-function procedure

An atomistic calibration must define the same collective $a(t)$ and $s(t)$,
hold the intended temperature/orientation/boundary conditions, and estimate

$$
\mathbf C(t)=\langle\delta\mathbf q(t)\delta\mathbf q(0)^T\rangle.
$$

Near a harmonic basin, normal-mode correlations follow

$$
C_\alpha(t)\simeq C_\alpha(0)e^{-t/\tau_\alpha}.
$$

The matrix fit supplies rates and mode vectors to combine with the same energy
Hessian. Atomic vibration periods are not substituted without demonstrating
the overdamped collective-coordinate mapping.

## 10. Diffusion/Einstein route

For the same collective coordinate, with a justified metric and Markovian
overdamped reduction,

$$
D_{q,\mathrm{phys}}=M_{q,\mathrm{phys}}k_BT,
\qquad M_{q,\mathrm{phys}}=D_{q,\mathrm{phys}}/(k_BT).
$$

Bulk Al self-diffusion is not a substitute: it describes a different variable
and mechanism.

## 11. Source status

No mapped $a/s$ relaxation trajectory, collective diffusivity, friction, or
matrix correlation function is present in the repository. The file
`data/aluminum_kinetic_calibration.json` contains known unit scales but null
mobility/time fields. Hypothetical values occur only in tests.

## 12. UI semantics

The default UI remains in **Model time** and displays cycles/model-time plus a
visible warning. **Physical time** is selectable only when a valid calibration
is loaded. It then accepts hertz, converts to the same internal model period,
uses seconds on time axes, converts rates to inverse seconds, and reports
$t_0$, both mobilities, relaxation times, and source. It is not a second
physics path.

Changing language or specimen areas never reruns the PDE. Changing physical
frequency does rerun it because $\omega\tau$ changes the dynamics.

## 13. First passage and frequency

Cumulative opening loss over a fixed number of cycles may depend on physical
frequency because each cycle provides different relaxation/escape time. The
code does not impose frequency-independent probability per cycle. Per-cycle
output records absorbed mass, minimum opening barrier, peak first-passage
flux, and end survival. A large first-cycle decrement can be initial Gibbs-tail
depletion, continuing escape, or a numerical moving-boundary artifact; the
table exposes but does not itself classify it.

## 14. Specimen probability remains separate

Time calibration neither calibrates $A_c$ nor proves independent regions. The
UI separately reports local PDE probability, an explicitly uncertified
mathematical independent-region extrapolation, and certified physical specimen
probability. The last requires an actual deterministic refinement study and a
local signal above its numerical floor. Preview alone never certifies it.

## 15. Current limitation

The mathematical dimensionalization is complete. Physical seconds and hertz
remain unavailable until $M_{a,\mathrm{phys}}$ and
$M_{s,\mathrm{phys}}$ are obtained from kinetic data mapping to these
collective coordinates. Static Al fitting, fatigue lifetime, and specimen
probability cannot supply them.

## 16. Evidence-bound import and physical-generator audit

`kinetic_calibration_workflow.py` and `import_kinetic_calibration.py` now accept
actual matrix correlation data, check the harmonic fluctuation/relaxation
relations, and bind the result to a static parameter fingerprint and coordinate
definition. The UI has a calibration-file loader; it no longer requires editing
the committed default. The default remains uncalibrated.

Physical mode uses the supplied temperature through kBT/E0 and both declared
model mobilities. The old axis conversion alone failed to propagate those
values if they differed from the reference. Model mode still uses exactly
M_a=1, M_s=.05, kT=.02. Keeping the old generator requires the old ratio;
independently measured different ratios now have an explicit same-PDE path.
No ratio or temperature is changed without selecting the supplied calibration.

See `KINETICS_LOADING_AND_STRESS_AUDIT.md` for matrix-log derivation, CSV units,
provenance requirements, limitations, and the rejected dislocation-mobility
shortcut. Validation of a file's internal consistency is not certification of
its experimental/MD provenance or quantitative Al fatigue predictions.
