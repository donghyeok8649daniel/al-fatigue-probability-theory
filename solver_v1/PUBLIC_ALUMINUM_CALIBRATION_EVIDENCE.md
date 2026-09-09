# Public Al evidence, coordinate normalization and kinetic limits (v15)

Actual primary-source data were downloaded and numerically inspected in this
work. This is not a list of hypothetical calibration routes. However, source
measurements, source-potential predictions and current LJ/Bessel predictions
are different records. No source coefficient is silently promoted into the
production probability solver. This note and the machine-readable artifacts
distinguish what was measured, derived, tested and remains uncalibrated.

## 1. Loading capabilities actually connected

| Path | Actual capability | Limitation |
|---|---|---|
| `interface_static_scenarios.resolved_tensor_tractions` | symmetric3x3 stress -> one normal/two shear tractions | projection, not spatial mechanics or a3D PDE |
| `vector_interface_reference` / static research | W_int(a,s1,s2), three-coordinate stationary solves | static research, not production |
| `run_low_stress_cyclic_diagnostic` | independent normal/shear mean, amplitude and phase in deterministic SG(a,s) | reflecting research problem; no validated opening probability |
| `app.desktop_ui -> solver_adapter -> build_energy_model -> run_probability_pde_2d` | scalar axial stress, fixed axial projection, P(a,s) | no independent UI shear or selectable physical loading direction |

For a common physical frame, t=sigma n, T_n=n.t and tau_i=m_i.t. Stress
components that do not change this interface traction can still change the
prestressed bulk environment; traction projection alone does not implement
that additional physics. A disabled UI direction placeholder is not a solver.
No production registry ID or default energy changed in this task.

## 2. The retrieved crystalline trajectory

[Fransson and Erhart (2023), Zenodo10014454](https://doi.org/10.5281/zenodo.10014454)
publishes crystalline Al trajectories used by the
[dynasor tutorial](https://dynasor.materialsmodeling.org/dev/tutorials/aluminum_solid_dynamic.html).
The actual `main.in` specifies LAMMPS metal units,6912 atoms, a_lat=4.065Å,
12³ conventional FCC cells, periodic boundaries,300K NVT,1ps thermostat
damping,5fs integration and25fs saved frames after25ps equilibration.
The source potential is the same Al99 release used for our static targets.

Verified provenance:

- main.in published MD5: `25f86114b8179f7fcc590b7545acdfa0`;
- Al99 SHA256: `60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`;
- log.lammps published MD5: `28a50d823517269ca0fef426c17a6874`.

The first completed analysis used1024 frames/25.575ps. Bounded HTTP ranges
read179,961,893 compressed bytes without saving the complete3.5GB archive.
Its whole-file MD5 is **not verified** by a prefix. The exact prefix/projection
hashes and source parameters are saved. The archive's code is not executed.
Longer-prefix studies, if completed, receive separate directories and status.

### The actual collective coordinate

For the cubic periodic torus with N=12 repeats, (111) sites have N plane
classes, not36 independent free layers. Each has Np=4N²=576 atoms. With fixed
site assignments, define mean displacements u_l and

    q_l = u_(l+1)-u_l,
    components = (normal[111],direct[1,-1,0],transverse[1,1,-2]).

Periodic minimum-image site displacement is allowed only while assignment is
unambiguous. Diffusion, moving cells or large displacements are refused, not
repaired. Plane classes are correlated; Np is geometric counting, NOT N_eff
or a statistical independence assumption. This coordinate lets all atoms move
and is not the production two-row cell or a rigid-half-crystal interface.

## 3. Equal-time energy/length/atom-count normalization

Let H(k) be the per-atom Bloch restoring matrix in physical J/m², and use a
unitary discrete Fourier transform over the N plane classes. The harmonic
energy restricted to uniform in-plane displacement is

    E_harm = (Np/2) sum_k u_k^dagger H(k) u_k.

Consequently,

    <u_k u_k^dagger> = kBT/Np H(k)^(-1),
    <q_k q_k^dagger> = 4 sin²(pi k/N) kBT/Np H(k)^(-1),
    <q_l q_l^T> = (1/N) sum_(k=1)^(N-1) <q_k q_k^dagger>.

The translation k=0 cancels exactly; it is not inverted or regularized.
Allowed cubic wavevectors are (2pi/(N a_lat))(k,k,k). Every conversion
eV/Å² -> J/m² is explicit. No atomic mass, frequency or mobility is used.

`SourceEAMBlochHessian` differentiates the actual published source EAM:

    H(q)=sum_R [1-cos(q.R)] Hess[phi+2F'_bulk rho](R)
         +F''_bulk g g^T,
    g=sum_R sin(q.R) grad rho(R).

The source's published cutoff defines that source only. Our LJ energy remains
an infinite Poisson/Bessel series. Independent sinusoidal-displacement
energies explicitly sum each atom's density before F, average site energies,
and refine displacement amplitude. This checks the Bloch factor-of-two and
coordinate conventions without assuming the covariance result.

Completed1024-frame comparison at the source's actual4.065Å/300K geometry:

| Plane-difference coordinate | Harmonic source RMS [m] | Actual MD RMS [m] | Variance difference |
|---|---:|---:|---:|
| normal |4.18743e-13|4.00226e-13|-8.65%|
| direct shear |8.90817e-13|1.00902e-12|+28.30%|
| transverse |8.90817e-13|9.64049e-13|+17.12%|

No scale factor was fitted. This supports the geometric/energy order of
magnitude, not exact harmonic agreement. Finite trajectory sampling,
anharmonicity and thermostat effects remain; individual mode/block results
are saved. Four-block variation is an empirical sensitivity, not a confidence
interval or box-size convergence proof. In particular this test does not
authorize treating a per-atom-cell energy as a finite nucleation event.

## 4. Why this trajectory does not supply an overdamped clock

For the reversible harmonic Smoluchowski model, whitened equilibrium
correlation is

    C0=kBT H^(-1),
    Cw(t)=C0^(-1/2) C(t) C0^(-1/2)=exp[-R_C t],
    R_C=H^(1/2) M H^(1/2)>0.

R_C has the same relaxation eigenvalues as sqrt(M) H sqrt(M), but the
eigenvectors/matrix entries need not coincide when H and M do not commute.
The covariance-whitening convention is not silently interchanged with the
mobility-whitening convention.

Every eigenvalue is positive. The measured1024-frame minimum is−0.594207
at0.150ps, versus an empirical four-block-plus-antisymmetry floor0.172295.
The negative mode remains in the first half, eight-block audit and sampling
stride2. This contradicts directly fitting these short-lag plane oscillations
as reversible harmonic overdamped a/s relaxation. Merely taking an absolute
correlation, discarding a negative eigenvalue, or fitting a preferred
exponential window would invent a clock. Longer-lag/coarse-grained validity
and coordinate equivalence need their own tests. This does not prove that no
future overdamped coarse-graining can work.

Actual source ps are physical simulation timestamps. They are **not** a
calibrated t0 for this solver. M_a,phys, M_s,phys and t0 remain null.

## 5. Actual dislocation-line drag coefficients

[Olmsted, Hector, Curtin and Clifton (2005)](https://doi.org/10.1088/0965-0393/13/3/007),
[author manuscript](https://arxiv.org/abs/cond-mat/0412324), Table3 reports
low-speed Al line drag from Ercolessi-Adams EAM MD, NOT our Mishin target:

| Character | B_line/T [Pa s/K], tau/T<.15MPa/K | B_line/T, tau/T<.30MPa/K |
|---|---:|---:|
| edge |3.7e-8|3.9e-8|
| screw |6.0e-8|7.5e-8|

At300K these give1.11/1.17e-5 Pa s and1.80/2.25e-5 Pa s. Different fit ranges
are alternatives, not confidence bounds. The API refuses extrapolation beyond
the stated temperature/stress ranges. The0K source lattice spacing is not
silently combined with finite-T velocity; v/b=tau/B is reported where useful.

[Gorman, Wood and Vreeland Jr (1969)](https://doi.org/10.1063/1.1657472),
[Caltech author copy](https://authors.library.caltech.edu/records/j5q43-k8847),
measures individual edge/mixed dislocation motion in99.999%Al[111] specimens
under microsecond torsional pulses. Published rounded endpoints are
1.5e-4 dyn s/cm² at123.15K and2.9e-4 at343.15K. Since1dyn/cm²=.1Pa,
these are1.5e-5 and2.9e-5 Pa s. The actual PDF checksum was verified.
No precise293K experimental value was invented from an un-digitized graph.

The coefficients and provenance are machine-readable in
`data/aluminum_line_drag_references.json`. This is meaningful physical kinetic
evidence. It still is not the missing collective opening mobility.

### Conditional metric, not automatic conversion

Line drag obeys B_line v=tau b with [B_line]=Pa s. For an explicitly derived
slip field and a HYPOTHETICAL uniform local friction density eta_s,

    Rayleigh dissipation/line length = (1/2) int eta_s (s_dot)^2 dx.

A rigidly translating core s(x-X(t)) would give

    I_profile=int (s')² dx [m],
    B_line=eta_s I_profile,
    zeta_cell=eta_s A_atomic_cell,
    M_s,cell=I_profile/(B_line A_atomic_cell) [m²/(J s)].

This derivation is unit-consistent and tested with known profiles. It requires
a condition-matched resolved core and validation that measured line drag can
be represented by this local dissipation. Far-field phonon dissipation is not
automatically local core friction. The helper marks its assumptions unvalidated
and does not return t0 or a normal mobility. No arbitrary line length or A_c
is inserted. Thus it opens a testable route, not a fabricated calibration.

## 6. Actual MPa stress-relaxation observations

[Verheyden, Deillon and Mortensen (2018)](https://doi.org/10.1016/j.dib.2018.11.047),
Data in Brief21,2134–2141, supplies raw microwire relaxation data and final
processed results. The83,134,534-byte supplementary archive was downloaded
and hashed. We read numeric text and cached XLSX values, never execute its
Mathematica notebooks/macros/formulas. The protocol is monotonic tension with
separate60s holds, NOT fatigue cycles.4N/5N as-cast and500°C/2h annealed
specimens are kept separate.

The raw audit parsed692 traces for29 table-matched specimens.23 reside in
source `Not considered/` folders. Neither692 nor669 should be called the
published686 count; filenames, selection flags and count mismatch remain
visible. The final-selection files supply643 apparent-activation rows for30
labels, including `5N_100_r_3`, which has no matching raw specimen table.
Final selection is source-reported, not independently re-derived here.

The `Numbers` sheet supplies analysis constants T=293K, b=.286nm,
kB=1.3806488e-23J/K. These are source analysis constants, not a new thermometer
or lattice measurement. Raw time is actual laboratory seconds, not PDE time.

For specimen Al_20_r_10: D14.7µm, axis[7,-1,-2], geometric maximum Schmid
factor0.4838498257. The raw axial stresses cover5.66888–18.92940MPa. Its eight
holds show measured shear relaxation; no empirical yield point is inferred
from this limited record. The exact unit audit is

    sigma_MPa=1000 F_mN/[pi D_um²/4].

This measured specimen cross section is used only to check the experimental
stress data. It is not A_c and does not enter atomistic work or strain.
The table/notebook diameter of5N_100_r_4 is114µm, whereas its raw force/stress
ratio implies112µm (a maximum1.00622MPa discrepancy using the table value).
It is not silently repaired. Duplicate DAQ timestamps and unnamed columns
are retained; no derivative is taken across a repeated timestamp.

## 7. Apparent activation sensitivity: a new scale constraint

Under the source's constant machine-stiffness reduction,

    a_app=(kBT/b) d ln(-tau_dot/M_machine)/d tau,
    V_app=b a_app=kBT d ln(rate)/d tau.

This is a stress derivative, NOT an invented characteristic volume, the local
opening coordinate a, geometric independent-region area or A_c. A varying
prefactor, mobile source population or internal stress can contribute to an
experimental apparent value; it is not automatically a single barrier's
activation volume. Negative signed inverse-error columns in the source are
preserved and not reported as positive standard deviations without explanation.

Five source-selected Al_20_r_10 values are1021.36–2988.90 b³. Their workbook
slopes5.90643–17.28456MPa^-1 are independently reproduced from T,b,kB with
relative error<4e-10. Source processed effective stress differs from workbook
tau0 by~0.420798MPa; that unverified correction is not silently identified with
the raw starting traction. Source rows remain validation references, not loss
terms or evidence that our model already predicts these experiments.

For our static interface-cell convention,

    G(q,T)=W(q)-A_atomic_cell L0 T.q,
    V_barrier,j=-d DeltaG/dT_j=A_atomic_cell L0(q_dagger,j-q_min,j).

The last equality follows from stationarity (envelope theorem), including
relaxed normal/transverse coordinates. It is not a finite-difference drift
law. The code differentiates positive verified minimum/index-one-saddle
branches and independently refines a traction finite difference.

Actual static mixed-load tests use Tn=−25/0/+25MPa, direct shear0/.5/1/2/4/6/
10MPa, source plusv14/v15 candidates:63 cases. At zero traction the SOURCE
per-cell derivative is0.321647 b³, v14 is0.307302 b³ and rejectedv15 is
0.411411 b³. Source barrier0.0762491eV/cell decreases only to0.0757781 at
10MPa. Its constant-prefactor-only rate slope at293K is~0.00187MPa^-1,
orders smaller than those five experimental apparent slopes.

This is a mechanism/normalization warning, not a direct293K-vs0K fit error.
Even the atomistically grounded source's uniform-interface cell is not a
finite collective dislocation event. The ratio V_app/b³ is not literally an
atom count; no factor of thousands is multiplied into W, kT or probability.
Resolved cores, finite-source geometry, obstacle/source populations and
condition-matched kinetics must generate the larger stress sensitivity.
The LJ/Bessel theory need not be discarded, but its collective state space
must represent that mechanism before claiming measured Al yield or fatigue.

## 8. Fixed-material validation at the ACTUAL MD box

Changing a lattice parameter during comparison must not reset the potential,
energy/length scales, density normalization or angular gauges. The MD source
box is4.065Angstrom, not the4.05Angstrom static-fit geometry. The new
`IsotropicBulkBasis` changes geometry by lambda=4.065/4.05 and retains every
material parameter and the original reference density.

At this cubic geometry x=rho_current/rho_ref need not equal1. For normalized
bond-density derivatives put R=sum(1-cos)Hess(x_bond), g=sum sin grad(x_bond).
The scalar embedding columns become

    H_A=-R/sqrt(x)+g g^T/(4 x^(3/2)),
    H_B=2R,
    H_C=4(x-1)R+2g g^T.

All full cubic STF moments vanish, but their displacement derivatives need
not. In particular the cross-invariant harmonic column is (x-1)H_D3, NOT
zero away from the fitted density. Its original-gauge reference limit is
exactly zero. Reciprocal density-tail errors and independent real-space
Hessian-tail bounds are propagated separately.

Four actual fixed candidates were checked at lambda=.99/1/4.065*4.05^-1/1.01,
radius12/16, MD plane modes and an off-axis static wedge. Independent actual
per-atom sinusoidal energies use every included neighbor BEFORE F and each
angular invariant. Refining amplitude4e-4/2e-4 and applying the stated
fourth-order Richardson comparison gives maximumH error1.304e-6eV/L0².
This is a derivative/finite-radius check, not a proof of whole-zone stability.

At300K and the actual MD box, R16 predicted RMS (m) is:

| Material | Normal | Either transverse polarization |
|---|---:|---:|
| source Al99 static |4.18743e-13|8.90817e-13|
| preserved v14 exact-bulk candidate |4.00332e-13|9.84199e-13|
| v15 quartic joint |3.82018e-13|9.47537e-13|
| v15 static-Bloch fit |4.12215e-13|1.15103e-12|
| v15 density/angular cross |4.15844e-13|8.88435e-13|
| actual MD |4.00226e-13|1.00902e-12 /9.64049e-13|

No atom count or temperature was fitted. The good v14 *total* variance is
not enough: individual normal modes differ by roughly-25% to+37%, with
large empirical block variations. Conjugate modes k and N-k are not treated
as independent observations. `fixed_material_modes/` saves every mode and
its block range, so aggregate cancellation cannot be mistaken for validation.
Finite-T anharmonic renormalization and longer sampling remain untested.

The periodic spatial sum of q_l vanishes identically, so a small spatially
averaged block-mean drift is NOT stationarity evidence. The repeated kinetic
audit now also stores each plane's mean drift. This changes no covariance or
mobility decision. The8192frame download twice failed; only the complete
1024frame study is reported. A128frame checkpoint is not a completed extension.

## 9. Independent core scales and a finite-source barrier

[Lu et al. (2000), PRB62,3099](https://doi.org/10.1103/PhysRevB.62.3099),
published TableII, was visually checked against the hashed publisher PDF.
Its Peierls stresses are **meV/Angstrom³**, so1 unit=160.2176634MPa.
DFT-input semidiscrete PN gives screw/edge256.35/3.204MPa; Ercolessi-Adams
input gives88.12/24.03MPa. These are static dislocation lattice resistance,
not measured yield, not a uniform interface spinodal, and not Mishin Al99.
The PN calculation specifies mu28.8GPa, nu.344 and b2.85Angstrom separately
from its LDA equilibrium geometry. Published width conventions are retained
without converting them into a fitted core cutoff or friction integral.
`aluminum_core_validation.json` records the8 character/method rows.

Our separate `FINITE_SOURCE_BARRIER_DERIVATION.md` derives minimum and
overhanging index-one saddle of the SAME anisotropic outer line energy.
For p=tau b, line stiffness T=gamma+gamma'', and
Q=gamma sin(theta)+gamma'cos(theta), the finite barrier is

    DeltaG=(2/p) int[theta0,pi-theta0] [Q(theta)-Q(theta0)] T sin(theta) dtheta,
    -dDeltaG/dtau=b Delta(swept area).

No guessed atomic multiplier, A_c or characteristic volume enters this
derivation.54 HYPOTHETICAL geometry/load cases were actually evaluated.
At r_core=2b,R=L, the shared0K bulk elasticity gives:

| Pin span [micrometre] | Critical shear [MPa] | Barrier at.95tau_c [eV] | Barrier at.999tau_c [eV] |
|---|---:|---:|---:|
|.2|59.1391|26.4882|.0700602|
|1|15.0788|168.843|.446583|
|5|3.66594|1026.22|2.71432|

Angular32/64 and quadrature refinement change these barriers by less than
2e-13 relative. Independent stress differentiation at two steps verifies the
envelope derivative; Richardson relative discrepancy is below5.2e-11 in the
recorded cases. Normal-variation stability is independently checked with
refined polygon energies. These are numerical checks of the OUTER line
model; core and nonlocal finite-part errors are not thereby bounded.

A micrometre source thus supplies an MPa force scale without reducing the
atomic potential. But at.999tau_c its whole-source stress derivative is
about3.03e5b³, much larger than the selected continuous-relaxation data.
The uniform cell (.3b³) and entire bowed source need not be that experiment's
operative activated segment. No source length was fitted to repair this
mismatch, and no barrier was converted into a probability or physical rate.

## 10. What is now known and what is not

Known: physical eV/Å/cell conversions; source plane atom counts; matched bulk
static calibration114/62/32GPa; published line-drag values with units;
actual experimental stress/relaxation and apparent derivative references.

Not established: a validated full interface/core energy, finite activated
source, production a/s mobility, physical PDE seconds/Hz, fatigue lifetime,
spatial initiation covariance or A_c. The numerical counterexamples and source
discrepancies are results, not missing numbers to be filled with guesses.
The original uncalibrated kinetic JSON remains unchanged. Existing production
UI still displays model time. Source data may legitimately show physical ps/s,
but their captions explicitly identify them as source observations.

Reproduction tools: `public_aluminum_kinetics`, `run_public_kinetic_audit`,
`public_wire_validation`, `public_activation_reference`,
`run_plane_covariance_validation`, `run_material_plane_covariance`,
`run_activation_scale_comparison`, `public_core_reference`,
`run_finite_source_barrier`.
Every completed study has its own source hashes, units, scope and actual output.

## 11. Lower-frequency mobility was tested, not ruled out from oscillation alone

`ZERO_FREQUENCY_MOBILITY_AUDIT.md` derives Gamma0=kBT C0^-1 K C0^-1,
K=integral C(t)dt, for a matched harmonic memory projection. A coupled
oscillatory control verifies that short-lag negative correlation does not
imply negative zero-frequency friction. Actual Al1024frame data, however,
have cutoff-dependent signed integrals below block/quadrature sensitivity.
At3.15ps the whitened integral eigenvalues are-.01078,-.00676,-.00419ps
against.03825ps empirical floor. No cutoff is chosen to manufacture M>0.
Longer source data and thermostat/coordinate mapping are still required;
production physical seconds/Hz remain unavailable.

## 12. Explicit dimensional ledger and final source-download status

`final_scaling_report/atomic_unit_ledger.csv` derives every conversion from
a_lat=4.05Angstrom and E0=1eV, without adjusting a strength or kinetic scale:

| Quantity | Value | Meaning |
|---|---:|---|
| L0=b |2.8637824638e-10m|fixed atomistic coordinate scale|
| h111 |2.3382685902e-10m|perfect plane spacing|
| A_atomic_cell |7.1024908428e-20m²|one crystallographic interface cell, not A_c|
| 1eV/cell |2.2557954237J/m²|interface-energy conversion|
| E0/(A_atomic_cell L0) |7.8769789682GPa|generalized-force unit, NOT yield strength|
| 1MPa traction |.0001269522242eV per reduced coordinate|conjugate interface force|
| E0/L0² |1.9535761427J/m²|physical coordinate Hessian unit|
| L0²/E0 |.51188176297m²/J|mobility conversion still divided by unknown t0|

The crystallographic per-atom volume A_atomic_cell*h111=a_lat³/4 is used to
convert a bulk energy curvature into GPa. It is neither a fitted activation
volume nor a specimen correlation parameter. A numerical force unit in GPa
does not establish a GPa experimental yield stress. Production's reduced
kappa*sigma/E bridge is unchanged; this table audits the separate physical
active-interface work convention.

There were three unsuccessful attempts to extend the MD sample: two8192frame
attempts and one2048frame attempt. The last failed before receiving data.
`trajectory_expansion_attempts.json` preserves these outcomes. Only the
completed1024frame/25.575ps projection is used for the reported kinetics.
Its small projected array and hashes are saved for offline replay; no partial
download is presented as a completed longer trajectory or a verified full
multi-gigabyte archive checksum.

The final direction-specific297-state static study and84-state saddle study
are recorded separately. Explicit tensor/normal/shear support in research
does not imply production tensor loading, a3D probability PDE, calibrated
experimental yield or physical PDE seconds. All kinetic, spatial-correlation
and production-adoption gates remain explicit.

## 13. Why a measured plane mobility cannot be fixed by one hidden scale factor

Even with a converged source-coordinate friction, the energy normalization
has to follow the same projection. For an explicitly coherent Np-atom mode,
suppose its derived energy is E_N(q)=Np E0 G*(q/L0). This is a conditional
normalization example, not an assertion that our full periodic plane process
has an independent two-coordinate per-cell PMF. If its physical coordinate
mobility is M_N, the projected physical Smoluchowski flux contains E_N and
the physical thermal energy kBT. After time/coordinate rescaling,

    M_star = t0 Np E0 M_N/L0²,
    kT_star_effective = kBT/(Np E0),
    M_star*kT_star_effective = t0 M_N kBT/L0².

Both drift and diffusion then transform consistently. Alternatively choose
the total-energy scale Np E0 rather than per-cell E0; the same physics must
result. Changing only mobility while leaving a mismatched per-cell thermal
energy cannot repair a wrong equilibrium covariance. The two are not
independently adjustable scale factors.

For the actual source plane count576 and300K, kBT/E0=.02585199979 and the
coherent per-cell thermal factor would be4.48819441e-5. This count follows
the source geometry, NOT A_c, a fitted area, or the number of independent
initiation regions. Neither number is inserted into the production PDE:
its state and energy are not that coherent periodic-plane mode. The actual
source energy also couples all plane classes through H(k), as derived in
section3. A matched PMF and memory/coordinate elimination remain necessary.

This identifies a concrete normalization obstacle rather than treating
unknown kinetics as a single missing arbitrary constant. The original
production dimensionalization and default kT are preserved, and no physical
temperature, atom count or mobility is silently retuned.

## 14. The source thermodynamic log was also checked

The checksum-verified published `log.lammps` distinguishes equilibration
from sampling. Its sampling header explicitly uses LAMMPS metal units:
time ps, temperature K, pressure bar (1bar=.1MPa), volume Angstrom³,
energies eV. Six logged snapshots fall within the1024frame coordinate
projection: T=293.51939–303.47586K, arithmetic snapshot mean299.10801K,
pressure mean+8.46413MPa (compression). These six snapshots are not claimed
to be the complete time-average temperature of1024 coordinate frames.

The whole500ps log contains101 such snapshots: mean299.93018K and
pressure+9.83225MPa. A complete scalar thermodynamic log is NOT a completed
500ps coordinate analysis. The public NVT box is not a zero-pressure
ensemble. Neither its pressure nor its source atomic mass is inserted into
the production stress bridge or overdamped time scale. Setpoint300K remains
explicitly distinguished from the logged temperatures; no fit or post-hoc
temperature adjustment is made. The data and source SHA are saved under
`source_thermodynamic_audit/`.
