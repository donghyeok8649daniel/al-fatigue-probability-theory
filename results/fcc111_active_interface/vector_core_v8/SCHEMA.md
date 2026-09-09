# Static vector-row results: interpretation contract

All material parameters are the unchanged `angular_monotone_opening` research
candidate, SHA256
`9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655`.
Neither this dataset nor its parent candidate is a calibrated Al fatigue model.

- `u,v,w_over_L0`: local displacement of an infinite e1 atomic row, with
  L0=2.863782463805517e-10m. Affine gamma*xz work is separate. All local
  components are relaxed, but transverse periodic cell vectors stay fixed.
- `energy_ev_per_line_repeat`: eV per b*L0 infinite-line repeat. Conversion
  to J/m divides by that repeat length. It is NOT a finite activation energy.
- `applied_shear_mpa`: imposed resolved shear traction, not nominal axial
  stress. The actual geometric volume in E-tau*V*gamma contains no A_c.
- `iteration`: numerical optimizer iteration, NOT model or physical time.
  `elapsed_seconds` measures computer execution only.
- `maximum_force_residual`: maximum component of actual three-direction
  energy gradient [eV/L0] at the declared transverse ring.
- `minimum_nongauge_curvature`: static curvature after translation removal,
  NOT a vibration frequency. Affine gamma uses h*sqrt(N) numerical scaling.
- `positive_cores/negative_cores/core_positions`: **x phase winding only**.
  Zero x winding does NOT imply a perfect three-component vector state.
- Raw case columns `registry_shear`, `intrabond_shear`, and their differences
  are the old x-projected kinematic decomposition. They are NOT a definition
  of vector plastic strain, notably at a vector minimum on x=b/2. Consolidated
  registry columns use `x_projected_` prefixes. No production PDE convention
  has been changed.
- `layer_registry`: actual vector additional shifts across adjacent layers,
  normal change, row dispersion, scalar partition margin and distances to
  0,tau,2tau modulo the triangular primitive lattice. These are geometry
  diagnostics, not automatic calibrated stacking-fault/core labels.
- `*_tail`: frozen-state changes with retained transverse neighborhood.
  LJ zero-mode bounds apply ONLY to that pair zero mode. Reference ring16/28
  is not an all-channel mathematical infinite-tail certificate.
- `perfect_vector_spectrum_tail`: independent 3x3 Fourier assembly checked
  against the actual Hessian product. Minimum eigenvalue stabilization alone
  can conceal polarization crossing; full matrix errors are also provided.
- `summary.completed=false`: original interrupted L-BFGS calculation, not an
  accepted equilibrium. Explicit `*_newton` continuations retain the source
  state and same energy/force tolerances. No incomplete case was discarded.
- The close-seed16x16 case is already perfect after the scalar seed relaxation;
  it is not the same initial core as the8x8 vector-fault state. Only the primary
  24/32 states share the original1.98408669nm core separation.

Static unloading is NOT a physical-time relaxation hold. Physical mobility,
seconds/Hz, finite-line/source activation, experimental yield, A_c, and
production/UI certification remain unavailable/unvalidated.
