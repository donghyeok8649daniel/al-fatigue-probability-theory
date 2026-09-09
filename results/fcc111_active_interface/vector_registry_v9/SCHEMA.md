# vector_registry_v9 — static full-vector registry audit

Energy: eV per primitive interface cell; surface energy: J/m2.
Coordinates: (a,ux,uy) / L0, L0=4.05/sqrt(2) angstrom.
Hessian: eV/L0^2; forces: eV/L0. MPa tractions are conjugate through the
actual atomic interface area, not through a statistical correlation area.
All calculations are STATIC. No probability/time/Hz/yield certification.

- `run_manifest.json`: input/source/code hashes and starting execution record.
  Its completed=false is the immutable START marker; consult scientific_status
  for final completion. parameter_metadata.path is the historical scalar
  constructor path; the evaluated registry space is fully two-dimensional.
- `scientific_status.json`: completed primary execution, unchanged parameter
  hash, numerical/scientific status, counts and actual wall time.
- `stationary_points.csv`: fixed-a special states, full stationary roots and
  the two downhill endpoint minima. `gradient_norm` is raw gradient, NOT force
  equilibrium residual under nonzero applied traction.
- `barrier_comparison.csv`: connected pristine->intrinsic-fault saddle and
  reverse barrier after opening relaxation. Not finite-loop activation energy.
- `relaxed_registry_branches.csv`: locally a/y-equilibrated fixed-x branch,
  41/81 samples per model. Schur curvature distinguishes stress stability;
  not a proof of the global minimum-energy path.
- `uniform_interface_folds.csv`: first bracketed/refined x-traction fold,
  full zero Hessian mode with positive eliminated subblock. NOT Al yield.
- `constraint_ablation.csv`: same energy with a/y fixed, only a free, both
  free. Nonzero holding tractions are recorded. Only the both-free row is a
  full-coordinate spinodal; the other maxima belong to constrained problems.
- `stress_and_unload.csv`: actual 32 low-MPa normal/two-shear static states.
  This is not cyclic time evolution. Tiny unload displacements are roundoff,
  not newly created/residual plasticity.
- `vector_energy_grid.csv`: full21x21 triangular primitive cell for each
  model, at its fixed bulk spacing. Negative curvature at a barrier does not
  imply instability of the pristine state. Maxima on this grid are samples,
  not a global-extremum proof.
- `direct_lattice_validation.csv`: independent finite real-space/layer
  convergence against the infinite model. Finite algebraic tails remain.
- `derivative_validation.csv`: all originally run stencil results, including
  the source spline-knot crossing error. Coarser failures are not deleted.
- `source_stencil_refinement.csv`: finer source/analytic checks showing the
  error collapse and eventual finite-difference cancellation floor.
- `reciprocal_neighborhood_refinement.csv`:2e-11/2e-13 comparison at stationary
  states. The small differences do not bound material/source error.
- `rechecks_status.json`: completed independent constraint/stencil experiments.
- `verification.json`: final tests, protected-file hashes, result/source binding.
  The original full-suite SOURCE NaN failure is retained. Its cause has NOT
  been identified despite successful reruns; it is not labeled a proved fix.
- `postguard_data_revalidation.json`: fresh calculation of58 special/load/fold
  states and882 grid points with final fail-fast guards, plus explicit rejection
  of nonfinite jets and flat zero-force states. Old/current differences were0.
- PNG files: plots of actual computed samples. No invented spatial fields.

`analytic_candidate` is the preserved angular_monotone_opening candidate,
not a newly calibrated Al surface. `Mishin_source_only` uses the checksum-
verified NIST Al99 file on matching rigid geometry, solely as an atomistic
benchmark. This is not experimental truth or a production potential option.

`A_c`, kinetic mobilities, PDE probability semantics and production UI are
unchanged. Read solver_v1/VECTOR_REGISTRY_AND_STRENGTH_AUDIT.md for equations,
the constrained-path caveat, numbers, and remaining scientific gates.
