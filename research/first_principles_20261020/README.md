# Isolated PART I–III verification

These tests concern mathematical identities, counterexamples and controlled
limits, using illustrative dimensionless parameters. They do not calibrate
Al or alter production dynamics. `foundations.py` derives row/plane kernels
independently; the direct sum reference does not call production forces.
One separately labelled row point checks the production reference implementation.

Requirements: Python 3.12, NumPy, SciPy. Executed versions are in the result.
From the repository root run:

```powershell
python -B -X utf8 research/first_principles_20261020/verify_foundations.py --output research/first_principles_20261020/results/my_new_run
```

The output directory must not already exist. Saved results are immutable
receipts, not overwritten baselines. `run_001/verification.json` records
8 passed/0 failed, actual truncations, tolerances, errors and limitations.
The owned process exited 0 and was reaped; no timeout occurred. No full
production solver or application regression suite was run for these new
documents and isolated controls.

The finite Hamiltonian bath uses the full nonlinear normal Bessel force,
deterministic initial conditions and normal cyclic forcing. FDT covariance
is computed analytically; no Monte Carlo fatigue counting is substituted.
Only the explicitly labelled OU limit control uses the local Bessel Hessian.
Its convergence illustration is not a nonlinear limit theorem or material
time-scale measurement. Real/reciprocal refinement at machine precision is
not a proof of a zero omitted tail. Plane analytic bounds apply to the
energy channel as stated in the output, not automatically to the Hessian.

See the [derivations and scope](../../docs/research/20261020_first_principles_i_iii/README.md).
