# Theory Core v1 solver

This package is the executable proof-of-principle solver for the paper's correlated configurational fatigue-initiation framework.

## What it solves

The state of one correlated representative region is

\[
q=(a_1,\ldots,a_N,s_1,\ldots,s_N),
\]

with interacting Lennard-Jones geometry, cyclic normal loading, finite configurational mobility, and thermal noise. The code integrates the coupled overdamped Langevin system, equivalent to the many-body Smoluchowski equation, without imposing a product closure \(P_N=\prod_i P_i\).

The user-facing outputs are macroscopic: axial strain and stress-strain hysteresis, configurational well-crossing activity, local opening-barrier evolution, survival, and first-passage fraction. The microscopic coordinates remain hidden internal states.

## Important scope

This is a **dimensionless mechanism solver**, not a calibrated aluminum fatigue-life predictor.

- The default LJ parameters are an analytical interaction prototype.
- `chi_axial_projection` is an explicit geometry/coarse-graining bridge and must ultimately be derived or calibrated from a physically selected configurational mode.
- The operational crack first-passage test uses the softest opening mode of the
  full many-body \(H_{aa}\) block at frozen \(\mathbf s\). A trajectory is
  absorbed after it enters the negative-curvature side with outward
  deterministic drift. The precomputed one-cell saddle/barrier remains a
  diagnostic only; it is no longer the absorbing boundary. A future full
  multidimensional minimum-energy-path calculation can refine this coupled
  dividing surface.
- Specimen-scale characteristic length/area aggregation is deliberately deferred.
- EAM/MEAM or validated Al energy data can replace the LJ layer without changing the stochastic/state architecture.

## Run

From the repository root:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m solver_v1.run_demo
```

Outputs are written to `solver_v1/output/`.

The default demo tests the required ordering: low-load elastic survival; intermediate configurational well crossing without crack initiation; larger-load opening-barrier reduction and nonzero first passage; and high-load first-passage accumulation.

## Theory-to-code map

- `model.py`: two-row LJ geometry, correlated interaction energy, macroscopic strain mapping, periodic well index \(s_i=b n_i+\xi_i\), and stable/opening-saddle lookup.
- `solver.py`: cyclic loading, Euler-Maruyama integration of the correlated state,
  strain-component output, and a full correlated opening-mode first-passage check.
- `run_demo.py`: reproducible screening run with CSV/JSON summaries and figures.

## Current screening parameters

- \(N=3\)
- \(b=1\)
- \(\epsilon_{\rm LJ}=1\)
- \(\sigma_{\rm LJ}=0.82\)
- \(k_BT=0.009\)
- \(M_a=1\)
- \(M_s=0.15\)
- \(\chi=0.40\)
- 32 trajectories
- 10 cycles
- \(\Delta t=0.02\)

These are **not** claimed as pure-Al material parameters.
