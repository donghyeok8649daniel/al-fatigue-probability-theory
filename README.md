# theory-core-v1-solver

Integration branch for the finalized Theory Core v1 manuscript and its executable proof-of-principle solver.

## Canonical manuscript

- Entry point: `paper/main.tex`
- Canonical body: `paper/fatigue_probability_main_v4.tex`
- Historical v1/v2/v3 manuscript bodies are preserved.
- v3 established the no-go result for strict pure-normal cycle accumulation.
- v4 adds the required slow correlated configurational state and promotes the theory to:
  1. interacting LJ/configurational energy;
  2. many-body correlated Smoluchowski dynamics without product closure;
  3. plastic memory through `s_i = b n_i + xi_i`;
  4. separate plastic and normal-opening stability criteria;
  5. mechanically derived opening saddle and first-passage crack initiation.

The macroscopic target remains stress/strain hysteresis, residual plastic response, irreversible evolution, and crack-initiation probability. Microscopic coordinates are hidden constitutive states rather than proposed experimental observables.

## Solver

`solver_v1/` contains the executable dimensionless mechanism solver.

- `model.py`: two-row LJ geometry, correlated interaction, strain bridge, well index, normal opening saddle/barrier.
- `solver.py`: cyclic load and Euler-Maruyama integration of the full correlated state.
- `run_demo.py`: reproducible four-load screening.
- `test_solver.py`: gradient, periodicity, and barrier-ordering checks.
- `output/summary.csv`: committed numerical screening provenance.

Run from repository root:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m solver_v1.run_demo
```

## Current mechanism-screening result

One dimensionless parameter set produces the required ordering:

| Fmax | first passage | max mean abs well index | min mean opening barrier |
|---:|---:|---:|---:|
| 2.5 | 0 | 0 | 0.3549 |
| 3.2 | 0 | 0.3021 | 0.1062 |
| 3.4 | 0.09375 | 3.1839 | 0.0609 |
| 3.6 | 0.84375 | 6.8000 | 0.0271 |

Interpretation: configurational well crossing can occur before crack initiation; further rearrangement lowers the normal opening barrier and increases first-passage probability.

## Scientific scope

This branch does **not** claim calibrated pure-Al fatigue life. The current LJ parameters, mobility values, thermal scale, and axial projection coefficient are mechanism-screening parameters. Quantitative aluminum work still requires an Al-specific EAM/MEAM or validated energy landscape, mobility calibration, a derived/configurational geometry bridge, and later characteristic correlation length/area for specimen-scale aggregation.

The fixed v1 architecture is therefore:

`interaction energy -> correlated state dynamics -> plastic/configurational memory -> opening stability/barrier -> first passage`.
