# Theta-based density reconstruction — numerical verification

This report checks the exact 1D smooth-moment shape identity against the nonlinear conservative layer-LJ chain.
No Boltzmann, Gaussian/Weibull physical PDF, Fokker–Planck, damping, Markov, or neighbor-independence assumption is used.
A Gaussian kernel is used only as a finite-sample estimator of the empirical density and conditional moments.

## Fixed numerical protocol

- atoms: 512
- dt: 0.01
- force amplitude: 0.03
- omega: 0.02
- time derivative half-window: 0.00025 cycle
- bandwidth: 1.5 × Silverman rule, fixed for all phases
- boundary spacings excluded from the bulk acceleration identity

## Results

| phase N | L1(P) | KS(P) | relative mean-energy error |
|---:|---:|---:|---:|
| 2.05 | 0.01131 | 0.00279 | 0.150% |
| 2.10 | 0.07626 | 0.03813 | 3.046% |
| 2.25 | 0.01890 | 0.00836 | 1.406% |
| 2.40 | 0.01757 | 0.00671 | 0.050% |
| 2.50 | 0.01858 | 0.00919 | 1.733% |
| 2.60 | 0.00567 | 0.00231 | 0.963% |
| 2.75 | 0.01073 | 0.00394 | 0.010% |
| 2.90 | 0.00322 | 0.00158 | 0.246% |

Maximum L1 density error: **0.07626**.
Maximum KS density error: **0.03813**.
Maximum relative error in mean intrinsic energy on the same smoothed support: **3.046%**.

## Interpretation

The 1D Theta-based shape identity is numerically consistent with the deterministic nonlinear LJ-chain data at the tested phases.
The test does **not** prove a new closure: Theta, conditional acceleration, and the material acceleration are measured from the chain state.
The next derivation step is the two-coordinate (a,s) conditional velocity-covariance tensor and its integrability/compatibility condition.

## Assumption / validity ledger

- Exact: finite-M mechanics and the interior nonlinear LJ spacing acceleration.
- Exact at smooth-moment level: continuity and first/second velocity moment identities.
- Numerical approximation only: KDE/regression used to estimate smooth conditional fields from finite M.
- Not used: Boltzmann equilibrium, Gaussian physical state, Fokker–Planck, Markov bath, independent spacings.
- Invalid regime for the divided shape law: Theta = 0.
