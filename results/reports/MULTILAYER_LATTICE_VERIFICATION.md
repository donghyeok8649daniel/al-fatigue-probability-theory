# Multilayer lattice verification

Verified locally on 2026-09-02 for the active intrinsic potential

\[
U_0(a,s)=\sum_{k=1}^{\infty}W(ka,s).
\]

The counting is deliberately unweighted: there is no multiplicity prefactor
`k`, and every layer uses the same collective slip coordinate `s`, never `ks`
or `js`. External tensile work is excluded from `U0` and enters only through
the generalized forces in the probability currents.

## Independent checks

The script `simulations/verify_multilayer_lattice.py` compares the direct
double sum with the Fourier--Bessel/Bessel--Lambert representation for
`q = 6, 12` at three registries and three normal spacings. It independently
generates the half-integer Bessel coefficients with the factorial formula,
then verifies the resulting 12--6 polylogarithm closure.

| Check | Maximum relative/absolute error |
|---|---:|
| Direct double sum versus Bessel--Lambert | `1.1111094994615028e-11` |
| Bessel--Lambert versus 12--6 polylog form | `0.0` in float evaluation |
| Individual Lambert sum versus polylog form | `4.327421687874723e-16` |
| Analytic `dU0/da` versus direct differentiated sum | `6.297176770730253e-11` |
| Analytic `dU0/da` versus finite difference | `1.2345500989242514e-10` |
| Analytic `dU0/ds` versus finite difference | `9.596594027689862e-09` |
| Periodicity, even symmetry, and reference slip identities | `0.0` |

The independently obtained polynomial coefficients are

- `K_(5/2)`: `(1, 3, 3)`;
- `K_(11/2)`: `(1, 15, 105, 420, 945, 945)`.

For the dimensionless force-balance example `Q_a = 0.1`, the root finder
identifies a stable point at `a = 0.9056850907613117` with positive curvature
`112.11153147324117`, and an outer metastable barrier at
`a = 2.0486466236634353` with negative curvature
`-0.2917046732896391`. These values verify classification by curvature; they
are not an aluminum calibration.

## Regression and document build

- Full test suite: `167 passed in 27.91 s` with `MPLBACKEND=Agg`.
- Canonical PDF: `output/pdf/slip_lattice_energy_derivation.pdf`.
- PDF SHA-256: `e433296316a04a1ae4239da4ce323476c81feb1e8559e304789fd00bdf9fe23a`.
- PDF length: 11 pages (`137039` bytes).
- Local Tectonic 0.17.0 log scan: no LaTeX error, missing-character,
  overfull/underfull-box, or unresolved-reference messages.
- `.github/workflows/**`: unchanged; no persistent automation was introduced.

Machine-readable values are stored in
`results/data/registry_plasticity/multilayer_verification_summary.json`.
