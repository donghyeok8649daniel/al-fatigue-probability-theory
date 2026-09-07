# Pure-aluminum calibration sources

This source ledger separates electronic/0 K calculations, finite-temperature
atomistic values, room-temperature experiments, and later compilations. No
fatigue-life or empirical yield-stress datum is used.

## Mishin et al. Al EAM reference

Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos,
"Interatomic potentials for monoatomic metals from experimental data and ab
initio calculations," *Physical Review B* **59**, 3393--3407 (1999),
[doi:10.1103/PhysRevB.59.3393](https://doi.org/10.1103/PhysRevB.59.3393).

The [NIST Interatomic Potentials Repository entry](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/)
provides the original functions and stable model identifier
`1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al`. Its release
check reproduces the 4.05 angstrom minimum and -3.36 eV/atom cohesive energy.
The 0 K EAM validation values used here are explicitly model-reference values,
not experimental truth.

## Lu et al. DFT generalized stacking-fault surface

G. Lu, N. Kioussis, V. V. Bulatov, and E. Kaxiras,
"Generalized-stacking-fault energy surface and dislocation properties of
aluminum," *Physical Review B* **62**, 3099--3108 (2000),
[doi:10.1103/PhysRevB.62.3099](https://doi.org/10.1103/PhysRevB.62.3099).

The LDA-DFT calculation reports 3.94 angstrom and 82.52 GPa, intrinsic fault
energy 0.164 J/m^2, the relaxed lowest-path unstable fault 0.224 J/m^2, and the
direct `a/4<101>` unstable point 0.250 J/m^2. The last is the closest direct
mapping to the present symmetric full-Burgers one-coordinate registry path.

## Wei et al. clean Al(111) surface

X. Wei, C. Dong, Z. Chen, K. Xiao, and X. Li, "A DFT study of the adsorption
of O2 and H2O on Al(111) surfaces," *RSC Advances* **6**, 56303--56312 (2016),
[doi:10.1039/C6RA08958E](https://doi.org/10.1039/C6RA08958E).

The clean six-layer Al(111) DFT surface energy is 1.06 J/m^2. The paper also
quotes a 1.14 J/m^2 experimental comparison. Only the electronic DFT value is
used for the selected 0 K-like static fitting set; the room-temperature number
is retained separately.

## Finite-temperature cohesive/elastic context

V. Yamakov, E. Saether, D. R. Phillips, and E. H. Glaessgen,
"Molecular-dynamics simulation-based cohesive zone representation of
intergranular fracture processes in aluminum," *Journal of the Mechanics and
Physics of Solids* **54**, 1899--1928 (2006),
[doi:10.1016/j.jmps.2006.03.004](https://doi.org/10.1016/j.jmps.2006.03.004).

This gives 100 K Mishin-EAM elastic constants and an orientation-specific
surface energy. They are held out because they are neither 0 K nor clean
Al(111) values.

## Modern comparison ranges

Y. Choi and T. Brink, "Faceting transition in aluminum as a grain boundary
phase transition," *Physical Review Materials* **9**, 083607 (2025),
[doi:10.1103/2dnf-zdz8](https://doi.org/10.1103/2dnf-zdz8).

Its potential-comparison table supplies a useful provenance trail and broad
reference ranges. These compiled values are retained as uncertainty context;
they are not inserted into the selected loss without tracing their original
temperature and method.

The machine-readable values, conditions, mappings, and fit/held-out flags are
in `data/aluminum_reference_targets.csv`.
