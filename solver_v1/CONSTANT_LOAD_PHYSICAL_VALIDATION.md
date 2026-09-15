# Constant load: model consistency is not experimental agreement

The 1000 MPa constant-load audit in results/constant_hold_audit is a production
TwoRowLJ **model-time** calculation. Its outward flux is compatible with an
absorbing Smoluchowski process. This does not establish compatibility with a
particular physical aluminum specimen. The computed signal is also below the
required conservative probability resolution criterion.

## Scope correction: single crystal first

The user explicitly confirms SINGLE-CRYSTAL scope. No polycrystal model,
grain-boundary law or grain-size calibration is introduced. The production
TwoRowLJ is an idealized reduced row/registry model, not the full FCC specimen
or a complete dislocation-bearing single crystal. The separate full-FCC and
dislocation research modules do not automatically become the production PDE.

Distinguish a perfect/defect-free single crystal from a real single crystal
containing dislocations and sources. Missing grain boundaries are NOT a defect
of a model deliberately intended for single-crystal conditions. The relevant
question is whether the required single-crystal defect mechanisms, loading
orientation and collective-coordinate kinetics are represented.

An explicitly single-crystal reference is A. J. Kennedy and N. C. McGill,
*Orientation-dependence of Fatigue-accelerated Creep in Aluminium Single
Crystals*, Nature 207, 1086–1087 (1965), DOI 10.1038/2071086a0:
https://www.nature.com/articles/2071086a0 . The publisher abstract confirms
single-crystal experiments and orientation dependence. The full quantitative
data are subscription content and were not accessed; no stress, temperature,
rate or fitting target is invented from the abstract. It is a source lead,
not a completed calibration. The background studies below are not substitutes.

## Background references ONLY (not direct single-crystal calibration targets)

1. J. Shen, S. Yamasaki, K.-i. Ikeda, S. Hata, H. Nakashima (2011),
   *Low-Temperature Creep at Ultra-Low Strain Rates in Pure Aluminum Studied
   by a Helicoid Spring Specimen Technique*, Materials Transactions 52,
   1381–1387. DOI: https://doi.org/10.2320/matertrans.M2010405
   Publisher: https://www.jstage.jst.go.jp/article/matertrans/52/7/52_M2010405/_article

   The experiment covers 0.32–0.43 Tm and strain rates below 1e-10 s^-1.
   Reported microstructure dependence includes 5N Al with grain sizes above
   1600 micrometers versus 24 micrometers, and 2N Al at 25 micrometers.
   The abstract reports different stress exponents and activation energies
   across these conditions. A single ideal registry coordinate lacks these
   explicit grain-boundary and dislocation-population variables. No reported
   exponent is inserted as an empirical production creep law.

2. H. Nakashima, H. Yoshinaga (1987), *Transient Creep Mechanism in Pure
   Aluminum at High-Temperature*, Transactions of the Japan Institute of
   Metals 28, 644–654. DOI: https://doi.org/10.2320/matertrans1960.28.644
   Publisher: https://www.jstage.jst.go.jp/article/matertrans1960/28/8/28_8_644/_article

   The reported tests use 623–823 K and 0.81–6.7 MPa, measuring transient
   plastic strain and strain rate. Their discussion distinguishes early glide
   contributions and later recovery. These high-temperature experiments are
   not a room-temperature 1000 MPa target and are not mobility data for the
   current local a/s coordinates.

Both publisher abstracts were read on 2026-09-15. Full time-series digitization,
uncertainty extraction and a condition-matched fit have NOT been performed.
We do not use the search engine's crawl date as the publication year.

## Required comparison, in order

Define single-crystal purity, loading/slip orientation, initial dislocation
structure, temperature, geometry, loading mode and measured quantity first. A spring
shear creep measurement, uniaxial creep strain, rupture time and local crack
initiation probability are different observables. A load-controlled test also
requires the distinction between engineering and true stress during changing
section area. No universal numerical yield strength is assumed for all Al.

For a mapped strain history, a possible residual is

    r_k = [epsilon_model(t_k) - epsilon_experiment(t_k)] / uncertainty_k.

This residual is currently **undefined**, not zero: t_model to t_seconds is
unavailable and the experimental state has not been mapped to the production
collective variables. Fitting a mobility to rupture lifetime would not supply
that missing mapping. Nor is a dimensionless probability comparable to strain.

1. Validate the same bulk/interface elastic and slip/opening energy landscape
   against atomistic targets without changing the LJ/Bessel premise.
2. Represent and independently constrain the defect/source and boundary
   mechanisms required by the chosen specimen; preserve local energy and
   nonlocal elastic consistency instead of scaling ideal strength by hand.
3. Validate collective-coordinate kinetics at the correct temperature and
   normalization; only then enable seconds and compare strain-rate histories.
4. Require spatial/time/domain convergence and unloading controls. Compare
   creep deformation separately from opening loss and specimen rupture.
5. Only after the above, validate independent held-out loads/temperatures.

## Present conclusion

The production result is **not experimentally validated creep**. A mechanically
nearly stationary ensemble with an extremely small cumulative boundary loss
is possible within the equations, but is not a prediction that a macroscopic
real Al specimen can sustain 1000 MPa with negligible plastic deformation.
The microscopic/mesoscopic-to-specimen bridge and kinetics are substantive
missing physics, not a plot-scale problem. No probability multiplier, arbitrary
stress rescaling, mobility change, or A_c adjustment is made to hide this.
