# Pure-normal frequency-scaling no-go

This report accompanies `simulations/audit_pure_normal_frequency_no_go.py`.

## Analytical result

For a phase-controlled periodic stress waveform in the strict quasistatic reduced model,

$$
\Lambda_f(t)=\Lambda_*(ft),
$$

and the instantaneous thermal escape rate has no explicit dependence on loading rate. Therefore

$$
\mathcal H_f
=\int_0^{1/f}k_c[\Lambda_f(t),T]dt
=\frac1f\int_0^1k_c[\Lambda_*(\theta),T]d\theta.
$$

Hence

$$
\boxed{\mathcal H_f\propto f^{-1}},
$$

which implies, in the rare-event narrow-state limit,

$$
\boxed{N_{50}\propto f}
$$

while median elapsed time is approximately frequency independent.

This scaling is not a fit and does not depend on the chosen value of $A_c$.

## Historical aluminum frequency interval

Daniels and Dorn reported that room-temperature fatigue strength of high-purity aluminum was insensitive to test frequency over 25--1440 cycles/min. Those endpoints are approximately 0.4167 Hz and 24 Hz.

The strict fast-equilibrium normal model predicts

$$
\frac{\mathcal H(0.4167\ \mathrm{Hz})}{\mathcal H(24\ \mathrm{Hz})}
=57.6.
$$

So, at equal stress waveform and equal number of cycles, it predicts about 57.6 times more local cumulative hazard at the lower frequency. That is a strong mechanism-level tension with a broad frequency-insensitive room-temperature fatigue response.

The historical experiment is not an exact single-crystal match, so this is a falsification warning, not a standalone proof that every normal-instability contribution is absent.

## Additional single-crystal mechanism evidence

Room-temperature aluminum single-crystal literature reports persistent-slip-band-related structure, secondary slip, net irreversible slip, subgrain development, and strong slip-system/orientation dependence of initiation life. These observations identify a cycle-evolving microstructural mechanism that the strict pure-normal state does not carry.

Relevant sources:

- N. H. G. Daniels and J. E. Dorn, “The Effect of Temperature, Frequency, and Grain Size on the Fatigue Properties of High-Purity Aluminum,” ASTM STP 196, p. 94 (1957), DOI `10.1520/STP19619570007`.
- T. Zhai, G. A. D. Briggs, and J. W. Martin, “Fatigue damage at room temperature in aluminium single crystals—IV. Secondary slip,” *Acta Materialia* 44(9), 3489--3496 (1996), DOI `10.1016/1359-6454(96)00025-0`.
- M. Hayashi, “Effect of crystal orientation on fatigue crack initiation life in pure aluminum single crystals,” *International Journal of Fatigue* 156, 106661 (2022), DOI `10.1016/j.ijfatigue.2021.106661`.

## Final numerical verdict

The final reduced thermal survivor equation remains a valid **mathematical normal-instability submodel under its assumptions**, but it cannot be promoted to a complete room-temperature pure-aluminum fatigue mechanism merely because it is closed.

Under the strict constraints

$$
\text{pure normal}+\text{quasistatic reversible mechanics}+\text{no slow state},
$$

deterministic first passage saturates after repeated identical cycles, whereas fast thermal renewal produces elapsed-time-controlled accumulation. A complete cycle-controlled material model therefore requires additional cycle-evolving physical structure.
