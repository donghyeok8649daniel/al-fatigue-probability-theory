# Literature checkpoint: room-temperature high-purity and single-crystal Al fatigue

Status: **EXTERNAL-EVIDENCE CHECKPOINT / PHYSICAL WARNING FOR THE REDUCED NORMAL CLOSURE.**

This note compares the active one-dimensional normal thermal-first-passage closure with published observations. The mathematical closure remains useful as a null/baseline model, but the literature evidence below prevents it from being presented as an established dominant mechanism of room-temperature pure-Al fatigue.

## 1. Frequency signature of the reduced normal closure

In the strict quasistatic, fast-intrawell-equilibration limit, the current reduced normal model predicts

$$
\mathcal H_c\propto\frac1f,
$$

and therefore, for a fixed local initiation probability,

$$
N_p\propto f,
\qquad
\frac{N_p}{f}\approx\mathrm{constant}.
$$

That is a **time-controlled** rare-event signature: increasing frequency gives proportionally more cycles before the same probability is reached.

## 2. High-purity Al frequency evidence

N. H. G. Daniels and J. E. Dorn, *The Effect of Temperature, Frequency, and Grain Size on the Fatigue Properties of High-Purity Aluminum*, ASTM STP 196, pp. 94--110, DOI 10.1520/STP19619570007, reported that at room temperature the fatigue strength of high-purity aluminum was insensitive to frequency over 25--1440 cycles per minute, approximately 0.42--24 Hz. The same study found a stronger time/temperature/frequency coupling at elevated temperature, above roughly 150 C.

Source links:
- ASTM: https://store.astm.org/stp19619570007.html
- ASTM volume metadata identifies STP 196 as a 1956 publication; later indexing also cites this Daniels--Dorn work as 1956/1957-era literature.

### Consequence

Frequency-insensitive fatigue strength at room temperature is **not the signature expected from the strict thermal-only normal closure**, which predicts cycles to a fixed probability to scale approximately with frequency in its quasistatic regime.

This is not a strict one-to-one falsification because the historical study reports fatigue strength of high-purity material rather than the exact local single-crystal first-passage observable defined in this project. It is nevertheless a strong physical warning and a direct motivation for a dedicated frequency sweep in the planned experiment.

## 3. Pure-Al single-crystal crack-initiation evidence

M. Hayashi et al., *Effect of crystal orientation on fatigue crack initiation life in pure aluminum single crystals*, International Journal of Fatigue 156 (2022) 106661, DOI 10.1016/j.ijfatigue.2021.106661, studied 99.99% pure Al single crystals.

The published abstract reports:

- crack-initiation life depends on multiple slip systems rather than only the primary slip system;
- secondary slip plays an important role;
- X-ray micro-beam observations show strongly developed subgrain structure;
- total misorientation at crack initiation is about

$$
5.5\times10^{-2}\ \mathrm{rad},
$$

approximately independent of orientation and shear-stress amplitude in the reported tests;
- the authors interpret this as suggesting crack initiation when dislocation density at the subgrain boundary reaches a critical state.

Source: https://doi.org/10.1016/j.ijfatigue.2021.106661

## 4. Persistent-slip-band and irreversible-slip evidence

T. Zhai, J. W. Martin, and G. A. D. Briggs reported a series of room-temperature Al single-crystal fatigue studies at 20 Hz and resolved shear-stress amplitude 4 MPa.

Relevant papers include:

1. *Fatigue damage in aluminum single crystals---I. On the surface containing the slip Burgers vector*, Acta Metallurgica et Materialia 43 (1995) 3813--3825, DOI 10.1016/0956-7151(95)90165-5.
   - Microcracks, microvoids, macrobands, extrusions, and intrusions were observed.
   - Net irreversible slip in one direction was reported in most persistent slip bands.

2. *Fatigue damage at room temperature in aluminium single crystals---II. TEM*, Acta Materialia 44 (1996) 1729--1739, DOI 10.1016/1359-6454(95)00330-4.
   - Dislocation cell/band structures, dislocation walls associated with persistent slip bands, and many dislocation loops were observed.
   - The loop population was interpreted as an important form of irreversible deformation.

3. *Fatigue damage at room temperature in aluminium single crystals---IV. Secondary slip*, Acta Materialia 44 (1996) 3489--3496, DOI 10.1016/1359-6454(96)00025-0.
   - Secondary slip, extrusions/intrusions, and short cracks were observed at later fatigue stages.
   - Internal stress associated with net irreversible secondary slip and lattice rotation was proposed as important to crack initiation.

These observations provide direct evidence for a slow cycle-evolving structural state in real room-temperature Al single-crystal fatigue.

## 5. Revised interpretation of the current normal thermal closure

The closed equation

$$
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P_b
\right]
=-k_cP_b
$$

remains mathematically valid **under its declared assumptions** and remains useful for:

- a normal-instability null model;
- quantifying how far a purely thermal normal mechanism can go without empirical fatigue laws;
- generating strong falsification signatures;
- testing whether a special normal-dominated regime exists.

However, current literature does **not** justify promoting it as the complete or dominant room-temperature pure-Al fatigue mechanism.

The strongest current physical reading is

$$
\boxed{
\text{fast/reversible normal spacing response}
+
\text{slow slip/dislocation structural evolution}
\to
\text{fatigue initiation}.
}
$$

The slow state is not yet given an active governing equation because introducing an empirical dislocation-density or damage law would violate the project rules.

## 6. Why the previous registry-kink no-go does not eliminate slip as the real mechanism

The previous reduced registry audit used an ideal finite row and found weak post-formation bulk trapping. That result only rejects **long-lived bulk kink-pair trapping in that ideal reduced surface** as the missing memory state.

It does not reject experimentally observed irreversible slip mechanisms involving:

- dislocation escape to a free surface;
- persistent-slip-band surface steps;
- extrusion/intrusion formation;
- interactions among multiple slip systems;
- dislocation walls, loops, and subgrain boundaries;
- residual internal stresses after slip.

These ingredients are absent from the ideal two-row bulk registry calculation.

## 7. New modeling target

Do not force permanent drift into $P_a$ directly. The next physics target is a **derived slow structural coordinate** $\chi$ that modifies the local normal spacing environment or its first-passage boundary:

$$
P_0
\to
\chi(t),\ P_a(a,t\mid\chi)
\to
F_{\mathrm{ci}}(t).
$$

Candidates for the physical meaning of $\chi$ include a topological slip count, residual internal-stress state, or dislocation/subgrain measure. No candidate is active until its evolution follows from mechanics or an independently justified reduced balance.

A particularly promising route is to revisit the spatial registry model with a **finite free surface and defect flux through the boundary**. A kink/dislocation need not remain trapped in the bulk to produce permanent plastic memory: it can migrate out of the modeled region and leave a net slip offset or surface step. This mechanism is qualitatively consistent with the published persistent-slip-band extrusion/intrusion observations, but it has not yet been derived or validated in the repository.

## 8. Current verdict

The literature checkpoint changes the confidence level, not the exact mathematics:

- exact finite-chain probability projection: **retained reference layer**;
- quasistatic structural $P_0$ transport: **retained normal baseline**;
- thermal normal first passage: **retained as a falsifiable null/baseline mechanism**;
- thermal normal first passage as the dominant room-temperature pure-Al fatigue clock: **not supported; strong warning from frequency and single-crystal microstructure evidence**;
- physically derived slow slip/dislocation state: **new primary mechanism target**.

This is a productive falsification result. The project now has a closed null model whose failure signature points directly to the missing physical state instead of to another arbitrary probability closure.
