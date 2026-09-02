# Registry symmetry stability under normal-only cyclic loading

## Scope

This audit uses only the active reduced row/layer potential and the existing
normal-only deterministic spatial chain.  No FCC geometry, imposed registry
force, damping, random noise, Boltzmann distribution, or Fokker--Planck closure
is introduced.

The question is whether normal cyclic loading can destabilize the symmetric
registry coordinate internally.

## 1. Static curvature ordering

For the normalized diagnostic parameters

$$
m=12,\quad n=6,\quad b=\sigma_{LJ}=\epsilon_{LJ}=1,
$$

with symmetric stable registry phase $s_0/b=0.5$, direct double-sum
curvatures give

$$
\boxed{a_a^*=1.130690887\quad(U_{aa}=0)}
$$

and

$$
\boxed{a_s^*=1.264187982\quad(U_{ss}=0).}
$$

Thus the static registry-curvature loss occurs about $11.81\%$ beyond the
normal-curvature loss.  At the normal marginal point,

$$
U_{ss}(a_a^*,s_0)\approx2.38763>0.
$$

Therefore static normal opening does **not** make the registry minimum unstable
before the normal opening instability in this parameter set.

The ordering converges from $(k_{max},p_{max})=(20,50)$ through $(200,500)$;
see `curvature_convergence.csv`.

## 2. Normal-only chain modulates registry stiffness without moving registry

In the existing five-cycle normal-only spatial-chain test, $s_i=s_0$ remains
invariant to numerical precision, but the boundary-generated $a_i(t)$ field
changes the linear registry stiffness

$$
K_{s,i}(t)=U_{ss}(a_i(t),s_0).
$$

Using the same numerical settings as the earlier spatial-chain test, the
observed $K_s$ range over five cycles is approximately

$$
15.47\lesssim K_s\lesssim34.99,
$$

so it remains positive throughout that run.

Hence there is no static registry instability in the tested trajectory.

## 3. The coefficient is not exactly period-one

Because the chain is conservative and undamped, boundary-generated waves do
not relax to a unique steady period-one state.  Comparing successive cycle
histories of $K_s$ gives relative RMS mismatches

- cycle 1 to 2: $5.21\%$;
- cycle 2 to 3: $5.69\%$;
- cycle 3 to 4: $5.92\%$;
- cycle 4 to 5: $6.25\%$.

Therefore calling the five-cycle coefficient exactly periodic and applying a
strict Floquet interpretation would be unjustified.  The numerical test below
is reported instead as a **finite-time variational amplification** calculation.

## 4. Finite-time registry perturbation equation

An infinitesimal registry perturbation $\xi_i$ around the exact symmetric state
obeys

$$
\boxed{
\mu_s\ddot\xi_i+K_{s,i}(t)\xi_i=0.
}
$$

No seed amplitude is selected.  The code integrates the two fundamental
solutions, builds the transfer matrix, and reports the largest singular value
in the initial quadratic-energy scaling

$$
\left(\sqrt{K_s(0)}\,\xi,\sqrt{\mu_s}\,\dot\xi\right).
$$

Thus the result measures how strongly *any existing infinitesimal physical
registry perturbation* could be amplified by the computed normal-only history.

## 5. Amplification depends strongly on the unresolved registry inertia

Selected five-cycle results are:

| $\mu_s$ | max scaled singular amplification | max transfer spectral radius |
|---:|---:|---:|
| 1 | 1.0569 | 1.0359 |
| 100 | 1.1071 | 1.0000 |
| 500 | 1.2432 | 1.0000 |
| 650 | 1.8284 | 1.0000 |
| 700 | 2.0326 | 1.0000 |
| 750 | 2.1589 | 1.8293 |
| 800 | 2.2082 | 2.2076 |
| 812.5 | **2.2101** | 2.1835 |
| 850 | 2.1949 | 1.8721 |
| 900 | 2.1371 | 1.0000 |
| 1000 | 1.9512 | 1.0000 |
| 2000 | 1.1167 | 1.0000 |

The strongest point in this coarse/refined diagnostic lies near
$\mu_s\approx8.1\times10^2$ in the present normalized units.  This agrees in
scale with the small-harmonic principal-resonance estimate

$$
\mu_s\sim\frac{U_{ss}(a_0,s_0)}{(\Omega/2)^2}
\approx8.40\times10^2.
$$

This is **not** a calibrated material prediction because $\mu_s$ is unresolved.
It demonstrates only that normal cyclic modulation of the existing $U_0$
contains a real parametric-amplification channel for suitable inertial scaling.

## 6. What this does and does not solve

The result eliminates one candidate and keeps another:

### Rejected for this normalized model

$$
\text{normal opening}\to U_{ss}<0\to\text{static slip instability before crack}.
$$

Normal curvature is lost first.

### Still viable internally

$$
\text{normal cyclic wave}\to U_{ss}(a_i(t),s_0)\text{ modulation}
\to\text{amplification of a nonzero registry seed}.
$$

However exact symmetry remains invariant:

$$
\xi_i(0)=\dot\xi_i(0)=0\Rightarrow\xi_i(t)=0.
$$

The reduced model therefore still needs a **physical source of the initial
registry seed** before it can predict a nontrivial $P(s,t)$.  The seed must not
be replaced by an arbitrary Gaussian/noise law merely to make the model work.
Possible sources must be treated explicitly later (finite-temperature
microstate, defect/asymmetry, previous-cycle residual registry motion, etc.).

## 7. Consequence for the current research mainline

For the ideal exact-symmetry $T=0$ reduced state, the safe statement is

$$
P(s,t)=\delta(s-s_0)
$$

under pure normal loading.

The normal probability density $P(a,t)$ remains mechanically generated by
spatial wave diversity.  A genuinely two-dimensional $P(a,s,t)$ becomes active
only after a physically justified registry seed is introduced and then
amplified/suppressed by the deterministic stability equation above.

This means the immediate theory should not invent $q_s$.  The next unresolved
quantity is the physical mapping of $\mu_s$ and the admissible source of the
registry seed.
