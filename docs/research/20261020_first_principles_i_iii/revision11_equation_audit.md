# Revision 11 equation audit: foundations and the two bridges

20 October 2026. The authority for this audit is the **downloaded actual
Papeeria Revision 11**, not an older Git branch's manuscript. Its exact
original source is preserved under `paper/rev11_original/`; SHA-256 and base
commit are in `source_provenance.json`. Numbering below follows the original
display-equation sequence. The LaTeX labels/expressions are the stable keys.
The full 55-display downstream audit requested by the long directive is
**not claimed complete**: this first delivery concentrates on PART I–III.

Statuses: **identity** = established under specified mathematical hypotheses;
**definition** = coordinate/model declaration; **conditional closure** =
consistent given additional physical assumptions; **open bridge** = those
assumptions are not established for the same material/coordinate. These are
not interchangeable with "incorrect."

| Original display | Expression / source label | Verdict and action |
|---|---|---|
| 1 | $S_p=\sum_j[a^2+((j+1/2)b-s)^2]^{-p}$ | Identity/definition for (a>0,p>1/2). Poisson transform with all modes independently derived and checked. It is a row sum. |
| 2 | $W_{\rm row}=4\epsilon(\sigma^{12}S_6-\sigma^6S_3)$ | Exact chosen pair model. A single atom against a row; not a bulk or counted-cut energy and not automatically a PMF. Retain. |
| 3 | $h=(a-a_0)+\chi s, e=h/a_0, G_{\rm row}=W_{\rm row}-fh$ | Definition of observable and projected work. It generates $(f,\chi f)$; pure normal traction requires zero tangential component, e.g. $\chi=0$. It supplies no kinetic derivation. |
| 4 | `eq:r10-bulk-count` | Correct pair counting with declared inversion symmetry, common geometry/gauge and self exclusion. Half the raw neighbor field is a pair energy, not an embedding density. |
| 5 | `eq:r10-cut-count`, $c_k=k$ | Correct for the declared infinite half-stack rigid cut; finite stacks need actual counts. Require the correct row/plane kernel, registries, zero mode and convergence. Independently verified for an explicit zigzag stack; not validation of every FCC orientation. |
| 6 | `eq:r11-selected-energy` | Correct dimensional map and work definition. **Open bridge:** $U_\Gamma^*=W_\Gamma/E_\Gamma$ is a rigid energetic input, whereas diffusion's thermodynamic potential is the PMF under the chosen measure. Replace only after (H4), or retain as an explicit approximation. |
| 7 | `eq:r11-counted-derivatives` | Correct fixed-geometry differentiation with local uniform convergence. No second $k$ for cut opening. Sitewise embedding derivatives are separate. Confirmed by direct/reciprocal and finite-difference checks in the declared example. |
| 8 | `eq:r11-chain-work` | Correct affine chain rule at fixed area/scales. These are derivatives of the adopted energy, not proof that its force equals the finite-temperature conditional mean force. |
| 9 | `eq:retained-fpe` | **Conditional closure.** Conservation and constant-mobility current are consistent with the assumed Itô diffusion. A material derivation additionally needs the PMF bridge, bath/Markov limit, overdamping, compatible load/temperature and domain/boundary contract. Diagonal mobility needs symmetry/evidence. |
| 10 | `eq:r11-cut-sde` | Equivalent interior Itô SDE to display 9 under constant positive mobility and temperature. For state-dependent $M$, add $k_BT\nabla\cdot M$ in Itô drift. Neither white noise nor inertia removal follows from the lattice identity. |
| 11 | `eq:r11-gap-density` | Correct affine pushforward with $\det D_\Gamma>0$ and corresponding times; densities are relative to the stated coordinate measure. This is unrelated to a missing coarea term for a different nonlinear microscopic coordinate. |

Rev.11 explicitly states after display 11 that the diffusion is assumed and
that inertia removal, memory loss and Brownian forcing lack an accepted
same-coordinate atomistic limit/calibration. The audit **retains that honest
qualification** and makes the two missing bridges constructive. The abstract
word "drives" should be qualified as a candidate input to an assumed closure.

The physical scaling written after display 10,

$$
M_{\rm phys}=D_\Gamma M^*D_\Gamma^T/(t_\Gamma E_\Gamma),
$$

is algebraically consistent. It declares the conversion if $t_\Gamma$ is
known; it does not identify $t_\Gamma$. Existing v31 finite-band Al
plane-coordinate estimates are relevant evidence but leave zero drag within
their error sets and do not establish the same local PMF/thermal map.

Downstream first-opening balances, SG fluxes, cell interactions, observables,
and dissipation claims keep their existing stated assumptions. They are
not re-certified by this audit. In particular, row-specific bounds do not
transfer to the counted cut, and finite-box normalization is not proof of a
global continuous equilibrium or physical absorbing event.

## Recommended manuscript action

1. Preserve Revision 11 exactly as a source backup.
2. Add a foundations audit notice and standalone PART I–III companion,
   retaining current equation labels and conditional model status.
3. Use $A=U-k_BT\log Z_\perp+C$ as the explicit first bridge and the exact
   memory reduction as the intermediate dynamical object.
4. Promote the existing FPE to a derived material closure only after the
   matching-coordinate validation gates pass. No production energy or
   mobility is changed by this delivery.
