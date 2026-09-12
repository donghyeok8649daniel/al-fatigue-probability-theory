# v31 material compatibility studies

The production LJ/Bessel energy and calibration files are unchanged. All these
fits are research candidates or diagnostic ablations, not adopted Al materials.
All previously inspected interface targets are development data.

1. `corrected`: three inherited radial shapes, coefficient minimax profiles
   for finite q, interface curvature, and both. `main` contains only the
   preserved definition of an early interrupted role-label attempt.
2. `shape_search`, `shape_followup`, `shape_continuation`: bounded existing
   five-positive-shape searches of 40, 160, and 320 objective evaluations.
3. `shape_validation`, `continuation_validation`: full-matrix replay and
   excluded-wavevector tests. Evaluation-budget termination is not optimizer
   convergence.
4. `target_conflicts`: seven fixed-shape target ablations and LP dual witnesses
   after the 320-evaluation continuation. A fixed-shape incompatibility is not
   a proof against the entire nonlinear family.
5. `screening_followup`, `screening_continuation`: additionally vary the
   **already existing** signed density-screening exponent in its inherited
   [-1,1] range; no new interaction term. The separate validation directories
   replay the final returned candidate and evaluate excluded q points.
6. `screening_continuation_conflicts`: final-shape LP ablations/duals. The
   quartic nonnegative bound is active here; do not reuse the earlier
   sign-free certificate without this changed condition.

Final 600-evaluation continuation: eta12.50738498, optimizer max-evaluation
termination, excluded-q maximum relative matrix error210.08%. All seven
exact anchors hold numerically, but the material fails the declared box.
The final candidate is not adopted. See the final section of the v31 theory
document for the actual parameter vector, competing curvatures and limitations.

`summary.json.completed` means the study ran, not that the material passed.
Inspect `optimizer_success`, minimax eta (required <= 1 for the declared
discrepancy box), positive LJ, target residuals, numerical tails and excluded-q
errors together. Positive eigenvalues at the eight tested wavevectors are not
an all-wavevector stability proof. Shape saturation near its numerical upper
bound is not a measured Al parameter. No fatigue, yield or desired event order
is in the objective. No A_c, kinetic mobility or atomic mass is used here.
