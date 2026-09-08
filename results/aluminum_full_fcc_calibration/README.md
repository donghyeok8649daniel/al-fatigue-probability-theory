# Calibration evidence: current versus historical

Use **audited_v2/** for the current independent-elastic-mode audit.
It actually reruns deterministic fits and records their failures as well as
their numerical residuals. No candidate is production-valid.

The CSVs directly in this parent directory are retained historical results
from an underidentified normal/biaxial experiment. Their rank-five and
bulk-exact labels are superseded: they did not constrain C11-C12. They must
not be cited as an independently validated aluminum calibration.

The original handoff is preserved in legacy_handoff_2026-09-09.md to make
that correction traceable. Read the repository CURRENT_WORK_HANDOFF.md for
the current state, not the archived handoff.

The archival replay writes to legacy_replay/ and cannot replace audited_v2.
See solver_v1/ALUMINUM_FULL_FCC_CALIBRATION.md for equations and interpretation.
