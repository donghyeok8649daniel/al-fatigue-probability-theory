# desktop-app

End-user desktop packaging and OS-integration branch.

Owned here:
- `app/`: Windows packaging, application metadata, installer/portable layout, file-association design.
- `docs/`: `.ftgsim`, desktop application, distribution/platform documentation.

Not owned here: FEM/simulation implementation, theory, manuscript sources, or fatigue-tester firmware/hardware.

The application consumes stable UI/runtime functionality from `numerical-fem`; combined validation happens in `integration`.
