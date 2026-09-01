# integration

Cross-module integration branch. Unlike the dedicated ownership branches, this branch intentionally contains the combined project tree needed to validate interfaces and end-to-end behavior.

Integrated modules:
- `theory/`, `libraries/`: theory-core
- `simulations/`, `fem1d/`, `tests/`, `results/`, `examples/`: numerical-fem
- `app/`: desktop-app packaging/OS integration
- `paper/`, `output/`, `research/`: paper-manuscript
- `fatigue_tester/`: physical tester implementation

This is the only development branch where those domains are expected to coexist. Stable, reviewed checkpoints may flow from here to `main`.
